"""Workbench Entry 170: can injected neutrons solve the lithium problem? A 2026 proposal tested on deuterium.

Arbuzova & Dolgov (2026), arXiv:2606.03316, propose that primordial black holes of about 10^13 g, evaporating at t ~ 1.5e11 s
(plasma at ~1.8 eV), release baryons whose neutrons turn 7Li into 8Li (7Li + n -> 8Li + gamma; 8Li -> 2 4He), cutting lithium to the
observed third. Their Eq. (30) puts the depletion exponent K = sigma n_B tau c (v_Lin / v_Hepbar) ~ 1.

Two questions here, both with no data, only nuclear physics:
 (1) Kinetics. A free neutron lives 879 s, not tau_cosm, and the 7Li(n,gamma) cross-section follows 1/v: 5.4 mb at 1.8 eV, not the
     paper's 0.2 b. The factor c in Eq. (30) comes from the paper's helium bound, n_B <= 1e3 c / v_Hepbar. What is K then?
 (2) Deuterium. Whatever the source, a neutron that is captured is captured mostly by protons (n + p -> D + gamma). For neutrons
     injected at any epoch, how much must be injected to cut A = 7 to a third, and what does that do to D/H, measured to 1.2%?
     (a) during and just after nucleosynthesis: the network of Entry 53 (bbn_net.py) with a neutron source added;
     (b) later, after the network's end: capture ratios from thermal cross-sections, which follow 1/v, so their ratios hold at any
         temperature below a few keV. Below about 25 keV 7Be(n,p)7Li only converts 7Be into 7Li, which then survives; at the paper's
         epoch A = 7 is 7Li and only 7Li(n,gamma)8Li removes it.
Thermal (0.0253 eV) cross-sections: n + p -> D + gamma 0.3326 b; 3He(n,p)3H 5333 b; 7Be(n,p)7Li 38,400 b; 7Li(n,gamma)8Li 0.0454 b;
D(n,gamma)3H 0.000508 b (ENDF/B-VIII.0 thermal values; Koehler et al. 1988 for 7Be).
Run from code/:  ../.venv/bin/python neutron_li.py      (about 10 minutes: about 60 network runs)
"""
import os, sys, json, inspect, textwrap, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import bbn_net as bn

# ---- the network with a neutron source: A neutrons per baryon, injected as a Gaussian in ln T9 (width 0.15) about T9c
src = textwrap.dedent(inspect.getsource(bn.run))
src = src.replace("def run(eta10=6.10,", "def run_inj(A_inj=0.0, T9c=0.3, wid=0.15, eta10=6.10,", 1)
src = src.replace("d[1+N]=-w-r1+r3-r5+r6-r10",
                  "d[1+N]=-w-r1+r3-r5+r6-r10+A_inj*np.exp(-0.5*(np.log(T9/T9c)/wid)**2)/(np.sqrt(2*np.pi)*wid)*(-dT/T9)", 1)
assert "A_inj*np.exp" in src
exec(src, bn.__dict__)
run_inj = bn.run_inj

ETA = 6.124                                   # eta10 from omega_b = 0.02237, as in Entry 113
T9F = 0.006
D_OBS, D_ERR = 2.527e-5, 0.030e-5             # Cooke et al. 2018
LI_OBS = 1.6e-10                              # Spite plateau, about a third of the prediction


def checked(r):
    assert r['t_end'] > 1e4, 'integration stopped early'
    return r


def solve_A(T9c, base, target=1/3.):
    """Neutrons per baryon that bring total A = 7 to target x standard, by secant/bisection in log A."""
    lo, hi = None, None
    A = 1e-5
    hist = []
    for it in range(14):
        r = checked(run_inj(A_inj=A, T9c=T9c, eta10=ETA, T9f=T9F))
        ratio = r['Li7'] / base['Li7']
        hist.append((A, ratio, r))
        if abs(np.log(ratio / target)) < 0.01:
            return A, r, hist
        if ratio > target:
            lo = (A, ratio)
        else:
            hi = (A, ratio)
        if lo and hi:                                                     # interpolate in log A against log ratio
            la, lb = np.log(lo[0]), np.log(hi[0]); ra, rb = np.log(lo[1]), np.log(hi[1])
            A = float(np.exp(la + (np.log(target) - ra) * (lb - la) / (rb - ra)))
        elif lo:
            A *= 10.0 if ratio > 0.9 else 3.0
        else:
            A /= 3.0
        if A > 0.01:                                                      # 1% of baryons as neutrons: far past anything allowed
            return None, r, hist
    return A, r, hist


def main():
    t0 = time.time()
    out = dict(eta10=ETA, T9_end=T9F, method='bbn_net.py (Entry 53) with a neutron source; Gaussian in ln T9, width 0.15')
    base = checked(run_inj(A_inj=0.0, eta10=ETA, T9f=T9F))
    out['standard'] = dict(Yp=base['Yp'], D=base['D'], He3=base['He3'], Li7=base['Li7'], Be_frac=base['Be_frac'])
    print('standard', {k: (round(v, 6) if k == 'Yp' else float('%.4g' % v)) for k, v in out['standard'].items()}, flush=True)
    rows = []
    for T9c in (1.0, 0.7, 0.5, 0.4, 0.35, 0.3, 0.25, 0.15, 0.05):
        A, r, hist = solve_A(T9c, base)
        if A is None:                                     # no dose reaches a third: report the run that came closest
            r = min(hist, key=lambda h: h[1])[2]
        row = dict(T9c=T9c, keV=round(T9c * 86.17, 2), solved=A is not None, A_per_baryon=A, Li7=r['Li7'], Li_ratio=r['Li7'] / base['Li7'], D=r['D'],
                   D_rise=r['D'] / base['D'] - 1.0, D_pull_vs_obs=(r['D'] - base['D']) / D_ERR, He3=r['He3'], He3_rise=r['He3'] / base['He3'] - 1.0,
                   Yp=r['Yp'], dYp=r['Yp'] - base['Yp'], runs=len(hist),
                   scan=[[float('%.4g' % h[0]), round(h[1], 4), round(h[2]['D'] / base['D'] - 1.0, 4)] for h in hist])
        rows.append(row)
        print('T9c %.2f (%5.1f keV): A %.3g per baryon, Li x%.3f, D +%.1f%% (%.0f sigma), 3He %+.1f%%, Yp %+.4f  [%d runs, %d s]' % (
            T9c, T9c * 86.17, A or float('nan'), row['Li_ratio'], 100 * row['D_rise'], row['D_pull_vs_obs'], 100 * row['He3_rise'], row['dYp'],
            len(hist), time.time() - t0), flush=True)
    out['network'] = rows

    # ---- (b) late injection, after the network: capture competition with thermal (1/v) cross-sections.
    # Once 7Li + p has stopped (below about 25 keV, as the network shows), 7Be(n,p)7Li only turns 7Be into 7Li and leaves A = 7 intact.
    # The paper's epoch (1.8 eV, z ~ 7,700) comes after 7Be has captured an electron and decayed (z ~ 3e4), so A = 7 is 7Li and the only
    # neutron channel that removes it is 7Li(n,gamma)8Li, the paper's own.
    s = dict(np_=0.3326, He3=5333.0, Li7=0.0454, D=0.000508)                          # barns
    ab = dict(He3=base['He3'], D=base['D'], A7=base['Li7'])                            # per H, after nucleosynthesis
    Sigma = s['np_'] + ab['He3'] * s['He3'] + ab['D'] * s['D'] + ab['A7'] * s['Li7']  # per H, barns
    N_abs = np.log(3.0) * Sigma / s['Li7']                                             # absorbed neutrons per H to cut A = 7 by 3
    out['late_thermal'] = dict(channel='7Li(n,gamma)8Li', neutrons_absorbed_per_H=N_abs, share_to_protons=s['np_'] / Sigma,
                               D_made_per_H_first_captures=N_abs * s['np_'] / Sigma,
                               note='Destroying two thirds of 7Li by neutron capture needs about this many neutrons absorbed per hydrogen atom; '
                                    'most are captured by protons, so the gas would be turned largely into deuterium. 7Be(n,p)7Li conserves A = 7.')
    # ---- (1) the paper's exponent, redone. K = <sigma v>(7Li) x (number density of neutrons ever present) x (time each one lives)
    v_th = 2200.0e2                                        # cm/s at 0.0253 eV; sigma v is constant for 1/v
    sv_li = s['Li7'] * 1e-24 * v_th                        # cm3/s
    n_B, tau_cosm, tau_n = 1e3, 1.54e11, 879.0
    K_paper = 2e-25 * n_B * tau_cosm * 3e10
    K_true = sv_li * n_B * tau_n
    # cosmic baryon density at 1.8 eV, and whether a neutron is captured by a proton before it decays
    Tcmb0 = 2.7255 * 8.617e-5                              # eV
    nb_now = 2.51e-7                                       # cm^-3 (omega_b 0.02237)
    nb = nb_now * (1.8 / Tcmb0) ** 3
    gp = 0.75 * nb * s["np_"] * 1e-24 * v_th               # capture rate on protons per neutron, per s (hydrogen nuclei: 0.75 per baryon)
    # Two readings of the paper's Eq. (25) and (30). Its helium bound is n_B <= 1e3 / (v_Hepbar / c) cm^-3, and Eq. (30) is
    # K = sigma x 1e3 cm^-3 x tau x c x (v_Lin / v_Hepbar), with "v_Lin / v_Hepbar ~ 1" (checked against the paper's text, 21 Nov 2026).
    #  A: v_Hepbar = c, so n_B = 1e3 cm^-3 and the lithium term moves at c too (K_paper as printed).
    #  B: both velocities thermal at 1.8 eV, ratio 1: n_B = 1e3 c / v_th(1.8 eV), which exceeds the cosmic baryon density.
    v18 = v_th * (1.8 / 0.0253) ** 0.5                     # cm/s, neutron at 1.8 eV
    nB_B = 1e3 * 3e10 / v18
    K_B = sv_li * nB_B * tau_n
    out['kinetics'] = dict(K_paper=K_paper, K_redone=K_true, factor=K_paper / K_true,
                           pieces=dict(paper_sigma_times_c_over_thermal_sigma_v=(2e-25 * 3e10) / sv_li, tau_cosm_over_tau_n=tau_cosm / tau_n),
                           reading_B=dict(v_thermal_1p8eV_cm_s=v18, n_B=nB_B, K=K_B, factor=K_paper / K_B,
                                          pieces=dict(paper_sigma_over_true_sigma_at_1p8eV=2e-25 / (s['Li7'] * 1e-24 * v_th / v18),
                                                      tau_cosm_over_tau_n=tau_cosm / tau_n),
                                          antibaryons_per_cosmic_baryon=nB_B / nb),
                           baryon_density_at_1p8eV=nb, neutron_capture_on_p_per_s=gp, fraction_captured_before_decay=gp * tau_n / (1 + gp * tau_n),
                           note='Corrected 21 Nov 2026 after re-reading the paper: the factor c in Eq. (30) comes from its helium bound, n_B <= 1e3 c / v_Hepbar, '
                                'not from a lithium velocity set to c. Reading A keeps n_B = 1e3 cm^-3 (K_redone); reading B takes both velocities thermal, '
                                'which needs more antibaryons than there are cosmic baryons. Either way K falls short by 6e9 to 1e14: the neutron '
                                'lifetime (1.75e8 in both) and the 7Li(n,gamma) cross-section, 0.2 b in the paper against 5.4 mb at 1.8 eV by the 1/v law. '
                                'All evaporated baryons are counted as neutrons (an upper bound).')
    # ---- (3) the energy of the antibaryons. Each evaporated antibaryon annihilates with a cosmic baryon (the paper's section 5),
    # releasing 2 x 0.938 GeV, at z ~ 7,660. Energy that late cannot thermalise: it leaves a Compton y distortion, y ~ (delta rho / rho) / 4.
    # The neutrinos from pion decays carry about half of it away; f_em = 0.5 is assumed. FIRAS: |y| < 1.5e-5 (95%, Fixsen et al. 1996).
    T = 1.8; rho_g = 2.7012 * T * 410.7 * (T / Tcmb0) ** 3       # eV cm^-3: mean photon energy 2.701 kT times photon density at T
    E_ann = n_B * 2 * 0.938e9 * 0.5                               # eV cm^-3, electromagnetic share
    out['cmb_distortion'] = dict(antibaryons_per_baryon=n_B / nb, delta_rho_over_rho=E_ann / rho_g, y=E_ann / rho_g / 4, firas_y_limit=1.5e-5,
                                 over_limit=E_ann / rho_g / 4 / 1.5e-5, f_em=0.5)
    print(json.dumps(out['late_thermal'], indent=1)); print(json.dumps(out['kinetics'], indent=1)); print(json.dumps(out['cmb_distortion'], indent=1))
    out['seconds'] = round(time.time() - t0)
    json.dump(out, open(os.path.join(HERE, 'neutron_li.json'), 'w'), indent=1)
    print('written neutron_li.json', out['seconds'], 's')


if __name__ == '__main__':
    main()
