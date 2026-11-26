"""Workbench Entry 170, addition of 22 November 2026, revised 25 November: can PROTONS from evaporating black holes deplete 7Li?

The version of record of Arbuzova & Dolgov (Eur. Phys. J. C 85, 1362, 2025; arXiv:2606.03316, posted 2 June 2026 after it, with the
same equation numbers) adds radiative proton capture, 7Li + p -> 8Be + gamma, beside neutron capture, "possibly less efficient due to
the Coulomb barrier", and gives no estimate for it.  Three cases, no data, only kinematics and cross-sections:

 (a) Thermal protons at T = 1.8 eV.  The Coulomb barrier: non-resonant rate for 7Li(p,alpha)4He, the strong branch of 7Li + p, with
     S(0) = 55 keV b (Lattuada et al. 2001, ApJ 562, 1076: 55 +- 3).  The radiative branch the paper names passes the same barrier.
     Evaporated protons that have thermalised are also indistinguishable from the ambient ones, which outnumber them.
 (b) Energetic protons below the break-up threshold of 4He (23.0 MeV in the laboratory, from Q = -18.35 MeV for 4He(p,d)3He).
     They slow down on the plasma electrons; the column they cross before stopping, times the cross-section, is the share of 7Li
     destroyed per energetic proton per ambient hydrogen atom.
     Stopping on free electrons: -dE/dx = (4 pi e^4 n_e / m_e v^2) ln(2 m_e v^2 gamma^2 / hbar omega_p).
     Cross-section: flat values, and a ceiling that needs no compilation, the strong-absorption limit of every reaction on 7Li,
     pi (R + lambda-bar)^2 (1 - V_C / E_cm), with a generous R = 1.4 (7^(1/3) + 1) fm.
 (c) Protons above that threshold.  The same fluence acts on 7Li and on 4He, so the destroyed shares are tied by the ratio r of the
     destruction cross-sections, whatever the number of protons: 1 - f_He = (1 - f_Li)^r.  For the GeV protons a 1 GeV black hole
     emits, r is taken from the Glauber-model inelastic cross-sections tabulated by the Particle Data Group (Atomic and Nuclear
     Properties of Materials, 2024: nuclear interaction length 71.0 g/cm^2 for helium, 71.3 g/cm^2 for lithium, A = 6.94; for
     200 GeV/c, "reasonably independent of momentum").  The share of a 1 GeV proton's path spent near the threshold is computed,
     since r is smaller there.  Deuterons per 4He broken up: only a floor of 1 in 100 is used.
 (d) For the photon route of the paper's Section 6 (cooler black holes, after recombination): photodisintegration thresholds from
     the AME2020 atomic masses (Wang et al. 2021, Chin. Phys. C 45, 030003).
Run from code/:  python3 neutron_li_protons.py      (a second)          Output: neutron_li_protons.json
"""
import json, math, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))

# ---- the plasma at T = 1.8 eV (Eq. 13 of the paper)
T_eV = 1.8; kB = 8.617333e-5; T_K = T_eV / kB
eta = 6.1e-10; Yp = 0.2467
n_gamma = 410.7 * (T_K / 2.7255) ** 3
n_b = eta * n_gamma; n_H = (1 - Yp) * n_b; n_He = Yp / 4 * n_b; n_e = n_H + 2 * n_He      # fully ionised at z = 7,700
tau_cosm = 1.54e11
out = dict(plasma=dict(T_eV=T_eV, z=T_K / 2.7255 - 1, n_b=n_b, n_H=n_H, n_He=n_He, n_e=n_e))

# ---- (a) thermal 7Li(p,alpha)4He
A1, A2 = 1.007276, 7.014358; mu = A1 * A2 / (A1 + A2); Z12 = 3
T6 = T_K / 1e6
tau = 42.46 * (Z12 ** 2 * mu / T6) ** (1 / 3)                      # 3 E_0 / kT
E0_keV = tau * T_eV / 3 / 1e3
S0 = 55.0                                                          # keV b
log10_sv = math.log10(7.20e-19 / (mu * Z12) * S0 * tau ** 2) - tau / math.log(10)     # cm^3/s
log10_K = log10_sv + math.log10(n_H * tau_cosm)
out['thermal'] = dict(gamow_peak_keV=E0_keV, tau=tau, log10_sigma_v_cm3_s=log10_sv, log10_K=log10_K)

# ---- (b) energetic protons below the 4He break-up threshold
MP, ME = 938.272, 0.510999                                         # MeV
K4 = 4 * math.pi * (2.8179403e-13) ** 2 * ME                       # 4 pi r_e^2 m_e c^2 = 5.10e-25 MeV cm^2
hbar_wp_eV = 6.582120e-16 * 5.64146e4 * math.sqrt(n_e)             # plasma frequency of the electrons
E_th = 18.353 * (A1 + 4.001506) / 4.001506                         # laboratory threshold of 4He(p,d)3He, MeV
def stopping(E, n=20000):
    """energies (MeV) and electron column per MeV, dN_e/dE (cm^-2 MeV^-1), on a grid up to E"""
    e = np.linspace(0.1, E, n); g = 1 + e / MP; b2 = 1 - 1 / g ** 2
    lnL = np.log(2 * ME * 1e6 * b2 * g ** 2 / hbar_wp_eV)
    return e, b2 / (K4 * lnL), lnL
def integ(y, x):
    return float(np.sum((y[1:] + y[:-1]) / 2 * np.diff(x)))
e, dNe, lnL = stopping(E_th)
Ne = integ(dNe, e); lnL_mid = float(lnL[len(e) // 2])
N_H = Ne * n_H / n_e
# strong-absorption ceiling on every reaction with 7Li (fm^2 -> cm^2), along the slowing-down path
R_fm = 1.4 * (A2 ** (1 / 3) + 1.0)
Ecm = e * A2 / (A1 + A2)
lam_fm = 197.327 / np.sqrt(2 * mu * 931.494 * Ecm)
VC = 1.44 * Z12 / R_fm
sig_max = np.pi * (R_fm + lam_fm) ** 2 * np.clip(1 - VC / Ecm, 0, None) * 1e-26
per_max = integ(sig_max * dNe, e) * n_H / n_e
ratioA = 1e3 / n_H                                                 # reading A of the paper's Eq. (25): 1e3 evaporated baryons per cm^3
sub = dict(threshold_MeV=E_th, coulomb_log_mid=lnL_mid, hbar_omega_p_eV=hbar_wp_eV, column_e_cm2=Ne, column_H_cm2=N_H,
           evaporated_per_H_reading_A=ratioA, cases=[])
for sig_b, label in ((0.1, 'flat 0.1 b'), (0.5, 'flat 0.5 b'), (1.0, 'flat 1 b')):
    per = N_H * sig_b * 1e-24                                      # share of 7Li destroyed per energetic proton per ambient H
    sub['cases'].append(dict(sigma_b=sig_b, label=label, destroyed_per_proton_per_H=per, destroyed_reading_A=per * ratioA,
                             protons_per_H_for_ln3=math.log(3) / per, shortfall_reading_A=math.log(3) / per / ratioA))
sub['cases'].append(dict(sigma_b=per_max / N_H * 1e24, label='strong-absorption ceiling, R = %.2f fm' % R_fm,
                         sigma_max_at_10_and_23_MeV_b=[float(np.interp(10, e, sig_max)) * 1e24, float(sig_max[-1]) * 1e24],
                         destroyed_per_proton_per_H=per_max, destroyed_reading_A=per_max * ratioA,
                         protons_per_H_for_ln3=math.log(3) / per_max, shortfall_reading_A=math.log(3) / per_max / ratioA))
out['below_threshold'] = sub

# ---- (c) protons above that threshold: lithium against helium
NA = 6.02214076e23
sig_He = 4.002602 / (NA * 71.0); sig_Li = 6.94 / (NA * 71.3)        # cm^2, PDG Glauber inelastic, 200 GeV/c
r_pdg = sig_He / sig_Li
eG, dNG, _ = stopping(1000.0, 200000)                              # a 1 GeV proton slowing on the electrons only
tot = integ(dNG[eG >= E_th], eG[eG >= E_th])
share = {str(E1): integ(dNG[(eG >= E_th) & (eG <= E1)], eG[(eG >= E_th) & (eG <= E1)]) / tot for E1 in (50, 100, 200)}
DH_obs, DH_err = 2.527e-5, 0.030e-5
He_per_H = n_He / n_H
above = dict(He_per_H=He_per_H, sigma_inel_He_mb=sig_He * 1e27, sigma_inel_Li_mb=sig_Li * 1e27, ratio_pdg=r_pdg,
             path_share_between_threshold_and_E_MeV=share, cases=[])
for r in (0.2, 0.4, r_pdg):                                        # sigma_He / sigma_Li
    fHe = 1 - 3 ** (-r)
    above['cases'].append(dict(sigma_He_over_sigma_Li=r, He_destroyed_for_Li_to_a_third=fHe,
                               D_over_H_added_one_deuteron_per_100=fHe * He_per_H * 0.01))
# the other way round: the most 7Li that can go if D/H may rise by 3 sigma, with 1 deuteron per 100 break-ups
fHe_max = 3 * DH_err / (He_per_H * 0.01)
above['Li_destroyed_at_most'] = dict(deuterons_per_breakup=0.01, He_destroyed_at_most=fHe_max,
                                     Li_at_ratio_pdg=1 - (1 - fHe_max) ** (1 / r_pdg), Li_at_ratio_0p2=1 - (1 - fHe_max) ** (1 / 0.2))
out['above_threshold'] = above

# ---- (d) photodisintegration thresholds, MeV (atomic masses in u, AME2020)
U = 931.49410242
m = {'n': 1.00866491595, 'H1': 1.00782503223, 'H2': 2.01410177812, 'H3': 3.01604928132, 'He3': 3.01602932197,
     'He4': 4.00260325413, 'Li6': 6.0151228874, 'Li7': 7.0160034366, 'Be7': 7.016928717}
out['photon_thresholds_MeV'] = {
    'D(g,n)p': (m['H1'] + m['n'] - m['H2']) * U,
    '7Li(g,t)4He': (m['He4'] + m['H3'] - m['Li7']) * U,
    '7Li(g,n)6Li': (m['Li6'] + m['n'] - m['Li7']) * U,
    '7Be(g,3He)4He': (m['He4'] + m['He3'] - m['Be7']) * U,
    '4He(g,p)3H': (m['H3'] + m['H1'] - m['He4']) * U}

json.dump(out, open(os.path.join(HERE, 'neutron_li_protons.json'), 'w'), indent=1)
p = out['plasma']; t = out['thermal']
print('plasma at 1.8 eV: z = %.0f, n_b = %.3g, n_H = %.3g, n_He = %.3g, n_e = %.3g cm^-3' % (p['z'], p['n_b'], p['n_H'], p['n_He'], p['n_e']))
print('(a) thermal: Gamow peak %.3f keV, 3E0/kT = %.1f, log10 <sigma v> = %.1f, log10 K = %.1f' % (t['gamow_peak_keV'], t['tau'], t['log10_sigma_v_cm3_s'], t['log10_K']))
print('(b) below %.1f MeV: Coulomb log %.1f, column %.3g electrons = %.3g H per cm^2; evaporated per H (reading A) %.4f' % (E_th, lnL_mid, Ne, N_H, ratioA))
for c in sub['cases']:
    print('    %s (mean %.2f b): 7Li destroyed %.2e per proton per H; reading A %.1e; protons per H for a third: %.0f (%.0f times reading A)' % (
        c['label'], c['sigma_b'], c['destroyed_per_proton_per_H'], c['destroyed_reading_A'], c['protons_per_H_for_ln3'], c['shortfall_reading_A']))
print('(c) PDG: sigma_inel 4He %.1f mb, Li %.1f mb, ratio %.3f; share of a 1 GeV proton path from threshold to 50/100/200 MeV: %s' % (
    above['sigma_inel_He_mb'], above['sigma_inel_Li_mb'], r_pdg, ', '.join('%.4f' % v for v in share.values())))
for c in above['cases']:
    print('    ratio %.2f: 4He destroyed %.1f%%; D/H added %.2e at 1 deuteron per 100 break-ups; measured D/H %.3e' % (
        c['sigma_He_over_sigma_Li'], 100 * c['He_destroyed_for_Li_to_a_third'], c['D_over_H_added_one_deuteron_per_100'], DH_obs))
v = above['Li_destroyed_at_most']
print('    D/H + 3 sigma, 1 deuteron per 100: 4He at most %.1e, 7Li at most %.1e (PDG ratio) or %.1e (ratio 0.2)' % (v['He_destroyed_at_most'], v['Li_at_ratio_pdg'], v['Li_at_ratio_0p2']))
print('(d) photon thresholds (MeV): ' + ', '.join('%s %.3f' % kv for kv in out['photon_thresholds_MeV'].items()))
