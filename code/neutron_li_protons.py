"""Workbench Entry 170, addition of 22 November 2026: can PROTONS from evaporating black holes deplete 7Li?

The journal version of Arbuzova & Dolgov (Eur. Phys. J. C 85, 1362, 2025; arXiv:2606.03316) names proton capture,
7Li + p -> 8Be -> 2 4He, beside neutron capture.  Three cases, no data, only kinematics and cross-sections:

 (a) Thermal protons at T = 1.8 eV.  The Coulomb barrier: non-resonant rate with S(0) = 55 keV b for 7Li(p,alpha)4He.
     Evaporated protons that have thermalised are also indistinguishable from the ambient ones, which outnumber them.
 (b) Energetic protons below the break-up threshold of 4He (22.9 MeV in the laboratory, from Q = -18.35 MeV for
     4He(p,d)3He).  They slow down on the plasma electrons; the column they cross before stopping, times the cross-section,
     is the chance that a given 7Li nucleus is destroyed per energetic proton per ambient hydrogen atom.
     Stopping on free electrons: -dE/dx = (4 pi e^4 n_e / m_e v^2) ln(2 m_e v^2 gamma^2 / hbar omega_p).
 (c) Protons above that threshold.  The same fluence acts on 7Li and on 4He, so the destroyed shares are tied by the ratio
     of cross-sections, whatever the number of protons: 1 - f_He = (1 - f_Li)^(sigma_He / sigma_Li).

Cross-sections (representative values, to be checked against the compilations before submission; see the draft's notes):
  7Li(p,alpha)4He below 23 MeV: 0.1 b taken as representative, 0.5 b as a ceiling that no single channel reaches;
  p + 4He inelastic above threshold: about 0.10 b;  p + 7Li reaction cross-section at 0.1-1 GeV: about 0.25 b (range 0.2-0.5).
Sources for the physics: Meyer (1972, A&AS 7, 417); Kawasaki, Kohri & Moroi (2005, PRD 71, 083502); Jedamzik (2006, PRD 74, 103509).
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
def column_e(E_hi, E_lo=0.1, n=20000):
    """electron column (cm^-2) crossed while slowing from E_hi to E_lo (MeV)"""
    E = np.linspace(E_lo, E_hi, n); g = 1 + E / MP; b2 = 1 - 1 / g ** 2
    lnL = np.log(2 * ME * 1e6 * b2 * g ** 2 / hbar_wp_eV)
    y = b2 / (K4 * lnL)
    return float(np.sum((y[1:] + y[:-1]) / 2 * np.diff(E))), float(lnL[n // 2])
Ne, lnL_mid = column_e(E_th)
N_H = Ne * n_H / n_e
ratioA = 1e3 / n_H                                                 # reading A of the paper's Eq. (25): 1e3 evaporated baryons per cm^3
sub = dict(threshold_MeV=E_th, coulomb_log_mid=lnL_mid, hbar_omega_p_eV=hbar_wp_eV, column_e_cm2=Ne, column_H_cm2=N_H, evaporated_per_H_reading_A=ratioA, cases=[])
for sig_b, label in ((0.1, 'representative'), (0.5, 'ceiling')):
    per = N_H * sig_b * 1e-24                                      # share of 7Li destroyed per energetic proton per ambient H
    sub['cases'].append(dict(sigma_b=sig_b, label=label, destroyed_per_proton_per_H=per, destroyed_reading_A=per * ratioA, protons_per_H_for_ln3=math.log(3) / per))
out['below_threshold'] = sub

# ---- (c) protons above the threshold: lithium against helium
DH_obs, DH_err = 2.527e-5, 0.030e-5
He_per_H = n_He / n_H
above = dict(He_per_H=He_per_H, cases=[])
for r in (0.2, 0.4, 0.5):                                          # sigma_He / sigma_Li
    fHe = 1 - 3 ** (-r)
    row = dict(sigma_He_over_sigma_Li=r, He_destroyed_for_Li_to_a_third=fHe, D_over_H_added={})
    for y in (0.01, 0.3):                                          # deuterons per 4He broken up
        row['D_over_H_added'][str(y)] = fHe * He_per_H * y
    above['cases'].append(row)
# the other way round: the most 7Li that can go if D/H may rise by 3 sigma
for y in (0.3, 0.1):
    fHe_max = 3 * DH_err / (He_per_H * y)
    above['Li_destroyed_at_most_y%s' % y] = dict(He_destroyed_at_most=fHe_max, Li_destroyed_at_most_ratio0p4=1 - (1 - fHe_max) ** (1 / 0.4), Li_destroyed_at_most_ratio0p2=1 - (1 - fHe_max) ** (1 / 0.2))
out['above_threshold'] = above

json.dump(out, open(os.path.join(HERE, 'neutron_li_protons.json'), 'w'), indent=1)
p = out['plasma']; t = out['thermal']
print('plasma at 1.8 eV: z = %.0f, n_b = %.3g, n_H = %.3g, n_He = %.3g, n_e = %.3g cm^-3' % (p['z'], p['n_b'], p['n_H'], p['n_He'], p['n_e']))
print('(a) thermal: Gamow peak %.3f keV, 3E0/kT = %.1f, log10 <sigma v> = %.1f, log10 K = %.1f' % (t['gamow_peak_keV'], t['tau'], t['log10_sigma_v_cm3_s'], t['log10_K']))
print('(b) below %.1f MeV: Coulomb log %.1f, column %.3g electrons = %.3g H per cm^2; evaporated per H (reading A) %.4f' % (E_th, lnL_mid, Ne, N_H, ratioA))
for c in sub['cases']:
    print('    sigma %.1f b (%s): 7Li destroyed %.2e per proton per H; reading A %.1e; protons per H needed for a third: %.0f' % (
        c['sigma_b'], c['label'], c['destroyed_per_proton_per_H'], c['destroyed_reading_A'], c['protons_per_H_for_ln3']))
for c in above['cases']:
    print('(c) sigma_He/sigma_Li %.1f: 4He destroyed %.1f%%; D/H added %.1e (1 deuteron per 100 break-ups), %.1e (30 per 100); measured D/H %.3e' % (
        c['sigma_He_over_sigma_Li'], 100 * c['He_destroyed_for_Li_to_a_third'], c['D_over_H_added']['0.01'], c['D_over_H_added']['0.3'], DH_obs))
for k in ('Li_destroyed_at_most_y0.3', 'Li_destroyed_at_most_y0.1'):
    v = above[k]; print('    %s: 4He at most %.1e, 7Li at most %.1e (ratio 0.4) or %.1e (ratio 0.2)' % (k, v['He_destroyed_at_most'], v['Li_destroyed_at_most_ratio0p4'], v['Li_destroyed_at_most_ratio0p2']))
