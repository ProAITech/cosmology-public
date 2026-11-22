"""helium.py - Workbench Entry 102: is primordial helium low?

Every route to the primordial helium mass fraction Y_p, set against the Big Bang prediction and against
Lerner's claim that the oldest stars hold "less than half" the predicted helium.

Inputs (all published; see the entry's references):
  Prediction, standard nucleosynthesis with Planck's baryon density: 0.2467 +- 0.0002 (as quoted by LBT IV)
  Routes, one value each (independent methods):
    H II regions   LBT Y_p Project IV (Aver et al. 2026), 15 targets with O/H < 4e-5: 0.2458 +- 0.0013
    CMB            ACT DR6 extended models (Calabrese et al. 2025): 0.2312 +- 0.0092
    oldest stars   57 globular clusters, R-parameter (Salaris et al. 2004): 0.250 +- 0.006
    pristine gas   absorber at z = 1.724, ~1/100 solar metals (Cooke & Fumagalli 2018): 0.250 +0.033 -0.025
  Same method, other teams (H II regions): EMPRESS XV 2025 0.2402 +- 0.0040, EMPRESS VIII 2022 0.2370 +0.0034 -0.0033,
    and the analyses tabulated by LBT IV: Izotov 2014 0.2551 +- 0.0022, Aver 2015 0.2449 +- 0.0040,
    Hsyu 2020 0.2436 +0.0039 -0.0040, Kurichin 2021 0.2462 +- 0.0022, Aver 2022 0.2448 +- 0.0033
  LBT IV regression of y = n(He)/n(H) on O/H, 41 targets: y_p = 0.08117 +- 0.00096, slope 14 +- 12;
    15-target weighted mean y_p = 0.08146 +- 0.00058 with error scale factor S = 1.15; lowest unflagged O/H 1.49e-5 (I Zw 18 SE)
  Sensitivities (Steigman 2012): dY = 0.013 dN_eff, dY = -0.23 xi_e
  N_eff and H0 with N_eff free, from helium + deuterium + CMB (Planck, ACT, SPT) + DESI BAO (Goldstein & Hill 2026):
    N_eff = 2.990 +- 0.070, H0 = 68.16 +- 0.50; direct one-step H0 average from Entry 56: 71.4 +- 1.2
"""
import json, math, os
from scipy.stats import chi2 as CHI2
HERE = os.path.dirname(os.path.abspath(__file__))
BBN, BBN_E = 0.2467, 0.0002
HALF = BBN/2

def side(y, em, ep, ref):          # the error on the side facing ref
    return ep if y < ref else em

ROUTES = [dict(key='hii', name='H II regions, LBT 2026', Y=0.2458, em=0.0013, ep=0.0013),
          dict(key='cmb', name='CMB damping tail, ACT DR6', Y=0.2312, em=0.0092, ep=0.0092),
          dict(key='gc', name='Oldest stars: 57 globular clusters', Y=0.250, em=0.006, ep=0.006),
          dict(key='igm', name='Gas cloud at z = 1.72', Y=0.250, em=0.025, ep=0.033)]
OTHER = [dict(name='EMPRESS 2025', Y=0.2402, em=0.0040, ep=0.0040),
         dict(name='Izotov 2014', Y=0.2551, em=0.0022, ep=0.0022),
         dict(name='Aver 2015', Y=0.2449, em=0.0040, ep=0.0040),
         dict(name='Hsyu 2020', Y=0.2436, em=0.0040, ep=0.0039),
         dict(name='Kurichin 2021', Y=0.2462, em=0.0022, ep=0.0022),
         dict(name='EMPRESS 2022', Y=0.2370, em=0.0033, ep=0.0034),
         dict(name='Aver 2022', Y=0.2448, em=0.0033, ep=0.0033)]
out = {'BBN': [BBN, BBN_E], 'half_claim_Y': round(HALF, 4)}

# 1. each route against the prediction and against "less than half"
for r in ROUTES + OTHER:
    r['sig_BBN'] = round((r['Y'] - BBN)/math.hypot(side(r['Y'], r['em'], r['ep'], BBN), BBN_E), 2)
    r['sig_above_half'] = round((r['Y'] - HALF)/r['em'], 1)
out['routes'], out['other_HII'] = ROUTES, OTHER

# 2. the four independent routes combined (iterate so each asymmetric error faces the mean)
mean = BBN
for _ in range(20):
    w = [1/side(r['Y'], r['em'], r['ep'], mean)**2 for r in ROUTES]
    mean = sum(wi*r['Y'] for wi, r in zip(w, ROUTES))/sum(w)
err = sum(w)**-0.5
chi = sum(wi*(r['Y'] - mean)**2 for wi, r in zip(w, ROUTES))
out['combined'] = dict(Y=round(mean, 5), err=round(err, 5), sig_BBN=round((mean - BBN)/math.hypot(err, BBN_E), 2),
                       chi2=round(chi, 2), dof=len(ROUTES) - 1, p=round(float(CHI2.sf(chi, len(ROUTES) - 1)), 2))
# without the H II regions: does the rest of the sky agree on its own?
R3 = [r for r in ROUTES if r['key'] != 'hii']; m3 = BBN
for _ in range(20):
    w3 = [1/side(r['Y'], r['em'], r['ep'], m3)**2 for r in R3]; m3 = sum(wi*r['Y'] for wi, r in zip(w3, R3))/sum(w3)
out['combined_without_HII'] = dict(Y=round(m3, 4), err=round(sum(w3)**-0.5, 4), sig_BBN=round((m3 - BBN)/sum(w3)**-0.5, 2))

# 3. the same method, other teams: how far apart are the H II analyses?
ys = [0.2458] + [o['Y'] for o in OTHER if o['name'] != 'EMPRESS 2022']
out['HII_spread'] = dict(min=min(ys), max=max(ys), range=round(max(ys) - min(ys), 4))
out['EMPRESS25_vs_LBT'] = round((0.2458 - 0.2402)/math.hypot(0.0013, 0.0040), 2)
out['Izotov14_vs_LBT'] = round((0.2551 - 0.2458)/math.hypot(0.0013, 0.0022), 2)

# 4. what helium says about new physics, at Planck's baryon density
imp = {}
for lab, y, e in (('LBT', 0.2458, 0.0013), ('EMPRESS 2025', 0.2402, 0.0040), ('combined', mean, err)):
    dY, eY = y - BBN, math.hypot(e, BBN_E)
    imp[lab] = dict(dN_eff=round(dY/0.013, 3), dN_eff_err=round(eY/0.013, 3), xi_e=round(-dY/0.23, 4), xi_e_err=round(eY/0.23, 4))
out['implied'] = imp
out['xi_005_vs_LBT_sigma'] = round((0.05 - imp['LBT']['xi_e'])/imp['LBT']['xi_e_err'], 1)   # EMPRESS's lepton asymmetry, read at standard N_eff
H0_bbn, H0_bbn_e, H0_dir, H0_dir_e = 68.16, 0.50, 71.4, 1.2
out['H0'] = dict(N_eff=[2.990, 0.070], H0=[H0_bbn, H0_bbn_e], direct=[H0_dir, H0_dir_e],
                 sigma=round((H0_dir - H0_bbn)/math.hypot(H0_dir_e, H0_bbn_e), 2))

# 5. the floor in the most metal-poor galaxies
yp, S, n = 0.08146, 1.15, 15
dYdy = 4/(1 + 4*yp)**2                                           # Y = 4y(1-Z)/(1+4y), Z negligible here
sig_i = (0.00058/S)*math.sqrt(n)                                 # per-target error in y if all were equal
chi_obs = S**2*(n - 1)
q05 = float(CHI2.ppf(0.05, n - 1))
sint_best = sig_i*math.sqrt(max(S**2 - 1, 0)); sint_95 = sig_i*math.sqrt(chi_obs/q05 - 1)
out['floor'] = dict(per_target_err_Y=round(sig_i*dYdy, 4), intrinsic_scatter_Y=round(sint_best*dYdy, 4),
                    intrinsic_scatter_Y_95=round(sint_95*dYdy, 4), chi2=round(chi_obs, 1), dof=n - 1)
XH = 1 - BBN
dXO = 16*XH                                                      # oxygen mass fraction per unit O/H (number)
slope_Y = dYdy*14; slope_Y_e = dYdy*12
out['helium_per_oxygen'] = dict(measured=round(slope_Y/dXO, 2), err=round(slope_Y_e/dXO, 2),
                                needed_at_IZw18=round(0.2458/(dXO*1.49e-5)),
                                ratio=round(0.2458/(dXO*1.49e-5)/(slope_Y/dXO)),
                                ratio_2sigma=round(0.2458/(dXO*1.49e-5)/((slope_Y + 2*slope_Y_e)/dXO)),
                                dY_across_sample=round(slope_Y*(15e-5 - 1.49e-5), 4), dY_across_sample_err=round(slope_Y_e*(15e-5 - 1.49e-5), 4))
json.dump(out, open(os.path.join(HERE, 'helium.json'), 'w'), indent=1, ensure_ascii=False)
for k in ('combined', 'combined_without_HII', 'HII_spread', 'EMPRESS25_vs_LBT', 'Izotov14_vs_LBT', 'implied', 'xi_005_vs_LBT_sigma', 'H0', 'floor', 'helium_per_oxygen'):
    print(k, out[k])
for r in ROUTES + OTHER: print(r['name'], r['sig_BBN'], r['sig_above_half'])
