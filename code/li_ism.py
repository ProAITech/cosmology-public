"""li_ism.py - Workbench Entry 98
Lithium in gas that halo stars never processed: the Small Magellanic Cloud cloud toward Sk 143.
Measured (Molaro et al. 2024, ESPRESSO): log N(7Li I)/N(K I) = -2.31 +- 0.02; with the Li-K ionization correction +0.54 +- 0.08,
log(7Li/K) = -1.77 +- 0.08, i.e. [Li/K] = +0.07 +- 0.08 against the solar -1.84 (Lodders 2021). Howk et al. (2012) measured the same
columns; the two papers differ only in the cloud's metallicity, which sets K/H:  A(Li) = A(Li)sun + [K/H] + [Li/K].
Metallicity tracers of this cloud (each + [K/Fe] = +0.1): Zn -1.28 +- 0.09; the full depletion pattern (De Cia et al. 2024) -1.09 +- 0.24;
S -0.72 +- 0.11; Howk et al. used SMC stars, [Fe/H] = -0.59 +- 0.06 with [K/Fe] = 0.
Two floors: standard nucleosynthesis with Planck's baryons, Li/H = 5.464 +- 0.220 e-10 (Fields et al. 2020, A = 2.737);
the halo-star plateau, A = 2.20 +- 0.09 (Sbordone et al. 2010). Galactic production is bracketed: none, or a linear rise in metallicity
to the solar value, or a delayed (quadratic) rise.
"""
import json, math, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
A_SUN, LIK, LIK_E = 3.27, 0.07, 0.08
BBN, PLAT = 2.737, 2.20
MET = {'Zn (Molaro et al. 2024)': (-1.18, 0.10), 'depletion pattern (De Cia et al. 2024)': (-0.99, 0.25), 'S (Molaro et al. 2024)': (-0.62, 0.12), 'SMC stars (Howk et al. 2012)': (-0.59, 0.06)}
def A_meas(KH): return A_SUN + KH + LIK
def A_model(KH, floor, gce):
    Z = 10**KH; fl = 10**(floor - 12); sun = 10**(A_SUN - 12)
    f = {'none': 0.0, 'linear': Z, 'delayed': Z*Z}[gce]
    return 12 + math.log10(fl + max(0.0, sun - fl)*f)
out = {'Li_over_K': [LIK, LIK_E], 'A_sun': A_SUN, 'BBN': BBN, 'plateau': PLAT, 'by_metallicity': {}}
for k, (kh, e) in MET.items():
    am = A_meas(kh); err = math.hypot(LIK_E, e, 0.03)
    row = dict(K_H=kh, err=e, A=round(am, 3), A_err=round(err, 3))
    for fl_name, fl in (('BBN', BBN), ('plateau', PLAT)):
        for g in ('none', 'delayed', 'linear'):
            m = A_model(kh, fl, g); row['%s_%s' % (fl_name, g)] = dict(A=round(m, 3), sigma=round((am - m)/math.hypot(err, 0.09 if fl_name == 'plateau' else 0.017), 2))
    out['by_metallicity'][k] = row
    print('%-42s [K/H] %+.2f  A(Li) %.2f +- %.2f  | BBN: %s  plateau: %s' % (k, kh, am, err,
          ' '.join('%s %+.1fs' % (g, row['BBN_' + g]['sigma']) for g in ('none', 'delayed', 'linear')),
          ' '.join('%s %+.1fs' % (g, row['plateau_' + g]['sigma']) for g in ('none', 'delayed', 'linear'))))
# the metallicity each floor needs for the measured Li/K
from scipy.optimize import brentq
need = {}
for fl_name, fl in (('BBN', BBN), ('plateau', PLAT)):
    for g in ('none', 'delayed', 'linear'):
        need['%s_%s' % (fl_name, g)] = round(brentq(lambda kh: A_meas(kh) - A_model(kh, fl, g), -3.0, -0.001), 2)
out['K_H_needed'] = need; print('metallicity needed:', need)
# what [Li/K] the BBN floor predicts at each tracer, without galactic production (a lower bound)
out['LiK_pred_BBN_none'] = {k: round(BBN - A_SUN - kh, 2) for k, (kh, e) in MET.items()}
out['LiK_pred_plateau_none'] = {k: round(PLAT - A_SUN - kh, 2) for k, (kh, e) in MET.items()}
print('BBN floor, no production, predicts [Li/K]:', out['LiK_pred_BBN_none'])
# curves for the figure
kh = np.round(np.linspace(-2.0, 0.2, 111), 3)
out['curves'] = {'kh': kh.tolist(), 'meas': [round(A_meas(x), 3) for x in kh]}
for fl_name, fl in (('BBN', BBN), ('plateau', PLAT)):
    for g in ('none', 'linear'):
        out['curves']['%s_%s' % (fl_name, g)] = [round(A_model(x, fl, g), 3) for x in kh]
json.dump(out, open(os.path.join(HERE, 'li_ism.json'), 'w'), indent=1)
