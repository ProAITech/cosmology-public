"""Entry 91: is matter less clumped than the CMB predicts? The S8 measurements of 2026.

S8 = sigma8 (Omega_m/0.3)^0.5, all in flat LCDM, as published. Reference: the primary CMB of
Planck + ACT DR6 + SPT-3G, 0.836 (+0.012 / -0.013). Tension uses the error on the side facing
the other value. The weighted mean treats the low-redshift measurements as independent, which
they are not quite (shared lensing maps, overlapping sky), so its error is a lower bound.
"""
import json, math, os, numpy as np
os.chdir(os.path.dirname(os.path.abspath(__file__)))
CMB = (0.836, 0.012, 0.013)                    # value, +err, -err
M = [  # name, group, value, +err, -err, independent-set flag
    ('KiDS-Legacy cosmic shear (2025)', 'shear', 0.815, 0.016, 0.021, True),
    ('DES Y6 cosmic shear, NLA (2026)', 'shear', 0.798, 0.014, 0.015, False),
    ('DES Y6 clustering + lensing, 3×2pt (2026)', 'shear', 0.789, 0.012, 0.012, True),
    ('HSC Y3 cosmic shear (2023)', 'shear', 0.776, 0.032, 0.033, True),
    ('HSC lensing × BOSS galaxies (2025)', 'galaxy lensing', 0.8294, 0.011, 0.011, True),
    ('unWISE × ACT+Planck CMB lensing (2024)', 'CMB lensing', 0.810, 0.015, 0.015, True),
    ('ACT+Planck+SPT CMB lensing, with DESI BAO (2026)', 'CMB lensing', 0.828, 0.012, 0.012, True),
    ('eROSITA clusters (2024)', 'clusters', 0.86, 0.01, 0.01, True),
    ('SPT clusters (2024)', 'clusters', 0.795, 0.029, 0.029, True),
    ('DES Y3 shear with kSZ-measured feedback (2024)', 'feedback', 0.823, 0.019, 0.020, False),
]
def tension(v, ep, em):
    c, cp, cm = CMB
    if v < c: return (c - v)/math.hypot(ep, cm)
    return -(v - c)/math.hypot(em, cp)
rows = []
for n, g, v, ep, em, ind in M:
    rows.append(dict(name=n, group=g, S8=v, ep=ep, em=em, independent=ind, tension_sigma=round(tension(v, ep, em), 2),
                     ratio=round(v/CMB[0], 3)))
ind = [r for r in rows if r['independent']]
v = np.array([r['S8'] for r in ind]); e = np.array([(r['ep']+r['em'])/2 for r in ind])
w = 1/e**2; mean = float(np.sum(w*v)/np.sum(w)); err = float(1/math.sqrt(np.sum(w)))
chi2 = float(np.sum((v-mean)**2*w)); dof = len(v)-1
from scipy import stats
p = float(stats.chi2.sf(chi2, dof))
# scale factor for inconsistency (PDG style)
sf = math.sqrt(max(chi2/dof, 1.0))
res = dict(cmb=CMB, rows=rows, low_z_mean=round(mean, 4), low_z_err=round(err, 4), chi2=round(chi2, 2), dof=dof, p=float('%.2g' % p),
           scale_factor=round(sf, 2), low_z_err_scaled=round(err*sf, 4),
           tension_mean=round((CMB[0]-mean)/math.hypot(err, CMB[2]), 2), tension_mean_scaled=round((CMB[0]-mean)/math.hypot(err*sf, CMB[2]), 2),
           deficit_percent=round(100*(1-mean/CMB[0]), 2))
# leave-one-out: how much does each measurement move the mean / chi2
loo = {}
for i, r in enumerate(ind):
    m = np.ones(len(ind), bool); m[i] = False
    mm = float(np.sum(w[m]*v[m])/np.sum(w[m])); c2 = float(np.sum((v[m]-mm)**2*w[m]))
    loo[r['name']] = dict(mean=round(mm, 4), chi2=round(c2, 2), dof=int(m.sum()-1))
res['leave_one_out'] = loo
# without the two extremes (eROSITA high, DES Y6 low)
m = np.array([not (r['name'].startswith('eROSITA') or r['name'].startswith('DES Y6')) for r in ind])
mm = float(np.sum(w[m]*v[m])/np.sum(w[m])); ee = float(1/math.sqrt(np.sum(w[m])))
res['without_extremes'] = dict(mean=round(mm, 4), err=round(ee, 4), chi2=round(float(np.sum((v[m]-mm)**2*w[m])), 2), dof=int(m.sum()-1),
                               tension=round((CMB[0]-mm)/math.hypot(ee, CMB[2]), 2))
json.dump(res, open('s8_board.json', 'w'), indent=1)
for r in rows: print(r['name'], r['S8'], r['tension_sigma'])
for k, val in res.items():
    if k != 'rows': print(k, val)
