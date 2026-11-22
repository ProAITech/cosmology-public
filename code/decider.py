"""decider.py - Workbench Entry 108: the numbers that decide MOND against a dark halo at the Milky Way, compiled.

Entry 104 left two deciders: the mass inside 50 kpc and the escape speed at the Sun. Here every published value with the
Large Magellanic Cloud modelled (for the masses) is compiled, tested for consistency, and set against the two pictures:
  dark halo of Entry 103 (NFW M200 = 6.0e11, c = 16.9, fitted to Gaia's curve alone)
  MOND at the external field its surroundings supply (Entry 104, e = 0.00138 a0), stellar mass refitted to Gaia + tracers.
"""
import json, os, math, sys, io, contextlib
import numpy as np
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
with contextlib.redirect_stdout(io.StringIO()):
    import mw_mass as M
    import efe as F
G = M.G
M50 = [('Outer halo stars, H3+SEGUE+MagE (Brooks et al. 2026)', 3.36, 0.15, 0.15),
       ('Orphan stream (Erkal et al. 2019)', 3.80, 0.11, 0.14),
       ('Orphan-Chenab stream, S5 (Koposov et al. 2023)', 3.77, 0.19, 0.26),
       ('87 streams, STREAMFINDER (Ibata et al. 2024)', 4.6, 0.3, 0.3)]
M100 = [('Sagittarius stream (Vasiliev et al. 2021)', 5.6, 0.4, 0.4),
        ('Globular clusters and satellites (Correa Magnus & Vasiliev 2022)', 7.5, 1.0, 1.0),
        ('Globular clusters (Eadie & Juric 2019)', 5.3, 1.2, 2.1)]
VESC = [('LAMOST K giants (Wu et al. 2025)', 523.7, 13.5, 12.8),
        ('Gaia DR2 (Deason et al. 2019)', 528.0, 25.0, 24.0),
        ('Gaia DR2 (Necib & Lin 2022)', 484.6, 7.4, 17.8),
        ('Gaia EDR3 (Koppelman & Helmi 2021)', 497.0, 8.0, 8.0),
        ('RAVE (Williams et al. 2017)', 521.0, 30.0, 46.0)]
def combine(rows):
    m = np.mean([r[1] for r in rows])
    for _ in range(50):
        w = np.array([1/(r[3] if r[1] < m else r[2])**2 for r in rows]); m = float(np.sum(w*np.array([r[1] for r in rows]))/w.sum())
    e = w.sum()**-0.5; chi = float(np.sum(w*(np.array([r[1] for r in rows]) - m)**2)); dof = len(rows) - 1
    birge = math.sqrt(chi/dof) if chi > dof else 1.0
    return dict(mean=round(m, 3), err=round(e, 3), chi2=round(chi, 2), dof=dof, p=float('%.3g' % stats.chi2.sf(chi, dof)), birge=round(birge, 2), err_scaled=round(e*birge, 3),
                spread=[min(r[1] for r in rows), max(r[1] for r in rows)])
out = {'M50': combine(M50), 'M100': combine(M100), 'vesc': combine(VESC)}
out['vesc_gaia_only'] = combine(VESC[2:4]); out['vesc_other'] = combine([VESC[0], VESC[1], VESC[4]])
# predictions
nc = M.out['nfw']['curve_only']
pred = {'halo': dict(M50=nc['M50']/1e11, M100=nc['M100']/1e11, vesc_2r200=F.out['halo']['vesc_2r200'], vesc_inf=F.out['halo']['vesc_inf'])}
def mond_at(e):
    u, c2 = F.fit_curve(e)
    v50 = float(np.sqrt(M.gmond(M.gbar(np.array([50.0]), u), e)*50)[0]); v100 = float(np.sqrt(M.gmond(M.gbar(np.array([100.0]), u), e)*100)[0])
    return dict(e=e, M50=v50**2*50/G/1e11, M100=v100**2*100/G/1e11, vesc=F.vesc(u, e), chi2_curve_tracers=c2)
pred['mond_pull'] = mond_at(F.out['pull']['e_total'])
pred['mond_joint'] = mond_at(F.out['joint']['e'])
for k in pred:
    for q in pred[k]: pred[k][q] = (round(pred[k][q], 5) if q == 'e' else round(pred[k][q], 3)) if isinstance(pred[k][q], float) else pred[k][q]
out['pred'] = pred
# score each picture against each compiled value, with the Birge-scaled error
def sig(val, comb):
    e = comb['err_scaled']; return round((val - comb['mean'])/e, 2)
score = {}
for k, p in pred.items():
    s = {'M50': sig(p['M50'], out['M50']), 'M100': sig(p['M100'], out['M100'])}
    vv = p['vesc'] if 'vesc' in p else p['vesc_2r200']
    s['vesc_all'] = sig(vv, out['vesc']); s['vesc_gaia'] = round((vv - out['vesc_gaia_only']['mean'])/out['vesc_gaia_only']['err_scaled'], 2)
    s['vesc_other'] = round((vv - out['vesc_other']['mean'])/out['vesc_other']['err_scaled'], 2)
    score[k] = s
out['score'] = score
out['rows'] = dict(M50=M50, M100=M100, vesc=VESC)
json.dump(out, open(os.path.join(HERE, 'decider.json'), 'w'), indent=1)
for k in ('M50', 'M100', 'vesc', 'vesc_gaia_only', 'vesc_other', 'pred', 'score'): print(k, json.dumps(out[k]))
