"""Entry 175, exploratory: how firmly does the f = 1 result sit below 2 sigma?

Made after the frozen result; nothing here decides anything.  Two things are varied: the aperture over which the dust map is
averaged (0.25, 0.5, 1 and 2 degrees; the rule used 0.5), and the conversion ratio P_353 / p_V within its 10% uncertainty
(equivalently f = 0.9 and 1.1).  10,000 shuffles each, one seed for all so that the rows differ only by the template.
Run from the repository root:  .venv/bin/python -I entries/175-quasar-polarization-dust/variants.py     (about 4 minutes)
"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import qso_dust as q

Q = q.load_quasars()
nvs = q.CFG['nvs']
prep = q.prepare(Q['ra'], Q['dec'], Q['z'], max(nvs))
SEED = q.CFG['seed'] + 4242
out = dict(aperture=[], ratio=[])
maps = None
for disc in (0.25, 0.5, 1.0, 2.0):
    q.CFG['disc_deg'] = disc
    T, maps = q.template_at(Q['ra'], Q['dec'], maps)
    fs, cs = q.f_star(Q['p'], Q['pa'], T)
    p2, pa2 = q.correct(Q['p'], Q['pa'], T['q'], T['u'], 1.0)
    r = q.align_test(prep, pa2 * q.D2R, nvs, 10000, SEED)
    out['aperture'].append(dict(disc_deg=disc, template_median=float(np.median(T['p'])), link0=cs[0], f_star=fs, sigma=r['sigma'],
                                p_global=r['p_global'], S=r['sigma_global']))
    print('aperture %.2f deg: template median %.3f%%, link %+.3f, f* %.2f, sigma %s, global p %.4f (%.2f sigma)' % (
        disc, np.median(T['p']), cs[0], fs, [round(x, 2) for x in r['sigma']], r['p_global'], r['sigma_global']), flush=True)
    if disc == 0.5:
        T05 = T
for f in (0.8, 0.9, 1.0, 1.1, 1.2):
    p2, pa2 = q.correct(Q['p'], Q['pa'], T05['q'], T05['u'], f)
    r = q.align_test(prep, pa2 * q.D2R, nvs, 10000, SEED)
    out['ratio'].append(dict(f=f, ratio_MJy_sr=round(5.4 / f, 2), sigma=r['sigma'], p_global=r['p_global'], S=r['sigma_global']))
    print('f %.1f (ratio %.2f MJy/sr): sigma %s, global p %.4f (%.2f sigma)' % (f, 5.4 / f, [round(x, 2) for x in r['sigma']], r['p_global'], r['sigma_global']), flush=True)
json.dump(out, open(os.path.join(HERE, 'variants.json'), 'w'), indent=1)
