"""Entry 175, exploratory: can a foreground of the predicted size, with the sample's cut at p >= 0.6%, make what is observed?

Made after the frozen result; nothing here decides anything.  At the 355 real positions, with the real predicted foreground
vectors (work/template.npz) and the real measurement errors, quasars are given random intrinsic angles and an intrinsic
polarization drawn from a log-normal parent; the foreground times g is added, measurement noise is added, and a quasar is redrawn
until its observed polarization passes 0.6%.  g = 1 is Planck's prediction; g = 2 is a foreground twice as strong.
Then the same quantities as for the data: the link before correction and after removing the predicted foreground (f = 1), the
level f* that ends the link, and the alignment (30 neighbours) before and after.
Run from the repository root:  .venv/bin/python -I entries/175-quasar-polarization-dust/simulate.py     (about 10 minutes)
"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import qso_dust as q

z = np.load(os.path.join(HERE, 'work', 'template.npz'))
ra, dec, zz, tq, tu, tpa = z['ra'], z['dec'], z['z'], z['q'], z['u'], z['pa']
Q = q.load_quasars()
ep = Q['ep']
n = len(ra)
nvs = [30]
prep = q.prepare(ra, dec, zz, 30)
T = dict(q=tq, u=tu, pa=tpa)
obs = json.load(open(os.path.join(HERE, 'results.json')))
data = dict(link0=obs['link']['mean_cos2delta'][0], link1=obs['link']['mean_cos2delta'][4], f_star=obs['link']['f_star'],
            sigma30_0=obs['alignment']['f0']['sigma'][2], sigma30_1=obs['alignment']['f1']['sigma'][2])


def draw(g, med, sig, gain):
    """One synthetic sample: observed p and angle at each position."""
    p_out, pa_out = np.empty(n), np.empty(n)
    todo = np.arange(n)
    while len(todo):
        m = len(todo)
        pi = g.lognormal(np.log(med), sig, m); th = g.uniform(0, np.pi, m)
        qo = pi * np.cos(2 * th) + gain * tq[todo] + g.normal(0, 1, m) * ep[todo]
        uo = pi * np.sin(2 * th) + gain * tu[todo] + g.normal(0, 1, m) * ep[todo]
        po = np.hypot(qo, uo)
        ok = po >= q.CFG['pcut']
        p_out[todo[ok]] = po[ok]; pa_out[todo[ok]] = np.mod(0.5 * np.arctan2(uo[ok], qo[ok]) / q.D2R, 180)
        todo = todo[~ok]
    return p_out, pa_out


out = dict(data=data, cases=[])
for med, sig, gain in ((0.5, 0.9, 1.0), (0.3, 0.9, 1.0), (0.8, 0.9, 1.0), (0.5, 0.9, 2.0), (0.3, 0.9, 2.0), (0.5, 0.9, 0.0)):
    g = np.random.default_rng(q.CFG['seed'] + int(1000 * med) + int(10 * gain))
    rec = {k: [] for k in ('link0', 'link1', 'f_star', 'sigma30_0', 'sigma30_1', 'pmed')}
    for it in range(120):
        p, pa = draw(g, med, sig, gain)
        fs, cs = q.f_star(p, pa, T)
        rec['link0'].append(cs[0]); rec['link1'].append(cs[4]); rec['f_star'].append(np.nan if fs is None else fs)
        rec['sigma30_0'].append(q.align_test(prep, pa * q.D2R, nvs, 200, 1000 + it)['sigma'][0])
        p1, pa1 = q.correct(p, pa, tq, tu, 1.0)
        rec['sigma30_1'].append(q.align_test(prep, pa1 * q.D2R, nvs, 200, 5000 + it)['sigma'][0])
        rec['pmed'].append(float(np.median(p)))
    s = dict(parent_median=med, parent_sigma_ln=sig, foreground_gain=gain, realisations=120,
             share_of_parent_passing=None)
    for k, v in rec.items():
        v = np.array(v, float); ok = np.isfinite(v)
        s[k] = dict(mean=float(v[ok].mean()), sd=float(v[ok].std(ddof=1)), q16=float(np.percentile(v[ok], 16)), q84=float(np.percentile(v[ok], 84)))
    s['f_star']['undefined'] = int(np.isnan(np.array(rec['f_star'], float)).sum())
    for k in ('link0', 'sigma30_0', 'sigma30_1'):
        v = np.array(rec[k]); s[k]['share_at_or_above_data'] = float((v >= data[k]).mean())
    out['cases'].append(s)
    print('parent median %.1f%%, foreground x%.0f: link %.3f +- %.3f (data %.3f), after f=1 %.3f (data %.3f), f* %.2f +- %.2f (data %.2f), '
          'alignment n=30 %.2f +- %.2f (data %.2f; share >= data %.3f), after f=1 %.2f +- %.2f (data %.2f), median p %.2f' % (
              med, gain, s['link0']['mean'], s['link0']['sd'], data['link0'], s['link1']['mean'], data['link1'], s['f_star']['mean'], s['f_star']['sd'],
              data['f_star'], s['sigma30_0']['mean'], s['sigma30_0']['sd'], data['sigma30_0'], s['sigma30_0']['share_at_or_above_data'],
              s['sigma30_1']['mean'], s['sigma30_1']['sd'], data['sigma30_1'], s['pmed']['mean']), flush=True)
json.dump(out, open(os.path.join(HERE, 'simulate.json'), 'w'), indent=1)
