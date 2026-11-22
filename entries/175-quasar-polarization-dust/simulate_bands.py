"""Entry 175, exploratory: does the simulation (random angles, predicted foreground, the 0.6% cut) reproduce the link by polarization band?

The data keep a link of +0.19 +- 0.08 for quasars above 3%, more than an added 0.16% should give.  The simulation of simulate.py is
repeated and the link computed in the same four bands of observed polarization, for the parents with median 0.5% and 0.8%, and
for a two-part parent: 80% log-normal with median 0.4% and 20% with median 2.5% (a tail of strongly polarized objects).
Run from the repository root:  .venv/bin/python -I entries/175-quasar-polarization-dust/simulate_bands.py     (about 2 minutes)
"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import qso_dust as q

z = np.load(os.path.join(HERE, 'work', 'template.npz'))
tq, tu, tpa = z['q'], z['u'], z['pa']
Q = q.load_quasars(); ep = Q['ep']; n = len(tq)
bands = [(0.6, 1.0), (1.0, 1.6), (1.6, 3.0), (3.0, 99.0)]
data = json.load(open(os.path.join(HERE, 'results.json')))['reported_only']['link_by_polarization']
out = dict(data=data, cases=[])
def parent(g, m, kind):
    if kind == 'two':
        tail = g.uniform(size=m) < 0.2
        return np.where(tail, g.lognormal(np.log(2.5), 0.7, m), g.lognormal(np.log(0.4), 0.7, m))
    return g.lognormal(np.log(kind), 0.9, m)
for kind in (0.5, 0.8, 'two'):
    g = np.random.default_rng(q.CFG['seed'] + 31 + (7 if kind == 'two' else int(100 * kind)))
    acc = {b: [] for b in bands}; pm = []; allc = []
    for it in range(400):
        p_out, pa_out = np.empty(n), np.empty(n); todo = np.arange(n)
        while len(todo):
            m = len(todo); pi = parent(g, m, kind); th = g.uniform(0, np.pi, m)
            qo = pi * np.cos(2 * th) + tq[todo] + g.normal(0, 1, m) * ep[todo]; uo = pi * np.sin(2 * th) + tu[todo] + g.normal(0, 1, m) * ep[todo]
            po = np.hypot(qo, uo); okk = po >= 0.6
            p_out[todo[okk]] = po[okk]; pa_out[todo[okk]] = np.mod(0.5 * np.arctan2(uo[okk], qo[okk]) / q.D2R, 180); todo = todo[~okk]
        c = np.cos(2 * (pa_out - tpa) * q.D2R)
        allc.append(c.mean()); pm.append(np.median(p_out))
        for b in bands:
            mm = (p_out >= b[0]) & (p_out < b[1])
            if mm.sum() >= 5:
                acc[b].append((c[mm].mean(), mm.sum()))
    row = dict(parent=str(kind), median_observed_p=float(np.mean(pm)), link_all=float(np.mean(allc)), bands=[])
    for b, dd in zip(bands, data):
        v = np.array([a for a, _ in acc[b]]); cnt = np.mean([k for _, k in acc[b]])
        row['bands'].append(dict(p_lo=b[0], p_hi=b[1], simulated_link=float(v.mean()), sd_between_samples=float(v.std(ddof=1)), mean_count=float(cnt),
                                 data_link=dd['mean_cos2delta'], data_error=dd['error'], data_count=dd['n'],
                                 share_of_samples_at_or_above_data=float((v >= dd['mean_cos2delta']).mean())))
    out['cases'].append(row)
    print('parent %s: median observed p %.2f, link all %.3f | ' % (kind, row['median_observed_p'], row['link_all']) + '; '.join(
        '%.1f-%.1f%%: sim %.3f +- %.3f (n %.0f) data %.3f (n %d), P(sim>=data) %.3f' % (x['p_lo'], min(x['p_hi'], 9.9), x['simulated_link'], x['sd_between_samples'],
                                                                                   x['mean_count'], x['data_link'], x['data_count'], x['share_of_samples_at_or_above_data']) for x in row['bands']), flush=True)
json.dump(out, open(os.path.join(HERE, 'simulate_bands.json'), 'w'), indent=1)
