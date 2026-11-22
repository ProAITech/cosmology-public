"""Entry 176, exploratory: how the excess inside R500 scales with cluster mass, and what the tails do.

Made after the frozen result was seen; nothing here decides anything.
  mass      the excess variance inside R500, and at 1-2 R500, for four quarters of the clusters by M500 (200 rotations each)
  inner     the same for x < 0.5 and 0.5 <= x < 1 by mass half
  tails     the share of sources with |dRM| above 30 rad/m^2 in each range, against the rotated clusters
Run from the repository root:  .venv/bin/python -I entries/176-cluster-faraday-profile/explore.py     (about 8 minutes)
"""
import json, os, sys
import numpy as np
from scipy.spatial import cKDTree
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cluster_rm as c

log = lambda s: print(s, flush=True)
rows, src, _ = c.load_spice()
cl, _ = c.load_clusters()
S, _ = c.sample(src, c.CFG['fg_radius'], c.CFG['min_neigh'])
out = dict(mass=[], tails={})
q = np.nanpercentile(cl['m500'], [0, 25, 50, 75, 100])
for i in range(4):
    sel = (cl['m500'] >= q[i]) & (cl['m500'] <= q[i + 1] if i == 3 else cl['m500'] < q[i + 1])
    sub = {k: v[sel] for k, v in cl.items()}
    r, _, _ = c.profile_test(S, sub, [[0.0, 0.5], [0.5, 1.0], [0.0, 1.0], [1.0, 2.0]], 200, c.CFG['seed'] + 60 + i, log=log, label='mass quarter %d' % (i + 1))
    out['mass'].append(dict(quarter=i + 1, m500_lo_1e13=float(q[i]), m500_hi_1e13=float(q[i + 1]), m500_median_1e13=float(np.nanmedian(sub['m500'])),
                            z_median=float(np.median(sub['z'])), clusters=int(sel.sum()), ranges=r))
# power law: excess variance inside R500 against median mass, weighted least squares in the linear quantity
m = np.array([x['m500_median_1e13'] for x in out['mass']]); ex = np.array([x['ranges'][2]['excess_variance'] for x in out['mass']])
sd = np.array([x['ranges'][2]['null_sd'] for x in out['mass']])
best = None
for a in np.linspace(0.0, 3.0, 301):
    f = m ** a
    A = (ex * f / sd ** 2).sum() / (f * f / sd ** 2).sum()
    chi = (((ex - A * f) / sd) ** 2).sum()
    if best is None or chi < best[0]:
        best = (chi, a, A)
chis = []
for a in np.linspace(0.0, 3.0, 301):
    f = m ** a; A = (ex * f / sd ** 2).sum() / (f * f / sd ** 2).sum(); chis.append((((ex - A * f) / sd) ** 2).sum())
chis = np.array(chis); grid = np.linspace(0.0, 3.0, 301); ok = grid[chis <= best[0] + 1]
out['mass_scaling'] = dict(exponent=float(best[1]), exponent_lo=float(ok.min()), exponent_hi=float(ok.max()), chi2=float(best[0]), dof=2,
                           note='excess variance inside R500 proportional to M500^exponent')
log('mass scaling: excess variance ~ M500^%.2f (%.2f to %.2f), chi2 %.1f for 2' % (best[1], ok.min(), ok.max(), best[0]))
# tails
tree = cKDTree(c.gal_unit(S['l'], S['b']))
n = len(S['d'])
x, _ = c.scaled_distance(tree, n, cl)
big = np.abs(S['d']) > 30.0
g = np.random.default_rng(c.CFG['seed'] + 99)
names = list(c.CFG['ranges'])
real = {k: float(big[(x >= c.CFG['ranges'][k][0]) & (x < c.CFG['ranges'][k][1])].mean()) for k in names}
nul = {k: [] for k in names}
for i in range(300):
    xr, _ = c.scaled_distance(tree, n, cl, shift=g.uniform(10.0, 350.0))
    for k in names:
        nul[k].append(big[(xr >= c.CFG['ranges'][k][0]) & (xr < c.CFG['ranges'][k][1])].mean())
for k in names:
    v = np.array(nul[k])
    out['tails'][k] = dict(share_above_30=real[k], null_mean=float(v.mean()), null_sd=float(v.std(ddof=1)), t=float((real[k] - v.mean()) / v.std(ddof=1)))
    log('tails %s: %.3f against %.3f +- %.3f, t %+.2f' % (k, real[k], v.mean(), v.std(ddof=1), out['tails'][k]['t']))
json.dump(out, open(os.path.join(HERE, 'explore.json'), 'w'), indent=1, default=float)
