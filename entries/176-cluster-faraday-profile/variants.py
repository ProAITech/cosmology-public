"""Entry 176, exploratory: catalogue cuts, and the fractional polarization of sources against distance from clusters.

Made after the frozen result was seen; nothing here decides anything.
  cuts      the three ranges with a stricter contamination cut (PCONT < 0.1), with redshift floors of 0.03 and 0.1, and with
            the 2,000 clusters of largest angular R500 removed (so that no cluster spans more than a few arcminutes)
  fracpol   the median fractional polarization (SPICE-RACS column fracpol, first component of each source) in the fine ranges,
            against the rotated clusters.  Plasma that scrambles rotation across the beam lowers it.  Needs Entry 171's
            work/extra.npz (the extra columns of the same rows as work/spice_rm.npz).
Run from the repository root:  .venv/bin/python -I entries/176-cluster-faraday-profile/variants.py     (about 6 minutes)
"""
import json, os, sys
import numpy as np
from scipy.spatial import cKDTree
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cluster_rm as c
from astropy.io import fits

log = lambda s: print(s, flush=True)
rows, src, _ = c.load_spice()
S, _ = c.sample(src, c.CFG['fg_radius'], c.CFG['min_neigh'])
names = list(c.CFG['ranges']); R3 = [c.CFG['ranges'][k] for k in names]
out = dict(cuts={}, fracpol={})


def clusters(pcont, zmin):
    c.CFG['pcont_max'], c.CFG['zmin'] = pcont, zmin
    cl, info = c.load_clusters()
    return cl, info


for name, pc, zm in (('PCONT below 0.1', 0.1, 0.05), ('z from 0.03', 0.5, 0.03), ('z from 0.1', 0.5, 0.1)):
    cl, info = clusters(pc, zm)
    r, _, _ = c.profile_test(S, cl, R3, 200, c.CFG['seed'] + 70, log=log, label=name)
    out['cuts'][name] = dict(clusters=info['kept'], ranges=dict(zip(names, r)))
cl, info = clusters(0.5, 0.05)
keep = np.argsort(cl['th'])[:-2000]
sub = {k: v[keep] for k, v in cl.items()}
r, _, _ = c.profile_test(S, sub, R3, 200, c.CFG['seed'] + 71, log=log, label='largest 2,000 on the sky removed')
out['cuts']['largest 2,000 on the sky removed'] = dict(clusters=int(len(keep)), theta500_arcmin_max=float(sub['th'].max() * 60), ranges=dict(zip(names, r)))

# ---- fractional polarization
zc = np.load(os.path.join(c.ROOT, 'entries', '171-spice-faraday-redshift', 'work', 'spice_rm.npz'))
ze = np.load(os.path.join(c.ROOT, 'entries', '171-spice-faraday-redshift', 'work', 'extra.npz'))
sid = zc['source_id']
u, first = np.unique(sid, return_index=True)               # first row of each source, as in merge_sources
fp_src = ze['fracpol'][first].astype(float); flux_src = ze['total_I_flux'][first].astype(float)
assert len(fp_src) == len(src['rm'])
idx = np.nonzero(np.abs(src['b']) >= c.CFG['bcut'])[0]     # sample() keeps every high-latitude source here (none lacks neighbours)
assert len(idx) == len(S['d'])
fp, fl = fp_src[idx], flux_src[idx]
tree = cKDTree(c.gal_unit(S['l'], S['b']))
n = len(fp)
edges = c.CFG['profile_edges']; E = [[edges[i], edges[i + 1]] for i in range(len(edges) - 1)]
x, _ = c.scaled_distance(tree, n, cl)
ok = np.isfinite(fp) & (fp > 0)
def med(xx, arr):
    return [float(np.median(arr[ok & (xx >= lo) & (xx < hi)])) for lo, hi in E]
# brighter sources can be detected at lower fractional polarization, and sources near cluster centres are brighter; so also
# divide each source's fracpol by the median fracpol of field sources (beyond 8 R500 of any cluster) in its decile of flux
field = ok & ~np.isfinite(x)
dec_edges = np.quantile(fl[ok], np.linspace(0, 1, 11)); dec_edges[0], dec_edges[-1] = -np.inf, np.inf
which = np.clip(np.searchsorted(dec_edges, fl, side='right') - 1, 0, 9)
expect = np.array([np.median(fp[field & (which == k)]) for k in range(10)])
fpn = fp / expect[which]
real_fp, real_fl, real_fn = med(x, fp), med(x, fl), med(x, fpn)
g = np.random.default_rng(c.CFG['seed'] + 88)
nf, nl, nn = [], [], []
for i in range(200):
    xr, _ = c.scaled_distance(tree, n, cl, shift=g.uniform(10.0, 350.0))
    nf.append(med(xr, fp)); nl.append(med(xr, fl)); nn.append(med(xr, fpn))
nf, nl, nn = np.array(nf), np.array(nl), np.array(nn)
out['fracpol'] = dict(units=str(ze['fracpol'].dtype), field_median=float(np.median(fp[ok & ~np.isfinite(x)])), ranges=[])
for j, (lo, hi) in enumerate(E):
    d = dict(lo=lo, hi=hi, median_fracpol=real_fp[j], null_mean=float(nf[:, j].mean()), null_sd=float(nf[:, j].std(ddof=1)),
             t=float((real_fp[j] - nf[:, j].mean()) / nf[:, j].std(ddof=1)),
             median_flux=real_fl[j], flux_null_mean=float(nl[:, j].mean()), flux_t=float((real_fl[j] - nl[:, j].mean()) / nl[:, j].std(ddof=1)),
             flux_matched_ratio=real_fn[j], ratio_null_mean=float(nn[:, j].mean()), ratio_null_sd=float(nn[:, j].std(ddof=1)),
             ratio_t=float((real_fn[j] - nn[:, j].mean()) / nn[:, j].std(ddof=1)))
    out['fracpol']['ranges'].append(d)
    log('fracpol x %.1f-%.1f: median %.4f against %.4f +- %.4f, t %+.2f; median flux %.4g against %.4g, t %+.2f; flux-matched ratio %.3f against %.3f +- %.3f, t %+.2f' % (
        lo, hi, d['median_fracpol'], d['null_mean'], d['null_sd'], d['t'], d['median_flux'], d['flux_null_mean'], d['flux_t'],
        d['flux_matched_ratio'], d['ratio_null_mean'], d['ratio_null_sd'], d['ratio_t']))
json.dump(out, open(os.path.join(HERE, 'variants.json'), 'w'), indent=1, default=float)
