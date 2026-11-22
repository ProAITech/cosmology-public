"""Entry 174, exploratory: the union of all sight lines (Entry 171's sample, the DESI-only sight lines, the new Quaia ones).

Made after the frozen results of Entries 171, 173 and 174 were seen.  Nothing here decides anything.  The union itself was a
reported-only item of Entry 174's rule; this script asks whether its growth survives the checks Entry 171's rule used, and whether
mixing three samples with different skies and redshift ranges could have made it.
  cells_by_sample   redshifts shuffled within sample x latitude band x e_RM tercile (36 cells), so that no sample's redshift
                    distribution is carried to another sample's sky
  checks            each Galactic hemisphere; the smaller-e_RM half; |b| >= 40 deg; four longitude sectors in the calibration
  variance          per bin, sigma^2 minus the mean e_RM^2
Run from the repository root:  .venv/bin/python -I entries/174-quaia-replication/explore_union.py     (about 15 minutes)
"""
import json, os, sys, time
import numpy as np
from scipy.spatial import cKDTree
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import quaia_growth as q

log = lambda s: print(s, flush=True)
rows, src, _ = q.load_spice()
quaia, _ = q.load_quaia()
mq = q.load_milliquas()
desi = q.load_desi()
tree = cKDTree(q.unit(src['ra'], src['dec']))
used = np.zeros(len(src['rm']), bool)
m_mq = q.match(mq, rows); used[m_mq['src']] = True
m_de = q.match(desi, rows); used[m_de['src']] = True
S, _ = q.build(q.match(quaia, rows), src, quaia, q.CFG['fg_radius'], q.CFG['min_neigh'], tree, exclude=used)
isQ = mq['cls'] == 'Q'
mq_Q = {k: v[isQ] for k, v in mq.items()}
M, _ = q.build(q.match(mq_Q, rows), src, mq_Q, q.CFG['fg_radius'], q.CFG['min_neigh'], tree)
Dd, _ = q.build(m_de, src, desi, q.CFG['fg_radius'], q.CFG['min_neigh'], tree)
donly = ~np.isin(Dd['src'], M['src'])
keys = ('z', 'd', 'e', 'b', 'l')
U = {k: np.r_[M[k], Dd[k][donly], S[k]] for k in keys}
U['sample'] = np.r_[np.zeros(len(M['z']), int), np.ones(int(donly.sum()), int), np.full(len(S['z']), 2)]
log('union %d = %d + %d + %d' % (len(U['z']), len(M['z']), donly.sum(), len(S['z'])))


def analyse_cells(Sx, cell, nshuf, label, seed=q.CFG['seed']):
    """q.analyse with the calibration cells given."""
    z, d = Sx['z'], Sx['d']
    real = q.bins_fit(z, d, np.random.default_rng(seed), Sx['e'])
    groups = [np.nonzero(cell == c)[0] for c in np.unique(cell)]
    R = np.random.default_rng(seed + 7)
    nulls = np.empty(nshuf)
    for s in range(nshuf):
        z2 = z.copy()
        for g in groups:
            z2[g] = z[g][R.permutation(len(g))]
        nulls[s] = q.bins_fit(z2, d, np.random.default_rng(seed + 1000 + s))['b']
    mu, sd = float(nulls.mean()), float(nulls.std(ddof=1))
    out = dict(label=label, n=int(len(z)), cells=len(groups), b=real['b'], null_mean=mu, null_sd=sd, t=float((real['b'] - mu) / sd),
               net=float(real['b'] - mu), p_upper=float((1 + (nulls >= real['b']).sum()) / (1 + nshuf)),
               bins=dict(zmed=real['zmed'], sigma=real['sigma'], err=real['err'], e2mean=real['e2mean'], n=real['n']))
    log('%s: n %d, cells %d, b %+.3f, null %+.3f +- %.3f, t %+.2f' % (label, out['n'], out['cells'], out['b'], mu, sd, out['t']))
    return out


out = {}
base = q.cells_for(U)
out['union_12_cells'] = analyse_cells(U, base, 1000, 'union, 12 cells (as reported in Entry 174)')
out['union_cells_by_sample'] = analyse_cells(U, base + 12 * U['sample'], 1000, 'union, cells also split by sample (36)')
sect = np.minimum(3, np.floor(U['l'] / 90).astype(int))
out['union_longitude_sectors'] = analyse_cells(U, base * 4 + sect, 300, 'union, longitude sectors (48 cells)')
for name, sel in (('north_galactic', U['b'] > 0), ('south_galactic', U['b'] < 0), ('smaller_eRM_half', U['e'] <= np.median(U['e'])),
                  ('larger_eRM_half', U['e'] > np.median(U['e'])), ('b_ge_40', np.abs(U['b']) >= 40), ('b_lt_40', np.abs(U['b']) < 40)):
    Sx = q.sub(U, sel)
    out[name] = analyse_cells(Sx, q.cells_for(Sx), 300, 'union, ' + name)
# independent samples, each with its own calibration (1,000 shuffles): the pieces of the union
parts = {}
for name, k in (('entry171', 0), ('desi_only', 1), ('quaia_new', 2)):
    Sx = q.sub(U, U['sample'] == k)
    parts[name] = analyse_cells(Sx, q.cells_for(Sx), 1000, name)
out['parts'] = parts
beta = 0.8346420362367484                                  # Entry 174's dilution factor for Quaia redshifts
nets = np.array([parts['entry171']['net'], parts['desi_only']['net'], parts['quaia_new']['net'] / beta])
sds = np.array([parts['entry171']['null_sd'], parts['desi_only']['null_sd'], parts['quaia_new']['null_sd'] / beta])
w = 1 / sds ** 2
out['combined_independent'] = dict(net_slopes=nets.tolist(), sds=sds.tolist(), weighted_mean=float((w * nets).sum() / w.sum()),
                                   sd=float(1 / np.sqrt(w.sum())), z=float((w * nets).sum() / np.sqrt(w.sum())),
                                   chi2_consistency=float((w * (nets - (w * nets).sum() / w.sum()) ** 2).sum()),
                                   stouffer=float(sum(p['t'] for p in parts.values()) / np.sqrt(3)),
                                   replications_only=dict(weighted_mean=float((w[1:] * nets[1:]).sum() / w[1:].sum()),
                                                          sd=float(1 / np.sqrt(w[1:].sum()))))
log(json.dumps(out['combined_independent']))
b = out['union_12_cells']['bins']
var = [s * s - e2 for s, e2 in zip(b['sigma'], b['e2mean'])]
verr = [2 * s * e for s, e in zip(b['sigma'], b['err'])]
out['variance'] = dict(zmed=b['zmed'], var_noise_removed=var, err=verr,
                       growth_first2_to_last2=float(np.mean(var[-2:]) - np.mean(var[:2])),
                       growth_err=float(0.5 * np.sqrt(verr[0] ** 2 + verr[1] ** 2 + verr[-1] ** 2 + verr[-2] ** 2)))
log(json.dumps(out['variance']))
json.dump(out, open(os.path.join(HERE, 'explore_union.json'), 'w'), indent=1, default=float)
os.makedirs(q.WORK, exist_ok=True)
np.savez(os.path.join(q.WORK, 'union.npz'), **U)
