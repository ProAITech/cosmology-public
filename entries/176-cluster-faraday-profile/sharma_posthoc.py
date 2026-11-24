"""Entry 176 follow-up, post hoc (made after sharma.py's frozen result was seen; nothing here changes its reading).

sharma.py's check (d) showed the mean dRM^2 around the real clusters at 0.6-10 Mpc to be about 80 rad^2/m^4 below that around the
clusters turned in Galactic longitude.  The eRASS1 clusters all lie at l >= 179.5 deg; turning them carries their low-latitude
sight lines into the other half, where the residual RM tail is heavier at |b| = 20-30 deg (mean dRM^2 712 against 438).  So the
random term of the frozen S1 is too high in every bin.  Two other random terms, each with sharma.py's jackknife and shuffles:
  P1  the clusters turned within the eRASS1 half: l' = 179.5 + ((l - 179.5 + s) mod 180.5), s uniform in [10, 170.5] deg,
      200 lists (seed 20261124); latitudes and redshifts unchanged.
  P2  a local background: R_k = the mean dRM^2 over the same clusters' pairs at 2,512-10,000 kpc (one value for every bin).
  tail  the mean dRM^2 at 2,512-10,000 kpc around the real clusters and around the 200 lists of sharma.py (turned over the whole
      sky), with all sources, with |dRM| < 100 and with |dRM| < 50, and (1.4826 MAD)^2: is the offset in the tail?
And the 2 x 2 of sharma.py again, with P1's turned lists as the null.
Run from the repository root:  .venv/bin/python -I entries/176-cluster-faraday-profile/sharma_posthoc.py   -> sharma_posthoc.json
"""
import json, os, sys, time
import numpy as np
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cluster_rm as c
import sharma as s

L0, HALF = 179.5, 180.5
log = s.log


def turn_within(l, shift):
    return L0 + np.mod(l - L0 + shift, HALF)


def s1_with(S, cl, tree, layouts, random_lists=None, local=None, label=''):
    """S1 with a given random term: a list of turned cluster lists, or a local background (edges in kpc)."""
    nb = len(s.EDGES) - 1
    DA = c.ang_diam_Mpc(cl['z'])
    w = S['d'] ** 2
    sD, kD, _ = s.pairs(tree, cl, DA, s.EDGES)
    if random_lists is not None:
        rs, rk, rr = [], [], []
        for i, clr in enumerate(random_lists):
            a, b, _ = s.pairs(tree, clr, DA, s.EDGES)
            rs.append(a); rk.append(b); rr.append(np.full(len(a), i))
        sR, kR, rR = np.concatenate(rs), np.concatenate(rk), np.concatenate(rr)
    else:
        a, _, _ = s.pairs(tree, cl, DA, np.asarray(local, float))
        sR = np.tile(a, nb); kR = np.repeat(np.arange(nb), len(a)); rR = np.zeros(len(sR), int)
    D, R = (sD, kD), (sR, kR)
    Dk, Rk, xi, nD, nR = s.xi_of(D, R, w, nb)
    keep = nD >= s.MINPAIRS
    C = sum(s.jackknife_cov(D, R, w, lab, nb) for lab in layouts) / len(layouts)
    chi2, sn, f, nk = s.snr(xi, C, keep, s.NPATCH)
    g = np.random.default_rng(s.SEED)
    Ck = np.linalg.inv(C[np.ix_(keep, keep)])
    chis = np.empty(s.NSHUF)
    for i in range(s.NSHUF):
        wp = w[g.permutation(len(w))]
        xs = s.xi_of(D, R, wp, nb)[2]
        chis[i] = f * xs[keep] @ Ck @ xs[keep]
    sig = np.sqrt(np.diag(C))
    out = dict(label=label, bins=[dict(lo_kpc=float(s.EDGES[j]), hi_kpc=float(s.EDGES[j + 1]), pairs=int(nD[j]), kept=bool(keep[j]), D=float(Dk[j]),
                                       R=float(Rk[j]), xi=float(xi[j]), sigma_jk=float(sig[j])) for j in range(nb)],
               chi2=chi2, snr=sn, bins_kept=nk, shuffle=dict(n=s.NSHUF, n_at_or_above=int((chis >= chi2).sum()), p=float((chis >= chi2).mean()),
                                                             chi2_mean=float(chis.mean()), chi2_max=float(chis.max())))
    log('%s: R %s; xi %s; sigma %s; chi2 %.2f on %d bins, S/N %.2f; shuffles at or above: %d of %d' % (
        label, np.round(Rk, 1).tolist(), np.round(xi, 1).tolist(), np.round(sig, 1).tolist(), chi2, nk, sn, (chis >= chi2).sum(), s.NSHUF))
    return out, dict(D=(sD, kD), R=(sR, kR, rR), keep=keep)


def main():
    t0 = time.time()
    _, src, _ = c.load_spice(); cl, _ = c.load_clusters()
    S, cnt = c.sample(src, c.CFG['fg_radius'], c.CFG['min_neigh'])
    assert cnt['kept'] == 92169 and len(cl['z']) == 9050 and (cl['l'] >= L0).all()
    tree = cKDTree(c.gal_unit(S['l'], S['b']))
    shifts = np.random.default_rng(s.SEED).uniform(10.0, HALF - 10.0, s.NROT)
    lists = [dict(cl, l=turn_within(cl['l'], sh)) for sh in shifts]
    layouts = s.layouts_for(S, s.NLAYOUT)
    out = dict(generated=time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime()), note='post hoc; sharma.py holds the frozen result')
    d2 = S['d'] ** 2; west = S['l'] >= L0
    out['halves'] = {}
    for name, m in (('eRASS1_half', west), ('other_half', ~west)):
        out['halves'][name] = dict(n=int(m.sum()), mean_d2=float(d2[m].mean()), by_abs_b={})
        for lo, hi in ((20, 30), (30, 45), (45, 60), (60, 90)):
            mm = m & (np.abs(S['b']) >= lo) & (np.abs(S['b']) < hi)
            out['halves'][name]['by_abs_b']['%d-%d' % (lo, hi)] = dict(n=int(mm.sum()), mean_d2=float(d2[mm].mean()))
    DA = c.ang_diam_Mpc(cl['z']); E = np.array([2512.0, 10000.0])
    full = np.random.default_rng(s.SEED).uniform(10.0, 350.0, s.NROT)              # sharma.py's angles
    si = s.pairs(tree, cl, DA, E)[0]; rr = [s.pairs(tree, cl, DA, E, shift=a)[0] for a in full]
    d = S['d']; out['tail'] = {}
    for name, m in (('all', np.ones(len(d), bool)), ('abs_dRM_below_100', np.abs(d) < 100), ('abs_dRM_below_50', np.abs(d) < 50)):
        rot = np.array([np.mean(d[r][m[r]] ** 2) for r in rr])
        out['tail'][name] = dict(real=float(np.mean(d[si][m[si]] ** 2)), turned_mean=float(rot.mean()), turned_sd=float(rot.std(ddof=1)),
                                 z=float((np.mean(d[si][m[si]] ** 2) - rot.mean()) / rot.std(ddof=1)))
    rotm = np.array([s.mad2(d[r]) for r in rr])
    out['tail']['mad2'] = dict(real=float(s.mad2(d[si])), turned_mean=float(rotm.mean()), turned_sd=float(rotm.std(ddof=1)))
    out['tail']['share_abs_dRM_ge_100'] = dict(real=float(np.mean(np.abs(d[si]) >= 100)), turned=float(np.mean(np.concatenate([np.abs(d[r]) >= 100 for r in rr]))))
    out['tail']['top_0p1pct_share_of_sum_dRM2'] = float(np.sort(d ** 2)[::-1][:int(round(0.001 * len(d)))].sum() / (d ** 2).sum())
    log('tail at 2.5-10 Mpc: %s' % json.dumps({k: (round(v['real'], 1), round(v['turned_mean'], 1), round(v['turned_sd'], 1)) for k, v in out['tail'].items() if isinstance(v, dict) and 'turned_mean' in v}))
    out['P1'], P = s1_with(S, cl, tree, layouts, random_lists=lists, label='P1 turned within the eRASS1 half')
    out['P2'], _ = s1_with(S, cl, tree, layouts, local=[2512.0, 10000.0], label='P2 local background 2.5-10 Mpc')
    # the 2 x 2 with P1's lists as the null
    nb = len(s.EDGES) - 1
    d, e = S['d'], S['e']
    (sD, kD), (sR, kR, rR), keep = P['D'], P['R'], P['keep']
    nrot = len(lists)
    rot_n = np.zeros((nrot, nb), int); rot_u = np.full((nrot, nb), np.nan); rot_m = np.full((nrot, nb), np.nan)
    for i in range(nrot):
        for j in range(nb):
            ss = sR[(rR == i) & (kR == j)]
            rot_n[i, j] = len(ss)
            if len(ss):
                rot_u[i, j] = np.mean(d[ss] ** 2) - np.mean(e[ss] ** 2)
            if len(ss) >= s.MIN_ROT_MAD:
                rot_m[i, j] = s.mad2(d[ss]) - np.mean(e[ss] ** 2)
    use = keep & (rot_n.min(0) >= s.MIN_ROT_MAD)
    tb = dict(rotated_pairs_min_per_bin=rot_n.min(0).tolist())
    for name, fn in (('unclipped', lambda x, ee: np.mean(x ** 2) - np.mean(ee ** 2)), ('mad', lambda x, ee: s.mad2(x) - np.mean(ee ** 2))):
        real = np.array([fn(d[sD[kD == j]], e[sD[kD == j]]) if keep[j] else np.nan for j in range(nb)])
        rot = rot_u if name == 'unclipped' else rot_m
        mu = rot[:, use].mean(0); Cr = np.atleast_2d(np.cov(rot[:, use], rowvar=False, ddof=1)); sd = np.sqrt(np.diag(Cr))
        delta = real[use] - mu; nbk = int(use.sum()); f = (nrot - nbk - 2) / (nrot - 1)
        chi2 = float(f * delta @ np.linalg.solve(Cr, delta))
        tb['physical_' + name] = dict(bins_used=np.nonzero(use)[0].tolist(), excess=delta.tolist(), t_per_bin=(delta / sd).tolist(),
                                      snr_rot=float(np.sqrt(max(chi2 - nbk, 0))))
        log('2x2 (P1 null) physical %-9s bins %s: excess %s, t %s, S/N_rot %.2f' % (name, np.nonzero(use)[0].tolist(), np.round(delta, 1).tolist(),
                                                                                 np.round(delta / sd, 2).tolist(), tb['physical_' + name]['snr_rot']))
    c.CFG['xmax'] = 1.0
    n = len(d)
    x, _ = c.scaled_distance(tree, n, cl); A = x < 1.0
    vals = {k: np.empty(nrot) for k in ('unclipped', 'mad')}
    for i, clr in enumerate(lists):
        xr, _ = c.scaled_distance(tree, n, clr); Ar = xr < 1.0
        vals['unclipped'][i] = np.mean(d[Ar] ** 2) - np.mean(e[Ar] ** 2); vals['mad'][i] = s.mad2(d[Ar]) - np.mean(e[Ar] ** 2)
    real = dict(unclipped=float(np.mean(d[A] ** 2) - np.mean(e[A] ** 2)), mad=float(s.mad2(d[A]) - np.mean(e[A] ** 2)))
    for k in ('unclipped', 'mad'):
        v = vals[k]
        tb['rangeA_' + k] = dict(value=real[k], null_mean=float(v.mean()), null_sd=float(v.std(ddof=1)), excess=float(real[k] - v.mean()),
                                 t=float((real[k] - v.mean()) / v.std(ddof=1)))
        log('2x2 (P1 null) range A %-9s: %.1f against %.1f +- %.1f, excess %+.1f, t %+.2f' % (k, real[k], v.mean(), v.std(ddof=1), real[k] - v.mean(),
                                                                                         (real[k] - v.mean()) / v.std(ddof=1)))
    tb['amplitude_ratio_rangeA'] = tb['rangeA_unclipped']['excess'] / tb['rangeA_mad']['excess']
    out['two_by_two_P1'] = tb
    out['seconds'] = round(time.time() - t0)
    json.dump(s.clean(out), open(os.path.join(HERE, 'sharma_posthoc.json'), 'w'), indent=1, allow_nan=False)
    log('sharma_posthoc.json written (%d s)' % out['seconds'])


if __name__ == '__main__':
    main()
