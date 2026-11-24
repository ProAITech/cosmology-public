"""Entry 176 follow-up: the estimator of Sharma, Ravi & Krause (2026, arXiv:2610.37646, Sections II.1, III.1, Appendix A) on
Entry 176's sample.

Rule: Workbench section 176, "The estimator of Sharma, Ravi and Krause on this sample", fixed 24 November 2026, 17:29 UTC
(commit 11b158a) and published as Workbench version 203 before anything here was computed.  Asked for by the review of the Entry 176
manuscript (issue #2).

  S1      every source-cluster pair at projected physical separation R = theta D_A(z), 10-1,000 kpc, five logarithmic bins;
          D_k = mean dRM^2 over the pairs in bin k (no clipping, no noise term); random positions = the cluster list turned in Galactic
          longitude by 200 angles in [10, 350) deg (seed 20261124); xi_k = D_k - R_k.
  errors  jackknife over 100 k-means patches of the sources (a pair belongs to its source's patch), averaged over 100 layouts
          (seeds 1-100); S/N = sqrt(f_H xi C^-1 xi - N_bin), f_H = (100 - N_bin - 2)/99; bins with < 30 real pairs left out.
  null    dRM shuffled among the sources 1,000 times; chi2 = f_H xi C^-1 xi with the data's C.
  2 x 2   unclipped mean of dRM^2 or (1.4826 MAD)^2, each minus mean e^2; Sharma's bins over pairs or Entry 176's range A (x < 1);
          all against the same 200 turned cluster lists.
  checks  (a) |b| > 25 deg for sources and clusters; (b) pair counts and mean e^2; (c) share of xi_k carried by the largest 1% of
          dRM^2 in bin k; (d) the profile to 10 Mpc in ten logarithmic bins, 20 turned lists, one patch layout.

Data: as cluster_rm.py (SPICE-RACS DR2 through Entry 171's cache; eRASS1 primary catalogue v3.2).  Reads no other file.
Run from the repository root:
  .venv/bin/python -I entries/176-cluster-faraday-profile/sharma.py synthetic   # random dRM at the real positions; checks the code
  .venv/bin/python -I entries/176-cluster-faraday-profile/sharma.py run         # the estimator; writes sharma.json
"""
import json, os, sys, time
from itertools import chain
import numpy as np
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cluster_rm as c

SEED = 20261124
EDGES = np.logspace(1.0, 3.0, 6)            # 10, 25.1, 63.1, 158.5, 398.1, 1000 kpc
EDGES_WIDE = np.logspace(1.0, 4.0, 11)      # check (d): 10 kpc to 10 Mpc
NROT, NSHUF, NPATCH, NLAYOUT, MINPAIRS = 200, 1000, 100, 100, 30
MIN_ROT_MAD = 10                            # implementation choice: a turned list gives a MAD in a bin only with >= 10 pairs there
KMEANS_ITERS = 25
log = lambda s: print(s, flush=True)


# ---------------------------------------------------------------- pairs
def pairs(tree, cl, DA, edges, shift=0.0):
    """All source-cluster pairs with edges[0] <= R < edges[-1] (kpc, physical): source index, bin index, cluster index."""
    v = c.gal_unit(np.mod(cl['l'] + shift, 360.0), cl['b'])
    th = edges[-1] / 1000.0 / DA / c.D2R                                 # degrees
    lists = tree.query_ball_point(v, c.chord(th))
    n = np.fromiter((len(L) for L in lists), int, len(lists))
    if n.sum() == 0:
        return np.zeros(0, int), np.zeros(0, int), np.zeros(0, int)
    si = np.fromiter(chain.from_iterable(lists), int, int(n.sum()))
    ci = np.repeat(np.arange(len(lists)), n)
    sep = 2 * np.arcsin(np.clip(np.linalg.norm(tree.data[si] - v[ci], axis=1) / 2, 0, 1))
    R = sep * DA[ci] * 1000.0
    k = np.searchsorted(edges, R, side='right') - 1
    ok = (k >= 0) & (k < len(edges) - 1)
    return si[ok], k[ok], ci[ok]


def kmeans_labels(X, k, seed, iters=KMEANS_ITERS):
    """Spherical k-means on unit vectors (Lloyd); returns the patch of each point."""
    g = np.random.default_rng(seed)
    C = X[g.choice(len(X), k, replace=False)].copy()
    for _ in range(iters):
        lab = np.argmax(X @ C.T, axis=1)
        M = np.stack([np.bincount(lab, weights=X[:, j], minlength=k) for j in range(3)], axis=1)
        nrm = np.linalg.norm(M, axis=1)
        good = nrm > 0
        C[good] = M[good] / nrm[good, None]
    return np.argmax(X @ C.T, axis=1)


def patch_sums(si, k, w, lab, nb):
    idx = lab[si] * nb + k
    S = np.bincount(idx, weights=w[si], minlength=NPATCH * nb).reshape(NPATCH, nb)
    N = np.bincount(idx, minlength=NPATCH * nb).reshape(NPATCH, nb).astype(float)
    return S, N


def jackknife_cov(D, R, w, lab, nb):
    """Jackknife covariance of xi = D - R for one patch layout.  D, R = (si, k) pair lists."""
    SD, ND = patch_sums(D[0], D[1], w, lab, nb)
    SR, NR = patch_sums(R[0], R[1], w, lab, nb)
    with np.errstate(invalid='ignore', divide='ignore'):
        xi = (SD.sum(0) - SD) / (ND.sum(0) - ND) - (SR.sum(0) - SR) / (NR.sum(0) - NR)
    m = xi.mean(0)
    return (NPATCH - 1) / NPATCH * (xi - m).T @ (xi - m)


def xi_of(D, R, w, nb):
    nD = np.bincount(D[1], minlength=nb); nR = np.bincount(R[1], minlength=nb)
    with np.errstate(invalid='ignore', divide='ignore'):
        Dk = np.bincount(D[1], weights=w[D[0]], minlength=nb) / nD
        Rk = np.bincount(R[1], weights=w[R[0]], minlength=nb) / nR
    return Dk, Rk, Dk - Rk, nD, nR


def snr(xi, C, keep, nsamp):
    nb = int(keep.sum())
    f = (nsamp - nb - 2) / (nsamp - 1)
    x = xi[keep]
    chi2 = float(f * x @ np.linalg.solve(C[np.ix_(keep, keep)], x))
    return chi2, float(np.sqrt(max(chi2 - nb, 0.0))), f, nb


# ---------------------------------------------------------------- the estimator
def estimator(S, cl, edges, angles, layouts, tree, label, shuffles=0, keep_rot_pairs=False):
    nb = len(edges) - 1
    DA = c.ang_diam_Mpc(cl['z'])
    w = S['d'] ** 2
    sD, kD, cD = pairs(tree, cl, DA, edges)
    rs, rk, rr = [], [], []
    for i, a in enumerate(angles):
        s_, k_, _ = pairs(tree, cl, DA, edges, shift=a)
        rs.append(s_); rk.append(k_); rr.append(np.full(len(s_), i))
    sR, kR, rR = np.concatenate(rs), np.concatenate(rk), np.concatenate(rr)
    D, R = (sD, kD), (sR, kR)
    Dk, Rk, xi, nD, nR = xi_of(D, R, w, nb)
    keep = nD >= MINPAIRS
    C = np.zeros((nb, nb)); C_layouts = []
    for lab in layouts:
        Ci = jackknife_cov(D, R, w, lab, nb)
        C += Ci / len(layouts)
        C_layouts.append(Ci)
    chi2, sn, f, nk = snr(xi, C, keep, NPATCH)
    eD = np.bincount(kD, weights=S['e'][sD] ** 2, minlength=nb) / np.maximum(nD, 1)
    eR = np.bincount(kR, weights=S['e'][sR] ** 2, minlength=nb) / np.maximum(nR, 1)
    sig = np.sqrt(np.diag(C))
    bins = []
    for j in range(nb):
        m = kD == j
        top = None
        if m.sum() >= MINPAIRS:
            v = np.sort(w[sD[m]])[::-1]; ntop = int(np.ceil(0.01 * len(v)))
            Dwo = v[ntop:].mean()
            top = dict(n_top=ntop, share=float((Dk[j] - Dwo) / xi[j]) if xi[j] != 0 else None, D_without_top=float(Dwo))
        bins.append(dict(lo_kpc=float(edges[j]), hi_kpc=float(edges[j + 1]), pairs=int(nD[j]), random_pairs_per_list=float(nR[j] / len(angles)),
                         kept=bool(keep[j]), D=float(Dk[j]), R=float(Rk[j]), xi=float(xi[j]), sigma_jk=float(sig[j]),
                         xi_over_sigma=float(xi[j] / sig[j]) if sig[j] > 0 else None, mean_e2_real=float(eD[j]), mean_e2_random=float(eR[j]),
                         clusters_with_pairs=int(len(np.unique(cD[m]))), top1pct=top))
    out = dict(label=label, sources=int(len(w)), clusters=int(len(DA)), bins=bins, bins_kept=nk, hartlap=f, chi2=chi2, snr=sn,
               snr_by_layout=None)
    # spread of S/N over single layouts, for the record
    sl = [snr(xi, Ci, keep, NPATCH)[1] for Ci in C_layouts]
    out['snr_by_layout'] = dict(min=float(min(sl)), median=float(np.median(sl)), max=float(max(sl)))
    log('%s: %d sources, %d clusters; pairs %s; xi %s; sigma %s; chi2 %.2f on %d bins, S/N %.2f (layouts %.2f-%.2f)' % (
        label, len(w), len(DA), nD.tolist(), np.round(xi, 1).tolist(), np.round(sig, 1).tolist(), chi2, nk, sn, min(sl), max(sl)))
    if shuffles:
        g = np.random.default_rng(SEED)
        Ck = np.linalg.inv(C[np.ix_(keep, keep)])
        chis, xis = np.empty(shuffles), np.empty((shuffles, nb))
        for i in range(shuffles):
            wp = w[g.permutation(len(w))]
            _, _, xs, _, _ = xi_of(D, R, wp, nb)
            xis[i] = xs
            chis[i] = f * xs[keep] @ Ck @ xs[keep]
        out['shuffle'] = dict(n=shuffles, chi2_mean=float(chis.mean()), chi2_sd=float(chis.std(ddof=1)), chi2_max=float(chis.max()),
                              n_at_or_above=int((chis >= chi2).sum()), p=float((chis >= chi2).mean()),
                              xi_mean=[float(v) for v in xis.mean(0)], xi_sd=[float(v) for v in xis.std(0, ddof=1)])
        log('  shuffles: chi2 %.2f +- %.2f (max %.2f); %d of %d at or above %.2f' % (chis.mean(), chis.std(ddof=1), chis.max(),
                                                                                  (chis >= chi2).sum(), shuffles, chi2))
    if keep_rot_pairs:
        out['_pairs'] = dict(D=(sD, kD), R=(sR, kR, rR), keep=keep)
    return out


def reading_S1(r):
    sn, p, nab = r['snr'], r['shuffle']['p'], r['shuffle']['n_at_or_above']
    if sn >= 5 and nab == 0:
        return 'detects the cluster signal on this sample'
    if sn >= 3 and p < 0.01:
        return 'sees it weakly'
    return 'does not detect it'


# ---------------------------------------------------------------- the 2 x 2
def mad2(d):
    m = np.median(d)
    return (1.4826 * np.median(np.abs(d - m))) ** 2


def two_by_two(S, cl, angles, tree, P):
    nb = len(EDGES) - 1
    d, e = S['d'], S['e']
    (sD, kD), (sR, kR, rR), keep = P['D'], P['R'], P['keep']
    nrot = len(angles)
    real_u, real_m = np.full(nb, np.nan), np.full(nb, np.nan)
    rot_u, rot_m = np.full((nrot, nb), np.nan), np.full((nrot, nb), np.nan)
    rot_n = np.zeros((nrot, nb), int)
    for j in range(nb):
        m = kD == j
        if m.sum() >= MINPAIRS:
            real_u[j] = np.mean(d[sD[m]] ** 2) - np.mean(e[sD[m]] ** 2)
            real_m[j] = mad2(d[sD[m]]) - np.mean(e[sD[m]] ** 2)
    order = np.lexsort((kR, rR))
    sR, kR, rR = sR[order], kR[order], rR[order]
    key = rR * nb + kR
    starts = np.searchsorted(key, np.arange(nrot * nb)); ends = np.searchsorted(key, np.arange(nrot * nb), side='right')
    for i in range(nrot):
        for j in range(nb):
            a, b = starts[i * nb + j], ends[i * nb + j]
            rot_n[i, j] = b - a
            if b - a >= 1:
                ss = sR[a:b]
                rot_u[i, j] = np.mean(d[ss] ** 2) - np.mean(e[ss] ** 2)
            if b - a >= MIN_ROT_MAD:
                ss = sR[a:b]
                rot_m[i, j] = mad2(d[ss]) - np.mean(e[ss] ** 2)
    use = keep & (rot_n.min(0) >= MIN_ROT_MAD)                           # same bins for both rows of this column
    res = {}
    if not use.any():
        raise SystemExit('no bin has enough pairs in every turned list')
    for name, real, rot in (('unclipped', real_u, rot_u), ('mad', real_m, rot_m)):
        mu = rot[:, use].mean(0); Cr = np.cov(rot[:, use], rowvar=False, ddof=1)
        delta = real[use] - mu
        nbk = int(use.sum()); f = (nrot - nbk - 2) / (nrot - 1)
        chi2 = float(f * delta @ np.linalg.solve(np.atleast_2d(Cr), delta))
        sd = np.sqrt(np.diag(np.atleast_2d(Cr)))
        res['physical_' + name] = dict(bins_used=[int(j) for j in np.nonzero(use)[0]], value=[float(v) for v in real[use]],
                                       null_mean=[float(v) for v in mu], null_sd=[float(v) for v in sd], excess=[float(v) for v in delta],
                                       t_per_bin=[float(v) for v in delta / sd], chi2=chi2, snr_rot=float(np.sqrt(max(chi2 - nbk, 0))), hartlap=f)
        log('2x2 physical %-9s bins %s: excess %s, t %s, S/N_rot %.2f' % (name, np.nonzero(use)[0].tolist(), np.round(delta, 1).tolist(),
                                                                     np.round(delta / sd, 2).tolist(), res['physical_' + name]['snr_rot']))
    res['rotated_pairs_min_per_bin'] = rot_n.min(0).tolist()
    # Entry 176's range A, nearest cluster, each source once.  x < 1 needs only clusters within one theta500, so xmax = 1 gives
    # the same x for every source with x < 1 (and is faster).
    c.CFG['xmax'] = 1.0
    n = len(d)
    x, _ = c.scaled_distance(tree, n, cl)
    A = x < 1.0
    real = dict(unclipped=float(np.mean(d[A] ** 2) - np.mean(e[A] ** 2)), mad=float(mad2(d[A]) - np.mean(e[A] ** 2)))
    check = dict(n=int(A.sum()), V=c.V_of(d[A], e[A]), entry176_V=179.28, entry176_n=589)
    log('range A: %d sources, V %.2f (Entry 176: 179.28 on 589)' % (A.sum(), check['V']))
    rots = dict(unclipped=np.empty(nrot), mad=np.empty(nrot)); nA = np.empty(nrot, int)
    for i, a in enumerate(angles):
        xr, _ = c.scaled_distance(tree, n, cl, shift=a)
        Ar = xr < 1.0; nA[i] = Ar.sum()
        rots['unclipped'][i] = np.mean(d[Ar] ** 2) - np.mean(e[Ar] ** 2)
        rots['mad'][i] = mad2(d[Ar]) - np.mean(e[Ar] ** 2)
    for name in ('unclipped', 'mad'):
        v = rots[name]
        res['rangeA_' + name] = dict(value=real[name], null_mean=float(v.mean()), null_sd=float(v.std(ddof=1)),
                                     excess=float(real[name] - v.mean()), t=float((real[name] - v.mean()) / v.std(ddof=1)), n=int(A.sum()),
                                     null_n_mean=float(nA.mean()))
        log('2x2 range A %-9s: %.1f against %.1f +- %.1f, excess %+.1f, t %+.2f' % (name, real[name], v.mean(), v.std(ddof=1),
                                                                                  real[name] - v.mean(), (real[name] - v.mean()) / v.std(ddof=1)))
    res['rangeA_check'] = check
    res['amplitude_ratio_rangeA'] = res['rangeA_unclipped']['excess'] / res['rangeA_mad']['excess']
    c.CFG['xmax'] = 8.0
    return res


# ---------------------------------------------------------------- main
def clean(o):
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


def load():
    _, src, sinfo = c.load_spice()
    cl, cinfo = c.load_clusters()
    S, cnt = c.sample(src, c.CFG['fg_radius'], c.CFG['min_neigh'])
    return S, cl, dict(spice=sinfo, clusters=cinfo, counts=cnt)


def layouts_for(S, n):
    X = c.gal_unit(S['l'], S['b'])
    t0 = time.time()
    L = [kmeans_labels(X, NPATCH, seed) for seed in range(1, n + 1)]
    log('%d patch layouts in %.0f s; patch sizes (first layout) %d-%d sources' % (n, time.time() - t0, np.bincount(L[0]).min(), np.bincount(L[0]).max()))
    return L


def synthetic():
    """Random dRM at the real positions: Gaussian, sd 8.5 rad/m^2, plus sd 10 inside 500 kpc of a real cluster.  No real dRM is used."""
    S, cl, info = load()
    g = np.random.default_rng(7)
    tree = cKDTree(c.gal_unit(S['l'], S['b']))
    DA = c.ang_diam_Mpc(cl['z'])
    sD, kD, _ = pairs(tree, cl, DA, np.array([10.0, 500.0]))
    near = np.zeros(len(S['d']), bool); near[np.unique(sD)] = True
    S = dict(S); S['d'] = g.normal(0, 8.5, len(S['d'])) + np.where(near, g.normal(0, 10.0, len(S['d'])), 0.0)
    log('synthetic: %d sources inside 500 kpc of a cluster carry +100 rad2/m4' % near.sum())
    angles = np.random.default_rng(SEED).uniform(10.0, 350.0, 40)
    L = layouts_for(S, 5)
    r = estimator(S, cl, EDGES, angles, L, tree, 'synthetic, signal', shuffles=200)
    S['d'] = g.normal(0, 8.5, len(S['d']))
    r0 = estimator(S, cl, EDGES, angles, L, tree, 'synthetic, no signal', shuffles=200)
    log('expected: xi about +100 in bins inside 500 kpc (diluted in the 398-1000 kpc bin), and S/N near 0 with no signal')


def run():
    t0 = time.time()
    S, cl, info = load()
    log('sample %s; clusters kept %d' % (json.dumps(info['counts']), info['clusters']['kept']))
    assert info['counts']['kept'] == 92169 and info['clusters']['kept'] == 9050
    tree = cKDTree(c.gal_unit(S['l'], S['b']))
    angles = np.random.default_rng(SEED).uniform(10.0, 350.0, NROT)
    L = layouts_for(S, NLAYOUT)
    out = dict(generated=time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime()), rule='Workbench section 176, follow-up rule of 24 November 2026, 17:29 UTC (commit 11b158a, version 203)',
               config=dict(seed=SEED, edges_kpc=EDGES.tolist(), nrot=NROT, nshuffle=NSHUF, npatch=NPATCH, nlayout=NLAYOUT, min_pairs=MINPAIRS,
                           min_rot_pairs_mad=MIN_ROT_MAD, kmeans_iters=KMEANS_ITERS), info=info)
    S1 = estimator(S, cl, EDGES, angles, L, tree, 'S1', shuffles=NSHUF, keep_rot_pairs=True)
    P = S1.pop('_pairs')
    S1['reading'] = reading_S1(S1)
    log('S1 reading: their estimator %s' % S1['reading'])
    out['S1'] = S1
    out['two_by_two'] = two_by_two(S, cl, angles, tree, P)
    # (a) their latitude cut
    sb = np.abs(S['b']) > 25.0; cb = np.abs(cl['b']) > 25.0
    Sa = {k: v[sb] for k, v in S.items()}; cla = {k: v[cb] for k, v in cl.items()}
    treea = cKDTree(c.gal_unit(Sa['l'], Sa['b']))
    La = [lab[sb] for lab in L]
    out['check_a_b25'] = estimator(Sa, cla, EDGES, angles, La, treea, 'check (a) |b| > 25')
    # (d) out to 10 Mpc
    out['check_d_10Mpc'] = estimator(S, cl, EDGES_WIDE, angles[:20], L[:1], tree, 'check (d) to 10 Mpc')
    out['seconds'] = round(time.time() - t0)
    json.dump(clean(out), open(os.path.join(HERE, 'sharma.json'), 'w'), indent=1, allow_nan=False)
    log('sharma.json written (%d s)' % out['seconds'])


if __name__ == '__main__':
    {'synthetic': synthetic, 'run': run}[sys.argv[1]]()
