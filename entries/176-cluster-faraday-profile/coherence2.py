"""Entry 176, analysis after the frozen result: the reversal length from source sizes and from pairs of components.

coherence.py showed that the reversal length Lambda is at least about the projected size of a source's polarized emission.
This script (1) reads the sizes, and (2) bounds Lambda from above with pairs of components of the same source.

  extract   one pass over the SPICE-RACS table; keeps, for rows passing Entry 171's cuts, the component identifier, position,
            RM, deconvolved major axis (dc_maj_axis), fitted major axis, beam, fractional polarization and flux.
  sizes     the deconvolved major axes of sources at 0.5-1 R500 and in the field; the bound of coherence.py in kpc, source by
            source, using each source's own size and its cluster's redshift; and the fractional polarization (flux-matched)
            of the smaller and larger halves, inside clusters against rotated clusters.
  pairs     for sources with two components at least 30 arcsec apart (so that the 25 arcsec beam does not blend them), the
            difference of the two RMs.  The Galactic foreground cancels in the difference.  If the cluster's rotation were
            uncorrelated between the two sight lines, the variance of the difference would rise by twice the single-sight-line
            excess; if the screen were coherent across the pair, it would not rise at all.  Tested against rotated clusters.
Run from the repository root:
  .venv/bin/python -I entries/176-cluster-faraday-profile/coherence2.py extract     # about 1 minute
  .venv/bin/python -I entries/176-cluster-faraday-profile/coherence2.py run         # about 3 minutes
"""
import gzip, json, os, sys, time
import numpy as np
from scipy.spatial import cKDTree
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cluster_rm as c

COMP = os.path.join(c.WORK, 'components.npz')
log = lambda s: print(s, flush=True)


def extract():
    t0 = time.time()
    keep = {k: [] for k in ('source_id', 'cat_id', 'ra', 'dec', 'rm', 'rm_err', 'dc_maj', 'maj', 'beam_maj', 'fracpol', 'flux', 'blended')}
    with gzip.open(c.SPICE, 'rb') as f:
        c.read_header(f); h = c.read_header(f)
        dt, rowbytes = c.table_dtype(h); nrows = int(c.card(h, 'NAXIS2')); done = 0
        while done < nrows:
            n = min(100_000, nrows - done)
            a = np.frombuffer(f.read(n * rowbytes), dtype=dt)
            clear = np.ones(n, bool)
            for k in c.FLAGS:
                clear &= a[k] == b'F'
            rm, e = a['rm'].astype('<f8'), a['rm_err'].astype('<f8')
            ok = clear & np.isfinite(rm) & np.isfinite(e) & (e > 0)
            keep['source_id'].append(a['source_id'][ok].copy()); keep['cat_id'].append(a['cat_id'][ok].copy())
            for k, col in (('ra', 'ra'), ('dec', 'dec'), ('dc_maj', 'dc_maj_axis'), ('maj', 'maj_axis'), ('beam_maj', 'beam_maj'), ('fracpol', 'fracpol'), ('flux', 'total_I_flux')):
                keep[k].append(a[col][ok].astype('<f8'))
            keep['rm'].append(rm[ok]); keep['rm_err'].append(e[ok]); keep['blended'].append(a['is_blended_flag'][ok] == b'T')
            done += n
    os.makedirs(c.WORK, exist_ok=True)
    np.savez(COMP, **{k: np.concatenate(v) for k, v in keep.items()})
    log('extract done, %d rows, %.0f s' % (sum(len(v) for v in keep['rm']), time.time() - t0))


def run():
    z = np.load(COMP)
    sid, cid = z['source_id'], z['cat_id']
    # ---- components: one per cat_id (repeat measurements of the same component averaged with inverse-variance weights)
    uc, cinv = np.unique(cid, return_inverse=True)
    w = 1 / z['rm_err'] ** 2; sw = np.bincount(cinv, weights=w)
    crm = np.bincount(cinv, weights=w * z['rm']) / sw; ce = 1 / np.sqrt(sw)
    first = np.zeros(len(uc), int); first[cinv[::-1]] = np.arange(len(cid))[::-1]
    cra, cdec, cdc, cmaj, cbeam, cfp, cfl, csid = (z[k][first] for k in ('ra', 'dec', 'dc_maj', 'maj', 'beam_maj', 'fracpol', 'flux', 'source_id'))
    log('components %d from %d rows; beam major axis median %.1f arcsec; dc_maj finite for %.3f, zero for %.3f' % (
        len(uc), len(cid), np.nanmedian(cbeam) * 3600, np.isfinite(cdc).mean(), (cdc == 0).mean()))
    # ---- sources as in Entry 176, with x
    rows, src, _ = c.load_spice()
    S, _ = c.sample(src, c.CFG['fg_radius'], c.CFG['min_neigh'])
    cl, _ = c.load_clusters()
    tree = cKDTree(c.gal_unit(S['l'], S['b']))
    n = len(S['d'])
    x, which = c.scaled_distance(tree, n, cl)
    zs = np.load(os.path.join(c.ROOT, 'entries', '171-spice-faraday-redshift', 'work', 'spice_rm.npz'))
    us = np.unique(zs['source_id'])                         # the order merge_sources uses
    idx_hi = np.nonzero(np.abs(src['b']) >= c.CFG['bcut'])[0]
    assert len(idx_hi) == n and len(us) == len(src['rm'])
    pos_in_S = np.full(len(us), -1); pos_in_S[idx_hi] = np.arange(n)
    comp_src = pos_in_S[np.searchsorted(us, csid)]           # component -> index in S (or -1)
    out = {}
    kpc_per_arcsec = c.ang_diam_Mpc(cl['z']) * 1000 * np.pi / 180 / 3600

    # ---- 1. sizes: the largest deconvolved major axis among a source's components
    okc = comp_src >= 0
    size = np.zeros(n); np.maximum.at(size, comp_src[okc], np.where(np.isfinite(cdc[okc]), cdc[okc], 0.0))
    ncomp = np.bincount(comp_src[okc], minlength=n)
    ring = (x >= 0.5) & (x < 1.0); field = ~np.isfinite(x)
    q = lambda v: [float(t) for t in np.percentile(v, [25, 50, 75, 90])]
    kpc = size[ring] * kpc_per_arcsec[which[ring]]
    co = json.load(open(os.path.join(HERE, 'coherence.json')))
    fac = (co['conservative']['lambda_over_diameter_at_least'], co['central']['lambda_over_diameter_at_least'])
    out['sizes'] = dict(ring_n=int(ring.sum()), ring_arcsec_quartiles=q(size[ring]), field_arcsec_quartiles=q(size[field]),
                        ring_share_unresolved=float((size[ring] == 0).mean()), field_share_unresolved=float((size[field] == 0).mean()),
                        ring_kpc_quartiles=q(kpc), lambda_at_least_kpc_median_source=[float(fac[0] * np.median(kpc)), float(fac[1] * np.median(kpc))],
                        lambda_at_least_kpc_upper_quartile_source=[float(fac[0] * np.percentile(kpc, 75)), float(fac[1] * np.percentile(kpc, 75))])
    log('sizes %s' % json.dumps(out['sizes']))
    # fractional polarization, flux-matched, by size half, inside R500 and in the ring, against rotated clusters
    zc = np.load(os.path.join(c.ROOT, 'entries', '171-spice-faraday-redshift', 'work', 'extra.npz'))
    u0, f0 = np.unique(zs['source_id'], return_index=True)
    fp = zc['fracpol'][f0].astype(float)[idx_hi]; fl = zc['total_I_flux'][f0].astype(float)[idx_hi]
    okf = np.isfinite(fp) & (fp > 0)
    dec_edges = np.quantile(fl[okf], np.linspace(0, 1, 11)); dec_edges[0], dec_edges[-1] = -np.inf, np.inf
    wh = np.clip(np.searchsorted(dec_edges, fl, side='right') - 1, 0, 9)
    med_size = float(np.median(size[okf]))
    small = size <= med_size
    out['depolarization_by_size'] = dict(size_split_arcsec=med_size, cases=[])
    g = np.random.default_rng(c.CFG['seed'] + 555)
    shifts = g.uniform(10.0, 350.0, 200)
    xrot = [c.scaled_distance(tree, n, cl, shift=sft)[0] for sft in shifts]
    for sname, ssel in (('smaller half', small), ('larger half', ~small)):
        expect = np.array([np.median(fp[okf & field & ssel & (wh == k)]) for k in range(10)])
        fpn = fp / expect[wh]
        for rname, lo, hi in (('0.5-1 R500', 0.5, 1.0), ('inside R500', 0.0, 1.0)):
            m = okf & ssel & (x >= lo) & (x < hi)
            real = float(np.median(fpn[m]))
            nul = np.array([np.median(fpn[okf & ssel & (xr >= lo) & (xr < hi)]) for xr in xrot])
            d = dict(size=sname, range=rname, n=int(m.sum()), ratio=real, null_mean=float(nul.mean()), null_sd=float(nul.std(ddof=1)),
                     t=float((real - nul.mean()) / nul.std(ddof=1)), median_size_arcsec=float(np.median(size[m])))
            out['depolarization_by_size']['cases'].append(d)
            log('fracpol %s, %s: n %d, ratio %.3f against %.3f +- %.3f, t %+.2f (median size %.1f arcsec)' % (
                sname, rname, d['n'], real, d['null_mean'], d['null_sd'], d['t'], d['median_size_arcsec']))

    # ---- 2. pairs of components at least 30 arcsec apart in the same source
    order = np.lexsort((ce, comp_src)); order = order[comp_src[order] >= 0]
    cs = comp_src[order]
    starts = np.r_[0, np.nonzero(cs[1:] != cs[:-1])[0] + 1]; ends = np.r_[starts[1:], len(cs)]
    P = dict(src=[], sep=[], d=[], e2=[], dev_a=[], dev_b=[], ea2=[], eb2=[])
    fgmed = src['rm'][idx_hi] - S['d']                      # the 3 degree foreground median of each source
    for a0, b0 in zip(starts, ends):
        if b0 - a0 < 2:
            continue
        ii = order[a0:b0]                                   # sorted by RM error: take the best component and its best partner beyond 30 arcsec
        i = ii[0]
        sep = np.hypot((cra[ii] - cra[i]) * np.cos(np.radians(cdec[i])), cdec[ii] - cdec[i]) * 3600
        far = np.nonzero(sep >= 30.0)[0]
        if len(far) == 0:
            continue
        j = ii[far[0]]
        s_ = comp_src[i]
        P['src'].append(s_); P['sep'].append(sep[far[0]]); P['d'].append(crm[i] - crm[j]); P['e2'].append(ce[i] ** 2 + ce[j] ** 2)
        P['dev_a'].append(crm[i] - fgmed[s_]); P['dev_b'].append(crm[j] - fgmed[s_]); P['ea2'].append(ce[i] ** 2); P['eb2'].append(ce[j] ** 2)
    P = {k: np.array(v) for k, v in P.items()}
    px = x[P['src']]
    log('pairs: %d sources with two components >= 30 arcsec apart; separation quartiles %s arcsec' % (len(px), q(P['sep'])))
    def Vp(m):
        return (1.4826 * np.median(np.abs(P['d'][m] - np.median(P['d'][m])))) ** 2 - float(np.mean(P['e2'][m]))
    def Vs(m):                                              # single sight lines: both components of the same pairs
        dd = np.r_[P['dev_a'][m], P['dev_b'][m]]; ee = np.r_[P['ea2'][m], P['eb2'][m]]
        return (1.4826 * np.median(np.abs(dd - np.median(dd)))) ** 2 - float(np.mean(ee))
    out['pairs'] = dict(n_pairs=int(len(px)), cases=[])
    pxr = [xr[P['src']] for xr in xrot]
    for rname, lo, hi in (('inside R500', 0.0, 1.0), ('0.5-1 R500', 0.5, 1.0), ('1-2 R500', 1.0, 2.0)):
        m = (px >= lo) & (px < hi)
        if m.sum() < 15:
            continue
        vp, vs = Vp(m), Vs(m)
        nulp = np.array([Vp((r >= lo) & (r < hi)) for r in pxr if ((r >= lo) & (r < hi)).sum() >= 15])
        nuls = np.array([Vs((r >= lo) & (r < hi)) for r in pxr if ((r >= lo) & (r < hi)).sum() >= 15])
        exp_, exs = vp - nulp.mean(), vs - nuls.mean()
        sepk = P['sep'][m] * kpc_per_arcsec[which[P['src'][m]]]
        d = dict(range=rname, n=int(m.sum()), separation_arcsec_median=float(np.median(P['sep'][m])), separation_kpc_quartiles=q(sepk)[:3],
                 V_pair=float(vp), pair_null_mean=float(nulp.mean()), pair_null_sd=float(nulp.std(ddof=1)), pair_excess=float(exp_),
                 t_pair=float(exp_ / nulp.std(ddof=1)), V_single=float(vs), single_null_mean=float(nuls.mean()), single_null_sd=float(nuls.std(ddof=1)),
                 single_excess=float(exs), t_single=float(exs / nuls.std(ddof=1)),
                 pair_excess_over_twice_single=float(exp_ / (2 * exs)) if exs > 0 else None)
        out['pairs']['cases'].append(d)
        log('pairs %s: n %d, separation %.0f arcsec (%.0f kpc median); pair variance %.1f against %.1f +- %.1f (excess %+.1f, t %+.2f); '
            'single sight lines %.1f against %.1f +- %.1f (excess %+.1f, t %+.2f); pair excess / (2 x single excess) = %s' % (
                rname, d['n'], d['separation_arcsec_median'], np.median(sepk), vp, nulp.mean(), nulp.std(ddof=1), exp_, d['t_pair'],
                vs, nuls.mean(), nuls.std(ddof=1), exs, d['t_single'], ('%.2f' % d['pair_excess_over_twice_single']) if d['pair_excess_over_twice_single'] is not None else 'undefined'))
    # ---- pairs: are the doubles cluster members?  errors by bootstrap; and the excess against separation
    foot = float(np.mean([(r >= 4.0).sum() + (~np.isfinite(r)).sum() for r in pxr]))
    foot = ((px >= 4.0).sum() + (~np.isfinite(px)).sum()) / foot          # real over rotated counts far from clusters
    gb = np.random.default_rng(c.CFG['seed'] + 556)
    extra = []
    for rname, lo, hi in (('inside R500', 0.0, 1.0), ('0.5-1 R500', 0.5, 1.0), ('1-2 R500', 1.0, 2.0)):
        m = np.nonzero((px >= lo) & (px < hi))[0]
        expect = foot * float(np.mean([((r >= lo) & (r < hi)).sum() for r in pxr]))
        nulp = np.array([Vp((r >= lo) & (r < hi)) for r in pxr if ((r >= lo) & (r < hi)).sum() >= 15])
        nuls = np.array([Vs((r >= lo) & (r < hi)) for r in pxr if ((r >= lo) & (r < hi)).sum() >= 15])
        bp, bs = [], []
        for _ in range(2000):
            k = gb.choice(m, len(m))
            mm = np.zeros(len(px), bool)                      # Vp, Vs take boolean masks; build the resample by index instead
            dd = P['d'][k]; vp_ = (1.4826 * np.median(np.abs(dd - np.median(dd)))) ** 2 - float(np.mean(P['e2'][k]))
            d2 = np.r_[P['dev_a'][k], P['dev_b'][k]]; vs_ = (1.4826 * np.median(np.abs(d2 - np.median(d2)))) ** 2 - float(np.mean(np.r_[P['ea2'][k], P['eb2'][k]]))
            bp.append(vp_ - nulp.mean()); bs.append(vs_ - nuls.mean())
        bp, bs = np.array(bp), np.array(bs)
        ratio = bp / (2 * np.where(bs > 0, bs, np.nan))
        sepk = P['sep'][m] * kpc_per_arcsec[which[P['src'][m]]]
        d = dict(range=rname, n=int(len(m)), expected_without_members=float(expect), pair_excess_bootstrap_sd=float(bp.std(ddof=1)),
                 pair_excess_q16_q84=[float(v) for v in np.percentile(bp, [16, 84])], single_excess_bootstrap_sd=float(bs.std(ddof=1)),
                 ratio_median=float(np.nanmedian(ratio)), ratio_q05=float(np.nanpercentile(ratio, 5)), ratio_q16_q84=[float(v) for v in np.nanpercentile(ratio, [16, 84])],
                 implied_single_variance=float(np.median(bp) / 2), implied_single_variance_sd=float(bp.std(ddof=1) / 2))
        if rname == 'inside R500':
            for label, sel2 in (('closer than the median separation', sepk < np.median(sepk)), ('wider than the median separation', sepk >= np.median(sepk))):
                k = m[sel2]; dd = P['d'][k]
                v = (1.4826 * np.median(np.abs(dd - np.median(dd)))) ** 2 - float(np.mean(P['e2'][k]))
                bb = []
                for _ in range(1000):
                    kk = gb.choice(k, len(k)); d3 = P['d'][kk]
                    bb.append((1.4826 * np.median(np.abs(d3 - np.median(d3)))) ** 2 - float(np.mean(P['e2'][kk])))
                d[label] = dict(n=int(len(k)), separation_kpc_median=float(np.median(sepk[sel2])), pair_variance=float(v), bootstrap_sd=float(np.std(bb, ddof=1)))
        extra.append(d)
        log('pairs %s: %d sources, %.0f expected without members; pair excess %+.0f +- %.0f (bootstrap); pair excess / (2 x single) median %.2f, 5th percentile %.2f; '
            'single-sight-line variance implied %.0f +- %.0f' % (rname, len(m), expect, np.median(bp), bp.std(ddof=1), d['ratio_median'], d['ratio_q05'],
                                                           d['implied_single_variance'], d['implied_single_variance_sd']))
        if rname == 'inside R500':
            log('   by separation: %s' % json.dumps({k2: v2 for k2, v2 in d.items() if isinstance(v2, dict)}))
    out['pairs']['bootstrap'] = extra
    out['pairs']['field_pair_variance'] = float(np.mean([Vp(~np.isfinite(r)) for r in pxr[:20]]))
    # upper bound on Lambda: correlation rho = 1 - ratio at the median separation; Gaussian correlation exp(-s^2 / 2 Lambda^2)
    ins = extra[0]; sep_med = [d for d in out['pairs']['cases'] if d['range'] == 'inside R500'][0]['separation_kpc_quartiles'][1]
    rho_max = max(1e-6, 1 - ins['ratio_q05'])
    out['pairs']['lambda_upper_kpc_95'] = float(sep_med / np.sqrt(2 * np.log(1 / rho_max))) if rho_max < 1 else None
    out['pairs']['rho_upper_95'] = float(rho_max); out['pairs']['separation_kpc_median_inside_R500'] = float(sep_med)
    log('correlation between sight lines %.0f kpc apart: below %.2f at 95%%; reversal length below %s kpc' % (
        sep_med, rho_max, ('%.0f' % out['pairs']['lambda_upper_kpc_95']) if out['pairs']['lambda_upper_kpc_95'] else 'unbounded'))
    json.dump(out, open(os.path.join(HERE, 'coherence2.json'), 'w'), indent=1, default=float)


if __name__ == '__main__':
    {'extract': extract, 'run': run}.get(sys.argv[1] if len(sys.argv) > 1 else '', lambda: print(__doc__))()
