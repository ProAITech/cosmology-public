"""Entry 171, exploratory checks made AFTER the frozen result was seen (not part of the rule; reported as exploratory).

The frozen test found the RM scatter growing with quasar redshift (t = +3.81).  These checks ask whether something other than the
sightline makes it: how many components each source merges, the RM error, the radio flux and polarized fraction, Faraday
complexity, and where on the sky the quasars of each redshift lie.  Same data and statistic as spice_faraday_z.py.

  .venv/bin/python -I entries/171-spice-faraday-redshift/explore.py extra     # re-stream the table for extra columns (~30 s)
  .venv/bin/python -I entries/171-spice-faraday-redshift/explore.py run       # writes explore.json
"""
import gzip, json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import spice_faraday_z as S  # noqa: E402

EXTRA = ('total_I_flux', 'polint', 'fracpol', 'snr_polint', 'Ncomp', 'sigma_add', 'rm_width', 'separation_tile_centre',
         'spectral_index', 's_code', 'is_blended_flag', 'local_rm_flag', 'goodRM_flag', 'complex_M2_CC_flag', 'complex_sigma_add_flag')


def extra():
    keep = {k: [] for k in EXTRA}
    with gzip.open(S.SPICE, 'rb') as f:
        S.read_header(f)
        h = S.read_header(f)
        dt, rowbytes = S.table_dtype(h)
        nrows, done = int(S.card(h, 'NAXIS2')), 0
        while done < nrows:
            n = min(100_000, nrows - done)
            a = np.frombuffer(f.read(n * rowbytes), dtype=dt)
            ok = np.ones(n, bool)
            for k in S.FLAGS:
                ok &= a[k] == b'F'
            rm, e = a['rm'].astype('<f8'), a['rm_err'].astype('<f8')
            ok &= np.isfinite(rm) & np.isfinite(e) & (e > 0)
            for k in EXTRA:
                v = a[k][ok]
                if v.dtype.kind == 'S':
                    keep[k].append(v.copy())
                else:
                    keep[k].append(v.astype('<f8'))
            done += n
    np.savez(os.path.join(S.WORK, 'extra.npz'), **{k: np.concatenate(v) for k, v in keep.items()})
    print('extra written')


def run():
    t0 = time.time()
    z = np.load(os.path.join(S.WORK, 'spice_rm.npz'))
    X = np.load(os.path.join(S.WORK, 'extra.npz'))
    rows = {k: z[k] for k in z.files}
    src, inv = S.merge_sources(rows)
    nsrc = len(src['rm'])
    # per-source extras: from the source's first row, plus counts
    first = np.zeros(nsrc, dtype=int); first[inv[::-1]] = np.arange(len(inv))[::-1]
    ex = {}
    for k in ('total_I_flux', 'polint', 'fracpol', 'snr_polint', 'Ncomp', 'sigma_add', 'rm_width', 'separation_tile_centre', 'spectral_index'):
        ex[k] = X[k][first]
    ex['s_code'] = X['s_code'][first]
    ex['goodRM_any'] = np.bincount(inv, weights=(X['goodRM_flag'] == b'T').astype(float)) > 0
    rows_ex = dict(row_ra=rows['ra'], row_dec=rows['dec'], row_src=inv)
    mq, _ = S.load_milliquas()
    m = S.match(mq, rows_ex, nsrc, ['Q'])
    SQ, cQ = S.build(m, src, mq, S.CFG['fg_radius'], S.CFG['min_neigh'])
    # recover the source index of each kept quasar (build keeps order of m['src'] filtered)
    lat = np.abs(src['b'][m['src']]) >= S.CFG['bcut']
    med, nn = S.foreground(src, m['src'][lat], S.CFG['fg_radius'], S.CFG['min_neigh'])
    kidx = m['src'][lat][np.isfinite(med)]
    assert len(kidx) == len(SQ['z'])
    for k, v in ex.items():
        SQ[k] = v[kidx]
    SQ['nrow'] = src['nrow'][kidx]
    SQ['ra'], SQ['dec'] = src['ra'][kidx], src['dec'][kidx]
    out = {'n': int(len(SQ['z']))}
    # 1. how properties change across the 8 redshift bins
    idx = np.argsort(SQ['z'], kind='stable'); n = len(idx)
    prof = {k: [] for k in ('z_med', 'e_med', 'nrow_multi_frac', 'I_med_mJy', 'fracpol_med', 'snr_med', 'Ncomp_gt1_frac', 'sigma_add_med',
                            'sep_tile_med', 'alpha_med', 'abs_b_med', 'e2_mean', 'sigma_minus_noise')}
    for k in range(8):
        q = idx[S.jsround(k * n / 8):S.jsround((k + 1) * n / 8)]
        prof['z_med'].append(float(np.median(SQ['z'][q])))
        prof['e_med'].append(float(np.median(SQ['e'][q])))
        prof['nrow_multi_frac'].append(float((SQ['nrow'][q] > 1).mean()))
        prof['I_med_mJy'].append(float(np.nanmedian(SQ['total_I_flux'][q]) * 1e3))
        prof['fracpol_med'].append(float(np.nanmedian(SQ['fracpol'][q])))
        prof['snr_med'].append(float(np.nanmedian(SQ['snr_polint'][q])))
        prof['Ncomp_gt1_frac'].append(float(np.nanmean(SQ['Ncomp'][q] > 1)))
        prof['sigma_add_med'].append(float(np.nanmedian(SQ['sigma_add'][q])))
        prof['sep_tile_med'].append(float(np.nanmedian(SQ['separation_tile_centre'][q])))
        prof['alpha_med'].append(float(np.nanmedian(SQ['spectral_index'][q])))
        prof['abs_b_med'].append(float(np.median(np.abs(SQ['b'][q]))))
        s = float(S.sig_mad(SQ['d'][q]))
        e2 = float(np.mean(SQ['e'][q] ** 2))
        prof['e2_mean'].append(e2)
        prof['sigma_minus_noise'].append(float(np.sqrt(max(0.0, s * s - e2))))
    out['profile'] = {k: [round(v, 3) for v in vals] for k, vals in prof.items()}
    print(json.dumps(out['profile']), flush=True)
    sub = S.sub

    def an(SS, label, cells=None, nsh=300):
        if cells is None:
            return S.analyse(SS, nsh, label)
        # custom cells: shuffle within the given labels
        Sx = dict(SS)
        return analyse_cells(Sx, cells, nsh, label)

    res = {}
    res['single_row'] = an(sub(SQ, SQ['nrow'] == 1), 'single-row sources')
    res['multi_row'] = an(sub(SQ, SQ['nrow'] > 1), 'multi-row sources')
    res['bright_I'] = an(sub(SQ, SQ['total_I_flux'] >= np.nanmedian(SQ['total_I_flux'])), 'Stokes I above median')
    res['faint_I'] = an(sub(SQ, SQ['total_I_flux'] < np.nanmedian(SQ['total_I_flux'])), 'Stokes I below median')
    res['high_fracpol'] = an(sub(SQ, SQ['fracpol'] >= np.nanmedian(SQ['fracpol'])), 'fracpol above median')
    res['low_fracpol'] = an(sub(SQ, SQ['fracpol'] < np.nanmedian(SQ['fracpol'])), 'fracpol below median')
    res['abs_b_ge_40'] = an(sub(SQ, np.abs(SQ['b']) >= 40), '|b| >= 40')
    res['goodRM'] = an(sub(SQ, SQ['goodRM_any']), 'goodRM_flag set')
    # finer sky cells: HEALPix nside 4 (about 15 deg) x e tercile
    import healpy as hp
    pix = hp.ang2pix(4, SQ['ra'], SQ['dec'], lonlat=True)
    es = np.sort(SQ['e']); nq = len(es); t1, t2 = es[nq // 3], es[2 * nq // 3]
    ter = np.where(SQ['e'] < t1, 0, np.where(SQ['e'] < t2, 1, 2))
    res['sky_cells_nside4'] = analyse_cells(SQ, pix * 3 + ter, 300, 'HEALPix nside 4 x e tercile')
    # matched cells on flux and fracpol and e (terciles each) x |b| band
    def terc(v):
        qs = np.nanquantile(v, [1 / 3, 2 / 3]); return np.digitize(v, qs)
    band = np.digitize(np.abs(SQ['b']), [35, 50, 65])
    res['cells_flux_fracpol_e_band'] = analyse_cells(SQ, ((terc(SQ['total_I_flux']) * 3 + terc(SQ['fracpol'])) * 3 + ter) * 4 + band, 300,
                                                     'I x fracpol x e terciles x |b| band')
    out['checks'] = {k: dict(n=v['n'], b=round(v['b'], 3), null_mean=round(v['null']['mean'], 3), null_sd=round(v['null']['sd'], 3),
                             t=round(v['t'], 2), sigma=[round(s, 2) for s in v['bins']['sigma']], x=[round(s, 3) for s in v['bins']['x']])
                     for k, v in res.items()}
    out['seconds'] = round(time.time() - t0)
    json.dump(out, open(os.path.join(HERE, 'explore.json'), 'w'), indent=1)
    print('written explore.json %d s' % out['seconds'], flush=True)


def analyse_cells(SS, cells, nshuf, label):
    seed = S.CFG['seed']
    z, d = SS['z'], SS['d']
    real = S.bins_fit(z, d, np.random.default_rng(seed))
    groups = [np.nonzero(cells == c)[0] for c in np.unique(cells)]
    R = np.random.default_rng(seed + 7)
    nulls = np.empty(nshuf)
    for s in range(nshuf):
        z2 = z.copy()
        for g in groups:
            z2[g] = z[g][R.permutation(len(g))]
        nulls[s] = S.bins_fit(z2, d, np.random.default_rng(seed + 1000 + s))['b']
    mu, sd = float(nulls.mean()), float(nulls.std(ddof=1))
    print('%s: n %d, cells %d, b %+.3f, null %+.3f +- %.3f, t %+.2f' % (label, len(z), len(groups), real['b'], mu, sd, (real['b'] - mu) / sd), flush=True)
    return dict(n=int(len(z)), b=real['b'], bins=dict(sigma=real['sigma'], x=real['x']), null=dict(mean=mu, sd=sd), t=(real['b'] - mu) / sd)


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else ''
    {'extra': extra, 'run': run}.get(cmd, lambda: print(__doc__))()
