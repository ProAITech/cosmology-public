"""Entry 177: does the rotation-measure difference between two components of one source grow with redshift?

Entries 171-174 found the scatter of rotation measure (RM) growing with quasar redshift at about 3.5 sigma on 10,449 single sight
lines, cause unknown.  A source with two components far enough apart gives two sight lines through the same Galactic foreground
and the same large-scale structure, a few hundred kpc apart at the distance of anything in between.  Rotation from screens smaller
than that (the gas around intervening galaxies, cluster cores) differs between the two and adds twice its variance to the
difference; rotation from anything coherent across both (the Galaxy, large-scale intergalactic structure) cancels.  So if the growth
of Entries 171-174 comes from small intervening screens, the variance of the difference must grow with redshift twice as fast.
The rule was frozen in the Workbench (section 177, "the rule, frozen first") before any difference was set against redshift.

Data (raw files are git-ignored under data/; URLs, sizes and checksums as in Entries 171, 173 and 174):
  data/spice-racs-dr2/spice-racs.dr2.fits.gz            SPICE-RACS DR2; columns source_id, cat_id, ra, dec, rm, rm_err, the four flags
  data/milliquas-v8/milliquas.zip                       Milliquas v8, classes Q, A, B, K, N with z > 0
  data/desi-dr1-qso/QSO_cat_iron_cumulative_v0.fits     DESI DR1 quasars
  data/quaia/quaia_G20.5.fits                           Quaia

Rule, in brief (section 177 has the full text):
  Components one RM per cat_id (inverse-variance mean of repeated rows), Entry 171's quality cuts.
  Pairs      per source: the component with the smallest RM error and the smallest-error component at least 30 arcsec from it.
             dRM = RM_a - RM_b, error^2 = e_a^2 + e_b^2.  Galactic |b| >= 20 deg.
  Redshift   an object within 10 arcsec of any component of the source: Milliquas first, then DESI, then Quaia.
  Statistic  K bins of equal count in z (6 if 600 pairs or more, else 4); V = (1.4826 MAD)^2 - mean error^2; weighted least-squares
             slope b_pair of V against mean log10(1+z), bootstrap errors; 1,000 shuffles of z within 9 cells (terciles of error
             by terciles of separation); t = (b_pair - mean) / sd.
  Reference  b_single, the same slope for single sight lines, from the eight bins of the union of Entries 171-174 (constants below).
  Decision   decide_177() below; the frozen text is section 177, item 6.

Run from the repository root:
  .venv/bin/python -I entries/177-pair-differences-redshift/pair_z.py synthetic   # synthetic pairs; reads no data
  .venv/bin/python -I entries/177-pair-differences-redshift/pair_z.py extract     # one pass over the catalogue, about 1 minute
  .venv/bin/python -I entries/177-pair-differences-redshift/pair_z.py run         # the test; writes results.json
"""
import gzip, json, os, sys, time, zipfile
import numpy as np
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
SPICE = os.path.join(ROOT, 'data', 'spice-racs-dr2', 'spice-racs.dr2.fits.gz')
QSO = os.path.join(ROOT, 'data', 'desi-dr1-qso', 'QSO_cat_iron_cumulative_v0.fits')
QUAIA = os.path.join(ROOT, 'data', 'quaia', 'quaia_G20.5.fits')
MQZIP = os.path.join(ROOT, 'data', 'milliquas-v8', 'milliquas.zip')
WORK = os.path.join(HERE, 'work')
D2R = np.pi / 180

# the union of Entries 171, 173 and 174 (entries/174-quaia-replication/explore_union.json): median z, noise-removed variance, error
UNION = dict(zmed=[0.518, 0.801, 1.017, 1.230, 1.460, 1.715, 2.062, 2.653],
             V=[64.06, 75.78, 74.02, 77.03, 78.09, 82.09, 88.03, 89.52], err=[5.27, 6.14, 6.24, 5.64, 6.01, 6.09, 6.91, 6.41])
CFG = dict(match_arcsec=10.0, min_sep_arcsec=30.0, bcut=20.0, nboot=500, nshuf=1000, seed=20261121, classes=['Q', 'A', 'B', 'K', 'N'],
           bins_many=6, bins_few=4, many=600, min_pairs=300, beta_quaia=0.835, share=0.5)

MG = np.array([[-0.0548755604, -0.8734370902, -0.4838350155],
               [0.4941094279, -0.4448296300, 0.7469822445],
               [-0.8676661490, -0.1980763734, 0.4559837762]])


def unit(ra, dec):
    a, d = np.asarray(ra) * D2R, np.asarray(dec) * D2R
    return np.c_[np.cos(d) * np.cos(a), np.cos(d) * np.sin(a), np.sin(d)]


def galactic(ra, dec):
    v = unit(ra, dec) @ MG.T
    return np.mod(np.arctan2(v[:, 1], v[:, 0]) / D2R, 360.0), np.arcsin(np.clip(v[:, 2], -1, 1)) / D2R


def chord(deg):
    return 2 * np.sin(0.5 * deg * D2R)


# ---------------------------------------------------------------- FITS tables, read by hand
def read_header(f):
    cards = []
    while True:
        blk = f.read(2880)
        if len(blk) < 2880:
            raise SystemExit('short read in header')
        for i in range(0, 2880, 80):
            c = blk[i:i + 80].decode('ascii', errors='replace')
            cards.append(c)
            if c.startswith('END     '):
                return cards


def card(cards, key):
    for c in cards:
        if c[:8].strip() == key:
            return c[10:].split('/')[0].strip().strip("'").strip()
    return None


def table_dtype(h):
    nf = int(card(h, 'TFIELDS'))
    fields, off = {}, 0
    for i in range(1, nf + 1):
        n, t = card(h, 'TTYPE%d' % i), card(h, 'TFORM%d' % i)
        if t.endswith('A'):
            code, size = 'S%d' % int(t[:-1] or 1), int(t[:-1] or 1)
        else:
            t = t.lstrip('1')
            code, size = {'D': ('>f8', 8), 'E': ('>f4', 4), 'K': ('>i8', 8), 'J': ('>i4', 4), 'I': ('>i2', 2), 'L': ('S1', 1)}[t]
        fields[n] = (code, off)
        off += size
    rowbytes = int(card(h, 'NAXIS1'))
    assert off == rowbytes, (off, rowbytes)
    return np.dtype({'names': list(fields), 'formats': [v[0] for v in fields.values()], 'offsets': [v[1] for v in fields.values()],
                     'itemsize': rowbytes}), rowbytes


def read_table(path, cols):
    """First binary table of an uncompressed FITS file; returns the named columns as native arrays, and the row count."""
    with open(path, 'rb') as f:
        read_header(f)
        h = read_header(f)
        dt, rowbytes = table_dtype(h)
        nrows = int(card(h, 'NAXIS2'))
        a = np.fromfile(f, dtype=dt, count=nrows)
    if len(a) != nrows:
        raise SystemExit('short table in %s' % path)
    out = {}
    for c in cols:
        v = a[c]
        out[c] = v.astype(v.dtype.newbyteorder('=')) if v.dtype.kind in 'fi' else v.copy()
    return out, nrows


FLAGS = ('snr_flag', 'leakage_flag', 'channel_flag', 'stokesI_fit_flag')


def load_desi():
    t, nrows = read_table(QSO, ['TARGETID', 'Z', 'TARGET_RA', 'TARGET_DEC', 'TSNR2_QSO'])
    tid, ts = t['TARGETID'], t['TSNR2_QSO'].astype(float)
    order = np.lexsort((-np.where(np.isfinite(ts), ts, -np.inf), tid))
    k = order[np.r_[True, tid[order][1:] != tid[order][:-1]]]
    z, ra, dec = t['Z'][k].astype(float), t['TARGET_RA'][k].astype(float), t['TARGET_DEC'][k].astype(float)
    ok = np.isfinite(z) & (z > 0) & np.isfinite(ra) & np.isfinite(dec)
    return dict(z=z[ok], ra=ra[ok], dec=dec[ok])


def load_milliquas():
    """Classes Q, A, B, K, N with z > 0 (Entry 171's 'all classes')."""
    with zipfile.ZipFile(MQZIP) as zf:
        raw = zf.read('milliquas.txt')
    ra, dec, z, cls = [], [], [], []
    for line in raw.split(b'\n'):
        line = line.rstrip(b'\r')
        if len(line) < 82:
            continue
        c = line[51:52].decode('ascii', 'replace')
        if c not in CFG['classes']:
            continue
        try:
            r, d, zz = float(line[0:11]), float(line[12:23]), float(line[76:82])
        except ValueError:
            continue
        if zz > 0:
            ra.append(r); dec.append(d); z.append(zz); cls.append(c)
    return dict(ra=np.array(ra), dec=np.array(dec), z=np.array(z), cls=np.array(cls))


def load_quaia():
    t, nrows = read_table(QUAIA, ['ra', 'dec', 'redshift_quaia', 'redshift_quaia_err', 'phot_g_mean_mag'])
    z = t['redshift_quaia'].astype(float)
    ok = np.isfinite(z) & (z > 0) & np.isfinite(t['ra']) & np.isfinite(t['dec'])
    q = dict(ra=t['ra'][ok].astype(float), dec=t['dec'][ok].astype(float), z=z[ok], zerr=t['redshift_quaia_err'][ok].astype(float),
             g=t['phot_g_mean_mag'][ok].astype(float))
    return q, {'rows': int(nrows), 'kept_z_positive': int(ok.sum())}


def jsround(x):
    return int(np.floor(x + 0.5))


def sig_mad(X):
    m = np.median(X, axis=-1, keepdims=True)
    return 1.4826 * np.median(np.abs(X - m), axis=-1)


def wls(x, y, e):
    x, y, w = np.asarray(x), np.asarray(y), 1 / np.asarray(e) ** 2
    S, Sx, Sy, Sxx, Sxy = w.sum(), (w * x).sum(), (w * y).sum(), (w * x * x).sum(), (w * x * y).sum()
    D = S * Sxx - Sx * Sx
    b = (S * Sxy - Sx * Sy) / D
    a = (Sy - b * Sx) / S
    return dict(a=a, b=b, eb=np.sqrt(S / D), chi2=float((((y - a - b * x) / np.asarray(e)) ** 2).sum()))


# ---------------------------------------------------------------- new in Entry 177
def extract():
    t0 = time.time()
    keep = {k: [] for k in ('source_id', 'cat_id', 'ra', 'dec', 'rm', 'rm_err')}
    with gzip.open(SPICE, 'rb') as f:
        read_header(f); h = read_header(f)
        dt, rowbytes = table_dtype(h); nrows = int(card(h, 'NAXIS2')); done = 0
        while done < nrows:
            n = min(100_000, nrows - done)
            a = np.frombuffer(f.read(n * rowbytes), dtype=dt)
            clear = np.ones(n, bool)
            for k in FLAGS:
                clear &= a[k] == b'F'
            rm, e = a['rm'].astype('<f8'), a['rm_err'].astype('<f8')
            ok = clear & np.isfinite(rm) & np.isfinite(e) & (e > 0)
            keep['source_id'].append(a['source_id'][ok].copy()); keep['cat_id'].append(a['cat_id'][ok].copy())
            keep['ra'].append(a['ra'][ok].astype('<f8')); keep['dec'].append(a['dec'][ok].astype('<f8'))
            keep['rm'].append(rm[ok]); keep['rm_err'].append(e[ok])
            done += n
    os.makedirs(WORK, exist_ok=True)
    np.savez(os.path.join(WORK, 'components.npz'), **{k: np.concatenate(v) for k, v in keep.items()})
    print('extract done, %.0f s' % (time.time() - t0), flush=True)


def components():
    p = os.path.join(WORK, 'components.npz')
    if not os.path.exists(p):                                # Entry 176 kept the same rows
        p = os.path.join(ROOT, 'entries', '176-cluster-faraday-profile', 'work', 'components.npz')
    z = np.load(p)
    uc, cinv = np.unique(z['cat_id'], return_inverse=True)
    w = 1 / z['rm_err'] ** 2; sw = np.bincount(cinv, weights=w)
    first = np.zeros(len(uc), int); first[cinv[::-1]] = np.arange(len(cinv))[::-1]
    us, sinv = np.unique(z['source_id'][first], return_inverse=True)
    return dict(ra=z['ra'][first], dec=z['dec'][first], rm=np.bincount(cinv, weights=w * z['rm']) / sw, e=1 / np.sqrt(sw), src=sinv, nsrc=len(us))


def source_redshifts(C, cats):
    """cats: list of (name, dict(ra, dec, z)) in order of priority.  Returns z and the catalogue index per source (-1: none)."""
    tree = cKDTree(unit(C['ra'], C['dec']))
    zs, which, sepbest = np.full(C['nsrc'], np.nan), np.full(C['nsrc'], -1), np.full(C['nsrc'], np.inf)
    for k, (name, cat) in enumerate(cats):
        d, j = tree.query(unit(cat['ra'], cat['dec']), k=1, distance_upper_bound=chord(CFG['match_arcsec'] / 3600.0))
        hit = np.nonzero(np.isfinite(d))[0]
        s = C['src'][j[hit]]
        for i, si, di in zip(hit, s, d[hit]):
            if which[si] == -1 or (which[si] == k and di < sepbest[si]):
                zs[si], which[si], sepbest[si] = cat['z'][i], k, di
    return zs, which


def build_pairs(C):
    order = np.lexsort((C['e'], C['src']))
    cs = C['src'][order]
    starts = np.r_[0, np.nonzero(cs[1:] != cs[:-1])[0] + 1]; ends = np.r_[starts[1:], len(cs)]
    P = dict(src=[], sep=[], d=[], e2=[], ra=[], dec=[])
    for a0, b0 in zip(starts, ends):
        if b0 - a0 < 2:
            continue
        ii = order[a0:b0]; i = ii[0]
        sep = np.hypot((C['ra'][ii] - C['ra'][i]) * np.cos(C['dec'][i] * D2R), C['dec'][ii] - C['dec'][i]) * 3600
        far = np.nonzero(sep >= CFG['min_sep_arcsec'])[0]
        if len(far) == 0:
            continue
        j = ii[far[0]]
        P['src'].append(C['src'][i]); P['sep'].append(sep[far[0]]); P['d'].append(C['rm'][i] - C['rm'][j]); P['e2'].append(C['e'][i] ** 2 + C['e'][j] ** 2)
        P['ra'].append(C['ra'][i]); P['dec'].append(C['dec'][i])
    return {k: np.array(v) for k, v in P.items()}


def V_of(d, e2):
    return float((1.4826 * np.median(np.abs(d - np.median(d)))) ** 2 - np.mean(e2))


def v_slope(z, d, e2, K, rng):
    n = len(z); idx = np.argsort(z, kind='stable')
    x, V, er, cnt, zm = [], [], [], [], []
    for k in range(K):
        q = idx[jsround(k * n / K):jsround((k + 1) * n / K)]; m = len(q)
        bs = [V_of(d[q][r], e2[q][r]) for r in rng.integers(0, m, size=(CFG['nboot'], m))]
        x.append(float(np.log10(1 + z[q]).mean())); V.append(V_of(d[q], e2[q])); er.append(float(np.std(bs, ddof=1))); cnt.append(m); zm.append(float(np.median(z[q])))
    f = wls(x, V, er)
    return dict(b=float(f['b']), a=float(f['a']), eb=float(f['eb']), chi2=f['chi2'], x=x, V=V, err=er, n=cnt, zmed=zm)


def analyse(z, d, e2, sep, K, nshuf, seed, log=print, label='pairs'):
    real = v_slope(z, d, e2, K, np.random.default_rng(seed))
    te = np.searchsorted(np.quantile(e2, [1 / 3, 2 / 3]), e2); ts = np.searchsorted(np.quantile(sep, [1 / 3, 2 / 3]), sep)
    cell = te * 3 + ts
    groups = [np.nonzero(cell == c_)[0] for c_ in np.unique(cell)]
    R = np.random.default_rng(seed + 7)
    nulls = np.empty(nshuf)
    for s in range(nshuf):
        z2 = z.copy()
        for g in groups:
            z2[g] = z[g][R.permutation(len(g))]
        nulls[s] = v_slope(z2, d, e2, K, np.random.default_rng(seed + 1000 + s))['b']
    mu, sd = float(nulls.mean()), float(nulls.std(ddof=1))
    out = dict(label=label, n=int(len(z)), bins=K, b=real['b'], eb_fit=real['eb'], chi2=real['chi2'], x=real['x'], V=real['V'], err=real['err'],
               n_bin=real['n'], zmed=real['zmed'], null_mean=mu, null_sd=sd, t=float((real['b'] - mu) / sd), net=float(real['b'] - mu),
               p_upper=float((1 + (nulls >= real['b']).sum()) / (1 + nshuf)))
    log('%s: n %d, slope of V %+.1f, null %+.1f +- %.1f, t %+.2f; V by bin %s' % (label, out['n'], out['b'], mu, sd, out['t'], [round(v, 1) for v in real['V']]))
    return out


def single_reference():
    """Slope of the single-sight-line variance against log10(1+z) in the union of Entries 171-174, with its fit error."""
    x = np.log10(1 + np.array(UNION['zmed']))
    f = wls(x, UNION['V'], UNION['err'])
    return float(f['b']), float(f['eb'])


def decide_177(P, share_quaia):
    """Section 177, item 6.  P: output of analyse() for the pairs."""
    bs, es = single_reference()
    beta = 1.0 - (1.0 - CFG['beta_quaia']) * share_quaia      # Quaia's estimated redshifts flatten a slope
    full = 2.0 * bs * beta                                    # what the pairs would show if the single-sight-line growth were uncorrelated between them
    # rotation made at the source differs between the components and is diluted by (1+z)^-2, which lowers V with redshift:
    # at most the whole of the lowest bin's V, falling as (1+z)^-4 from the first bin to the last
    mask = P['V'][0] * (1 - ((1 + P['zmed'][0]) / (1 + P['zmed'][-1])) ** 4) / (P['x'][-1] - P['x'][0]) if P['V'][0] > 0 else 0.0
    u = {f: float((P['net'] + mask - f * full) / np.sqrt(P['null_sd'] ** 2 + (f * 2.0 * es * beta) ** 2)) for f in (1.0, CFG['share'])}
    out = dict(t=P['t'], slope_pairs=P['net'], slope_pairs_sd=P['null_sd'], slope_single=bs, slope_single_err=es, beta=beta, slope_if_uncorrelated=full,
               masking_allowance=float(mask), ratio=float(P['net'] / full), ratio_upper95=float((P['net'] + mask + 1.645 * P['null_sd']) / full),
               u_all=u[1.0], u_half=u[CFG['share']])
    if P['t'] >= 3:
        o = 'the difference between components grows with redshift (t >= +3): screens smaller than the separation carry growth'
    elif u[CFG['share']] <= -1.645:
        o = ('screens smaller than the separation carry less than half of the single-sight-line growth (95%): '
             'if that growth is real, it comes from something coherent across the two components')
    elif u[1.0] <= -1.645:
        o = 'screens smaller than the separation do not carry all of the single-sight-line growth (95%); whether they carry half is not decided'
    else:
        o = 'unresolved: no growth in the difference is established, and the whole of the single-sight-line growth is not excluded'
    out['outcome'] = o
    return out


def synthetic():
    out = {}
    bs, es = single_reference()
    out['single_reference'] = dict(slope=bs, fit_error=es)
    for name, grow, intrinsic, seed in (('no growth', 0.0, 0.0, 41), ('no growth, second sky', 0.0, 0.0, 42), ('uncorrelated screens, full growth', 2 * bs, 0.0, 43),
                                        ('half the growth', bs, 0.0, 44), ('rotation at the source only, diluted', 0.0, 12.0, 45)):
        g = np.random.default_rng(seed)
        n = 850
        z = np.clip(g.gamma(3.0, 0.4, n), 0.05, 5)
        x = np.log10(1 + z)
        e2 = g.lognormal(np.log(4.0), 0.6, n); sep = g.uniform(30, 120, n)
        var = 12.0 + grow * (x - x.min()) + intrinsic * ((1 + z) ** -4 - 1)
        d = g.normal(0, 1, n) * np.sqrt(np.maximum(var, 0.5)) + g.normal(0, 1, n) * np.sqrt(e2)
        P = analyse(z, d, e2, sep, CFG['bins_many'], 200, seed, label=name)
        dec = decide_177(P, 0.3)
        out[name] = dict(t=round(P['t'], 2), slope=round(P['net'], 1), sd=round(P['null_sd'], 1), ratio=round(dec['ratio'], 2), u_all=round(dec['u_all'], 2), u_half=round(dec['u_half'], 2),
                         masking_allowance=round(dec['masking_allowance'], 1), outcome=dec['outcome'])
    print(json.dumps(out, indent=1))
    json.dump(out, open(os.path.join(HERE, 'synthetic.json'), 'w'), indent=1)


def run():
    t0 = time.time()
    log = lambda s: print(s, flush=True)
    C = components()
    mq = load_milliquas(); desi = load_desi(); quaia, _ = load_quaia()
    zs, which = source_redshifts(C, [('milliquas', mq), ('desi', desi), ('quaia', quaia)])
    P = build_pairs(C)
    l, b = galactic(P['ra'], P['dec'])
    zz = zs[P['src']]; wh = which[P['src']]
    keep = np.isfinite(zz) & (zz > 0) & (np.abs(b) >= CFG['bcut'])
    res = dict(generated=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), config=CFG, union_reference=UNION,
               counts=dict(components=int(len(C['rm'])), sources=int(C['nsrc']), pairs=int(len(zz)), pairs_high_latitude=int((np.abs(b) >= CFG['bcut']).sum()),
                           with_redshift=int(keep.sum()), by_catalogue={n_: int((wh[keep] == k).sum()) for k, n_ in enumerate(('milliquas', 'desi', 'quaia'))}))
    log('counts %s' % json.dumps(res['counts']))
    if keep.sum() < CFG['min_pairs']:
        res['decision'] = dict(outcome='not tested (fewer than 300 pairs with a redshift)')
    else:
        z, d, e2, sep = zz[keep], P['d'][keep], P['e2'][keep], P['sep'][keep]
        K = CFG['bins_many'] if keep.sum() >= CFG['many'] else CFG['bins_few']
        res['pairs'] = analyse(z, d, e2, sep, K, CFG['nshuf'], CFG['seed'], log=log)
        share_q = float((wh[keep] == 2).mean())
        res['decision'] = decide_177(res['pairs'], share_q)
        log('decision %s' % json.dumps(res['decision']))
        rep = {}
        spec = wh[keep] < 2
        if spec.sum() >= CFG['min_pairs']:
            rep['spectroscopic_only'] = analyse(z[spec], d[spec], e2[spec], sep[spec], CFG['bins_few'], 300, CFG['seed'] + 1, log=log, label='Milliquas and DESI redshifts only')
            rep['spectroscopic_only']['decision'] = decide_177(rep['spectroscopic_only'], 0.0)
        for name, m in (('separation below the median', sep < np.median(sep)), ('separation above the median', sep >= np.median(sep)),
                        ('smaller-error half', e2 <= np.median(e2))):
            rep[name] = analyse(z[m], d[m], e2[m], sep[m], CFG['bins_few'], 300, CFG['seed'] + 2, log=log, label=name)
        # the projected separation of the two sight lines at half the source's redshift, flat LCDM (H0 70, Om 0.3)
        def da(zv, n=200):
            zv = np.atleast_1d(zv); t = (np.arange(n) + 0.5)[None, :] * zv[:, None] / n
            return 299792.458 / 70 * (1 / np.sqrt(0.3 * (1 + t) ** 3 + 0.7)).sum(1) * zv / n / (1 + zv)
        rep['separation'] = dict(arcsec_quartiles=[float(v) for v in np.percentile(sep, [25, 50, 75])],
                                 kpc_at_half_redshift_quartiles=[float(v) for v in np.percentile(sep / 3600 * D2R * da(z / 2) * 1000, [25, 50, 75])],
                                 z_quartiles=[float(v) for v in np.percentile(z, [25, 50, 75])])
        rep['whole_sample'] = dict(V=V_of(d, e2), mean_error2=float(np.mean(e2)), raw_sigma=float(1.4826 * np.median(np.abs(d - np.median(d)))))
        res['reported_only'] = rep
        os.makedirs(WORK, exist_ok=True)
        np.savez(os.path.join(WORK, 'pairs.npz'), z=z, d=d, e2=e2, sep=sep, which=wh[keep])
    res['seconds'] = round(time.time() - t0)
    json.dump(res, open(os.path.join(HERE, 'results.json'), 'w'), indent=1, default=float)
    log('written results.json %d s' % res['seconds'])


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else ''
    {'extract': extract, 'synthetic': synthetic, 'run': run}.get(cmd, lambda: print(__doc__))()
