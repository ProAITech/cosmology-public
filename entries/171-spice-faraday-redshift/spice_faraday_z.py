"""Entry 171: does the scatter of Faraday rotation grow with quasar redshift?  Entry 164's rule, repeated on SPICE-RACS DR2.

The rule was frozen in the Workbench (section 171, "the rule, frozen first") before any SPICE-RACS rotation measure was matched to a quasar.
It is Entry 164's rule (section 164) with the NVSS rotation measures replaced by SPICE-RACS DR2 and Entry 166's quality selection.

Data (raw files are git-ignored under data/):
  data/spice-racs-dr2/spice-racs.dr2.fits.gz   SPICE-RACS DR2 (Thomson et al. 2026, arXiv:2605.16917), CSIRO DAP csiro:64891,
      https://data.csiro.au/collection/csiro:64891, 4,968,301,262 bytes, sha256 a1a2f3c1d12a8a8b3aa8526d74b8e0115054f610e2e4186108f0f068005d0053.
      One BINTABLE, 9,294,225 rows x 1030 bytes, 125 columns, read in chunks from the gzip stream.  Columns used: ra, dec, rm, rm_err,
      source_id, snr_flag, leakage_flag, channel_flag, stokesI_fit_flag (Entry 166's four quality flags) and complex_flag.
  data/milliquas-v8/milliquas.zip              Milliquas v8 (Flesch 2023, arXiv:2308.01505), https://quasars.org/milliquas.zip,
      32,701,849 bytes, sha256 9128d6aa91354d0ee1c83dde3ab8764a20d542aa82327e4a969a4357ba9ac984; holds milliquas.txt (1,021,800 records
      of 188 bytes).  The same catalogue Entry 164 read from VizieR (VII/294).  Fixed-width layout from Milliquas-ReadMe.txt:
      RA bytes 1-11, Dec 13-23 (J2000 degrees), Type 52-55, z 77-82.  Read straight from the zip; nothing is unpacked to disk.

Rule, in brief (section 171 has the full text):
  Sample     Milliquas Type beginning with Q, z > 0; matched to the nearest SPICE-RACS component within 10 arcsec that has the four flags
             clear and a finite RM with rm_err > 0; the component's source (source_id) carries the RM, its rows merged by inverse-variance
             weighted mean, at the position of its first row; one quasar per source (the closest); Galactic |b| >= 20 deg; at least
             10 other sources within 3 deg.
  Foreground dRM = RM - median RM of all other sources within 3 deg.
  Statistic  8 equal-count bins in z; sigma = 1.4826 MAD of dRM, error from 500 bootstraps; weighted least-squares slope b of sigma
             against the bin's mean log10(1+z).
  Calibration 1,000 shuffles of z within 12 cells (|b| 20-35-50-65-90 x tercile of e_RM); t = (b - mean)/sd.
  Decision   t >= +3 grows; t <= -3 falls; |t| < 2 flat (95% limit b + 1.645 sd); else unresolved.
  Checks     decisive for a growth claim (t >= 1.5 each, as in Entry 164): each Galactic hemisphere, the smaller-e_RM half, all classes
             (Q, A, B, K, N), longitude-sector calibration; 300 shuffles each.  Reported only: a 1 deg foreground (>= 5 others),
             Faraday-simple sources (complex_flag clear on every row), and chance matches with every quasar moved 5 arcmin north.

Run from the repository root:
  .venv/bin/python -I entries/171-spice-faraday-redshift/spice_faraday_z.py synthetic   # tests on synthetic skies; reads no data
  .venv/bin/python -I entries/171-spice-faraday-redshift/spice_faraday_z.py extract     # about 1 minute; writes work/spice_rm.npz
  .venv/bin/python -I entries/171-spice-faraday-redshift/spice_faraday_z.py run         # the test; writes results.json
"""
import gzip, json, os, sys, time, zipfile
import numpy as np
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
SPICE = os.path.join(ROOT, 'data', 'spice-racs-dr2', 'spice-racs.dr2.fits.gz')
MQZIP = os.path.join(ROOT, 'data', 'milliquas-v8', 'milliquas.zip')
WORK = os.path.join(HERE, 'work')
D2R = np.pi / 180

CFG = dict(match_arcsec=10.0, fg_radius=3.0, min_neigh=10, bcut=20.0, nbins=8, nboot=500, nshuf=1000, nshuf_check=300,
           seed=20261121, bands=[20, 35, 50, 65, 90], classes=['Q', 'A', 'B', 'K', 'N'],
           fg_radius_check=1.0, min_neigh_check=5, chance_shift_arcmin=5.0)

# J2000 equatorial -> Galactic rotation (the matrix Entry 164 used)
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


# ---------------------------------------------------------------- data
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
            code, size = {'D': ('>f8', 8), 'E': ('>f4', 4), 'K': ('>i8', 8), 'J': ('>i4', 4), 'L': ('S1', 1)}[t]
        fields[n] = (code, off)
        off += size
    rowbytes = int(card(h, 'NAXIS1'))
    assert off == rowbytes, (off, rowbytes)
    return np.dtype({'names': list(fields), 'formats': [v[0] for v in fields.values()], 'offsets': [v[1] for v in fields.values()],
                     'itemsize': rowbytes}), rowbytes


FLAGS = ('snr_flag', 'leakage_flag', 'channel_flag', 'stokesI_fit_flag')


def extract():
    """Stream the SPICE-RACS table once; keep rows with the four flags clear (False) and finite rm, rm_err > 0."""
    t0 = time.time()
    keep = {k: [] for k in ('ra', 'dec', 'l', 'b', 'rm', 'rm_err', 'source_id', 'complex')}
    counts = {'rows': 0, 'flags_clear': 0, 'usable': 0, 'complex_flag': {'T': 0, 'F': 0, 'other': 0}}
    with gzip.open(SPICE, 'rb') as f:
        read_header(f)
        h = read_header(f)
        dt, rowbytes = table_dtype(h)
        nrows = int(card(h, 'NAXIS2'))
        done = 0
        while done < nrows:
            n = min(100_000, nrows - done)
            buf = f.read(n * rowbytes)
            if len(buf) != n * rowbytes:
                raise SystemExit('short read at row %d' % done)
            a = np.frombuffer(buf, dtype=dt)
            clear = np.ones(n, bool)
            for k in FLAGS:
                clear &= a[k] == b'F'
            rm, e = a['rm'].astype('<f8'), a['rm_err'].astype('<f8')
            ok = clear & np.isfinite(rm) & np.isfinite(e) & (e > 0)
            counts['rows'] += n; counts['flags_clear'] += int(clear.sum()); counts['usable'] += int(ok.sum())
            cf = a['complex_flag'][ok]
            for sym in ('T', 'F'):
                counts['complex_flag'][sym] += int((cf == sym.encode()).sum())
            counts['complex_flag']['other'] += int(((cf != b'T') & (cf != b'F')).sum())
            for k in ('ra', 'dec', 'l', 'b'):
                keep[k].append(a[k][ok].astype('<f8'))
            keep['rm'].append(rm[ok]); keep['rm_err'].append(e[ok])
            keep['source_id'].append(a['source_id'][ok].copy())
            keep['complex'].append(cf == b'T')
            done += n
    os.makedirs(WORK, exist_ok=True)
    out = {k: np.concatenate(v) for k, v in keep.items()}
    np.savez(os.path.join(WORK, 'spice_rm.npz'), **out)
    counts['seconds'] = round(time.time() - t0)
    json.dump(counts, open(os.path.join(WORK, 'extract_counts.json'), 'w'), indent=1)
    print('extract', counts, flush=True)


def merge_sources(rows):
    """Entry 166's merge: per source_id, inverse-variance weighted mean RM, error 1/sqrt(sum w), position of the first row."""
    sid = rows['source_id']
    u, inv = np.unique(sid, return_inverse=True)
    w = 1.0 / rows['rm_err'] ** 2
    sw = np.bincount(inv, weights=w)
    rm = np.bincount(inv, weights=w * rows['rm']) / sw
    first = np.zeros(len(u), dtype=int); first[inv[::-1]] = np.arange(len(sid))[::-1]
    nrow = np.bincount(inv)
    cplx = np.bincount(inv, weights=rows['complex'].astype(float)) > 0
    ra, dec = rows['ra'][first], rows['dec'][first]
    l, b = galactic(ra, dec)
    return dict(ra=ra, dec=dec, l=l, b=b, rm=rm, e=1.0 / np.sqrt(sw), simple=~cplx, nrow=nrow), inv


def load_spice():
    z = np.load(os.path.join(WORK, 'spice_rm.npz'))
    rows = {k: z[k] for k in z.files}
    src, inv = merge_sources(rows)
    l_cat, b_cat = rows['l'], rows['b']
    l_me, b_me = galactic(rows['ra'], rows['dec'])
    dl = np.abs((l_me - l_cat + 180) % 360 - 180) * np.cos(b_cat * D2R)
    info = {'rows_usable': int(len(inv)), 'sources': int(len(src['rm'])), 'multi_row_sources': int((src['nrow'] > 1).sum()),
            'max_offset_galactic_vs_catalogue_deg': float(max(np.nanmax(dl), np.nanmax(np.abs(b_me - b_cat)))),
            'dec_range': [float(rows['dec'].min()), float(rows['dec'].max())],
            'median_e_source': float(np.median(src['e']))}
    return dict(row_ra=rows['ra'], row_dec=rows['dec'], row_src=inv), src, info


def load_milliquas():
    with zipfile.ZipFile(MQZIP) as zf:
        raw = zf.read('milliquas.txt')
    ra, dec, z, cls, bad = [], [], [], [], 0
    for line in raw.split(b'\n'):
        line = line.rstrip(b'\r')
        if len(line) < 82:
            if line.strip():
                bad += 1
            continue
        try:
            r, d = float(line[0:11]), float(line[12:23])
        except ValueError:
            bad += 1
            continue
        ty = line[51:55].decode('ascii', 'replace').strip()
        zs = line[76:82].strip()
        try:
            zz = float(zs) if zs else np.nan
        except ValueError:
            zz = np.nan
        ra.append(r); dec.append(d); z.append(zz); cls.append(ty[:1] if ty else '')
    cls = np.array(cls)
    mq = dict(ra=np.array(ra), dec=np.array(dec), z=np.array(z), cls=cls)
    info = {'records': int(len(cls)), 'unparsed_lines': bad,
            'with_z_by_class': {c: int(((cls == c) & (mq['z'] > 0)).sum()) for c in CFG['classes']}}
    return mq, info


# ---------------------------------------------------------------- matching and sample
def match(mq, rows, nsrc, classes, shift_arcmin=0.0):
    """Each quasar to the nearest usable component within 10"; per source keep the closest quasar. Returns {source: (sep, quasar)}."""
    sel = np.isin(mq['cls'], list(classes)) & (mq['z'] > 0)
    qi = np.nonzero(sel)[0]
    qdec = mq['dec'][qi] + shift_arcmin / 60.0
    tree = cKDTree(unit(rows['row_ra'], rows['row_dec']))
    d, j = tree.query(unit(mq['ra'][qi], qdec), k=1, distance_upper_bound=chord(CFG['match_arcsec'] / 3600.0))
    hit = np.isfinite(d)
    sep = 2 * np.arcsin(d[hit] / 2) / D2R * 3600
    src = rows['row_src'][j[hit]]
    q = qi[hit]
    order = np.lexsort((sep, src))                     # by source, then separation: first of each source is the closest
    src, sep, q = src[order], sep[order], q[order]
    first = np.r_[True, src[1:] != src[:-1]]
    return dict(src=src[first], sep=sep[first], q=q[first], n_quasars_matched=int(hit.sum()))


def foreground(src, idx, radius, min_neigh, tree=None):
    tree = tree or cKDTree(unit(src['ra'], src['dec']))
    v = unit(src['ra'][idx], src['dec'][idx])
    lists = tree.query_ball_point(v, chord(radius))
    med, nn = np.full(len(idx), np.nan), np.zeros(len(idx), int)
    for k, (i, L) in enumerate(zip(idx, lists)):
        L = [m for m in L if m != i]
        nn[k] = len(L)
        if len(L) >= min_neigh:
            med[k] = np.median(src['rm'][L])
    return med, nn


def build(m, src, mq, radius, min_neigh, tree=None):
    idx = m['src']
    lat = np.abs(src['b'][idx]) >= CFG['bcut']
    med, nn = foreground(src, idx[lat], radius, min_neigh, tree)
    okn = np.isfinite(med)
    k = idx[lat][okn]
    q = m['q'][lat][okn]
    S = dict(z=mq['z'][q], d=src['rm'][k] - med[okn], e=src['e'][k], b=src['b'][k], l=src['l'][k],
             simple=src['simple'][k], cls=mq['cls'][q], sep=m['sep'][lat][okn], nneigh=nn[okn])
    cnt = dict(matched=int(len(idx)), dropped_latitude=int((~lat).sum()), dropped_neighbours=int((~okn).sum()), kept=int(len(k)))
    return S, cnt


# ---------------------------------------------------------------- statistic, calibration, decision (Entry 164's, ported)
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


def bins_fit(z, d, rng):
    n, K = len(z), CFG['nbins']
    idx = np.argsort(z, kind='stable')
    x, sg, er, cn = [], [], [], []
    for k in range(K):
        lo, hi = jsround(k * n / K), jsround((k + 1) * n / K)
        q = idx[lo:hi]; m = len(q)
        dv = d[q]
        s0 = float(sig_mad(dv))
        bs = sig_mad(dv[rng.integers(0, m, size=(CFG['nboot'], m))])
        x.append(float(np.log10(1 + z[q]).mean())); sg.append(s0); er.append(float(bs.std(ddof=1))); cn.append(m)
    f = wls(x, sg, er)
    return dict(b=float(f['b']), a=float(f['a']), eb=float(f['eb']), chi2=f['chi2'], x=x, sigma=sg, err=er, n=cn)


def cells_for(S, sectors=False):
    n = len(S['z'])
    es = np.sort(S['e'])
    t1, t2 = es[int(np.floor(n / 3))], es[int(np.floor(2 * n / 3))]
    ab = np.abs(S['b'])
    band = np.zeros(n, int)
    for j in range(3):
        band += ab >= CFG['bands'][j + 1]
    ter = np.where(S['e'] < t1, 0, np.where(S['e'] < t2, 1, 2))
    c = band * 3 + ter
    if sectors:
        c = c * 4 + np.minimum(3, np.floor(S['l'] / 90).astype(int))
    return c


def analyse(S, nshuf, label, sectors=False, seed=None, log=print):
    seed = CFG['seed'] if seed is None else seed
    t0 = time.time()
    z, d = S['z'], S['d']
    real = bins_fit(z, d, np.random.default_rng(seed))
    cell = cells_for(S, sectors)
    groups = [np.nonzero(cell == c)[0] for c in np.unique(cell)]
    R = np.random.default_rng(seed + 7)
    nulls = np.empty(nshuf)
    for s in range(nshuf):
        z2 = z.copy()
        for g in groups:
            z2[g] = z[g][R.permutation(len(g))]
        nulls[s] = bins_fit(z2, d, np.random.default_rng(seed + 1000 + s))['b']
    mu, sd = float(nulls.mean()), float(nulls.std(ddof=1))
    t = (real['b'] - mu) / sd
    out = dict(label=label, n=int(len(z)), cells=len(groups), b=real['b'], eb_fit=real['eb'], chi2=real['chi2'],
               bins=dict(x=real['x'], sigma=real['sigma'], err=real['err'], n=real['n']),
               null=dict(mean=mu, sd=sd, n=nshuf, q=[float(np.quantile(nulls, q, method='nearest')) for q in (0.0013, 0.025, 0.5, 0.975, 0.9987)],
                         hist=hist(nulls)),
               t=float(t), p_upper=float((1 + (nulls >= real['b']).sum()) / (1 + nshuf)),
               p_lower=float((1 + (nulls <= real['b']).sum()) / (1 + nshuf)), limit95=float(real['b'] + 1.645 * sd),
               added_at_z1_limit=float(0.30103 * (real['b'] + 1.645 * sd)),
               z_median=float(np.median(z)), z_range=[float(z.min()), float(z.max())], e_median=float(np.median(S['e'])))
    log('%s: n %d, b %+.3f, null %+.3f +- %.3f, t %+.2f, %.0f s' % (label, out['n'], out['b'], mu, sd, t, time.time() - t0))
    return out


def hist(v):
    lo, hi = int(np.floor(v.min())), int(np.ceil(v.max()))
    nb = max(1, hi - lo)
    c = np.bincount(np.minimum(nb - 1, np.floor(v - lo).astype(int)), minlength=nb)
    return dict(lo=lo, width=1, counts=c.tolist())


def decide(t):
    if t >= 3:
        return 'grows with redshift (t >= +3)'
    if t <= -3:
        return 'falls with redshift (t <= -3)'
    if abs(t) < 2:
        return 'flat, the standard picture (|t| < 2)'
    return 'unresolved (2 <= |t| < 3)'


def sub(S, keep):
    return {k: v[keep] for k, v in S.items()}


def full_test(rows, src, mq, nshuf, nshuf_check, log=print):
    tree = cKDTree(unit(src['ra'], src['dec']))
    mQ = match(mq, rows, len(src['rm']), ['Q'])
    mA = match(mq, rows, len(src['rm']), CFG['classes'])
    SQ, cQ = build(mQ, src, mq, CFG['fg_radius'], CFG['min_neigh'], tree)
    SA, cA = build(mA, src, mq, CFG['fg_radius'], CFG['min_neigh'], tree)
    mC = match(mq, rows, len(src['rm']), ['Q'], shift_arcmin=CFG['chance_shift_arcmin'])
    res = dict(counts=dict(Q=dict(cQ, quasars_within_10arcsec=mQ['n_quasars_matched']),
                           all_classes=dict(cA, quasars_within_10arcsec=mA['n_quasars_matched']),
                           chance_matches_Q_shifted_5arcmin=int(len(mC['src'])),
                           match_sep_arcsec_median=float(np.median(SQ['sep'])) if len(SQ['z']) else None,
                           neighbours_median=float(np.median(SQ['nneigh'])) if len(SQ['z']) else None))
    log('counts %s' % json.dumps(res['counts']))
    if cQ['kept'] < 1000:
        res['decision'] = dict(outcome='not tested (fewer than 1,000 quasars kept)')
        return res
    res['primary'] = analyse(SQ, nshuf, 'primary', log=log)
    res['decision'] = dict(t=res['primary']['t'], outcome=decide(res['primary']['t']))
    ch = {}
    ch['north'] = analyse(sub(SQ, SQ['b'] > 0), nshuf_check, 'north', log=log)
    ch['south'] = analyse(sub(SQ, SQ['b'] < 0), nshuf_check, 'south', log=log)
    ch['low_eRM'] = analyse(sub(SQ, SQ['e'] <= np.median(SQ['e'])), nshuf_check, 'low e_RM', log=log)
    ch['all_classes'] = analyse(SA, nshuf_check, 'all classes', log=log)
    ch['longitude_sectors'] = analyse(SQ, nshuf_check, 'longitude sectors', sectors=True, log=log)
    res['checks'] = ch
    ts = [c['t'] for c in ch.values()]
    res['decision']['checks_t'] = {k: round(c['t'], 2) for k, c in ch.items()}
    res['decision']['growth_claim'] = (('claimed: every check reaches t >= 1.5' if all(t >= 1.5 for t in ts) else
                                        'unresolved, systematic suspected: not every check reaches t >= 1.5')
                                       if res['primary']['t'] >= 3 else 'not applicable (primary t < +3)')
    rep = {}
    S1, c1 = build(mQ, src, mq, CFG['fg_radius_check'], CFG['min_neigh_check'], tree)
    rep['foreground_1deg'] = analyse(S1, nshuf_check, 'foreground 1 deg', log=log); rep['foreground_1deg']['counts'] = c1
    rep['faraday_simple'] = analyse(sub(SQ, SQ['simple']), nshuf_check, 'Faraday-simple', log=log)
    res['reported_only'] = rep
    return res


# ---------------------------------------------------------------- synthetic skies (no data read)
def synthetic_sky(seed, n_src=60000, n_q=4000, grow=0.0):
    """Sources uniform on Dec < +49; RM = smooth 'Galactic' field (a few harmonics, stronger at low |b|) + noise.
    A random subset of sources are quasars with z drawn from a broad distribution; 'grow' adds Gaussian scatter
    of grow * log10(1+z) rad/m^2 to the quasars' RMs (an intergalactic term)."""
    g = np.random.default_rng(seed)
    u = g.uniform(np.sin(-90 * D2R), np.sin(49 * D2R), n_src)
    dec = np.arcsin(u) / D2R
    ra = g.uniform(0, 360, n_src)
    l, b = galactic(ra, dec)
    fg = 25 / (np.abs(np.sin(b * D2R)) + 0.3) * np.sin(l * D2R + 0.7) + 8 * np.cos(3 * l * D2R) * np.sin(2 * b * D2R)
    e = g.lognormal(np.log(2.0), 0.5, n_src)
    rm = fg + g.normal(0, 6, n_src) + g.normal(0, 1, n_src) * e
    qsrc = g.choice(n_src, n_q, replace=False)
    zq = np.clip(g.gamma(4.0, 0.35, n_q), 0.05, 6)
    rm[qsrc] += g.normal(0, 1, n_q) * grow * np.log10(1 + zq)
    rows = dict(row_ra=ra, row_dec=dec, row_src=np.arange(n_src))
    src = dict(ra=ra, dec=dec, l=l, b=b, rm=rm, e=e, simple=g.uniform(size=n_src) > 0.2, nrow=np.ones(n_src, int))
    off = g.uniform(0, 4, n_q) / 3600                  # quasar a few arcsec from its radio source
    pa = g.uniform(0, 2 * np.pi, n_q)
    qdec = dec[qsrc] + off * np.cos(pa)
    qra = ra[qsrc] + off * np.sin(pa) / np.cos(dec[qsrc] * D2R)
    nf = 3 * n_q                                        # unmatched field quasars
    fra, fdec = g.uniform(0, 360, nf), np.arcsin(g.uniform(-1, 1, nf)) / D2R
    mq = dict(ra=np.r_[qra, fra], dec=np.r_[qdec, fdec], z=np.r_[zq, g.gamma(4, 0.35, nf)],
              cls=np.r_[np.full(n_q, 'Q'), g.choice(['Q', 'A', 'N', ''], nf)])
    return rows, src, mq


def synthetic():
    out = {}
    for name, grow, seed in (('null_1', 0.0, 1), ('null_2', 0.0, 2), ('null_3', 0.0, 3), ('growth_20', 20.0, 4), ('growth_40', 40.0, 5)):
        rows, src, mq = synthetic_sky(seed, grow=grow)
        tree = cKDTree(unit(src['ra'], src['dec']))
        m = match(mq, rows, len(src['rm']), ['Q'])
        S, c = build(m, src, mq, CFG['fg_radius'], CFG['min_neigh'], tree)
        r = analyse(S, 200, name)
        out[name] = dict(counts=c, t=round(r['t'], 2), b=round(r['b'], 2), null_sd=round(r['null']['sd'], 2), outcome=decide(r['t']),
                         sigma=[round(s, 1) for s in r['bins']['sigma']])
    # parser check on two lines written to the ReadMe's layout
    line = ('%11.7f %11.7f  %-25s %-4s %5.2f %5.2f %-3s %1s %1s %6.3f %-6s %-6s' %
            (123.4567891, -12.3456789, 'SDSS J081349.62-122044.4', 'QR', 19.12, 19.80, 'p', '+', '+', 1.234, 'DR16Q', 'DR16Q')).ljust(188)
    b = line.encode()
    out['parser'] = dict(ra=float(b[0:11]), dec=float(b[12:23]), type=b[51:55].decode().strip(), z=float(b[76:82]))
    print(json.dumps(out, indent=1))
    json.dump(out, open(os.path.join(HERE, 'synthetic.json'), 'w'), indent=1)


def run():
    t0 = time.time()
    rows, src, sinfo = load_spice()
    print('spice', sinfo, flush=True)
    mq, minfo = load_milliquas()
    print('milliquas', minfo, flush=True)
    res = dict(generated=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), config=CFG, spice=sinfo, milliquas=minfo)
    res.update(full_test(rows, src, mq, CFG['nshuf'], CFG['nshuf_check'], log=lambda s: print(s, flush=True)))
    res['seconds'] = round(time.time() - t0)
    js = json.dumps(res, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
    open(os.path.join(HERE, 'results.json'), 'w').write(js)
    print('decision', json.dumps(res['decision']), flush=True)
    print('written results.json %d s' % res['seconds'], flush=True)


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else ''
    {'extract': extract, 'synthetic': synthetic, 'run': run}.get(cmd, lambda: print(__doc__))()
