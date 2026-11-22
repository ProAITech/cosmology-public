"""Entry 173: is the growth of Faraday scatter with quasar redshift (Entry 171) carried by sight lines with Mg II absorbers?

The rule was frozen in the Workbench (section 173, "the rule, frozen first") before any DESI quasar was matched to a rotation measure.
Entry 171's statistic is applied to SPICE-RACS sources behind DESI DR1 quasars (sample D), to the sight lines among them with no
Mg II absorber in DESI's catalogue (sample N), and the slope of N is set against random subsets of D with N's size and redshifts.

Data (raw files are git-ignored under data/):
  data/spice-racs-dr2/spice-racs.dr2.fits.gz   SPICE-RACS DR2 (Thomson et al. 2026, arXiv:2605.16917), CSIRO DAP csiro:64891,
      https://data.csiro.au/collection/csiro:64891, 4,968,301,262 bytes, sha256 a1a2f3c1d12a8a8b3aa8526d74b8e0115054f610e2e4186108f0f068005d0053.
      One BINTABLE, 9,294,225 rows.  Columns used: ra, dec, rm, rm_err, source_id, the four quality flags of Entry 166, complex_flag.
  data/desi-dr1-qso/QSO_cat_iron_cumulative_v0.fits   DESI DR1 quasar catalogue (DESI Collaboration 2025, arXiv:2503.14745),
      https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/QSO/iron/QSO_cat_iron_cumulative_v0.fits, 816,203,520 bytes,
      sha256 a65604303cc29af8865db9193a23934b2540591b16c5afd2a4bea2d61464aef4.  HDU QSO_CAT, 2,182,309 rows.
      Columns used: TARGETID, Z, TARGET_RA, TARGET_DEC, TSNR2_QSO.
  data/desi-dr1-mgii/MgII-Absorbers-DR1.fits   DESI DR1 Mg II absorber catalogue v1.0 (method: Napolitano et al. 2023, arXiv:2305.20016),
      https://data.desi.lbl.gov/public/dr1/vac/dr1/mgii-absorber/v1.0/MgII-Absorbers-DR1.fits, 65,750,400 bytes (checksum file beside it).
      HDU MGII_ABSORBERS, 270,529 rows.  Columns used: TARGETID, Z_QSO, Z_MGII, EW_2796.
  data/milliquas-v8/milliquas.zip              Milliquas v8 (Flesch 2023), used only to count how many sight lines Entry 171 also had.

Rule, in brief (section 173 has the full text):
  Sample D   DESI quasar rows, one per TARGETID (the row with the largest TSNR2_QSO), Z > 0; each matched to the nearest SPICE-RACS
             component within 10 arcsec (Entry 171's quality cuts and merge); one quasar per source (the closest); |b| >= 20 deg;
             at least 10 other sources within 3 deg.  dRM = RM - median RM of the other sources within 3 deg.
  A and N    A: the quasar's TARGETID has at least one row in the absorber catalogue.  N: it has none.
  Statistics Entry 171's slope b (8 equal-count z bins, sigma = 1.4826 MAD, bootstrap errors, weighted least squares against mean
             log10(1+z)) and its shuffle calibration (1,000 shuffles of z within 12 cells) for D and for N: t_D, t_N.
             t_share: b_N against 1,000 random subsets of D that take, in each of 16 equal-count redshift strata of D, as many
             sight lines as N has there; t_share = (b_N - mean)/sd.  drop = mean - b_N.
  Decision   see decide_173() below; it is the frozen text of section 173, item 6.

Run from the repository root:
  .venv/bin/python -I entries/173-spice-desi-absorbers/absorber_split.py synthetic   # synthetic skies; reads no data
  .venv/bin/python -I entries/173-spice-desi-absorbers/absorber_split.py extract     # about 1 minute; writes work/spice_rm.npz
  .venv/bin/python -I entries/173-spice-desi-absorbers/absorber_split.py run         # the test; writes results.json
The functions down to analyse() are Entry 171's (entries/171-spice-faraday-redshift/spice_faraday_z.py), copied unchanged so that
this folder stands alone; only the 'I' column type was added to the FITS reader.
"""
import gzip, json, os, sys, time, zipfile
import numpy as np
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
SPICE = os.path.join(ROOT, 'data', 'spice-racs-dr2', 'spice-racs.dr2.fits.gz')
QSO = os.path.join(ROOT, 'data', 'desi-dr1-qso', 'QSO_cat_iron_cumulative_v0.fits')
MGII = os.path.join(ROOT, 'data', 'desi-dr1-mgii', 'MgII-Absorbers-DR1.fits')
MQZIP = os.path.join(ROOT, 'data', 'milliquas-v8', 'milliquas.zip')
WORK = os.path.join(HERE, 'work')
D2R = np.pi / 180

CFG = dict(match_arcsec=10.0, fg_radius=3.0, min_neigh=10, bcut=20.0, nbins=8, nboot=500, nshuf=1000, nshuf_check=300,
           seed=20261121, bands=[20, 35, 50, 65, 90], fg_radius_check=1.0, min_neigh_check=5,
           nstrata=16, ndraw=1000, ndraw_check=300, min_D=1500, min_N=1000, strong_ew=1.0, zmin_search=0.4,
           growth_171=5.083 - 0.451)

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


def extract():
    """Stream the SPICE-RACS table once; keep rows with the four flags clear (False) and finite rm, rm_err > 0.  Entry 171's."""
    t0 = time.time()
    keep = {k: [] for k in ('ra', 'dec', 'l', 'b', 'rm', 'rm_err', 'source_id', 'complex')}
    counts = {'rows': 0, 'flags_clear': 0, 'usable': 0}
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
            for k in ('ra', 'dec', 'l', 'b'):
                keep[k].append(a[k][ok].astype('<f8'))
            keep['rm'].append(rm[ok]); keep['rm_err'].append(e[ok])
            keep['source_id'].append(a['source_id'][ok].copy())
            keep['complex'].append(a['complex_flag'][ok] == b'T')
            done += n
    os.makedirs(WORK, exist_ok=True)
    np.savez(os.path.join(WORK, 'spice_rm.npz'), **{k: np.concatenate(v) for k, v in keep.items()})
    counts['seconds'] = round(time.time() - t0)
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
    p = os.path.join(WORK, 'spice_rm.npz')
    if not os.path.exists(p):                              # Entry 171's cache holds the same rows
        p = os.path.join(ROOT, 'entries', '171-spice-faraday-redshift', 'work', 'spice_rm.npz')
    z = np.load(p)
    rows = {k: z[k] for k in z.files}
    src, inv = merge_sources(rows)
    info = {'rows_usable': int(len(inv)), 'sources': int(len(src['rm'])), 'cache': os.path.relpath(p, ROOT),
            'dec_range': [float(rows['dec'].min()), float(rows['dec'].max())], 'median_e_source': float(np.median(src['e']))}
    return dict(row_ra=rows['ra'], row_dec=rows['dec'], row_src=inv), src, info


def load_desi():
    """One row per TARGETID (largest TSNR2_QSO), Z > 0 and finite position."""
    t, nrows = read_table(QSO, ['TARGETID', 'Z', 'TARGET_RA', 'TARGET_DEC', 'TSNR2_QSO'])
    tid, ts = t['TARGETID'], t['TSNR2_QSO'].astype(float)
    tsf = np.where(np.isfinite(ts), ts, -np.inf)
    order = np.lexsort((-tsf, tid))
    first = np.r_[True, tid[order][1:] != tid[order][:-1]]
    k = order[first]
    z, ra, dec = t['Z'][k].astype(float), t['TARGET_RA'][k].astype(float), t['TARGET_DEC'][k].astype(float)
    ok = np.isfinite(z) & (z > 0) & np.isfinite(ra) & np.isfinite(dec)
    q = dict(tid=tid[k][ok], z=z[ok], ra=ra[ok], dec=dec[ok], tsnr2=ts[k][ok])
    info = {'rows': int(nrows), 'unique_targetid': int(len(k)), 'kept_z_positive': int(ok.sum())}
    return q, info


def load_absorbers():
    t, nrows = read_table(MGII, ['TARGETID', 'Z_QSO', 'Z_MGII', 'EW_2796'])
    info = {'rows': int(nrows), 'sightlines': int(len(np.unique(t['TARGETID'])))}
    return dict(tid=t['TARGETID'], zq=t['Z_QSO'].astype(float), zabs=t['Z_MGII'].astype(float), ew=t['EW_2796'].astype(float)), info


def load_milliquas():
    """Type Q, z > 0 only: positions and redshifts, to count the overlap with Entry 171's sample."""
    with zipfile.ZipFile(MQZIP) as zf:
        raw = zf.read('milliquas.txt')
    ra, dec, z = [], [], []
    for line in raw.split(b'\n'):
        line = line.rstrip(b'\r')
        if len(line) < 82 or line[51:52] != b'Q':
            continue
        try:
            r, d, zz = float(line[0:11]), float(line[12:23]), float(line[76:82])
        except ValueError:
            continue
        if zz > 0:
            ra.append(r); dec.append(d); z.append(zz)
    return dict(ra=np.array(ra), dec=np.array(dec), z=np.array(z))


# ---------------------------------------------------------------- matching and sample (Entry 171's, quasar classes removed)
def match(q, rows):
    """Each quasar to the nearest usable component within 10"; per source keep the closest quasar."""
    tree = cKDTree(unit(rows['row_ra'], rows['row_dec']))
    d, j = tree.query(unit(q['ra'], q['dec']), k=1, distance_upper_bound=chord(CFG['match_arcsec'] / 3600.0))
    hit = np.isfinite(d)
    sep = 2 * np.arcsin(d[hit] / 2) / D2R * 3600
    src = rows['row_src'][j[hit]]
    qi = np.nonzero(hit)[0]
    order = np.lexsort((sep, src))
    src, sep, qi = src[order], sep[order], qi[order]
    first = np.r_[True, src[1:] != src[:-1]]
    return dict(src=src[first], sep=sep[first], q=qi[first], n_quasars_matched=int(hit.sum()))


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


def build(m, src, q, ab, radius, min_neigh, tree=None):
    idx = m['src']
    lat = np.abs(src['b'][idx]) >= CFG['bcut']
    med, nn = foreground(src, idx[lat], radius, min_neigh, tree)
    okn = np.isfinite(med)
    k = idx[lat][okn]
    qq = m['q'][lat][okn]
    tid = q['tid'][qq]
    # absorber summary per kept sight line
    u, inv = np.unique(ab['tid'], return_inverse=True)
    nab_all = np.bincount(inv)
    ewmax_all = np.full(len(u), -np.inf); np.maximum.at(ewmax_all, inv, np.where(np.isfinite(ab['ew']), ab['ew'], -np.inf))
    pos = np.searchsorted(u, tid)
    pos = np.minimum(pos, len(u) - 1)
    has = u[pos] == tid
    nab = np.where(has, nab_all[pos], 0)
    ewmax = np.where(has, ewmax_all[pos], np.nan)
    S = dict(z=q['z'][qq], d=src['rm'][k] - med[okn], e=src['e'][k], b=src['b'][k], l=src['l'][k], simple=src['simple'][k],
             sep=m['sep'][lat][okn], nneigh=nn[okn], tid=tid, tsnr2=q['tsnr2'][qq], A=has, nabs=nab, ewmax=ewmax, src=k)
    cnt = dict(matched=int(len(idx)), dropped_latitude=int((~lat).sum()), dropped_neighbours=int((~okn).sum()), kept=int(len(k)),
               with_absorber=int(has.sum()), without_absorber=int((~has).sum()))
    return S, cnt


# ---------------------------------------------------------------- statistic and calibration (Entry 171's, unchanged)
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


def bins_fit(z, d, rng, e=None):
    n, K = len(z), CFG['nbins']
    idx = np.argsort(z, kind='stable')
    x, sg, er, cn, zm, em = [], [], [], [], [], []
    for k in range(K):
        lo, hi = jsround(k * n / K), jsround((k + 1) * n / K)
        q = idx[lo:hi]; m = len(q)
        dv = d[q]
        s0 = float(sig_mad(dv))
        bs = sig_mad(dv[rng.integers(0, m, size=(CFG['nboot'], m))])
        x.append(float(np.log10(1 + z[q]).mean())); sg.append(s0); er.append(float(bs.std(ddof=1))); cn.append(m)
        zm.append(float(np.median(z[q])))
        if e is not None:
            em.append(float(np.mean(e[q] ** 2)))
    f = wls(x, sg, er)
    return dict(b=float(f['b']), a=float(f['a']), eb=float(f['eb']), chi2=f['chi2'], x=x, sigma=sg, err=er, n=cn, zmed=zm, e2mean=em)


def cells_for(S):
    n = len(S['z'])
    es = np.sort(S['e'])
    t1, t2 = es[int(np.floor(n / 3))], es[int(np.floor(2 * n / 3))]
    ab = np.abs(S['b'])
    band = np.zeros(n, int)
    for j in range(3):
        band += ab >= CFG['bands'][j + 1]
    ter = np.where(S['e'] < t1, 0, np.where(S['e'] < t2, 1, 2))
    return band * 3 + ter


def analyse(S, nshuf, label, seed=None, log=print):
    seed = CFG['seed'] if seed is None else seed
    t0 = time.time()
    z, d = S['z'], S['d']
    real = bins_fit(z, d, np.random.default_rng(seed), S['e'])
    cell = cells_for(S)
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
               bins=dict(x=real['x'], sigma=real['sigma'], err=real['err'], n=real['n'], zmed=real['zmed'], e2mean=real['e2mean']),
               null=dict(mean=mu, sd=sd, n=nshuf, q=[float(np.quantile(nulls, q, method='nearest')) for q in (0.0013, 0.025, 0.5, 0.975, 0.9987)]),
               t=float(t), p_upper=float((1 + (nulls >= real['b']).sum()) / (1 + nshuf)), net=float(real['b'] - mu),
               z_median=float(np.median(z)), z_range=[float(z.min()), float(z.max())], e_median=float(np.median(S['e'])))
    log('%s: n %d, b %+.3f, null %+.3f +- %.3f, t %+.2f, %.0f s' % (label, out['n'], out['b'], mu, sd, t, time.time() - t0))
    return out


def sub(S, keep):
    return {k: v[keep] for k, v in S.items()}


# ---------------------------------------------------------------- new in Entry 173
def strata_of(z, K):
    r = np.empty(len(z), int)
    r[np.argsort(z, kind='stable')] = np.arange(len(z))
    return np.minimum(K - 1, r * K // len(z))


def share_test(S, inN, ndraw, label, seed=None, log=print):
    """b of the subset inN against random subsets of S with the same count in each of 16 equal-count redshift strata of S."""
    seed = CFG['seed'] if seed is None else seed
    t0 = time.time()
    z, d = S['z'], S['d']
    st = strata_of(z, CFG['nstrata'])
    members = [np.nonzero(st == s)[0] for s in range(CFG['nstrata'])]
    take = [int(inN[m].sum()) for m in members]
    bN = bins_fit(z[inN], d[inN], np.random.default_rng(seed))['b']
    R = np.random.default_rng(seed + 31)
    bs = np.empty(ndraw)
    for i in range(ndraw):
        pick = np.concatenate([R.choice(m, k, replace=False) for m, k in zip(members, take)])
        bs[i] = bins_fit(z[pick], d[pick], np.random.default_rng(seed + 5000 + i))['b']
    mu, sd = float(bs.mean()), float(bs.std(ddof=1))
    out = dict(label=label, n_parent=int(len(z)), n_subset=int(inN.sum()), b_subset=float(bN), random_mean=mu, random_sd=sd, n_draws=ndraw,
               t_share=float((bN - mu) / sd), drop=float(mu - bN), drop_upper95=float(mu - bN + 1.645 * sd),
               p_lower=float((1 + (bs <= bN).sum()) / (1 + ndraw)),
               strata=dict(z_lo=[float(z[m].min()) for m in members], z_hi=[float(z[m].max()) for m in members],
                           n=[int(len(m)) for m in members], n_subset=take))
    log('%s: subset %d of %d, b %+.3f, random subsets %+.3f +- %.3f, t_share %+.2f, drop %+.2f, %.0f s'
        % (label, out['n_subset'], out['n_parent'], bN, mu, sd, out['t_share'], out['drop'], time.time() - t0))
    return out


def excess_test(S, inA, ndraw, label, seed=None, log=print):
    """Scatter of the subset inA against (a) the rest drawn to inA's redshift strata and (b) random relabelling within strata."""
    seed = CFG['seed'] if seed is None else seed
    z, d = S['z'], S['d']
    st = strata_of(z, CFG['nstrata'])
    K = CFG['nstrata']
    allm = [np.nonzero(st == s)[0] for s in range(K)]
    restm = [np.nonzero((st == s) & ~inA)[0] for s in range(K)]
    kA = [int(inA[m].sum()) for m in allm]
    sA = float(sig_mad(d[inA]))
    R = np.random.default_rng(seed + 77)
    sN, sR = np.empty(ndraw), np.empty(ndraw)
    for i in range(ndraw):
        pn = np.concatenate([R.choice(m, k, replace=len(m) < k) for m, k in zip(restm, kA) if k and len(m)])
        pr = np.concatenate([R.choice(m, k, replace=False) for m, k in zip(allm, kA) if k])
        sN[i], sR[i] = sig_mad(d[pn]), sig_mad(d[pr])
    bsA = sig_mad(d[inA][np.random.default_rng(seed + 78).integers(0, int(inA.sum()), size=(CFG['nboot'], int(inA.sum())))])
    q2 = sA ** 2 - float(sN.mean()) ** 2
    e_q2 = float(np.sqrt((2 * sA * bsA.std(ddof=1)) ** 2 + (2 * sN.mean() * sN.std(ddof=1)) ** 2))
    out = dict(label=label, n_A=int(inA.sum()), n_rest=int((~inA).sum()), sigma_A=sA, sigma_A_err=float(bsA.std(ddof=1)),
               sigma_rest_matched=float(sN.mean()), sigma_rest_matched_sd=float(sN.std(ddof=1)),
               excess_variance=float(q2), excess_variance_err=e_q2, excess_quadrature=float(np.sign(q2) * np.sqrt(abs(q2))),
               relabel_mean=float(sR.mean()), relabel_sd=float(sR.std(ddof=1)), t_excess=float((sA - sR.mean()) / sR.std(ddof=1)),
               p_upper=float((1 + (sR >= sA).sum()) / (1 + ndraw)))
    log('%s: A %d, sigma_A %.2f, rest at A\'s redshifts %.2f, excess variance %+.1f +- %.1f, t_excess %+.2f'
        % (label, out['n_A'], sA, out['sigma_rest_matched'], q2, e_q2, out['t_excess']))
    return out


def decide_173(D, N, sh):
    """Section 173, item 6.  D, N: outputs of analyse(); sh: output of share_test()."""
    tD, tN, ts = D['t'], N['t'], sh['t_share']
    netD = D['net']
    u = (sh['drop'] - 0.5 * netD) / np.sqrt(sh['random_sd'] ** 2 + 0.25 * D['null']['sd'] ** 2)
    out = dict(t_D=tD, t_N=tN, t_share=ts, drop=sh['drop'], drop_upper95=sh['drop_upper95'], net_growth_D=netD,
               share_of_D_growth=(sh['drop'] / netD if tD >= 2.5 else None),
               share_of_entry171_growth=sh['drop'] / CFG['growth_171'],
               share_upper95_entry171=sh['drop_upper95'] / CFG['growth_171'], u_half=float(u))
    if ts <= -3:
        if tN < 2:
            o = 'absorber sight lines carry the growth: removing them lowers the slope (t_share <= -3) and the rest is flat within errors (t_N < 2)'
        elif tN >= 3:
            o = 'absorber sight lines carry part of the growth: removing them lowers the slope (t_share <= -3) and growth persists without them (t_N >= +3)'
        else:
            o = 'absorber sight lines carry growth (t_share <= -3); whether any remains without them is unresolved (2 <= t_N < 3)'
    elif ts <= -2:
        o = 'unresolved, leaning to absorbers (-3 < t_share <= -2)'
    elif tN >= 3:
        o = 'detected absorbers are not the cause: growth persists without absorber sight lines (t_N >= +3) and removing them does not lower the slope (t_share > -2)'
    elif tD >= 2.5 and u <= -1.645:
        o = 'detected absorbers carry less than half of the growth (95%): t_share > -2, t_D >= 2.5 and u <= -1.645'
    else:
        o = 'unresolved: removing absorber sight lines does not lower the slope beyond chance (t_share > -2), but growth without them is not established (t_N < 3) and the half-share bound is not reached'
    out['outcome'] = o
    return out


def bin_table(S):
    """For the 8 equal-count redshift bins of D: share with absorbers and the scatter of each class (descriptive)."""
    z, d, A = S['z'], S['d'], S['A']
    n, K = len(z), CFG['nbins']
    idx = np.argsort(z, kind='stable')
    rows = []
    g = np.random.default_rng(CFG['seed'] + 99)
    for k in range(K):
        q = idx[jsround(k * n / K):jsround((k + 1) * n / K)]
        r = dict(z_lo=float(z[q].min()), z_hi=float(z[q].max()), z_med=float(np.median(z[q])), n=int(len(q)), n_A=int(A[q].sum()),
                 share_A=float(A[q].mean()), mean_absorbers=float(S['nabs'][q].mean()), e2mean=float(np.mean(S['e'][q] ** 2)))
        for name, sel in (('all', np.ones(len(q), bool)), ('N', ~A[q]), ('A', A[q])):
            v = d[q][sel]
            if len(v) >= 8:
                r['sigma_' + name] = float(sig_mad(v))
                r['sigma_' + name + '_err'] = float(sig_mad(v[g.integers(0, len(v), size=(CFG['nboot'], len(v)))]).std(ddof=1))
            else:
                r['sigma_' + name] = None; r['sigma_' + name + '_err'] = None
        rows.append(r)
    return rows


def full_test(rows, src, q, ab, mq=None, nshuf=None, ndraw=None, nshuf_check=None, ndraw_check=None, log=print):
    nshuf = nshuf or CFG['nshuf']; ndraw = ndraw or CFG['ndraw']
    nshuf_check = nshuf_check or CFG['nshuf_check']; ndraw_check = ndraw_check or CFG['ndraw_check']
    tree = cKDTree(unit(src['ra'], src['dec']))
    m = match(q, rows)
    S, cnt = build(m, src, q, ab, CFG['fg_radius'], CFG['min_neigh'], tree)
    cnt['quasars_within_10arcsec'] = m['n_quasars_matched']
    res = dict(counts=cnt)
    log('counts %s' % json.dumps(cnt))
    if cnt['kept'] < CFG['min_D'] or cnt['without_absorber'] < CFG['min_N']:
        res['decision'] = dict(outcome='not tested (fewer than 1,500 sight lines, or fewer than 1,000 without absorbers)')
        return res, S
    A = S['A']
    res['D'] = analyse(S, nshuf, 'D (all DESI sight lines)', log=log)
    res['N'] = analyse(sub(S, ~A), nshuf, 'N (no absorber)', log=log)
    res['share'] = share_test(S, ~A, ndraw, 'share', log=log)
    res['decision'] = decide_173(res['D'], res['N'], res['share'])
    log('decision %s' % json.dumps(res['decision']))
    rep = {}
    rep['A_only'] = analyse(sub(S, A), nshuf_check, 'A (with absorber)', log=log) if A.sum() >= 200 else None
    rep['excess'] = excess_test(S, A, ndraw, 'A against N', log=log)
    rep['bins'] = bin_table(S)
    # (a) spectra good enough to show an absorber: z >= 0.4 and TSNR2_QSO at or above the lower quartile of A's
    if A.sum() >= 20:
        cut = float(np.quantile(S['tsnr2'][A], 0.25))
        good = (S['z'] >= CFG['zmin_search']) & (S['tsnr2'] >= cut)
        Sg = sub(S, good)
        rep['good_spectra'] = dict(tsnr2_cut=cut, n=int(good.sum()), n_A=int(Sg['A'].sum()),
                                   D=analyse(Sg, nshuf_check, 'good spectra, D', log=log),
                                   N=analyse(sub(Sg, ~Sg['A']), nshuf_check, 'good spectra, N', log=log),
                                   share=share_test(Sg, ~Sg['A'], ndraw_check, 'good spectra, share', log=log))
    # (b) strong absorbers only
    strong = A & (S['ewmax'] >= CFG['strong_ew'])
    if strong.sum() >= 50:
        rep['strong'] = dict(n_strong=int(strong.sum()),
                             N=analyse(sub(S, ~strong), nshuf_check, 'no strong absorber', log=log),
                             share=share_test(S, ~strong, ndraw_check, 'strong, share', log=log),
                             excess=excess_test(S, strong, ndraw_check, 'strong against rest', log=log))
    # (c) 1 deg foreground
    S1, c1 = build(m, src, q, ab, CFG['fg_radius_check'], CFG['min_neigh_check'], tree)
    rep['foreground_1deg'] = dict(counts=c1, D=analyse(S1, nshuf_check, '1 deg foreground, D', log=log),
                                  N=analyse(sub(S1, ~S1['A']), nshuf_check, '1 deg foreground, N', log=log),
                                  share=share_test(S1, ~S1['A'], ndraw_check, '1 deg foreground, share', log=log))
    # (d) scatter by number of absorbers
    rep['by_count'] = []
    for name, sel in (('0', S['nabs'] == 0), ('1', S['nabs'] == 1), ('2 or more', S['nabs'] >= 2)):
        if sel.sum() >= 8:
            v = S['d'][sel]
            bs = sig_mad(v[np.random.default_rng(CFG['seed'] + 5).integers(0, len(v), size=(CFG['nboot'], len(v)))])
            rep['by_count'].append(dict(absorbers=name, n=int(sel.sum()), z_median=float(np.median(S['z'][sel])),
                                        sigma=float(sig_mad(v)), err=float(bs.std(ddof=1))))
    # (e) overlap with Entry 171's sample
    if mq is not None:
        mm = match(mq, rows)
        common = np.isin(S['src'], mm['src'])
        zmq = dict(zip(mm['src'].tolist(), mq['z'][mm['q']].tolist()))
        dz = np.array([abs(zmq[s] - zz) for s, zz in zip(S['src'][common].tolist(), S['z'][common].tolist())])
        rep['entry171_overlap'] = dict(in_entry171_match=int(common.sum()), share=float(common.mean()),
                                       redshifts_agree_within_0p05=float((dz < 0.05).mean()) if len(dz) else None)
    rep['z_A_median'] = float(np.median(S['z'][A])) if A.sum() else None
    rep['z_N_median'] = float(np.median(S['z'][~A]))
    res['reported_only'] = rep
    return res, S


# ---------------------------------------------------------------- synthetic skies (no data read)
def synthetic_sky(seed, n_src=60000, n_q=2400, grow=0.0, abs_sigma=0.0):
    """Entry 171's synthetic sky, with a DESI-like quasar list and absorbers.  A quasar has an absorber with probability
    min(0.6, 0.17 z), so high-redshift quasars have more.  'grow' adds Gaussian scatter grow * log10(1+z) rad/m^2 to every
    quasar (a path-length term); 'abs_sigma' adds Gaussian scatter abs_sigma rad/m^2 to absorber sight lines only."""
    g = np.random.default_rng(seed)
    u = g.uniform(np.sin(-90 * D2R), np.sin(49 * D2R), n_src)
    dec = np.arcsin(u) / D2R
    ra = g.uniform(0, 360, n_src)
    l, b = galactic(ra, dec)
    fg = 25 / (np.abs(np.sin(b * D2R)) + 0.3) * np.sin(l * D2R + 0.7) + 8 * np.cos(3 * l * D2R) * np.sin(2 * b * D2R)
    e = g.lognormal(np.log(1.2), 0.5, n_src)
    rm = fg + g.normal(0, 7.5, n_src) + g.normal(0, 1, n_src) * e
    hi = np.nonzero(np.abs(b) >= 20)[0]
    qsrc = g.choice(hi, n_q, replace=False)
    zq = np.clip(g.gamma(4.0, 0.35, n_q), 0.05, 6)
    hasA = g.uniform(size=n_q) < np.minimum(0.6, 0.17 * zq)
    rm[qsrc] += g.normal(0, 1, n_q) * grow * np.log10(1 + zq) + g.normal(0, 1, n_q) * abs_sigma * hasA
    rows = dict(row_ra=ra, row_dec=dec, row_src=np.arange(n_src))
    src = dict(ra=ra, dec=dec, l=l, b=b, rm=rm, e=e, simple=np.ones(n_src, bool), nrow=np.ones(n_src, int))
    off = g.uniform(0, 4, n_q) / 3600
    pa = g.uniform(0, 2 * np.pi, n_q)
    q = dict(tid=np.arange(n_q, dtype=np.int64) + 1000, z=zq, ra=ra[qsrc] + off * np.sin(pa) / np.cos(dec[qsrc] * D2R),
             dec=dec[qsrc] + off * np.cos(pa), tsnr2=g.lognormal(3, 0.6, n_q))
    at = q['tid'][hasA]
    two = g.uniform(size=len(at)) < 0.3
    tid = np.r_[at, at[two]]
    ab = dict(tid=tid, zq=np.zeros(len(tid)), zabs=np.zeros(len(tid)), ew=g.lognormal(0, 0.6, len(tid)))
    return rows, src, q, ab


def synthetic():
    out = {}
    cases = (('flat', 0.0, 0.0, 11), ('growth_all_sightlines', 22.0, 0.0, 12), ('growth_from_absorbers', 0.0, 9.0, 13),
             ('half_and_half', 15.0, 6.5, 14))
    for name, grow, asig, seed in cases:
        rows, src, q, ab = synthetic_sky(seed, grow=grow, abs_sigma=asig)
        tree = cKDTree(unit(src['ra'], src['dec']))
        m = match(q, rows)
        S, c = build(m, src, q, ab, CFG['fg_radius'], CFG['min_neigh'], tree)
        D = analyse(S, 200, name + ' D')
        N = analyse(sub(S, ~S['A']), 200, name + ' N')
        sh = share_test(S, ~S['A'], 300, name + ' share')
        ex = excess_test(S, S['A'], 300, name + ' excess')
        dec = decide_173(D, N, sh)
        out[name] = dict(counts=c, b_D=round(D['b'], 2), t_D=round(D['t'], 2), b_N=round(N['b'], 2), t_N=round(N['t'], 2),
                         t_share=round(sh['t_share'], 2), drop=round(sh['drop'], 2), random_sd=round(sh['random_sd'], 2),
                         excess_quadrature=round(ex['excess_quadrature'], 2), t_excess=round(ex['t_excess'], 2), outcome=dec['outcome'])
    # calibration of t_share and t_N over repeated skies
    for name, grow, asig in (('growth_all_sightlines', 22.0, 0.0), ('growth_from_absorbers', 0.0, 9.0)):
        ts, tn = [], []
        for seed in range(100, 124):
            rows, src, q, ab = synthetic_sky(seed, grow=grow, abs_sigma=asig)
            m = match(q, rows)
            S, c = build(m, src, q, ab, CFG['fg_radius'], CFG['min_neigh'])
            sh = share_test(S, ~S['A'], 200, '%s sky %d' % (name, seed), log=lambda s: None)
            N = analyse(sub(S, ~S['A']), 100, 'N', log=lambda s: None)
            ts.append(sh['t_share']); tn.append(N['t'])
        out['repeat_' + name] = dict(skies=len(ts), t_share_mean=round(float(np.mean(ts)), 2), t_share_sd=round(float(np.std(ts, ddof=1)), 2),
                                     t_share_min=round(float(np.min(ts)), 2), t_share_max=round(float(np.max(ts)), 2),
                                     t_N_mean=round(float(np.mean(tn)), 2), t_N_sd=round(float(np.std(tn, ddof=1)), 2),
                                     share_le_minus3=int(np.sum(np.array(ts) <= -3)), share_le_minus2=int(np.sum(np.array(ts) <= -2)))
        print(name, out['repeat_' + name], flush=True)
    print(json.dumps(out, indent=1))
    json.dump(out, open(os.path.join(HERE, 'synthetic.json'), 'w'), indent=1)


def run():
    t0 = time.time()
    log = lambda s: print(s, flush=True)
    rows, src, sinfo = load_spice()
    log('spice %s' % sinfo)
    q, qinfo = load_desi()
    log('desi quasars %s' % qinfo)
    ab, ainfo = load_absorbers()
    log('absorbers %s' % ainfo)
    mq = load_milliquas() if os.path.exists(MQZIP) else None
    res = dict(generated=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), config=CFG, spice=sinfo, desi_quasars=qinfo, absorbers=ainfo)
    r, S = full_test(rows, src, q, ab, mq, log=log)
    res.update(r)
    res['seconds'] = round(time.time() - t0)
    js = json.dumps(res, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
    open(os.path.join(HERE, 'results.json'), 'w').write(js)
    os.makedirs(WORK, exist_ok=True)
    np.savez(os.path.join(WORK, 'sample.npz'), **S)
    log('decision %s' % json.dumps(res['decision'], default=float))
    log('written results.json %d s' % res['seconds'])


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else ''
    {'extract': extract, 'synthetic': synthetic, 'run': run}.get(cmd, lambda: print(__doc__))()
