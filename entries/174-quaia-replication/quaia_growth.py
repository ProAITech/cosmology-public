"""Entry 174: does the growth of Faraday scatter with redshift (Entry 171) appear on sight lines no earlier entry used?

The rule was frozen in the Workbench (section 174, "the rule, frozen first") before any Quaia quasar was matched to a rotation measure.
Entry 171's statistic, unchanged, on SPICE-RACS DR2 sources behind Quaia quasars, leaving out every source that carries a
Milliquas object (Entry 171, all classes) or a DESI DR1 quasar (Entry 173).

Data (raw files are git-ignored under data/):
  data/spice-racs-dr2/spice-racs.dr2.fits.gz   SPICE-RACS DR2 (Thomson et al. 2026, arXiv:2605.16917), CSIRO DAP csiro:64891,
      https://data.csiro.au/collection/csiro:64891, 4,968,301,262 bytes, sha256 a1a2f3c1d12a8a8b3aa8526d74b8e0115054f610e2e4186108f0f068005d0053.
  data/quaia/quaia_G20.5.fits                  Quaia, the Gaia-unWISE quasar catalogue (Storey-Fisher et al. 2024, arXiv:2306.17749),
      Zenodo doi:10.5281/zenodo.10403370, https://zenodo.org/api/records/10403370/files/quaia_G20.5.fits/content, 171,020,160 bytes,
      md5 98659ac4bd8a09da2c4ce653690d53df, sha256 9e357a760f2ff2db9417a6791fb71dc421ee30324301f22db898038c9e648584.
      One BINTABLE, 1,295,502 rows.  Columns used: ra, dec, redshift_quaia, redshift_quaia_err, phot_g_mean_mag.
  data/milliquas-v8/milliquas.zip              Milliquas v8 (Flesch 2023), https://quasars.org/milliquas.zip, 32,701,849 bytes,
      sha256 9128d6aa91354d0ee1c83dde3ab8764a20d542aa82327e4a969a4357ba9ac984.  Used to exclude Entry 171's sight lines.
  data/desi-dr1-qso/QSO_cat_iron_cumulative_v0.fits   DESI DR1 quasar catalogue, 816,203,520 bytes,
      sha256 a65604303cc29af8865db9193a23934b2540591b16c5afd2a4bea2d61464aef4.  Used to exclude Entry 173's sight lines.

Rule, in brief (section 174 has the full text):
  Sample     Quaia quasars with redshift_quaia > 0, matched as in Entry 171 (nearest qualifying SPICE-RACS component within 10 arcsec,
             one quasar per source, |b| >= 20 deg, at least 10 other sources within 3 deg), after removing every source that has a
             Milliquas object of class Q, A, B, K or N with z > 0, or a DESI DR1 quasar with Z > 0, within 10 arcsec of a component.
  Statistic  Entry 171's slope b and shuffle calibration (1,000 shuffles): t.
  Dilution   beta = least-squares slope of log10(1+z_Milliquas) on log10(1+z_Quaia) over Entry 171's kept sight lines that Quaia
             also has, clipped to 0.5-1.  Entry 171's net slope seen through Quaia redshifts is 4.63 beta.
  Decision   decide_174() below; the frozen text is section 174, item 5.

Run from the repository root:
  .venv/bin/python -I entries/174-quaia-replication/quaia_growth.py synthetic   # synthetic skies; reads no data
  .venv/bin/python -I entries/174-quaia-replication/quaia_growth.py run         # the test; writes results.json
The functions down to sub() are Entry 171's, copied from Entry 173's script so that this folder stands alone.
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

CFG = dict(match_arcsec=10.0, fg_radius=3.0, min_neigh=10, bcut=20.0, nbins=8, nboot=500, nshuf=1000, nshuf_check=300,
           seed=20261121, bands=[20, 35, 50, 65, 90], fg_radius_check=1.0, min_neigh_check=5, chance_shift_arcmin=5.0,
           classes=['Q', 'A', 'B', 'K', 'N'], min_sample=1000, growth_171=5.083 - 0.451, t_171=3.81, beta_clip=[0.5, 1.0])

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


# ---------------------------------------------------------------- new in Entry 174
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


def build(m, src, q, radius, min_neigh, tree, exclude=None):
    """Entry 171's sample from a match; 'exclude' is a boolean array over sources that may not be used."""
    idx, qq, sep = m['src'], m['q'], m['sep']
    cnt = dict(matched=int(len(idx)))
    if exclude is not None:
        keep = ~exclude[idx]
        cnt['dropped_used_before'] = int((~keep).sum())
        idx, qq, sep = idx[keep], qq[keep], sep[keep]
    lat = np.abs(src['b'][idx]) >= CFG['bcut']
    med, nn = foreground(src, idx[lat], radius, min_neigh, tree)
    okn = np.isfinite(med)
    k = idx[lat][okn]
    qk = qq[lat][okn]
    S = dict(z=q['z'][qk], d=src['rm'][k] - med[okn], e=src['e'][k], b=src['b'][k], l=src['l'][k], dec=src['dec'][k],
             sep=sep[lat][okn], src=k, qi=qk)
    for extra in ('zerr', 'g'):
        if extra in q:
            S[extra] = q[extra][qk]
    cnt.update(dropped_latitude=int((~lat).sum()), dropped_neighbours=int((~okn).sum()), kept=int(len(k)))
    return S, cnt


def beta_of(z_spec, z_phot):
    """Least-squares slope of log10(1+z_spec) on log10(1+z_phot): how much a slope against z_phot is diluted."""
    xs, xp = np.log10(1 + z_spec), np.log10(1 + z_phot)
    c = np.cov(xp, xs)
    return float(c[0, 1] / c[0, 0])


def decide_174(P, beta):
    """Section 174, item 5.  P: output of analyse() for the new sight lines."""
    t, net, sd = P['t'], P['net'], P['null']['sd']
    ref = CFG['growth_171'] * beta
    out = dict(t=t, net_slope=net, null_sd=sd, upper95_net=net + 1.645 * sd, beta=beta, entry171_net_through_quaia=ref,
               stouffer_with_entry171=(CFG['t_171'] + t) / np.sqrt(2))
    if t >= 3:
        o = 'growth confirmed on independent sight lines (t >= +3)'
    elif t >= 2:
        o = 'growth supported on independent sight lines (2 <= t < 3)'
    elif t <= -3:
        o = 'the scatter falls with redshift on the new sight lines (t <= -3)'
    elif net + 1.645 * sd < ref:
        o = "not replicated: Entry 171's slope, seen through Quaia redshifts, is excluded at 95%"
    else:
        o = "not replicated and not excluded (t < 2, but the 95% limit lies above Entry 171's slope): unresolved"
    out['outcome'] = o
    return out


def full_test(rows, src, quaia, mq, desi, nshuf=None, nshuf_check=None, log=print):
    nshuf = nshuf or CFG['nshuf']; nshuf_check = nshuf_check or CFG['nshuf_check']
    tree = cKDTree(unit(src['ra'], src['dec']))
    nsrc = len(src['rm'])
    used = np.zeros(nsrc, bool)
    m_mq = match(mq, rows); used[m_mq['src']] = True
    n_mq = int(used.sum())
    m_de = match(desi, rows); used[m_de['src']] = True
    res = dict(excluded=dict(sources_with_milliquas=n_mq, sources_with_desi=int(len(m_de['src'])), sources_with_either=int(used.sum())))
    mQ = match(quaia, rows)
    S, cnt = build(mQ, src, quaia, CFG['fg_radius'], CFG['min_neigh'], tree, exclude=used)
    cnt['quaia_within_10arcsec'] = mQ['n_quasars_matched']
    res['counts'] = cnt
    log('excluded %s' % json.dumps(res['excluded'])); log('counts %s' % json.dumps(cnt))
    # dilution: Entry 171's kept sight lines (class Q) that Quaia also has
    isQ = mq['cls'] == 'Q'
    mq_Q = {k: v[isQ] for k, v in mq.items()}
    mM = match(mq_Q, rows)
    M, cM = build(mM, src, mq_Q, CFG['fg_radius'], CFG['min_neigh'], tree)
    zq_by_src = dict(zip(mQ['src'].tolist(), quaia['z'][mQ['q']].tolist()))
    common = np.array([s in zq_by_src for s in M['src'].tolist()])
    zq_common = np.array([zq_by_src[s] for s in M['src'][common].tolist()])
    beta_raw = beta_of(M['z'][common], zq_common)
    beta = float(np.clip(beta_raw, *CFG['beta_clip']))
    dz = (zq_common - M['z'][common]) / (1 + M['z'][common])
    res['dilution'] = dict(entry171_kept=int(len(M['z'])), common=int(common.sum()), beta_raw=beta_raw, beta=beta,
                           share_within_0p01=float((np.abs(dz) < 0.01).mean()), share_within_0p2=float((np.abs(dz) < 0.2).mean()),
                           nmad=float(1.4826 * np.median(np.abs(dz - np.median(dz)))))
    log('dilution %s' % json.dumps(res['dilution']))
    if cnt['kept'] < CFG['min_sample']:
        res['decision'] = dict(outcome='not tested (fewer than 1,000 new sight lines)')
        return res, S
    res['primary'] = analyse(S, nshuf, 'new Quaia sight lines', log=log)
    res['decision'] = decide_174(res['primary'], beta)
    log('decision %s' % json.dumps(res['decision'], default=float))
    rep = {}
    rep['north_galactic'] = analyse(sub(S, S['b'] > 0), nshuf_check, 'north Galactic', log=log)
    rep['south_galactic'] = analyse(sub(S, S['b'] < 0), nshuf_check, 'south Galactic', log=log)
    for name, sel in (('dec_below_minus20', S['dec'] < -20), ('dec_above_minus20', S['dec'] >= -20),
                      ('smaller_redshift_error_half', S['zerr'] <= np.median(S['zerr'])), ('G_below_20', S['g'] < 20.0)):
        rep[name] = analyse(sub(S, sel), nshuf_check, name, log=log) if sel.sum() >= 400 else dict(n=int(sel.sum()), note='fewer than 400')
    S1, c1 = build(mQ, src, quaia, CFG['fg_radius_check'], CFG['min_neigh_check'], tree, exclude=used)
    rep['foreground_1deg'] = analyse(S1, nshuf_check, '1 deg foreground', log=log); rep['foreground_1deg']['counts'] = c1
    shifted = dict(quaia, dec=quaia['dec'] + CFG['chance_shift_arcmin'] / 60.0)
    mC = match(shifted, rows)
    rep['chance_matches_shifted_5arcmin'] = dict(all=int(len(mC['src'])), on_unused_sources=int((~used[mC['src']]).sum()))
    # Entry 171's sight lines that Quaia has: Milliquas redshifts against Quaia redshifts on the same sources
    Mc = sub(M, common)
    rep['entry171_common_milliquas_z'] = analyse(Mc, nshuf_check, "Entry 171's sight lines in Quaia, Milliquas z", log=log)
    Mq = dict(Mc, z=zq_common)
    rep['entry171_common_quaia_z'] = analyse(Mq, nshuf_check, "Entry 171's sight lines in Quaia, Quaia z", log=log)
    # everything together: Entry 171's sample, the DESI-only sight lines, and the new ones, each with its own redshift
    Dd, cD = build(m_de, src, desi, CFG['fg_radius'], CFG['min_neigh'], tree)
    donly = ~np.isin(Dd['src'], M['src'])
    U = {k: np.r_[M[k], Dd[k][donly], S[k]] for k in ('z', 'd', 'e', 'b', 'l')}
    rep['union_all'] = analyse(U, nshuf, 'union: Entry 171 + DESI-only + new', log=log)
    rep['union_all']['parts'] = dict(entry171=int(len(M['z'])), desi_only=int(donly.sum()), new=int(len(S['z'])))
    rep['z_median_new'] = float(np.median(S['z'])); rep['g_median_new'] = float(np.median(S['g']))
    rep['zerr_median_new'] = float(np.median(S['zerr']))
    res['reported_only'] = rep
    return res, S


# ---------------------------------------------------------------- synthetic skies (no data read)
def synthetic_sky(seed, n_src=60000, n_q=6000, grow=0.0, used_share=0.4):
    """Entry 171's synthetic sky.  Every quasar has a true redshift and a Quaia-like one: scatter 0.05 (1+z), and 6% thrown
    anywhere in 0.1-4.  A share 'used_share' of the quasars are also in a 'Milliquas' list with the true redshift."""
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
    zt = np.clip(g.gamma(4.0, 0.35, n_q), 0.05, 6)
    rm[qsrc] += g.normal(0, 1, n_q) * grow * np.log10(1 + zt)
    zp = zt + 0.05 * (1 + zt) * g.normal(0, 1, n_q)
    bad = g.uniform(size=n_q) < 0.06
    zp[bad] = g.uniform(0.1, 4, bad.sum())
    zp = np.clip(zp, 0.02, 6)
    rows = dict(row_ra=ra, row_dec=dec, row_src=np.arange(n_src))
    src = dict(ra=ra, dec=dec, l=l, b=b, rm=rm, e=e, simple=np.ones(n_src, bool), nrow=np.ones(n_src, int))
    off = g.uniform(0, 4, n_q) / 3600
    pa = g.uniform(0, 2 * np.pi, n_q)
    qra, qdec = ra[qsrc] + off * np.sin(pa) / np.cos(dec[qsrc] * D2R), dec[qsrc] + off * np.cos(pa)
    quaia = dict(ra=qra, dec=qdec, z=zp, zerr=0.05 * (1 + zt) * g.uniform(0.5, 1.5, n_q), g=g.uniform(17, 20.5, n_q))
    inmq = g.uniform(size=n_q) < used_share
    mq = dict(ra=qra[inmq], dec=qdec[inmq], z=zt[inmq], cls=np.full(int(inmq.sum()), 'Q'))
    desi = dict(ra=np.array([10.0]), dec=np.array([80.0]), z=np.array([1.0]))
    return rows, src, quaia, mq, desi, dict(n_new_true=int((~inmq).sum()), ztrue=zt, inmq=inmq)


def synthetic():
    out = {}
    for name, grow, seed in (('flat_1', 0.0, 21), ('flat_2', 0.0, 22), ('growth_10', 10.0, 23), ('growth_22', 22.0, 24)):
        rows, src, quaia, mq, desi, truth = synthetic_sky(seed, grow=grow)
        tree = cKDTree(unit(src['ra'], src['dec']))
        used = np.zeros(len(src['rm']), bool)
        used[match(mq, rows)['src']] = True
        mQ = match(quaia, rows)
        S, c = build(mQ, src, quaia, CFG['fg_radius'], CFG['min_neigh'], tree, exclude=used)
        P = analyse(S, 200, name)
        St = dict(S, z=truth['ztrue'][S['qi']])
        T = analyse(St, 200, name + ' (true redshifts)')
        beta = beta_of(truth['ztrue'][S['qi']], S['z'])
        out[name] = dict(counts=c, new_true=truth['n_new_true'], t=round(P['t'], 2), net=round(P['net'], 2), t_true_z=round(T['t'], 2),
                         net_true_z=round(T['net'], 2), beta=round(beta, 3), net_ratio=round(P['net'] / T['net'], 3) if grow else None,
                         outcome=decide_174(P, float(np.clip(beta, *CFG['beta_clip'])))['outcome'])
    print(json.dumps(out, indent=1))
    json.dump(out, open(os.path.join(HERE, 'synthetic.json'), 'w'), indent=1)


def run():
    t0 = time.time()
    log = lambda s: print(s, flush=True)
    rows, src, sinfo = load_spice()
    log('spice %s' % sinfo)
    quaia, qinfo = load_quaia()
    log('quaia %s' % qinfo)
    mq = load_milliquas()
    desi = load_desi()
    log('milliquas (Q, A, B, K, N; z > 0) %d, desi quasars %d' % (len(mq['z']), len(desi['z'])))
    res = dict(generated=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), config=CFG, spice=sinfo, quaia=qinfo,
               milliquas_objects=int(len(mq['z'])), desi_quasars=int(len(desi['z'])))
    r, S = full_test(rows, src, quaia, mq, desi, log=log)
    res.update(r)
    res['seconds'] = round(time.time() - t0)
    open(os.path.join(HERE, 'results.json'), 'w').write(json.dumps(res, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o)))
    os.makedirs(WORK, exist_ok=True)
    np.savez(os.path.join(WORK, 'sample.npz'), **S)
    log('decision %s' % json.dumps(res['decision'], default=float))
    log('written results.json %d s' % res['seconds'])


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else ''
    {'extract': extract, 'synthetic': synthetic, 'run': run}.get(cmd, lambda: print(__doc__))()
