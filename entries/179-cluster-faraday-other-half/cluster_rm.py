"""Entry 176: how far from galaxy clusters does Faraday rotation show their magnetised plasma?

The rule was frozen in the Workbench (section 176, "the rule, frozen first") before any eRASS1 cluster was matched to a rotation
measure.  SPICE-RACS DR2 sources are sorted by their projected distance from the nearest eRASS1 cluster in units of its R500, and
the scatter of the foreground-subtracted rotation measure in each range is set against the same clusters rotated in Galactic
longitude.

Data (raw files are git-ignored under data/):
  data/spice-racs-dr2/spice-racs.dr2.fits.gz   SPICE-RACS DR2 (Thomson et al. 2026, arXiv:2605.16917), CSIRO DAP csiro:64891,
      4,968,301,262 bytes, sha256 a1a2f3c1d12a8a8b3aa8526d74b8e0115054f610e2e4186108f0f068005d0053.
  data/erass1-clusters/erass1cl_primary_v3.2.fits   the eRASS1 catalogue of galaxy clusters and groups (Bulbul et al. 2024,
      arXiv:2402.08452), https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/BulbulE_DR1/erass1cl_primary_v3.2.fits.tgz,
      23,764,095 bytes, sha256 f9e9e8f750a2831ac68dad1fb24a028b024c00d3b5879f7a5760af81fc9a40ca; unpacked 50,636,160 bytes.
      HDU 1, 12,247 rows.  Columns used: RA, DEC, BEST_Z, PCONT, R500 (kpc), M500.

Rule, in brief (section 176 has the full text):
  Sources    SPICE-RACS DR2 with Entry 171's quality cuts and merge; |b| >= 20 deg; at least 10 other sources within 3 deg;
             dRM = RM - median RM of the other sources within 3 deg.
  Clusters   PCONT < 0.5, BEST_Z >= 0.05, R500 > 0; theta500 = R500 / D_A(z), flat LCDM with H0 = 67.74, Om = 0.3089.
  Distance   x = the smallest theta / theta500 over all clusters.
  Ranges     A: x < 1.  B: 1 <= x < 2.  C: 2 <= x < 4.
  Statistic  V = (1.4826 MAD of dRM)^2 - mean e_RM^2 in each range.
  Null       the cluster list rotated in Galactic longitude by 500 random angles between 10 and 350 deg; t = (V - mean) / sd.
  Decision   decide_176() below; the frozen text is section 176, item 6.

Run from the repository root:
  .venv/bin/python -I entries/176-cluster-faraday-profile/cluster_rm.py synthetic   # synthetic sky; reads no data
  .venv/bin/python -I entries/176-cluster-faraday-profile/cluster_rm.py run         # the test; writes results.json
The SPICE-RACS reader and merge are Entry 171's, copied so that this folder stands alone.
"""
import gzip, json, os, sys, time, zipfile
import numpy as np
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
SPICE = os.path.join(ROOT, 'data', 'spice-racs-dr2', 'spice-racs.dr2.fits.gz')
CLUSTERS = os.path.join(ROOT, 'data', 'erass1-clusters', 'erass1cl_primary_v3.2.fits')
WORK = os.path.join(HERE, 'work')
D2R = np.pi / 180

CFG = dict(fg_radius=3.0, min_neigh=10, bcut=20.0, seed=20261121, nrot=500, nperm=500, pcont_max=0.5, zmin=0.05, H0=67.74, Om=0.3089,
           ranges=dict(A=[0.0, 1.0], B=[1.0, 2.0], C=[2.0, 4.0]), xmax=8.0, min_sources=300,
           profile_edges=[0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0], fg_radius_check=1.0, min_neigh_check=5, nrot_check=200)

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



# ---------------------------------------------------------------- new in Entry 176
def ang_diam_Mpc(z, n=400):
    z = np.atleast_1d(np.asarray(z, float))
    zz = (np.arange(n) + 0.5)[None, :] * z[:, None] / n
    dc = 299792.458 / CFG['H0'] * (1 / np.sqrt(CFG['Om'] * (1 + zz) ** 3 + 1 - CFG['Om'])).sum(1) * z / n
    return dc / (1 + z)


def load_clusters():
    from astropy.io import fits
    with fits.open(CLUSTERS, memmap=True) as h:
        d = h[1].data
        ra, dec, z = np.asarray(d['RA'], float), np.asarray(d['DEC'], float), np.asarray(d['BEST_Z'], float)
        pc, r500, m500 = np.asarray(d['PCONT'], float), np.asarray(d['R500'], float), np.asarray(d['M500'], float)
    ok = np.isfinite(ra) & np.isfinite(dec) & np.isfinite(z) & np.isfinite(r500) & (pc < CFG['pcont_max']) & (z >= CFG['zmin']) & (r500 > 0)
    th = r500[ok] / 1000.0 / ang_diam_Mpc(z[ok]) / D2R       # theta500 in degrees
    l, b = galactic(ra[ok], dec[ok])
    info = dict(rows=int(len(ra)), kept=int(ok.sum()), z_median=float(np.median(z[ok])), theta500_arcmin_median=float(np.median(th) * 60),
                theta500_arcmin_max=float(th.max() * 60), M500_median_1e13=float(np.nanmedian(m500[ok])))
    return dict(l=l, b=b, th=th, z=z[ok], m500=m500[ok]), info


def gal_unit(l, b):
    return unit(l, b)


def scaled_distance(src_tree, nsrc, cl, shift=0.0):
    """For every source: the smallest theta / theta500 over the clusters (inf if none within xmax theta500), and that cluster."""
    x = np.full(nsrc, np.inf); which = np.full(nsrc, -1)
    v = gal_unit(np.mod(cl['l'] + shift, 360.0), cl['b'])
    lists = src_tree.query_ball_point(v, chord(np.minimum(CFG['xmax'] * cl['th'], 60.0)))
    data = src_tree.data
    for k, L in enumerate(lists):
        if not L:
            continue
        L = np.asarray(L)
        sep = 2 * np.arcsin(np.clip(np.linalg.norm(data[L] - v[k], axis=1) / 2, 0, 1)) / D2R
        xs = sep / cl['th'][k]
        better = xs < x[L]
        x[L[better]] = xs[better]; which[L[better]] = k
    return x, which


def V_of(d, e):
    m = np.median(d)
    return float((1.4826 * np.median(np.abs(d - m))) ** 2 - np.mean(e ** 2))


def stats_in(x, d, e, edges_list):
    out = []
    for lo, hi in edges_list:
        m = (x >= lo) & (x < hi)
        out.append((V_of(d[m], e[m]) if m.sum() >= 20 else np.nan, int(m.sum())))
    return out


def profile_test(S, cl, edges_list, nrot, seed, log=print, label=''):
    """V in each range for the real clusters and for nrot rotations of the cluster list in Galactic longitude."""
    tree = cKDTree(gal_unit(S['l'], S['b']))
    n = len(S['d'])
    x, which = scaled_distance(tree, n, cl)
    real = stats_in(x, S['d'], S['e'], edges_list)
    g = np.random.default_rng(seed)
    nulls = np.empty((nrot, len(edges_list))); counts = np.empty((nrot, len(edges_list)))
    for i in range(nrot):
        xr, _ = scaled_distance(tree, n, cl, shift=g.uniform(10.0, 350.0))
        r = stats_in(xr, S['d'], S['e'], edges_list)
        nulls[i] = [a for a, _ in r]; counts[i] = [b for _, b in r]
    res = []
    for j, (lo, hi) in enumerate(edges_list):
        v = nulls[:, j]; v = v[np.isfinite(v)]
        mu, sd = float(v.mean()), float(v.std(ddof=1))
        V, cnt = real[j]
        ex = V - mu
        res.append(dict(lo=lo, hi=hi, n=cnt, V=V, sigma=float(np.sqrt(max(V, 0))), null_mean=mu, null_sd=sd, null_n_mean=float(counts[:, j].mean()),
                        t=float(ex / sd), excess_variance=float(ex), excess_quadrature=float(np.sign(ex) * np.sqrt(abs(ex))),
                        upper95_excess_variance=float(ex + 1.645 * sd), p_upper=float((1 + (v >= V).sum()) / (1 + len(v)))))
        log('%s x %.1f-%.1f: n %d (rotated %.0f), V %.1f, null %.1f +- %.1f, t %+.2f, excess %+.1f rad2/m4' % (
            label, lo, hi, cnt, counts[:, j].mean(), V, mu, sd, ex / sd, ex))
    return res, x, which


def permutation_test(S, x, edges_list, nperm, seed):
    """Second null: dRM and its error shuffled together among sources within latitude band x e_RM tercile, positions fixed."""
    n = len(S['d'])
    es = np.sort(S['e']); t1, t2 = es[n // 3], es[2 * n // 3]
    ab = np.abs(S['b']); band = (ab >= 35).astype(int) + (ab >= 50) + (ab >= 65)
    cell = band * 3 + np.where(S['e'] < t1, 0, np.where(S['e'] < t2, 1, 2))
    groups = [np.nonzero(cell == c)[0] for c in np.unique(cell)]
    real = np.array([a for a, _ in stats_in(x, S['d'], S['e'], edges_list)])
    g = np.random.default_rng(seed)
    nulls = np.empty((nperm, len(edges_list)))
    for i in range(nperm):
        idx = np.arange(n)
        for gr in groups:
            idx[gr] = gr[g.permutation(len(gr))]
        nulls[i] = [a for a, _ in stats_in(x, S['d'][idx], S['e'][idx], edges_list)]
    return [dict(lo=lo, hi=hi, null_mean=float(nulls[:, j].mean()), null_sd=float(nulls[:, j].std(ddof=1)),
                 t=float((real[j] - nulls[:, j].mean()) / nulls[:, j].std(ddof=1))) for j, (lo, hi) in enumerate(edges_list)]


def decide_176(R):
    """Section 176, item 6.  R: dict of range name -> result of profile_test for A, B, C."""
    out = {}
    tA = R['A']['t']
    out['check_A'] = ('inside R500 the scatter is raised (t >= +3): the method sees cluster plasma' if tA >= 3 else
                      'inside R500 no excess is established (t < 3): the ranges beyond are reported, and no claim is made from them')
    for k, name in (('B', '1 to 2 R500'), ('C', '2 to 4 R500')):
        t = R[k]['t']
        if t >= 3:
            o = 'excess scatter at %s (t >= +3)' % name
        elif t <= -3:
            o = 'less scatter at %s than away from clusters (t <= -3)' % name
        elif abs(t) < 2:
            o = 'no excess at %s (|t| < 2)' % name
        else:
            o = 'unresolved at %s (2 <= |t| < 3)' % name
        out[k] = o
    out['claim'] = ('magnetised plasma is seen beyond R500' if (tA >= 3 and R['B']['t'] >= 3) else 'no claim beyond R500')
    return out


def sample(src, radius, min_neigh):
    idx = np.nonzero(np.abs(src['b']) >= CFG['bcut'])[0]
    med, nn = foreground(src, idx, radius, min_neigh)
    ok = np.isfinite(med)
    k = idx[ok]
    return dict(l=src['l'][k], b=src['b'][k], d=src['rm'][k] - med[ok], e=src['e'][k], simple=src['simple'][k]), dict(
        high_latitude=int(len(idx)), kept=int(len(k)))


def full_test(src, cl, log=print, nrot=None):
    nrot = nrot or CFG['nrot']
    S, cnt = sample(src, CFG['fg_radius'], CFG['min_neigh'])
    res = dict(counts=cnt)
    log('sources %s; field scatter V %.1f rad2/m4' % (json.dumps(cnt), V_of(S['d'], S['e'])))
    names = list(CFG['ranges'])
    rr, x, which = profile_test(S, cl, [CFG['ranges'][k] for k in names], nrot, CFG['seed'], log=log, label='primary')
    R = dict(zip(names, rr))
    res['ranges'] = R
    if min(R[k]['n'] for k in names) < CFG['min_sources']:
        res['decision'] = dict(outcome='not tested (a range holds fewer than 300 sources)')
        return res
    res['decision'] = decide_176(R)
    log('decision %s' % json.dumps(res['decision']))
    rep = {}
    rep['permutation_null'] = dict(zip(names, permutation_test(S, x, [CFG['ranges'][k] for k in names], CFG['nperm'], CFG['seed'] + 3)))
    log('permutation null t: %s' % {k: round(v['t'], 2) for k, v in rep['permutation_null'].items()})
    edges = CFG['profile_edges']
    fine, _, _ = profile_test(S, cl, [[edges[i], edges[i + 1]] for i in range(len(edges) - 1)], CFG['nrot_check'], CFG['seed'] + 11, log=log, label='profile')
    rep['profile'] = fine
    for name, sel in (('north_galactic', S['b'] > 0), ('south_galactic', S['b'] < 0), ('faraday_simple', S['simple'])):
        r, _, _ = profile_test({k: v[sel] for k, v in S.items()}, cl, [CFG['ranges'][k] for k in names], CFG['nrot_check'], CFG['seed'] + 21, log=log, label=name)
        rep[name] = dict(zip(names, r))
    med_m, med_z = np.nanmedian(cl['m500']), np.median(cl['z'])
    for name, sel in (('clusters_more_massive_half', cl['m500'] >= med_m), ('clusters_less_massive_half', cl['m500'] < med_m),
                      ('clusters_z_below_median', cl['z'] < med_z), ('clusters_z_above_median', cl['z'] >= med_z)):
        r, _, _ = profile_test(S, {k: v[sel] for k, v in cl.items()}, [CFG['ranges'][k] for k in names], CFG['nrot_check'], CFG['seed'] + 31, log=log, label=name)
        rep[name] = dict(zip(names, r))
    S1, c1 = sample(src, CFG['fg_radius_check'], CFG['min_neigh_check'])
    r, _, _ = profile_test(S1, cl, [CFG['ranges'][k] for k in names], CFG['nrot_check'], CFG['seed'] + 41, log=log, label='1 deg foreground')
    rep['foreground_1deg'] = dict(zip(names, r)); rep['foreground_1deg']['counts'] = c1
    # noise and counts by range, for the record
    rep['by_range'] = {}
    for k in names:
        lo, hi = CFG['ranges'][k]; m = (x >= lo) & (x < hi)
        rep['by_range'][k] = dict(n=int(m.sum()), mean_e2=float(np.mean(S['e'][m] ** 2)), median_e=float(np.median(S['e'][m])),
                                  raw_sigma=float(1.4826 * np.median(np.abs(S['d'][m] - np.median(S['d'][m])))),
                                  rms_clipped=float(np.std(S['d'][m][np.abs(S['d'][m]) < 100])))
    far = ~np.isfinite(x)
    rep['field'] = dict(n=int(far.sum()), V=V_of(S['d'][far], S['e'][far]), mean_e2=float(np.mean(S['e'][far] ** 2)))
    res['reported_only'] = rep
    return res


# ---------------------------------------------------------------- synthetic sky (no data read)
def synthetic():
    out = {}
    for name, amp, seed in (('no_cluster_effect', 0.0, 31), ('excess_inside_R500_only', 1.0, 32), ('excess_to_2_R500', 2.0, 33)):
        g = np.random.default_rng(seed)
        n_src = 120000
        dec = np.arcsin(g.uniform(np.sin(-90 * D2R), np.sin(49 * D2R), n_src)) / D2R; ra = g.uniform(0, 360, n_src)
        l, b = galactic(ra, dec)
        fg = 25 / (np.abs(np.sin(b * D2R)) + 0.3) * np.sin(l * D2R + 0.7) + 8 * np.cos(3 * l * D2R) * np.sin(2 * b * D2R)
        e = g.lognormal(np.log(1.5), 0.5, n_src)
        rm = fg + g.normal(0, 8, n_src) + g.normal(0, 1, n_src) * e
        ncl = 9000
        cl_l = g.uniform(180, 360, ncl); cl_b = np.arcsin(g.uniform(-1, 1, ncl)) / D2R
        keep = np.abs(cl_b) > 20
        cl = dict(l=cl_l[keep], b=cl_b[keep], th=g.lognormal(np.log(5 / 60), 0.4, int(keep.sum())), z=g.uniform(0.05, 0.6, int(keep.sum())),
                  m500=g.lognormal(3, 0.5, int(keep.sum())))
        tree = cKDTree(gal_unit(l, b))
        x, _ = scaled_distance(tree, n_src, cl)
        if amp > 0:
            add = np.where(x < 1, 20.0, 0.0) if amp == 1.0 else np.where(x < 1, 20.0, np.where(x < 2, 6.0, 0.0))
            rm = rm + g.normal(0, 1, n_src) * add
        src = dict(ra=ra, dec=dec, l=l, b=b, rm=rm, e=e, simple=np.ones(n_src, bool))
        S, cnt = sample(src, CFG['fg_radius'], CFG['min_neigh'])
        names = list(CFG['ranges'])
        rr, xx, _ = profile_test(S, cl, [CFG['ranges'][k] for k in names], 120, seed, label=name)
        R = dict(zip(names, rr))
        out[name] = dict(counts=cnt, n={k: R[k]['n'] for k in names}, t={k: round(R[k]['t'], 2) for k in names},
                         excess_quadrature={k: round(R[k]['excess_quadrature'], 2) for k in names}, decision=decide_176(R))
    print(json.dumps(out, indent=1))
    json.dump(out, open(os.path.join(HERE, 'synthetic.json'), 'w'), indent=1)


def run():
    t0 = time.time()
    log = lambda s: print(s, flush=True)
    rows, src, sinfo = load_spice()
    log('spice %s' % sinfo)
    cl, cinfo = load_clusters()
    log('clusters %s' % cinfo)
    res = dict(generated=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), config=CFG, spice=sinfo, clusters=cinfo)
    res.update(full_test(src, cl, log=log))
    res['seconds'] = round(time.time() - t0)
    open(os.path.join(HERE, 'results.json'), 'w').write(json.dumps(res, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o)))
    log('decision %s' % json.dumps(res['decision']))
    log('written results.json %d s' % res['seconds'])


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else ''
    {'extract': extract, 'synthetic': synthetic, 'run': run}.get(cmd, lambda: print(__doc__))()
