"""Entry 180: is the cluster Faraday signal there on sight lines to quasars known to lie behind the cluster?

Entries 176 and 179 found the rotation-measure scatter of SPICE-RACS sources raised inside R500 of clusters on both halves of the
sky. Most of those sources have no redshift. Inside 0.5 R500 about 40% are the clusters' own radio galaxies (Entry 176), and
between 0.5 and 1 R500 the absence of a source surplus argues that the excess is gas in front of unrelated sources. This entry
tests that directly. Only sources identified with a quasar whose redshift exceeds the cluster's by more than 0.1 are kept; such a
source cannot be a member and must lie behind the gas. The rule was frozen in the Workbench (section 180) before any quasar
redshift was set beside a cluster.

Data (git-ignored under data/), all used before:
  data/spice-racs-dr2/spice-racs.dr2.fits.gz    through Entry 171's cache (entries/171-spice-faraday-redshift/work/spice_rm.npz)
  data/milliquas-v8/milliquas.zip               Milliquas v8 (Flesch 2023), classes Q, A, B, K, N with z > 0 (Entry 171)
  data/desi-dr1-qso/QSO_cat_iron_cumulative_v0.fits   DESI DR1 quasars (Entry 173)
  data/quaia/quaia_G20.5.fits                   Quaia G < 20.5 (Entry 174)
  data/erass1-clusters/erass1cl_primary_v3.2.fits    eRASS1 clusters, as Entry 176
  data/mcxc/table.tsv, data/psz2/HFI_PCCS_SZ-union_R2.08.fits.gz   MCXC and PSZ2 at longitude < 179.5 deg, as Entry 179
cluster_rm.py (Entry 176) and other_half.py (Entry 179) are copied unchanged.

Rule (section 180):
  Sources    Entry 176's sample (|b| >= 20 deg, dRM against neighbours within 3 deg).  A source gets a redshift from an object
             within 10" of any of its catalogue rows: Milliquas first, then DESI DR1, then Quaia (Entry 177's order).
  Clusters   Entry 176's 9,050 eRASS1 clusters plus Entry 179's 1,101 MCXC and PSZ2 clusters.
  Background x = smallest theta / theta500 over all of them; a source is used only if it has a redshift above that cluster's + 0.1.
             Sources more than 8 theta500 from every cluster form the field.
  Statistic  V in A: x < 1, B: 1-2, C: 2-4, over background quasars only; 500 rotations of the whole cluster list in Galactic
             longitude, the background condition applied to each rotated nearest cluster; t = (V - mean) / sd.
  Decision   decide_180() below.

Run from the repository root:
  .venv/bin/python -I entries/180-cluster-faraday-background/background.py synthetic
  .venv/bin/python -I entries/180-cluster-faraday-background/background.py run
"""
import json, os, sys, time, zipfile
import numpy as np
from scipy.spatial import cKDTree
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cluster_rm as c
import other_half as oh

ROOT = c.ROOT
QSO = os.path.join(ROOT, 'data', 'desi-dr1-qso', 'QSO_cat_iron_cumulative_v0.fits')
QUAIA = os.path.join(ROOT, 'data', 'quaia', 'quaia_G20.5.fits')
MQZIP = os.path.join(ROOT, 'data', 'milliquas-v8', 'milliquas.zip')
CFG = dict(match_arcsec=10.0, classes=['Q', 'A', 'B', 'K', 'N'], dz=0.1, nrot=500, seed=20261180, min_A=30,
           ranges=[[0.0, 1.0], [1.0, 2.0], [2.0, 4.0]], names=['A', 'B', 'C'])


# ---------------------------------------------------------------- redshift catalogues (Entry 177's readers)
def load_desi():
    t, nrows = c.read_table(QSO, ['TARGETID', 'Z', 'TARGET_RA', 'TARGET_DEC', 'TSNR2_QSO'])
    tid, ts = t['TARGETID'], t['TSNR2_QSO'].astype(float)
    order = np.lexsort((-np.where(np.isfinite(ts), ts, -np.inf), tid))
    k = order[np.r_[True, tid[order][1:] != tid[order][:-1]]]
    z, ra, dec = t['Z'][k].astype(float), t['TARGET_RA'][k].astype(float), t['TARGET_DEC'][k].astype(float)
    ok = np.isfinite(z) & (z > 0) & np.isfinite(ra) & np.isfinite(dec)
    return dict(z=z[ok], ra=ra[ok], dec=dec[ok])


def load_milliquas():
    with zipfile.ZipFile(MQZIP) as zf:
        raw = zf.read('milliquas.txt')
    ra, dec, z = [], [], []
    for line in raw.split(b'\n'):
        line = line.rstrip(b'\r')
        if len(line) < 82:
            continue
        if line[51:52].decode('ascii', 'replace') not in CFG['classes']:
            continue
        try:
            r, d, zz = float(line[0:11]), float(line[12:23]), float(line[76:82])
        except ValueError:
            continue
        if zz > 0:
            ra.append(r); dec.append(d); z.append(zz)
    return dict(ra=np.array(ra), dec=np.array(dec), z=np.array(z))


def load_quaia():
    t, nrows = c.read_table(QUAIA, ['ra', 'dec', 'redshift_quaia'])
    z = t['redshift_quaia'].astype(float)
    ok = np.isfinite(z) & (z > 0) & np.isfinite(t['ra']) & np.isfinite(t['dec'])
    return dict(ra=t['ra'][ok].astype(float), dec=t['dec'][ok].astype(float), z=z[ok])


def source_redshifts(rows, nsrc, cats):
    """Redshift per source from the nearest catalogue object within 10" of any of its rows, catalogues in priority order."""
    tree = cKDTree(c.unit(rows['row_ra'], rows['row_dec']))
    zs, which, best = np.full(nsrc, np.nan), np.full(nsrc, -1), np.full(nsrc, np.inf)
    for k, cat in enumerate(cats):
        d, j = tree.query(c.unit(cat['ra'], cat['dec']), k=1, distance_upper_bound=c.chord(CFG['match_arcsec'] / 3600.0))
        hit = np.nonzero(np.isfinite(d))[0]
        for i, si, di in zip(hit, rows['row_src'][j[hit]], d[hit]):
            if which[si] == -1 or (which[si] == k and di < best[si]):
                zs[si], which[si], best[si] = cat['z'][i], k, di
    return zs, which


def sample_with_index(src):
    """Entry 176's sample, also returning the source indices kept."""
    idx = np.nonzero(np.abs(src['b']) >= c.CFG['bcut'])[0]
    med, nn = c.foreground(src, idx, c.CFG['fg_radius'], c.CFG['min_neigh'])
    ok = np.isfinite(med)
    k = idx[ok]
    return dict(l=src['l'][k], b=src['b'][k], d=src['rm'][k] - med[ok], e=src['e'][k]), k


def cells(x, which, zq, zcl, d, e):
    bg = np.isfinite(zq) & (which >= 0)
    bg[bg] = zq[bg] > zcl[which[bg]] + CFG['dz']
    out = []
    for lo, hi in CFG['ranges']:
        m = bg & (x >= lo) & (x < hi)
        out.append((c.V_of(d[m], e[m]) if m.sum() >= 15 else np.nan, int(m.sum())))
    return out


def analyse(S, zq, cl, nrot, seed, log=print, label='background quasars'):
    tree = cKDTree(c.gal_unit(S['l'], S['b']))
    n = len(S['d'])
    x, which = c.scaled_distance(tree, n, cl)
    real = cells(x, which, zq, cl['z'], S['d'], S['e'])
    g = np.random.default_rng(seed)
    nulls = np.empty((nrot, 3)); cnts = np.empty((nrot, 3))
    for i in range(nrot):
        xr, wr = c.scaled_distance(tree, n, cl, shift=g.uniform(10.0, 350.0))
        r = cells(xr, wr, zq, cl['z'], S['d'], S['e'])
        nulls[i] = [a for a, _ in r]; cnts[i] = [b for _, b in r]
    res = {}
    for j, name in enumerate(CFG['names']):
        v = nulls[:, j]; v = v[np.isfinite(v)]
        mu, sd = float(v.mean()), float(v.std(ddof=1))
        V, cnt = real[j]
        res[name] = dict(lo=CFG['ranges'][j][0], hi=CFG['ranges'][j][1], n=cnt, null_n_mean=float(cnts[:, j].mean()), V=V, null_mean=mu,
                         null_sd=sd, t=float((V - mu) / sd), excess_variance=float(V - mu))
        log('%s x %.0f-%.0f: n %d (rotated %.0f), V %.1f, null %.1f +- %.1f, t %+.2f, excess %+.1f' % (
            label, CFG['ranges'][j][0], CFG['ranges'][j][1], cnt, cnts[:, j].mean(), V, mu, sd, (V - mu) / sd, V - mu))
    field = np.isfinite(zq) & ~np.isfinite(x)
    res['field'] = dict(n=int(field.sum()), V=c.V_of(S['d'][field], S['e'][field]))
    return res, x, which


def decide_180(R):
    if R['A']['n'] < CFG['min_A']:
        return dict(outcome='not tested: fewer than %d background quasars inside R500' % CFG['min_A'])
    t = R['A']['t']
    out = dict(A=('seen: the scatter of background quasars is raised inside R500 (t >= +3); the gas is in front of the sources' if t >= 3 else
                  'not seen on background quasars inside R500 (t < +2)' if t < 2 else 'unresolved inside R500 (2 <= t < 3)'))
    for k, name in (('B', '1 to 2 R500'), ('C', '2 to 4 R500')):
        tk = R[k]['t']
        out[k] = ('excess at %s (t >= +3)' % name if tk >= 3 else 'no excess at %s (|t| < 2)' % name if abs(tk) < 2 else
                  'unresolved at %s' % name)
    return out


def all_clusters():
    e1, i1 = c.load_clusters()
    e2, i2 = oh.load_clusters()
    cl = {k: np.concatenate([e1[k], e2[k]]) for k in ('l', 'b', 'th', 'z', 'm500')}
    return cl, dict(erass1=i1['kept'], other_half=i2['clusters'], total=int(len(cl['z'])))


def synthetic():
    res = {}
    for name, mode, seed in (('gas_in_front_of_all', 'gas', 71), ('members_only', 'members', 72), ('no_effect', 'none', 73)):
        g = np.random.default_rng(seed)
        n_src = 120000
        dec = np.arcsin(g.uniform(np.sin(-90 * c.D2R), np.sin(49 * c.D2R), n_src)) / c.D2R; ra = g.uniform(0, 360, n_src)
        l, b = c.galactic(ra, dec)
        e = g.lognormal(np.log(1.7), 0.5, n_src)
        rm = 25 / (np.abs(np.sin(b * c.D2R)) + 0.3) * np.sin(l * c.D2R + 0.7) + g.normal(0, 8, n_src) + g.normal(0, 1, n_src) * e
        nk = 11000
        cl = dict(l=g.uniform(0, 360, nk), b=np.arcsin(g.uniform(-1, 1, nk)) / c.D2R, th=g.lognormal(np.log(3.5 / 60), 0.4, nk),
                  z=g.uniform(0.05, 0.6, nk), m500=g.lognormal(np.log(20.0), 0.5, nk))
        tree = cKDTree(c.gal_unit(l, b))
        x, which = c.scaled_distance(tree, n_src, cl)
        zq = np.where(g.uniform(size=n_src) < 0.11, g.uniform(0.3, 3.0, n_src), np.nan)        # 11% of sources are quasars
        inside = x < 1
        if mode == 'gas':                       # every source behind the cluster is rotated, quasar or not
            zsrc = np.where(np.isfinite(zq), zq, g.uniform(0.0, 2.0, n_src))
            behind = inside & (zsrc > cl['z'][np.maximum(which, 0)])
            rm = rm + g.normal(0, 1, n_src) * np.where(behind, 11.0, 0.0)
        elif mode == 'members':                 # only non-quasar sources inside R500 (members) are rotated
            mem = inside & ~np.isfinite(zq)
            rm = rm + g.normal(0, 1, n_src) * np.where(mem, 11.0, 0.0)
        src = dict(ra=ra, dec=dec, l=l, b=b, rm=rm, e=e, simple=np.ones(n_src, bool))
        S, k = sample_with_index(src)
        R, _, _ = analyse(S, zq[k], cl, 100, seed, log=lambda s: print(name, s, flush=True), label=name)
        res[name] = dict(n={a: R[a]['n'] for a in CFG['names']}, t={a: round(R[a]['t'], 2) for a in CFG['names']}, decision=decide_180(R))
    print(json.dumps(res, indent=1))
    json.dump(res, open(os.path.join(HERE, 'synthetic.json'), 'w'), indent=1)


def run():
    t0 = time.time()
    log = lambda s: print(s, flush=True)
    rows, src, sinfo = c.load_spice()
    S, k = sample_with_index(src)
    mq, de, qu = load_milliquas(), load_desi(), load_quaia()
    zs, cat = source_redshifts(rows, len(src['rm']), [mq, de, qu])
    zq = zs[k]
    cl, cinfo = all_clusters()
    info = dict(sources=int(len(k)), with_redshift=int(np.isfinite(zq).sum()),
                by_catalogue={n: int((cat[k] == i).sum()) for i, n in enumerate(['milliquas', 'desi', 'quaia'])}, clusters=cinfo)
    log('sample %s' % json.dumps(info))
    R, x, which = analyse(S, zq, cl, CFG['nrot'], CFG['seed'], log=log)
    res = dict(generated=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), config=CFG, info=info, ranges=R, decision=decide_180(R))
    log('decision %s' % json.dumps(res['decision']))
    # reported only: the same clusters with every source (quasar or not), and quasars in front of or near the cluster's redshift
    allz = np.zeros(len(zq)) + 99.0
    Ra, _, _ = analyse(S, allz, cl, 200, CFG['seed'] + 1, log=log, label='all sources')
    rep = dict(all_sources=Ra, ratio_A=float(R['A']['excess_variance'] / Ra['A']['excess_variance']))
    bgA = np.isfinite(zq) & (x < 1) & (which >= 0)
    if bgA.any():
        dzs = zq[bgA] - cl['z'][which[bgA]]
        rep['quasars_inside_R500'] = dict(total=int(bgA.sum()), behind_by_more_than_0p1=int((dzs > CFG['dz']).sum()),
                                          median_quasar_z=float(np.median(zq[bgA])), median_cluster_z=float(np.median(cl['z'][which[bgA]])))
    res['reported_only'] = rep
    log('reported: background-quasar excess / all-source excess inside R500 = %.2f' % rep['ratio_A'])
    res['seconds'] = round(time.time() - t0)
    json.dump(res, open(os.path.join(HERE, 'results.json'), 'w'), indent=1, default=float)
    log('written results.json %d s' % res['seconds'])


if __name__ == '__main__':
    {'synthetic': synthetic, 'run': run}[sys.argv[1]]()
