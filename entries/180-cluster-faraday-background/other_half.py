"""Entry 179: does the Faraday excess of Entry 176 appear on the other half of the sky, around clusters eRASS1 did not select?

Entry 176 used eRASS1 clusters, which cover only the western Galactic hemisphere (Galactic longitude 180-360 deg). SPICE-RACS
also covers much of the eastern hemisphere, where Entry 176 used no cluster. This entry repeats Entry 176's statistic there with
clusters from two other catalogues: MCXC, compiled from ROSAT X-ray surveys, and the Planck PSZ2 catalogue, selected by the SZ
effect. The rule was frozen in the Workbench (section 179) before this script read any cluster position or rotation measure.

Data (raw files git-ignored under data/):
  data/mcxc/table.tsv   MCXC, Piffaretti et al. 2011 (A&A 534, A109), VizieR J/A+A/534/A109/mcxc from vizier.cfa.harvard.edu
      (asu-tsv, all columns), fetched 21 November 2026, 551,285 bytes,
      sha256 9ff800289c01547bce757b48af21cffbf40585139ab3d9ce553962b70d2e4948.  Columns: RAJ2000, DEJ2000, z, M500 (1e14 Msun),
      R500 (Mpc).
  data/psz2/HFI_PCCS_SZ-union_R2.08.fits.gz   Planck 2015 SZ union catalogue (Planck Collaboration XXVII 2016, A&A 594, A27),
      https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/HFI_PCCS_SZ-union_R2.08.fits.gz, fetched 21 November 2026,
      155,324 bytes, sha256 afba0c3e2b7341579aa5b0212dd0ee4ae15ba19493329e2a4827d17db52053a5.  Columns: GLON, GLAT, REDSHIFT,
      MSZ (1e14 Msun), MCXC (name of the MCXC counterpart, blank if none).
  SPICE-RACS DR2 as in Entry 176 (through Entry 171's cache).

Rule, in brief (section 179 has the full text):
  Clusters  MCXC with z >= 0.05 and Galactic longitude below 179.5 deg; R500 as listed.  PSZ2 entries with REDSHIFT >= 0.05,
            MSZ > 0, longitude below 179.5 deg, no MCXC counterpart listed and none within 5' in the list above; R500 from MSZ and
            the redshift (critical density, flat LCDM H0 = 67.74, Om = 0.3089).  theta500 = R500 / D_A(z).
  Sources   Entry 176's 92,169 (|b| >= 20 deg, dRM against neighbours within 3 deg).
  Statistic x = smallest theta / theta500 over these clusters; V in A: x < 1, B: 1-2, C: 2-4; 500 rotations in Galactic longitude.
  Decision  decide_179() below.

Run from the repository root:
  .venv/bin/python -I entries/179-cluster-faraday-other-half/other_half.py synthetic   # synthetic sky; reads no data
  .venv/bin/python -I entries/179-cluster-faraday-other-half/other_half.py run         # the test; writes results.json
"""
import json, os, sys, time
import numpy as np
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cluster_rm as c

ROOT = c.ROOT
MCXC = os.path.join(ROOT, 'data', 'mcxc', 'table.tsv')
PSZ2 = os.path.join(ROOT, 'data', 'psz2', 'HFI_PCCS_SZ-union_R2.08.fits.gz')
CFG = dict(lmax=179.5, zmin=0.05, dedup_arcmin=5.0, min_sources_A=40, nrot=500, seed=20261179,
           ranges=dict(A=[0.0, 1.0], B=[1.0, 2.0], C=[2.0, 4.0]), nrot_check=200,
           ref_176=dict(excess_at_2p54e14=176.0, exponent=0.86, note='Entry 176: excess inside R500 against M500, mass quarters'))


def sexa(s, hours):
    p = s.replace(':', ' ').split()
    sign = -1.0 if p[0].startswith('-') else 1.0
    v = abs(float(p[0])) + float(p[1]) / 60 + (float(p[2]) if len(p) > 2 else 0.0) / 3600
    return sign * v * (15.0 if hours else 1.0)


def load_mcxc():
    lines = [ln.rstrip('\n') for ln in open(MCXC) if ln.strip() and not ln.startswith('#')]
    head = lines[0].split('\t')
    rows = [ln.split('\t') for ln in lines[3:]]                  # after the header, units and dashes rows
    col = {k: i for i, k in enumerate(head)}
    ra = np.array([sexa(r[col['RAJ2000']], True) for r in rows])
    dec = np.array([sexa(r[col['DEJ2000']], False) for r in rows])
    f = lambda k: np.array([float(r[col[k]]) if r[col[k]].strip() else np.nan for r in rows])
    return dict(ra=ra, dec=dec, z=f('z'), m500=f('M500'), r500=f('R500') * 1000.0)   # r500 in kpc, m500 in 1e14


def r500_from_m(m14, z):
    Ez2 = c.CFG['Om'] * (1 + z) ** 3 + 1 - c.CFG['Om']
    rho_c = 2.775e11 * (c.CFG['H0'] / 100) ** 2 * Ez2          # Msun / Mpc^3
    return (3 * m14 * 1e14 / (4 * np.pi * 500 * rho_c)) ** (1 / 3) * 1000.0


def load_clusters():
    from astropy.io import fits
    m = load_mcxc()
    ml, mb = c.galactic(m['ra'], m['dec'])
    okm = np.isfinite(m['z']) & (m['z'] >= CFG['zmin']) & (ml < CFG['lmax']) & np.isfinite(m['r500']) & (m['r500'] > 0)
    with fits.open(PSZ2) as h:
        d = h[1].data
        pl, pb = np.asarray(d['GLON'], float), np.asarray(d['GLAT'], float)
        pz, pm = np.asarray(d['REDSHIFT'], float), np.asarray(d['MSZ'], float)
        pmc = np.array([str(s).strip() for s in d['MCXC']])
    okp = np.isfinite(pz) & (pz >= CFG['zmin']) & np.isfinite(pm) & (pm > 0) & (pl < CFG['lmax']) & (pmc == '')
    # PSZ2 entries within 5' of any MCXC cluster (all MCXC, any longitude) are dropped as duplicates
    allm = cKDTree(c.gal_unit(ml, mb))
    near = np.array([len(L) > 0 for L in allm.query_ball_point(c.gal_unit(pl, pb), c.chord(CFG['dedup_arcmin'] / 60.0))])
    okp &= ~near
    l = np.concatenate([ml[okm], pl[okp]]); b = np.concatenate([mb[okm], pb[okp]])
    z = np.concatenate([m['z'][okm], pz[okp]]); m14 = np.concatenate([m['m500'][okm], pm[okp]])
    r500 = np.concatenate([m['r500'][okm], r500_from_m(pm[okp], pz[okp])])
    src = np.concatenate([np.zeros(okm.sum(), int), np.ones(okp.sum(), int)])
    th = r500 / 1000.0 / c.ang_diam_Mpc(z) / c.D2R
    info = dict(mcxc_rows=int(len(m['z'])), mcxc_kept=int(okm.sum()), psz2_rows=int(len(pz)), psz2_kept=int(okp.sum()),
                psz2_dropped_as_duplicates=int((np.isfinite(pz) & (pz >= CFG['zmin']) & (pm > 0) & (pl < CFG['lmax']) & (pmc == '') & near).sum()),
                clusters=int(len(z)), z_median=float(np.median(z)), m500_median_1e13=float(np.median(m14) * 10),
                theta500_arcmin_median=float(np.median(th) * 60), theta500_arcmin_max=float(th.max() * 60))
    return dict(l=l, b=b, th=th, z=z, m500=m14 * 10.0, cat=src), info      # m500 in 1e13 as in Entry 176


def decide_179(R, nA):
    if nA < CFG['min_sources_A']:
        return dict(outcome='not tested: fewer than %d sources inside R500' % CFG['min_sources_A'])
    tA = R['A']['t']
    out = dict(A=('replicated: the scatter inside R500 is raised (t >= +3)' if tA >= 3 else
                  'not replicated: no excess inside R500 (t < +2)' if tA < 2 else 'unresolved inside R500 (2 <= t < 3)'))
    for k, name in (('B', '1 to 2 R500'), ('C', '2 to 4 R500')):
        t = R[k]['t']
        out[k] = ('excess scatter at %s (t >= +3)' % name if t >= 3 else 'no excess at %s (|t| < 2)' % name if abs(t) < 2 else
                  'less scatter at %s (t <= -3)' % name if t <= -3 else 'unresolved at %s (2 <= |t| < 3)' % name)
    return out


def analyse(S, cl, nrot, seed, log=print, variants=True):
    names = list(CFG['ranges'])
    rr, x, which = c.profile_test(S, cl, [CFG['ranges'][k] for k in names], nrot, seed, log=log, label='primary')
    R = dict(zip(names, rr))
    out = dict(ranges=R, decision=decide_179(R, R['A']['n']))
    log('decision %s' % json.dumps(out['decision']))
    if not variants:
        return out
    rep = {}
    rep['permutation_null'] = dict(zip(names, c.permutation_test(S, x, [CFG['ranges'][k] for k in names], 500, seed + 3)))
    log('permutation null t: %s' % {k: round(v['t'], 2) for k, v in rep['permutation_null'].items()})
    # amplitude against Entry 176's mass relation, at the median mass of the clusters that hold the sources inside R500
    inA = x < 1
    mA = float(np.median(cl['m500'][which[inA]])) if inA.any() else np.nan
    ref = CFG['ref_176']
    pred = ref['excess_at_2p54e14'] * (mA / 25.4) ** ref['exponent']
    rep['amplitude'] = dict(median_m500_1e13_of_clusters_holding_A=mA, median_z_of_those=float(np.median(cl['z'][which[inA]])) if inA.any() else None,
                            excess_A=R['A']['excess_variance'], error_A=R['A']['null_sd'], entry176_relation_predicts=float(pred),
                            ratio=float(R['A']['excess_variance'] / pred))
    log('amplitude: excess inside R500 %+.1f +- %.1f; Entry 176 relation at M500 %.1fe13 predicts %.1f (ratio %.2f)' % (
        R['A']['excess_variance'], R['A']['null_sd'], mA, pred, R['A']['excess_variance'] / pred))
    for name, sel in (('mcxc_only', cl['cat'] == 0), ('psz2_only', cl['cat'] == 1)):
        if sel.sum() < 30:
            continue
        r, _, _ = c.profile_test(S, {k: v[sel] for k, v in cl.items()}, [CFG['ranges'][k] for k in names], CFG['nrot_check'], seed + 21, log=log, label=name)
        rep[name] = dict(zip(names, r))
    edges = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0]
    fine, _, _ = c.profile_test(S, cl, [[edges[i], edges[i + 1]] for i in range(len(edges) - 1)], CFG['nrot_check'], seed + 11, log=log, label='profile')
    rep['profile'] = fine
    out['reported_only'] = rep
    return out


def synthetic():
    res = {}
    for name, amp, seed in (('no_cluster_effect', 0.0, 51), ('excess_inside_R500', 1.0, 52)):
        g = np.random.default_rng(seed)
        n_src = 120000
        dec = np.arcsin(g.uniform(np.sin(-90 * c.D2R), np.sin(49 * c.D2R), n_src)) / c.D2R; ra = g.uniform(0, 360, n_src)
        l, b = c.galactic(ra, dec)
        fg = 25 / (np.abs(np.sin(b * c.D2R)) + 0.3) * np.sin(l * c.D2R + 0.7)
        e = g.lognormal(np.log(1.7), 0.5, n_src)
        rm = fg + g.normal(0, 8, n_src) + g.normal(0, 1, n_src) * e
        ncl = 1500
        cl_l = g.uniform(0, 179.5, ncl); cl_b = np.arcsin(g.uniform(-1, 1, ncl)) / c.D2R
        keep = np.abs(cl_b) > 15; nk = int(keep.sum())
        z = g.uniform(0.05, 0.3, nk)
        cl = dict(l=cl_l[keep], b=cl_b[keep], th=g.lognormal(np.log(6 / 60), 0.4, nk), z=z, m500=g.lognormal(np.log(30), 0.5, nk), cat=(g.uniform(size=nk) < 0.3).astype(int))
        tree = cKDTree(c.gal_unit(l, b))
        x, _ = c.scaled_distance(tree, n_src, cl)
        if amp > 0:
            rm = rm + g.normal(0, 1, n_src) * np.where(x < 1, 12.0, 0.0)
        src = dict(ra=ra, dec=dec, l=l, b=b, rm=rm, e=e, simple=np.ones(n_src, bool))
        S, cnt = c.sample(src, c.CFG['fg_radius'], c.CFG['min_neigh'])
        out = analyse(S, cl, 100, seed, log=lambda s: print(name, s, flush=True), variants=False)
        res[name] = dict(n={k: v['n'] for k, v in out['ranges'].items()}, t={k: round(v['t'], 2) for k, v in out['ranges'].items()},
                         excess={k: round(v['excess_variance'], 1) for k, v in out['ranges'].items()}, decision=out['decision'])
    print(json.dumps(res, indent=1))
    json.dump(res, open(os.path.join(HERE, 'synthetic.json'), 'w'), indent=1)


def run():
    t0 = time.time()
    log = lambda s: print(s, flush=True)
    rows, src, sinfo = c.load_spice()
    S, cnt = c.sample(src, c.CFG['fg_radius'], c.CFG['min_neigh'])
    cl, cinfo = load_clusters()
    log('sources %s; clusters %s' % (cnt, cinfo))
    res = dict(generated=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), config=CFG, spice=sinfo, sources=cnt, clusters=cinfo)
    res.update(analyse(S, cl, CFG['nrot'], CFG['seed'], log=log))
    res['seconds'] = round(time.time() - t0)
    json.dump(res, open(os.path.join(HERE, 'results.json'), 'w'), indent=1, default=float)
    log('written results.json %d s' % res['seconds'])


if __name__ == '__main__':
    {'synthetic': synthetic, 'run': run}[sys.argv[1]]()
