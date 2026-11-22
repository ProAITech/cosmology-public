"""Entry 167: the Giant Arc on DESI DR1 Mg II absorbers (Workbench section 167, rule frozen in commit 96074ac, 14 Nov 2026).

Data: data/desi-dr1-mgii/MgII-Absorbers-DR1.fits (DESI DR1 Mg II absorber VAC v1.0, https://data.desi.lbl.gov/doc/releases/dr1/vac/mgii-absorber/,
65,750,400 bytes, sha256 38c2e6281ce983a2b34b5810fdf0e3021cebaad3848cfc1cca974f185901bddf), table MGII_ABSORBERS, 270,529 rows.
Columns used: TARGETID, RA, DEC, ZWARN, Z_QSO, Z_MGII, LINE_SNR_MIN, TSNR2_QSO.
Run from the repository root:  .venv/bin/python -I entries/167-desi-arc/desi_arc.py        (about a minute; deterministic, no random numbers)

Rule as frozen:
 1 Sample: ZWARN = 0, 0.862 < Z_QSO <= 3.1, LINE_SNR_MIN >= 3; one row per absorber (same TARGETID within dz 0.001, or two TARGETIDs within 1 arcsec and dz 0.001: keep the
   larger LINE_SNR_MIN, ties the larger TSNR2_QSO).
 2 Slice 0.742 <= Z_MGII <= 0.862 (2,665 to 3,003 Mpc, flat, Om 0.27, H0 70); cluster RA 215.405, Dec +38.361; gnomonic offsets x (east), y (north) times 2,837 Mpc;
   field |x| <= 770.5, |y| <= 807.5.
 3 Band |x| <= 500, |y| <= 100 Mpc; N_band = sample absorbers in the slice inside it.
 4 N_oth = sample absorbers in the same band at comoving distance 1,815-2,515 or 3,153-3,853 Mpc; R = (slice absorbers) / (such outer absorbers) over the whole table with the
   Arc's field left out; N_exp = R N_oth; E = N_band / N_exp - 1.
 5 Null: same band in the same slice at lattice nodes Dec -10..+80 step 5, RA from 0 in steps of 22/cos(Dec); each of the band's five 200-Mpc cells along x must hold at least
   half the median outer count per cell; nodes whose band overlaps the Arc's dropped; >= 40 nodes; m, sd of E.
 6 Decision: supported if E >= +0.20 and E exceeds every null value; no excess if (E - m)/sd < 2 (95% limit E + 1.645 sd); else unresolved. z against Entry 66's +0.34 +- 0.19.
 7 Checks: (a) 300-Mpc-wide band; (b) the slices just below (0.630-0.742) and above (0.862-0.991); (c) linking at 95 (504/N_field)^(1/3) Mpc in the field and its rank among lattice fields.
Implementation choices the rule left open are listed in README.md.
"""
import json, os, time
import numpy as np
from astropy.io import fits
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
import warnings
warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
PATH = os.path.join(ROOT, 'data', 'desi-dr1-mgii', 'MgII-Absorbers-DR1.fits')
D2R = np.pi / 180.0
OM, H0, C = 0.27, 70.0, 299792.458
RA_C, DEC_C, SCALE = 215.405, 38.361, 2837.0
FIELD_HX, FIELD_HY = 770.5, 807.5
BAND_HX, BAND_HY = 500.0, 100.0
SLICE = (0.742, 0.862)
PRIOR_E, PRIOR_ERR = 0.34, 0.19                      # Entry 66


# ---------------------------------------------------------------- cosmology
_zg = np.linspace(0, 4, 400001)
_Ez = np.sqrt(OM * (1 + _zg) ** 3 + 1 - OM)
_dg = np.concatenate([[0], np.cumsum((1 / _Ez[1:] + 1 / _Ez[:-1]) / 2 * np.diff(_zg))]) * C / H0


def dist(z):
    return np.interp(z, _zg, _dg)


def zofd(d):
    return np.interp(d, _dg, _zg)


# ---------------------------------------------------------------- sample
def load():
    with fits.open(PATH, memmap=True) as h:
        d = h[1].data
        t = {k: np.array(d[k]) for k in ('TARGETID', 'RA', 'DEC', 'ZWARN', 'Z_QSO', 'Z_MGII', 'LINE_SNR_MIN', 'TSNR2_QSO')}
    n0 = len(t['RA'])
    ok = (t['ZWARN'] == 0) & (t['Z_QSO'] > 0.862) & (t['Z_QSO'] <= 3.1) & (t['LINE_SNR_MIN'] >= 3)
    t = {k: v[ok] for k, v in t.items()}
    n1 = len(t['RA'])
    # duplicates: better row = larger LINE_SNR_MIN, then larger TSNR2_QSO
    order = np.lexsort((-t['TSNR2_QSO'], -t['LINE_SNR_MIN']))             # best first
    t = {k: v[order] for k, v in t.items()}
    keep = np.ones(n1, bool)
    ra, dec = t['RA'] * D2R, t['DEC'] * D2R
    xyz = np.column_stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])
    pairs = cKDTree(xyz).query_pairs(r=2 * np.sin(0.5 * 1.0 / 3600 * D2R), output_type='ndarray')
    same_pos = pairs[np.abs(t['Z_MGII'][pairs[:, 0]] - t['Z_MGII'][pairs[:, 1]]) <= 0.001]
    # same TARGETID within dz 0.001 (positions identical, so caught above unless coordinates differ)
    idt = np.argsort(t['TARGETID'], kind='stable'); tid = t['TARGETID'][idt]
    same_t = []
    start = 0
    while start < n1:
        end = start
        while end + 1 < n1 and tid[end + 1] == tid[start]:
            end += 1
        if end > start:
            g = idt[start:end + 1]
            for a in range(len(g)):
                for b in range(a + 1, len(g)):
                    if abs(t['Z_MGII'][g[a]] - t['Z_MGII'][g[b]]) <= 0.001:
                        same_t.append((g[a], g[b]))
        start = end + 1
    allp = np.vstack([same_pos, np.array(same_t, dtype=same_pos.dtype).reshape(-1, 2)]) if len(same_t) else same_pos
    for a, b in allp:
        keep[max(a, b)] = False                                            # sorted best-first, so the later index is the worse row
    t = {k: v[keep] for k, v in t.items()}
    info = dict(rows=int(n0), after_cuts=int(n1), duplicates_removed=int(n1 - keep.sum()), sample=int(keep.sum()))
    t['dist'] = dist(t['Z_MGII'])
    t['sra'], t['cra'] = np.sin(t['RA'] * D2R), np.cos(t['RA'] * D2R)
    t['sdec'], t['cdec'] = np.sin(t['DEC'] * D2R), np.cos(t['DEC'] * D2R)
    return t, info


def gnomonic(t, ra0, dec0, idx=None):
    a0, d0 = ra0 * D2R, dec0 * D2R
    sa, ca, sd, cd = (t[k] if idx is None else t[k][idx] for k in ('sra', 'cra', 'sdec', 'cdec'))
    dra_s = sa * np.cos(a0) - ca * np.sin(a0)                              # sin(a - a0)
    dra_c = ca * np.cos(a0) + sa * np.sin(a0)                              # cos(a - a0)
    cosc = np.sin(d0) * sd + np.cos(d0) * cd * dra_c
    xi = cd * dra_s / cosc
    eta = (np.cos(d0) * sd - np.sin(d0) * cd * dra_c) / cosc
    return xi * SCALE, eta * SCALE, cosc


# ---------------------------------------------------------------- band counts
def windows(sl):
    dlo, dhi = float(dist(sl[0])), float(dist(sl[1]))
    return (dlo, dhi), [(dlo - 850.0, dlo - 150.0), (dhi + 150.0, dhi + 850.0)]


def band_stats(t, ra0, dec0, sl, hy=BAND_HY, cache=None):
    """Counts in the band at (ra0, dec0): slice count, outer count, and outer counts in the five 200-Mpc cells along x."""
    # restrict to a cap around the node first
    a0, d0 = ra0 * D2R, dec0 * D2R
    cosang = np.sin(d0) * t['sdec'] + np.cos(d0) * t['cdec'] * (t['cra'] * np.cos(a0) + t['sra'] * np.sin(a0))
    idx = np.nonzero(cosang > np.cos(25 * D2R))[0]
    x, y, cosc = gnomonic(t, ra0, dec0, idx)
    inb = (np.abs(x) <= BAND_HX) & (np.abs(y) <= hy) & (cosc > 0.5)
    idx, x, y = idx[inb], x[inb], y[inb]
    (dlo, dhi), outer = windows(sl)
    d = t['dist'][idx]
    in_slice = (d >= dlo) & (d <= dhi)
    in_out = ((d >= outer[0][0]) & (d <= outer[0][1])) | ((d >= outer[1][0]) & (d <= outer[1][1]))
    edges = np.array([-500, -300, -100, 100, 300, 500.0])
    cells = np.histogram(x[in_out], bins=edges)[0]
    return int(in_slice.sum()), int(in_out.sum()), cells.astype(int)


def global_ratio(t, sl):
    """R = slice absorbers / outer absorbers over the whole table with the Arc's field left out."""
    x, y, cosc = gnomonic(t, RA_C, DEC_C)
    infield = (np.abs(x) <= FIELD_HX) & (np.abs(y) <= FIELD_HY) & (cosc > 0.5)
    (dlo, dhi), outer = windows(sl)
    d = t['dist']
    sl_m = (d >= dlo) & (d <= dhi) & ~infield
    out_m = (((d >= outer[0][0]) & (d <= outer[0][1])) | ((d >= outer[1][0]) & (d <= outer[1][1]))) & ~infield
    return float(sl_m.sum() / out_m.sum()), int(sl_m.sum()), int(out_m.sum())


def lattice():
    nodes = []
    for dec in np.arange(-10, 80.001, 5.0):
        step = 22.0 / np.cos(dec * D2R)
        n = int(np.floor(360.0 / step))
        for k in range(n):
            nodes.append((k * step, float(dec)))
    return nodes


def excess_table(t, sl, hy, nodes, arc_xy, R):
    rows = []
    for ra0, dec0 in nodes:
        nb, no, cells = band_stats(t, ra0, dec0, sl, hy)
        rows.append((ra0, dec0, nb, no, cells))
    return rows


def evaluate(t, sl, hy, nodes, R, node_xy):
    """E at the Arc's band and at every lattice node for one slice and band width; applies the coverage rule and the overlap rule."""
    arc = band_stats(t, RA_C, DEC_C, sl, hy)
    rows = [band_stats(t, ra0, dec0, sl, hy) for ra0, dec0 in nodes]
    allcells = np.concatenate([r[2] for r in rows] + [arc[2]])
    med = float(np.median(allcells[allcells > 0])) if (allcells > 0).any() else 0.0      # median over cells holding at least one outer absorber
    thr = 0.5 * med
    def cov(cells):
        return bool((cells >= thr).all())
    out = dict(median_cell=med, threshold=thr, arc=dict(N_band=arc[0], N_oth=arc[1], cells=arc[2].tolist(), covered=cov(arc[2])))
    if arc[1] > 0:
        out['arc']['N_exp'] = R * arc[1]; out['arc']['E'] = arc[0] / (R * arc[1]) - 1.0
        out['arc']['E_err_poisson'] = (1 + out['arc']['E']) * float(np.sqrt(1.0 / max(arc[0], 1) + 1.0 / arc[1]))
    nodesE, kept = [], 0
    for (ra0, dec0), r in zip(nodes, rows):
        x, y = node_xy[(ra0, dec0)]
        if abs(x) < 1000.0 and abs(y) < 200.0:                                 # overlaps the Arc's band
            continue
        if r[1] == 0 or not cov(r[2]):
            continue
        nodesE.append(dict(ra=round(ra0, 3), dec=dec0, N_band=r[0], N_oth=r[1], E=r[0] / (R * r[1]) - 1.0))
    E = np.array([n['E'] for n in nodesE])
    out['null'] = dict(n_nodes=len(nodesE), mean=float(E.mean()) if len(E) else None, sd=float(E.std(ddof=1)) if len(E) > 1 else None,
                       quantiles=[float(v) for v in np.quantile(E, [0.025, 0.16, 0.5, 0.84, 0.975, 1.0])] if len(E) else None, nodes=nodesE)
    return out


# ---------------------------------------------------------------- linking (check c)
def field_members(t, ra0, dec0, sl):
    a0, d0 = ra0 * D2R, dec0 * D2R
    cosang = np.sin(d0) * t['sdec'] + np.cos(d0) * t['cdec'] * (t['cra'] * np.cos(a0) + t['sra'] * np.sin(a0))
    idx = np.nonzero(cosang > np.cos(30 * D2R))[0]
    x, y, cosc = gnomonic(t, ra0, dec0, idx)
    d = t['dist'][idx]; dlo, dhi = float(dist(sl[0])), float(dist(sl[1]))
    m = (np.abs(x) <= FIELD_HX) & (np.abs(y) <= FIELD_HY) & (cosc > 0.5) & (d >= dlo) & (d <= dhi)
    return x[m], y[m], d[m] - SCALE


def largest_group(x, y, zc):
    n = len(x)
    if n < 10:
        return None
    L = 95.0 * (504.0 / n) ** (1.0 / 3.0)
    P = np.column_stack([x, y, zc])
    pairs = cKDTree(P).query_pairs(r=L, output_type='ndarray')
    g = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    nc, lab = connected_components(g, directed=False)
    sizes = np.bincount(lab); k = int(np.argmax(sizes)); m = lab == k
    Q = P[m]
    c = Q.mean(axis=0); u, s, vt = np.linalg.svd(Q - c, full_matrices=False)
    proj = (Q - c) @ vt[0]
    ax = vt[0][:2]
    ang = float(np.degrees(np.arctan2(abs(ax[1]), abs(ax[0]))))            # angle of the long axis' sky projection from east-west
    return dict(N_field=int(n), link_mpc=float(L), members=int(m.sum()), span_mpc=float(proj.max() - proj.min()), centre_xy=[float(c[0]), float(c[1])],
                centre_offset_mpc=float(np.hypot(c[0], c[1])), axis_angle_from_EW_deg=ang, second_largest=int(np.sort(sizes)[-2]) if len(sizes) > 1 else 0)


def main():
    t0 = time.time()
    t, info = load()
    print('sample', info, flush=True)
    sl = SLICE
    R, nsl, nout = global_ratio(t, sl)
    nodes = lattice()
    node_xy = {}
    for ra0, dec0 in nodes:                                                # node position in the Arc's tangent plane, for the overlap rule
        a0, d0, ac, dc = ra0 * D2R, dec0 * D2R, RA_C * D2R, DEC_C * D2R
        cosc = np.sin(dc) * np.sin(d0) + np.cos(dc) * np.cos(d0) * np.cos(a0 - ac)
        xi = np.cos(d0) * np.sin(a0 - ac) / cosc; eta = (np.cos(dc) * np.sin(d0) - np.sin(dc) * np.cos(d0) * np.cos(a0 - ac)) / cosc
        node_xy[(ra0, dec0)] = (xi * SCALE if cosc > 0 else 1e9, eta * SCALE if cosc > 0 else 1e9)
    res = dict(generated=time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime()), rule='Workbench section 167, frozen in 96074ac', sample=info,
               distances_mpc=dict(z0742=float(dist(0.742)), z0802=float(dist(0.802)), z0862=float(dist(0.862))), windows_mpc=windows(sl)[1],
               R=dict(value=R, slice=nsl, outer=nout), lattice_nodes=len(nodes))
    main_ = evaluate(t, sl, BAND_HY, nodes, R, node_xy)
    res['primary'] = main_
    a, nl = main_['arc'], main_['null']
    if not a['covered'] or a.get('N_oth', 0) == 0 or nl['n_nodes'] < 40:
        res['decision'] = dict(outcome='not tested', why='Arc band covered=%s, null nodes=%d' % (a['covered'], nl['n_nodes']))
    else:
        E = a['E']; nullE = np.array([n['E'] for n in nl['nodes']]); m, sd = nl['mean'], nl['sd']
        zscore = (E - m) / sd; above_all = bool(E > nullE.max())
        if E >= 0.20 and above_all:
            out = 'excess supported (independence from the discovery sightlines not yet established)'
        elif zscore < 2:
            out = 'no excess'
        else:
            out = 'unresolved'
        res['decision'] = dict(outcome=out, E=E, E_err_poisson=a['E_err_poisson'], null_mean=m, null_sd=sd, z_vs_null=float(zscore), above_every_null=above_all,
                               rank_p_upper=float((1 + (nullE >= E).sum()) / (1 + len(nullE))), limit95=float(E + 1.645 * sd),
                               z_vs_entry66=float((E - PRIOR_E) / np.sqrt(sd ** 2 + PRIOR_ERR ** 2)))
    # ---- checks
    chk = {}
    wide = evaluate(t, sl, 150.0, nodes, R, node_xy)
    chk['wide_band_300'] = dict(arc=wide['arc'], null_n=wide['null']['n_nodes'], null_mean=wide['null']['mean'], null_sd=wide['null']['sd'],
                                rank_p_upper=(float((1 + (np.array([n['E'] for n in wide['null']['nodes']]) >= wide['arc']['E']).sum()) / (1 + wide['null']['n_nodes']))
                                              if 'E' in wide['arc'] and wide['null']['n_nodes'] else None), entry66='52 against 40.1 expected (+30%)')
    for name, s2 in (('slice_below_0630_0742', (0.630, 0.742)), ('slice_above_0862_0991', (0.862, 0.991))):
        R2, n2, o2 = global_ratio(t, s2)
        ev = evaluate(t, s2, BAND_HY, nodes, R2, node_xy)
        ne = np.array([n['E'] for n in ev['null']['nodes']])
        chk[name] = dict(R=R2, windows_mpc=windows(s2)[1], arc=ev['arc'], null_n=ev['null']['n_nodes'], null_mean=ev['null']['mean'], null_sd=ev['null']['sd'],
                         rank_p_upper=(float((1 + (ne >= ev['arc']['E']).sum()) / (1 + len(ne))) if 'E' in ev['arc'] and len(ne) else None))
    x, y, zc = field_members(t, RA_C, DEC_C, sl)
    lg = largest_group(x, y, zc)
    ranks = []
    if lg:
        for n in main_['null']['nodes']:
            xn, yn, zn = field_members(t, n['ra'], n['dec'], sl)
            g = largest_group(xn, yn, zn)
            if g:
                ranks.append(g['members'])
        ranks = np.array(ranks)
        lg['null_fields'] = int(len(ranks)); lg['null_members_quantiles'] = [float(v) for v in np.quantile(ranks, [0.16, 0.5, 0.84, 1.0])] if len(ranks) else None
        lg['rank_p_upper_members'] = float((1 + (ranks >= lg['members']).sum()) / (1 + len(ranks))) if len(ranks) else None
    chk['linking'] = lg
    res['checks'] = chk
    res['seconds'] = round(time.time() - t0)
    json.dump(res, open(os.path.join(HERE, 'results.json'), 'w'), indent=1)
    print(json.dumps(res['decision'], indent=1)); print('written results.json', res['seconds'], 's')


if __name__ == '__main__':
    main()
