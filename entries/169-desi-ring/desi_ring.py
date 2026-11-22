"""Entry 169: the Big Ring on DESI DR1 Mg II absorbers (Workbench section 169, rule frozen in commit f3374e3, 20 Nov 2026).

Data: data/desi-dr1-mgii/MgII-Absorbers-DR1.fits (DESI DR1 Mg II absorber VAC v1.0, https://data.desi.lbl.gov/doc/releases/dr1/vac/mgii-absorber/,
65,750,400 bytes, sha256 38c2e6281ce983a2b34b5810fdf0e3021cebaad3848cfc1cca974f185901bddf), table MGII_ABSORBERS, 270,529 rows.
Columns used: TARGETID, RA, DEC, ZWARN, Z_QSO, Z_MGII, LINE_SNR_MIN, TSNR2_QSO.
Claims tested: Lopez, Clowes & Williger (2024), A Big Ring on the Sky, JCAP 07, 055, arXiv:2402.07591 (centre read from its Figure 6);
Lopez & Clowes (2026), A Giant Ring on the Sky, arXiv:2604.17534 (ellipse a 588, b 442, theta 0.278 rad, centre (5.12, 52.2) Mpc).
Run from the repository root:  .venv/bin/python -I entries/169-desi-ring/desi_ring.py        (a few minutes; deterministic, no random numbers)

Rule as frozen:
 1 Sample, slice (0.742-0.862), scale (2,837 Mpc), outer windows (1,815-2,515 and 3,153-3,853 Mpc) and E = N / (R N_oth) - 1 as in Entry 167;
   R over the whole table without the Arc's field and without the sky within 10 deg of the ring's nominal centre.
 2 Nominal centre RA 216.0, Dec +50.0; annulus 140 <= r <= 260 Mpc; centre searched over 81 offsets (-100..+100 step 25 Mpc in x and y);
   statistic E_max over the 81; E_nom at offset (0, 0) reported.
 3 Coverage: eight 45-deg sectors of the nominal annulus each hold >= half the median outer count per sector (median over all sectors with >= 1).
 4 Null: lattice Dec -10..+80 step 10, RA in n = floor(360 cos(Dec) / 10.5) equal steps from 0; same search and coverage; nodes within 15 deg dropped;
   >= 30 nodes; m, sd of E_max; sd1 of E_nom.
 5 Decision: supported if E_max >= +0.20 and above every node; no excess if (E_max - m)/sd < 2, limit E_nom + 1.645 sd1; else unresolved.
 6 Checks: (a) annuli 160-240, 120-280; (b) slices 0.630-0.742 and 0.862-0.991; (c) disc r <= 260 at the best centre and its rank;
   (d) interior r < 140 and surround 260-380 at the best centre; (e) the Giant Ring elliptical shell, 40 and 120 Mpc thick, in a frame centred
   240 Mpc south of the nominal centre, rank among nodes whose 16-sector shell passes coverage if >= 10 do.
"""
import json, os, time
import numpy as np
from astropy.io import fits
from scipy.spatial import cKDTree
import warnings
warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
PATH = os.path.join(ROOT, 'data', 'desi-dr1-mgii', 'MgII-Absorbers-DR1.fits')
D2R = np.pi / 180.0
OM, H0, C = 0.27, 70.0, 299792.458
SCALE = 2837.0
ARC_RA, ARC_DEC, ARC_HX, ARC_HY = 215.405, 38.361, 770.5, 807.5        # Entry 167's Arc field
RING_RA, RING_DEC = 216.0, 50.0
ANN = (140.0, 260.0)
OFFS = np.arange(-100.0, 100.001, 25.0)
SLICE = (0.742, 0.862)
GR = dict(a=588.0, b=442.0, theta=0.278, xc=5.12, yc=52.2, frame_dy=-240.0)
CAP_DEG = 13.0                                                          # absorbers projected per node: 380 Mpc + 141 Mpc search is 10.5 deg
GR_CAP_DEG = 24.0

# ---------------------------------------------------------------- cosmology
_zg = np.linspace(0, 4, 400001)
_Ez = np.sqrt(OM * (1 + _zg) ** 3 + 1 - OM)
_dg = np.concatenate([[0], np.cumsum((1 / _Ez[1:] + 1 / _Ez[:-1]) / 2 * np.diff(_zg))]) * C / H0


def dist(z):
    return np.interp(z, _zg, _dg)


# ---------------------------------------------------------------- sample (Entry 167's, unchanged)
def load():
    with fits.open(PATH, memmap=True) as h:
        d = h[1].data
        t = {k: np.array(d[k]) for k in ('TARGETID', 'RA', 'DEC', 'ZWARN', 'Z_QSO', 'Z_MGII', 'LINE_SNR_MIN', 'TSNR2_QSO')}
    n0 = len(t['RA'])
    ok = (t['ZWARN'] == 0) & (t['Z_QSO'] > 0.862) & (t['Z_QSO'] <= 3.1) & (t['LINE_SNR_MIN'] >= 3)
    t = {k: v[ok] for k, v in t.items()}
    n1 = len(t['RA'])
    order = np.lexsort((-t['TSNR2_QSO'], -t['LINE_SNR_MIN']))             # best first
    t = {k: v[order] for k, v in t.items()}
    keep = np.ones(n1, bool)
    ra, dec = t['RA'] * D2R, t['DEC'] * D2R
    xyz = np.column_stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])
    pairs = cKDTree(xyz).query_pairs(r=2 * np.sin(0.5 * 1.0 / 3600 * D2R), output_type='ndarray')
    same_pos = pairs[np.abs(t['Z_MGII'][pairs[:, 0]] - t['Z_MGII'][pairs[:, 1]]) <= 0.001]
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
        keep[max(a, b)] = False
    t = {k: v[keep] for k, v in t.items()}
    info = dict(rows=int(n0), after_cuts=int(n1), duplicates_removed=int(n1 - keep.sum()), sample=int(keep.sum()))
    t['dist'] = dist(t['Z_MGII'])
    t['sra'], t['cra'] = np.sin(t['RA'] * D2R), np.cos(t['RA'] * D2R)
    t['sdec'], t['cdec'] = np.sin(t['DEC'] * D2R), np.cos(t['DEC'] * D2R)
    return t, info


# ---------------------------------------------------------------- geometry
def cosang(t, ra0, dec0):
    a0, d0 = ra0 * D2R, dec0 * D2R
    return np.sin(d0) * t['sdec'] + np.cos(d0) * t['cdec'] * (t['cra'] * np.cos(a0) + t['sra'] * np.sin(a0))


def project(t, ra0, dec0, cap_deg):
    """Indices, gnomonic x (east), y (north) in Mpc and comoving distance for absorbers within cap_deg of (ra0, dec0)."""
    idx = np.nonzero(cosang(t, ra0, dec0) > np.cos(cap_deg * D2R))[0]
    a0, d0 = ra0 * D2R, dec0 * D2R
    sa, ca, sd, cd = t['sra'][idx], t['cra'][idx], t['sdec'][idx], t['cdec'][idx]
    dra_s = sa * np.cos(a0) - ca * np.sin(a0)
    dra_c = ca * np.cos(a0) + sa * np.sin(a0)
    cosc = np.sin(d0) * sd + np.cos(d0) * cd * dra_c
    x = cd * dra_s / cosc * SCALE
    y = (np.cos(d0) * sd - np.sin(d0) * cd * dra_c) / cosc * SCALE
    return x, y, t['dist'][idx]


def inverse_gnomonic(ra0, dec0, x, y):
    xi, eta = x / SCALE, y / SCALE
    a0, d0 = ra0 * D2R, dec0 * D2R
    rho = np.hypot(xi, eta); c = np.arctan(rho)
    if rho == 0:
        return ra0, dec0
    dec = np.arcsin(np.cos(c) * np.sin(d0) + eta * np.sin(c) * np.cos(d0) / rho)
    ra = a0 + np.arctan2(xi * np.sin(c), rho * np.cos(d0) * np.cos(c) - eta * np.sin(d0) * np.sin(c))
    return float((ra / D2R) % 360.0), float(dec / D2R)


def sep_deg(ra1, dec1, ra2, dec2):
    a1, d1, a2, d2 = ra1 * D2R, dec1 * D2R, ra2 * D2R, dec2 * D2R
    return float(np.degrees(np.arccos(np.clip(np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(a1 - a2), -1, 1))))


def windows(sl):
    dlo, dhi = float(dist(sl[0])), float(dist(sl[1]))
    return (dlo, dhi), [(dlo - 850.0, dlo - 150.0), (dhi + 150.0, dhi + 850.0)]


def masks(d, sl):
    (dlo, dhi), outer = windows(sl)
    ins = (d >= dlo) & (d <= dhi)
    out = ((d >= outer[0][0]) & (d <= outer[0][1])) | ((d >= outer[1][0]) & (d <= outer[1][1]))
    return ins, out


def global_ratio(t, sl):
    """R over the whole table without the Arc's field and without the sky within 10 deg of the ring's nominal centre."""
    x, y, _ = project(t, ARC_RA, ARC_DEC, 89.0)
    idx_arc = np.nonzero(cosang(t, ARC_RA, ARC_DEC) > np.cos(89.0 * D2R))[0]
    infield = np.zeros(len(t['RA']), bool)
    infield[idx_arc[(np.abs(x) <= ARC_HX) & (np.abs(y) <= ARC_HY)]] = True
    near_ring = cosang(t, RING_RA, RING_DEC) > np.cos(10.0 * D2R)
    ins, out = masks(t['dist'], sl)
    keep = ~infield & ~near_ring
    return float((ins & keep).sum() / (out & keep).sum()), int((ins & keep).sum()), int((out & keep).sum())


# ---------------------------------------------------------------- ring statistics at one position
def ring_stats(t, ra0, dec0, sl, R, ann=ANN):
    """For the annulus about each of the 81 centres: slice and outer counts and E; sector outer counts at the nominal centre;
    and, at the best centre, the disc, interior and surround."""
    x, y, d = project(t, ra0, dec0, CAP_DEG)
    ins, out = masks(d, sl)
    xs, ys, xo, yo = x[ins], y[ins], x[out], y[out]
    grid = []
    for dx in OFFS:
        for dy in OFFS:
            rs = np.hypot(xs - dx, ys - dy); ro = np.hypot(xo - dx, yo - dy)
            ns = int(((rs >= ann[0]) & (rs <= ann[1])).sum()); no = int(((ro >= ann[0]) & (ro <= ann[1])).sum())
            grid.append((float(dx), float(dy), ns, no, ns / (R * no) - 1.0 if no > 0 else np.nan))
    G = np.array(grid)
    ro = np.hypot(xo, yo)
    m = (ro >= ann[0]) & (ro <= ann[1])
    sect = np.histogram(np.mod(np.arctan2(yo[m], xo[m]), 2 * np.pi), bins=np.linspace(0, 2 * np.pi, 9))[0]
    nom = G[(G[:, 0] == 0) & (G[:, 1] == 0)][0]
    valid = ~np.isnan(G[:, 4])
    res = dict(N=int(nom[2]), N_oth=int(nom[3]), E_nom=float(nom[4]) if nom[3] > 0 else None, sectors=sect.astype(int))
    if valid.any():
        k = int(np.nanargmax(G[:, 4]))
        bx, by = G[k, 0], G[k, 1]
        res.update(E_max=float(G[k, 4]), best_offset=[bx, by], N_best=int(G[k, 2]), N_oth_best=int(G[k, 3]))
        rs = np.hypot(xs - bx, ys - by); rob = np.hypot(xo - bx, yo - by)
        def e_in(lo, hi):
            ns = int(((rs >= lo) & (rs < hi)).sum()); no = int(((rob >= lo) & (rob < hi)).sum())
            return dict(N=ns, N_oth=no, E=(ns / (R * no) - 1.0) if no > 0 else None,
                        E_err=((ns / (R * no)) * np.sqrt(1 / max(ns, 1) + 1 / no)) if no > 0 else None)
        res['disc'] = e_in(0.0, ann[1] + 1e-9)
        res['interior'] = e_in(0.0, ann[0])
        res['surround'] = e_in(ann[1] + 1e-9, 380.0)
    return res


def lattice():
    nodes = []
    for dec in np.arange(-10.0, 80.001, 10.0):
        n = int(np.floor(360.0 * np.cos(dec * D2R) / 10.5))
        for k in range(n):
            nodes.append((k * 360.0 / n, float(dec)))
    return nodes


def evaluate(t, sl, R, nodes, ann=ANN):
    ring = ring_stats(t, RING_RA, RING_DEC, sl, R, ann)
    rows = [ring_stats(t, ra0, dec0, sl, R, ann) for ra0, dec0 in nodes]
    allsec = np.concatenate([r['sectors'] for r in rows] + [ring['sectors']])
    med = float(np.median(allsec[allsec > 0]))
    thr = 0.5 * med
    cov = lambda s: bool((s >= thr).all())
    ring['covered'] = cov(ring['sectors'])
    kept = []
    for (ra0, dec0), r in zip(nodes, rows):
        if sep_deg(ra0, dec0, RING_RA, RING_DEC) < 15.0:
            continue
        if r['N_oth'] == 0 or 'E_max' not in r or not cov(r['sectors']):
            continue
        kept.append(dict(ra=round(ra0, 3), dec=dec0, **r))
    Em = np.array([k['E_max'] for k in kept]); En = np.array([k['E_nom'] for k in kept])
    null = dict(n_nodes=len(kept), m=float(Em.mean()), sd=float(Em.std(ddof=1)), sd1=float(En.std(ddof=1)), m1=float(En.mean()),
                quantiles_Emax=[float(v) for v in np.quantile(Em, [0.025, 0.16, 0.5, 0.84, 0.975, 1.0])]) if len(kept) > 1 else dict(n_nodes=len(kept))
    ring['sectors'] = ring['sectors'].tolist()
    return dict(median_sector=med, threshold=thr, ring=ring, null=null, nodes=kept)


def rank_upper(val, arr):
    arr = np.asarray(arr, float)
    return float((1 + (arr >= val).sum()) / (1 + len(arr)))


# ---------------------------------------------------------------- Giant Ring shell (check e)
def shell_stats(t, ra0, dec0, sl, R, dr):
    x, y, d = project(t, ra0, dec0, GR_CAP_DEG)
    ins, out = masks(d, sl)
    th = GR['theta']
    u = (x - GR['xc']) * np.cos(th) + (y - GR['yc']) * np.sin(th)
    v = -(x - GR['xc']) * np.sin(th) + (y - GR['yc']) * np.cos(th)
    h = dr / 2.0
    inner = (u / (GR['a'] - h)) ** 2 + (v / (GR['b'] - h)) ** 2
    outer = (u / (GR['a'] + h)) ** 2 + (v / (GR['b'] + h)) ** 2
    sh = (outer <= 1.0) & (inner > 1.0)
    ns, no = int((sh & ins).sum()), int((sh & out).sum())
    ang = np.mod(np.arctan2(v[sh & out] / GR['b'], u[sh & out] / GR['a']), 2 * np.pi)
    sect = np.histogram(ang, bins=np.linspace(0, 2 * np.pi, 17))[0]
    E = ns / (R * no) - 1.0 if no > 0 else None
    return dict(N=ns, N_oth=no, E=E, E_err=((ns / (R * no)) * np.sqrt(1 / max(ns, 1) + 1 / no)) if no > 0 else None, sectors=sect.astype(int))


def giant_ring(t, sl, R, nodes):
    fra, fdec = inverse_gnomonic(RING_RA, RING_DEC, 0.0, GR['frame_dy'])
    out = dict(frame_centre=[round(fra, 3), round(fdec, 3)])
    for dr in (40.0, 120.0):
        g = shell_stats(t, fra, fdec, sl, R, dr)
        rows = []
        for ra0, dec0 in nodes:
            if sep_deg(ra0, dec0, fra, fdec) < 25.0:
                continue
            rows.append((ra0, dec0, shell_stats(t, ra0, dec0, sl, R, dr)))
        allsec = np.concatenate([r[2]['sectors'] for r in rows] + [g['sectors']])
        thr = 0.5 * float(np.median(allsec[allsec > 0]))
        g['covered'] = bool((g['sectors'] >= thr).all())
        kept = [dict(ra=round(a, 3), dec=b, E=r['E']) for a, b, r in rows if r['N_oth'] > 0 and (r['sectors'] >= thr).all()]
        g['sectors'] = g['sectors'].tolist()
        res = dict(shell=g, null_n=len(kept))
        if len(kept) >= 10 and g['E'] is not None:
            E = np.array([k['E'] for k in kept])
            res.update(null_mean=float(E.mean()), null_sd=float(E.std(ddof=1)), rank_p_upper=rank_upper(g['E'], E), null_nodes=kept)
        out['dr_%d' % dr] = res
    return out


def summarize(ev):
    """Keep results.json small: ring in full, nodes as (ra, dec, E_max, E_nom, disc E)."""
    return dict(median_sector=ev['median_sector'], threshold=ev['threshold'], ring=ev['ring'], null=ev['null'],
                nodes=[[k['ra'], k['dec'], round(k['E_max'], 4), round(k['E_nom'], 4), round(k['disc']['E'], 4) if k['disc']['E'] is not None else None]
                       for k in ev['nodes']])


def main():
    t0 = time.time()
    t, info = load()
    print('sample', info, flush=True)
    sl = SLICE
    R, nsl, nout = global_ratio(t, sl)
    nodes = lattice()
    res = dict(generated=time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime()), rule='Workbench section 169, frozen in f3374e3', sample=info,
               ring_nominal=[RING_RA, RING_DEC], annulus_mpc=list(ANN), search_offsets_mpc=OFFS.tolist(), windows_mpc=windows(sl)[1],
               R=dict(value=R, slice=nsl, outer=nout), lattice_nodes=len(nodes))
    ev = evaluate(t, sl, R, nodes)
    print('primary done', round(time.time() - t0), 's', flush=True)
    ring, nl = ev['ring'], ev['null']
    res['primary'] = summarize(ev)
    if not ring['covered'] or ring['N_oth'] == 0 or nl['n_nodes'] < 30:
        res['decision'] = dict(outcome='not tested', why='ring covered=%s, null nodes=%d' % (ring['covered'], nl['n_nodes']))
    else:
        Em = np.array([k['E_max'] for k in ev['nodes']])
        E, m, sd = ring['E_max'], nl['m'], nl['sd']
        z = (E - m) / sd; above_all = bool(E > Em.max())
        if E >= 0.20 and above_all:
            out = 'excess supported (independence from the discovery sightlines not yet established)'
        elif z < 2:
            out = 'no excess'
        else:
            out = 'unresolved'
        res['decision'] = dict(outcome=out, E_max=E, best_offset=ring['best_offset'],
                               best_centre=list(inverse_gnomonic(RING_RA, RING_DEC, *ring['best_offset'])),
                               E_nom=ring['E_nom'], E_nom_err_poisson=(1 + ring['E_nom']) * float(np.sqrt(1 / max(ring['N'], 1) + 1 / ring['N_oth'])),
                               null_m=m, null_sd=sd, z_vs_null=float(z), above_every_node=above_all, rank_p_upper=rank_upper(E, Em),
                               null_m1=nl['m1'], null_sd1=nl['sd1'], z_nom=float((ring['E_nom'] - nl['m1']) / nl['sd1']),
                               limit95_nom=float(ring['E_nom'] + 1.645 * nl['sd1']))
    print(json.dumps(res['decision'], indent=1), flush=True)
    chk = {}
    # (a) annulus widths
    for name, ann in (('annulus_160_240', (160.0, 240.0)), ('annulus_120_280', (120.0, 280.0))):
        e2 = evaluate(t, sl, R, nodes, ann)
        Em2 = [k['E_max'] for k in e2['nodes']]
        chk[name] = dict(ring=dict((k, e2['ring'][k]) for k in ('N', 'N_oth', 'E_nom', 'E_max', 'best_offset', 'covered')), null=e2['null'],
                         z_vs_null=(e2['ring']['E_max'] - e2['null']['m']) / e2['null']['sd'], rank_p_upper=rank_upper(e2['ring']['E_max'], Em2))
    print('check a done', round(time.time() - t0), 's', flush=True)
    # (b) neighbouring slices
    for name, s2 in (('slice_below_0630_0742', (0.630, 0.742)), ('slice_above_0862_0991', (0.862, 0.991))):
        R2, n2, o2 = global_ratio(t, s2)
        e2 = evaluate(t, s2, R2, nodes)
        Em2 = [k['E_max'] for k in e2['nodes']]
        chk[name] = dict(R=R2, windows_mpc=windows(s2)[1], ring=dict((k, e2['ring'][k]) for k in ('N', 'N_oth', 'E_nom', 'E_max', 'best_offset', 'covered')),
                         null=e2['null'], z_vs_null=(e2['ring']['E_max'] - e2['null']['m']) / e2['null']['sd'], rank_p_upper=rank_upper(e2['ring']['E_max'], Em2))
    print('check b done', round(time.time() - t0), 's', flush=True)
    # (c) disc and (d) shape, at the best centre
    discs = np.array([k['disc']['E'] for k in ev['nodes'] if k['disc']['E'] is not None])
    chk['disc_r260'] = dict(ring=ring['disc'], null_n=len(discs), null_mean=float(discs.mean()), null_sd=float(discs.std(ddof=1)),
                            z_vs_null=float((ring['disc']['E'] - discs.mean()) / discs.std(ddof=1)), rank_p_upper=rank_upper(ring['disc']['E'], discs))
    ins_ = np.array([k['interior']['E'] for k in ev['nodes'] if k['interior']['E'] is not None])
    sur_ = np.array([k['surround']['E'] for k in ev['nodes'] if k['surround']['E'] is not None])
    chk['shape'] = dict(annulus_E=ring['E_max'], interior=ring['interior'], surround=ring['surround'],
                        null_interior_mean_sd=[float(ins_.mean()), float(ins_.std(ddof=1))], null_surround_mean_sd=[float(sur_.mean()), float(sur_.std(ddof=1))])
    # (e) Giant Ring
    chk['giant_ring'] = giant_ring(t, sl, R, nodes)
    print('check e done', round(time.time() - t0), 's', flush=True)
    res['checks'] = chk
    res['seconds'] = round(time.time() - t0)

    def clean(o):
        if isinstance(o, dict):
            return {k: clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [clean(v) for v in o]
        if isinstance(o, np.ndarray):
            return clean(o.tolist())
        if isinstance(o, (np.floating, float)):
            return None if not np.isfinite(o) else float(o)
        if isinstance(o, np.integer):
            return int(o)
        return o
    json.dump(clean(res), open(os.path.join(HERE, 'results.json'), 'w'), indent=1)
    print('written results.json', res['seconds'], 's')


if __name__ == '__main__':
    main()
