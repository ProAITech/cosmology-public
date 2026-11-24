"""Entry 167, follow-up rule (Workbench section 167, "Two further counts", fixed 23 November 2026, 02:38 UTC, commit ba0886c,
published as Workbench version 201 before either count was made).  Asked for by the review of the Entry 167 manuscript (issue #2).

A. Strong absorbers.  Entry 167's sample with EW_2796 >= 0.6 A (deciding), the same band, outer windows, R (recomputed, field left out)
   and the 52 lattice nodes Entry 167 kept; item 6's decision.  Checks reported only: EW_2796 >= 1.0 A, and >= 0.6 A in the slice just
   below (0.630-0.742).
B. The published members.  The 44 main-group positions of ga_members.py against DESI DR1 quasars (ZWARN = 0, 0.862 < Z <= 3.1) within 1';
   chance matches from positions moved 0.3 deg in 100 random directions (seed 20261122); recovered if the Entry 167 sample (any strength)
   has an absorber at 0.742 <= Z_MGII <= 0.862 on that quasar (TARGETID, or within 1"); base rate p0 over the field's quasars.

Data:
  data/desi-dr1-mgii/MgII-Absorbers-DR1.fits (65,750,400 bytes, sha256 38c2e628...bddf; as desi_arc.py), adds column EW_2796.
  data/desi-dr1-qso/QSO_cat_iron_cumulative_v0.fits (DESI DR1 quasar catalogue, 816,203,520 bytes; DATA.md gives its source and
  checksum), table QSO_CAT: TARGETID, TARGET_RA, TARGET_DEC, Z, ZWARN.
Writes followup.json.  Run from the repository root:  .venv/bin/python -I entries/167-desi-arc/followup.py
"""
import ast, json, os, sys, time
import numpy as np
from astropy.io import fits
from scipy.spatial import cKDTree
from scipy.stats import binom
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import desi_arc as A

ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
QSO = os.path.join(ROOT, 'data', 'desi-dr1-qso', 'QSO_cat_iron_cumulative_v0.fits')
D2R = np.pi / 180.0
CUT, CUT_CHECK = 0.6, 1.0
BELOW = (0.630, 0.742)
SEED, NDIR, MOVE_DEG = 20261122, 100, 0.3
MATCH_ARCMIN, NEAR_ARCMIN, SAME_ARCSEC = 1.0, 10.0, 1.0


def chord(arcmin):
    return 2 * np.sin(0.5 * arcmin / 60 * D2R)


def unit(ra, dec):
    ra, dec = np.asarray(ra) * D2R, np.asarray(dec) * D2R
    return np.column_stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


# ---------------------------------------------------------------- the Entry 167 sample, with EW_2796 carried along
def load_ew():
    """desi_arc.load() line for line, with EW_2796 added; checked against desi_arc.load() below."""
    with fits.open(A.PATH, memmap=True) as h:
        d = h[1].data
        t = {k: np.array(d[k]) for k in ('TARGETID', 'RA', 'DEC', 'ZWARN', 'Z_QSO', 'Z_MGII', 'LINE_SNR_MIN', 'TSNR2_QSO', 'EW_2796')}
    ok = (t['ZWARN'] == 0) & (t['Z_QSO'] > 0.862) & (t['Z_QSO'] <= 3.1) & (t['LINE_SNR_MIN'] >= 3)
    t = {k: v[ok] for k, v in t.items()}
    n1 = len(t['RA'])
    order = np.lexsort((-t['TSNR2_QSO'], -t['LINE_SNR_MIN']))
    t = {k: v[order] for k, v in t.items()}
    keep = np.ones(n1, bool)
    xyz = unit(t['RA'], t['DEC'])
    pairs = cKDTree(xyz).query_pairs(r=2 * np.sin(0.5 * 1.0 / 3600 * D2R), output_type='ndarray')
    same_pos = pairs[np.abs(t['Z_MGII'][pairs[:, 0]] - t['Z_MGII'][pairs[:, 1]]) <= 0.001]
    idt = np.argsort(t['TARGETID'], kind='stable'); tid = t['TARGETID'][idt]
    same_t = []; start = 0
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
    t['dist'] = A.dist(t['Z_MGII'])
    t['sra'], t['cra'] = np.sin(t['RA'] * D2R), np.cos(t['RA'] * D2R)
    t['sdec'], t['cdec'] = np.sin(t['DEC'] * D2R), np.cos(t['DEC'] * D2R)
    return t


def subset(t, m):
    return {k: v[m] for k, v in t.items()}


# ---------------------------------------------------------------- part A
def kept_nodes():
    """The 52 lattice nodes Entry 167 kept (results.json), at the exact positions desi_arc.lattice() gives."""
    res = json.load(open(os.path.join(HERE, 'results.json')))
    want = {(n['ra'], n['dec']) for n in res['primary']['null']['nodes']}
    nodes = [(ra, dec) for ra, dec in A.lattice() if (round(ra, 3), dec) in want]
    assert len(nodes) == len(want) == 52, (len(nodes), len(want))
    return nodes, res


def count(t, sl, nodes):
    R, nsl, nout = A.global_ratio(t, sl)
    nb, no, _ = A.band_stats(t, A.RA_C, A.DEC_C, sl)
    rows = []
    for ra0, dec0 in nodes:
        b, o, _ = A.band_stats(t, ra0, dec0, sl)
        if o > 0:
            rows.append(dict(ra=round(ra0, 3), dec=dec0, N_band=b, N_oth=o, E=b / (R * o) - 1.0))
    E = nb / (R * no) - 1.0
    nullE = np.array([r['E'] for r in rows]); m, sd = float(nullE.mean()), float(nullE.std(ddof=1))
    out = dict(sample=int(len(t['RA'])), R=dict(value=R, slice=nsl, outer=nout), arc=dict(N_band=nb, N_oth=no, N_exp=R * no, E=E,
               E_err_poisson=(1 + E) * float(np.sqrt(1.0 / max(nb, 1) + 1.0 / no))),
               null=dict(n_nodes=len(rows), mean=m, sd=sd, max=float(nullE.max()), nodes=rows),
               z_vs_null=(E - m) / sd, above_every_null=bool(E > nullE.max()),
               rank_p_upper=float((1 + (nullE >= E).sum()) / (1 + len(nullE))), limit95=E + 1.645 * sd,
               z_vs_entry66=(E - A.PRIOR_E) / np.sqrt(sd ** 2 + A.PRIOR_ERR ** 2))
    return out


def decide(c):
    if c['null']['n_nodes'] < 40:
        return 'not tested (fewer than 40 null nodes)'
    if c['arc']['E'] >= 0.20 and c['above_every_null']:
        return 'excess supported'
    return 'no excess' if c['z_vs_null'] < 2 else 'unresolved'


# ---------------------------------------------------------------- part B
def members():
    src = open(os.path.join(HERE, 'ga_members.py'), encoding='utf-8').read()
    raw = next(ast.literal_eval(n.value) for n in ast.walk(ast.parse(src))
               if isinstance(n, ast.Assign) and getattr(n.targets[0], 'id', None) == 'RAW')
    rows = [r.split() for r in raw.split(';')]
    B = [(float(r[1]), float(r[2])) for r in rows if r[0] == 'B']
    assert len(B) == 44
    return np.array(B)


def load_qso():
    with fits.open(QSO, memmap=True) as h:
        d = h[1].data
        q = {k: np.array(d[k]) for k in ('TARGETID', 'TARGET_RA', 'TARGET_DEC', 'Z', 'ZWARN')}
    n0 = len(q['Z'])
    m = (q['ZWARN'] == 0) & (q['Z'] > 0.862) & (q['Z'] <= 3.1)
    q = {k: v[m] for k, v in q.items()}
    n1 = len(q['Z'])
    # one entry per position: quasars within 1" of an earlier one are dropped (the catalogue can list one object under two TARGETIDs)
    xyz = unit(q['TARGET_RA'], q['TARGET_DEC'])
    pairs = cKDTree(xyz).query_pairs(r=chord(SAME_ARCSEC / 60), output_type='ndarray')
    keep = np.ones(n1, bool); keep[pairs.max(axis=1)] = False
    q = {k: v[keep] for k, v in q.items()}
    return q, xyz[keep], dict(rows=int(n0), after_cuts=int(n1), same_position_dropped=int(n1 - keep.sum()), quasars=int(keep.sum()))


def move(ra, dec, deg, bearing):
    """Point at angular distance deg along the given bearing (radians, from north through east)."""
    d, a1, p1 = deg * D2R, ra * D2R, dec * D2R
    p2 = np.arcsin(np.sin(p1) * np.cos(d) + np.cos(p1) * np.sin(d) * np.cos(bearing))
    a2 = a1 + np.arctan2(np.sin(bearing) * np.sin(d) * np.cos(p1), np.cos(d) - np.sin(p1) * np.sin(p2))
    return (a2 / D2R) % 360, p2 / D2R


def part_b(t):
    M = members()
    q, qxyz, qinfo = load_qso()
    tree = cKDTree(qxyz)
    mx = unit(M[:, 0], M[:, 1])
    dd, ii = tree.query(mx)
    sep = 2 * np.arcsin(dd / 2) / D2R * 60                                         # arcmin
    near10 = np.array([len(tree.query_ball_point(v, chord(NEAR_ARCMIN))) for v in mx])
    rng = np.random.default_rng(SEED)
    pc = []
    for ra, dec in M:
        b = rng.uniform(0, 2 * np.pi, NDIR)
        r2, d2 = move(ra, dec, MOVE_DEG, b)
        dm, _ = tree.query(unit(r2, d2))
        pc.append(float((dm <= chord(MATCH_ARCMIN)).mean()))
    pc = np.array(pc); Nc = float(pc.sum())
    matched = sep <= MATCH_ARCMIN
    Nm = int(matched.sum())
    # slice absorbers of the Entry 167 sample (any strength)
    sl = (t['Z_MGII'] >= A.SLICE[0]) & (t['Z_MGII'] <= A.SLICE[1])
    s = subset(t, sl)
    stree = cKDTree(unit(s['RA'], s['DEC']))
    tid2rows = {}
    for k, tid in enumerate(s['TARGETID']):
        tid2rows.setdefault(int(tid), []).append(k)

    def absorbers_on(j):
        rows = set(tid2rows.get(int(q['TARGETID'][j]), [])) | set(stree.query_ball_point(qxyz[j], chord(SAME_ARCSEC / 60)))
        return sorted(rows)
    per = []
    for i in range(44):
        j = int(ii[i]); row = dict(ra=float(M[i, 0]), dec=float(M[i, 1]), sep_arcmin=round(float(sep[i]), 3), quasars_within_10arcmin=int(near10[i]),
                                    chance_match=round(float(pc[i]), 3), matched=bool(matched[i]))
        if matched[i]:
            ab = absorbers_on(j)
            row.update(qso_targetid=int(q['TARGETID'][j]), qso_z=round(float(q['Z'][j]), 4),
                       absorbers=[dict(z_mgii=round(float(s['Z_MGII'][k]), 4), ew_2796=round(float(s['EW_2796'][k]), 3),
                                       line_snr_min=round(float(s['LINE_SNR_MIN'][k]), 2)) for k in ab],
                       recovered=bool(len(ab) > 0))
        per.append(row)
    Nr = int(sum(r.get('recovered', False) for r in per))
    # descriptive, not part of the decision: where the members sit relative to the band (|x| <= 500, |y| <= 100 Mpc)
    mx_, my_, _ = A.gnomonic(dict(sra=np.sin(M[:, 0] * D2R), cra=np.cos(M[:, 0] * D2R), sdec=np.sin(M[:, 1] * D2R), cdec=np.cos(M[:, 1] * D2R)),
                             A.RA_C, A.DEC_C)
    for r, x_, y_ in zip(per, mx_, my_):
        r['x_mpc'], r['y_mpc'] = round(float(x_), 1), round(float(y_), 1)
        r['in_band'] = bool(abs(x_) <= A.BAND_HX and abs(y_) <= A.BAND_HY)
    ib = [r for r in per if r['in_band']]
    in_band = dict(members=len(ib), matched=sum(r['matched'] for r in ib), recovered=sum(bool(r.get('recovered')) for r in ib),
                   recovered_ew_ge_06=sum(any(a['ew_2796'] >= CUT for a in r.get('absorbers', [])) for r in ib))
    # base rate over the field's quasars
    qx, qy, qc = A.gnomonic(dict(sra=np.sin(q['TARGET_RA'] * D2R), cra=np.cos(q['TARGET_RA'] * D2R), sdec=np.sin(q['TARGET_DEC'] * D2R),
                                 cdec=np.cos(q['TARGET_DEC'] * D2R)), A.RA_C, A.DEC_C)
    infield = np.nonzero((np.abs(qx) <= A.FIELD_HX) & (np.abs(qy) <= A.FIELD_HY) & (qc > 0.5))[0]
    carry = np.array([len(absorbers_on(j)) > 0 for j in infield])
    p0 = float(carry.mean())
    pbin = float(binom.sf(Nr - 1, Nm, p0)) if Nm else None
    if Nm - Nc < 10:
        outcome = 'not decided (N_m - N_c < 10)'
    elif Nr / Nm >= 0.5 and pbin < 0.001:
        outcome = 'members recovered'
    elif pbin > 0.05:
        outcome = 'not recovered'
    else:
        outcome = 'partly recovered'
    return dict(qso=qinfo, N_members=44, N_m=Nm, N_c=round(Nc, 2), N_r=Nr, fraction=(Nr / Nm if Nm else None),
                field_quasars=int(len(infield)), field_quasars_with_slice_absorber=int(carry.sum()), p0=p0, binomial_p=pbin,
                chance_corrected_fraction=((Nr - Nc * p0) / (Nm - Nc) if Nm - Nc > 0 else None),
                members_without_quasar_within_10arcmin=int((near10 == 0).sum()), in_band=in_band,
                recovered_ew_2796_median=float(np.median([a['ew_2796'] for r in per for a in r.get('absorbers', [])])) if Nr else None,
                sep_quantiles_arcmin=[float(v) for v in np.quantile(sep, [0.1, 0.25, 0.5, 0.75, 0.9])],
                outcome=outcome, members=per)


def main():
    t0 = time.time()
    t = load_ew()
    ref, info = A.load()
    assert np.array_equal(ref['TARGETID'], t['TARGETID']) and np.array_equal(ref['Z_MGII'], t['Z_MGII']), 'sample differs from desi_arc.load()'
    nodes, res167 = kept_nodes()
    out = dict(generated=time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime()),
               rule='Workbench section 167, follow-up rule fixed 23 November 2026 02:38 UTC (commit ba0886c), published as version 201',
               sample=int(len(t['RA'])), ew_nonfinite=int((~np.isfinite(t['EW_2796'])).sum()))
    # check: no cut, Entry 167's nodes, must give Entry 167's numbers
    base = count(t, A.SLICE, nodes)
    p = res167['primary']
    out['reproduce_entry167'] = dict(E=base['arc']['E'], E_167=p['arc']['E'], null_mean=base['null']['mean'], null_mean_167=p['null']['mean'],
                                     null_sd=base['null']['sd'], null_sd_167=p['null']['sd'])
    assert abs(base['arc']['E'] - p['arc']['E']) < 1e-12 and abs(base['null']['sd'] - p['null']['sd']) < 1e-12, out['reproduce_entry167']
    print('Entry 167 reproduced: E %+.4f, null %+.4f +- %.4f' % (base['arc']['E'], base['null']['mean'], base['null']['sd']), flush=True)
    fin = np.isfinite(t['EW_2796'])
    ga = json.load(open(os.path.join(HERE, 'ga_members.json')))['bands']['200']['cases']
    for name, c, sl in (('A_ew06', CUT, A.SLICE), ('check_ew10', CUT_CHECK, A.SLICE), ('check_ew06_slice_below', CUT, BELOW)):
        tc = subset(t, fin & (t['EW_2796'] >= c))
        r = count(tc, sl, nodes)
        r['cut_ew2796'] = c; r['slice'] = list(sl)
        if sl == A.SLICE:
            r['below_published_members_sd'] = [dict(field=k['field'], delta=k['delta'], predicted_E=k['predicted_E'],
                                                    sd_below=(k['predicted_E'] - (r['arc']['E'] - r['null']['mean'])) / r['null']['sd']) for k in ga]
        if name == 'A_ew06':
            r['decision'] = decide(r)
        out[name] = r
        print('%-24s n=%6d  N_band %3d  N_exp %6.1f  E %+.3f  null %+.3f +- %.3f (%d nodes)  z %+.2f  rank p %.3f  95%% limit %+.3f%s' % (
            name, r['sample'], r['arc']['N_band'], r['arc']['N_exp'], r['arc']['E'], r['null']['mean'], r['null']['sd'], r['null']['n_nodes'],
            r['z_vs_null'], r['rank_p_upper'], r['limit95'], ('  -> ' + r['decision']) if 'decision' in r else ''), flush=True)
    out['B_members'] = b = part_b(t)
    print('B: matched %d of 44 within 1\' (chance %.1f); recovered %d; p0 %.3f over %d field quasars; binomial p %.2g; members with no quasar within 10\': %d -> %s' % (
        b['N_m'], b['N_c'], b['N_r'], b['p0'], b['field_quasars'], b['binomial_p'] if b['binomial_p'] is not None else float('nan'),
        b['members_without_quasar_within_10arcmin'], b['outcome']))
    out['seconds'] = round(time.time() - t0)
    json.dump(out, open(os.path.join(HERE, 'followup.json'), 'w'), indent=1)
    print('written followup.json', out['seconds'], 's')


if __name__ == '__main__':
    main()
