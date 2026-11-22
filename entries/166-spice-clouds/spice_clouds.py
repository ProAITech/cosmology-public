"""Entry 166, step 2: the frozen cloud-current test on SPICE-RACS DR2 (Workbench section 166, rule frozen in commit 1955eae, 14 Nov 2026).

Input: work/compact.npz from extract.py (columns of data/spice-racs-dr2/spice-racs.dr2.fits.gz, CSIRO DAP csiro:64891, the SPICE-RACS DR2 table:
l, b, rm, rm_err, source_id and the four quality flags snr_flag, leakage_flag, channel_flag, stokesI_fit_flag; True = flagged, so a source is kept
only when all four are False).  Sources with several rows are merged per source_id by weighted mean (weights 1/rm_err^2).
Output: results.json.   Run from the repository root:  .venv/bin/python -I entries/166-spice-clouds/spice_clouds.py

Rule, as frozen (the section numbers are those of the Workbench list):
 2  Spines: polylines in Galactic (l, b) for Orion A (420 pc) and Perseus (300 pc); flat frame x = (l - l0) cos b0, y = b - b0 at the arc-length midpoint.
 3  As amended (Workbench Amendment 1, commit c069c9f, before any analysis): chord-frame coordinates. t = distance along the chord from the first vertex to the
    last, n = offset perpendicular to it (positive to the left); the spine is n = h(t), piecewise linear through the vertices, continued parallel to the
    chord beyond each end.  s = t, d = n - h(t).  Kept: -1 <= s <= Lc + 1 (Lc = chord length), 0.05 < |d| <= 10 deg, |RM - median| < 5 MAD
    (MAD = median |RM - median|, unscaled); tested only if >= MIN_SOURCES kept (1000 as frozen; 400 after Amendment 2, commit 9a58d95: the windows hold 556 and 476 sources,
    see results_gate.json).  (The first wording, nearest-point distance to the polyline, degenerates far
    from a short kinked polyline: s clumps at a few vertices and one side of a segment can be empty.)
 4  WLS with weights 1/(rm_err^2 + sig^2), sig^2 re-estimated in 6 passes (sig^2 = max(0, sum r^2/(N - p) - mean rm_err^2)):
    RM = c0 + c1 d + c2 s + A_in sgn(d)[|d|<1.5] + A_mid sgn(d)[1.5<=|d|<5] + A_out sgn(d)[5<=|d|<=10];  s* = sign(A_in + A_mid); a = A s*.
 5  Null: the same shape and orientation, centre shifted by dl = +-10..60 (step 2.5), db = -4, 0, +4 deg; dropped if its midpoint is within 15 deg of either
    cloud's real midpoint or it holds fewer than half the sources of the real window; p_out = (1 + #{null a_out >= real}) / (1 + N); >= 40 positions.
 6  Pinned: four equal parts of the chord length, step position d0 scanned -6..+6 (0.25), sigma fixed from the fit without a step; pinned if |d0| <= 1;
    cloud pinned if >= 3 of 4 segments are pinned with the same sign of A and p_pin <= 0.05 (rank of the same-sign count among the null positions).
 7  Decision for each cloud separately.
 8  Checks: (a) replication of Entries 161-162 on straight chords; (b) curved against straight; (c) orientation scan; (d) quadratic background.
"""
import json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D2R = np.pi / 180.0
ZIN, ZMID, DMAX, DMIN = 1.5, 5.0, 10.0, 0.05
NPASS = 6
MIN_SOURCES = 400            # Amendment 2 (frozen value: 1000)
D0GRID = np.arange(-6.0, 6.0001, 0.25)

CLOUDS = {
    'ori': dict(name='Orion A', dist_pc=420, pc_per_deg=7.3,
                poly=[(208.60, -19.10), (208.99, -19.38), (210.06, -19.60), (210.93, -19.45), (212.63, -19.00), (214.50, -19.60)],
                nvss=dict(Ain=9.7, eAin=9.7, Aout=23.8, eAout=6.6, d0=0.25)),
    'per': dict(name='Perseus', dist_pc=300, pc_per_deg=5.2,
                poly=[(158.35, -20.50), (160.49, -17.80), (160.57, -16.80)],
                nvss=dict(Ain=20.2, eAin=8.6, Aout=34.7, eAout=9.5, d0=-1.0)),
}


# ---------------------------------------------------------------- data
def load():
    z = np.load(os.path.join(HERE, 'work', 'compact.npz'))
    fl = z['flags']
    ok = (fl == 0).all(axis=1) & np.isfinite(z['rm']) & np.isfinite(z['rm_err']) & (z['rm_err'] > 0) & np.isfinite(z['l']) & np.isfinite(z['b'])
    out = {'rows': int(len(fl)), 'rows_flags_clear': int(((fl == 0).all(axis=1)).sum()), 'rows_usable': int(ok.sum())}
    l, b, rm, e, sid = z['l'][ok], z['b'][ok], z['rm'][ok], z['rm_err'][ok], z['source_id'][ok]
    u, inv, cnt = np.unique(sid, return_inverse=True, return_counts=True)
    out['unique_sources'] = int(len(u)); out['multi_row_sources'] = int((cnt > 1).sum())
    if len(u) < len(sid):                                     # merge components per source by weighted mean
        w = 1.0 / e ** 2
        sw = np.bincount(inv, weights=w)
        rm = np.bincount(inv, weights=w * rm) / sw
        first = np.zeros(len(u), dtype=int); first[inv[::-1]] = np.arange(len(sid))[::-1]
        l, b = l[first], b[first]
        e = 1.0 / np.sqrt(sw)
    out['sources'] = int(len(l))
    return dict(l=np.mod(l, 360.0), b=b, rm=rm, e=e), out


# ---------------------------------------------------------------- geometry
def make_frame(poly_lb):
    """Local flat frame at the arc-length midpoint. Returns (l0, b0, V) with V the polyline in the frame (degrees)."""
    P = np.array(poly_lb, float)
    lm, bm = P[:, 0].mean(), P[:, 1].mean()
    V0 = np.c_[(P[:, 0] - lm) * np.cos(bm * D2R), P[:, 1] - bm]
    seg = np.hypot(*np.diff(V0, axis=0).T); cum = np.r_[0, np.cumsum(seg)]
    half = cum[-1] / 2; i = np.searchsorted(cum, half) - 1; t = (half - cum[i]) / seg[i]
    mid = V0[i] + t * (V0[i + 1] - V0[i])
    b0 = bm + mid[1]; l0 = lm + mid[0] / np.cos(bm * D2R)
    V = np.c_[((P[:, 0] - l0) * np.cos(b0 * D2R)), P[:, 1] - b0]
    return l0, b0, V


def polyline_length(V):
    return float(np.hypot(*np.diff(V, axis=0).T).sum())


def rotate(V, theta_deg):
    c, s = np.cos(theta_deg * D2R), np.sin(theta_deg * D2R)
    return V @ np.array([[c, s], [-s, c]])


def project(V, l0, b0, D, ext=1.0, box=26.0):
    """Chord-frame coordinates (Amendment 1). Sources near (l0, b0): returns s = distance along the chord from the first vertex, d = n - h(t) with
    n the offset to the left of the chord and h the piecewise-linear spine through the vertices (continued parallel to the chord beyond the ends),
    and the index into D.  Sources with s outside [-ext, Lc + ext] or |d| outside (DMIN, DMAX] are dropped."""
    dl = (D['l'] - l0 + 180.0) % 360.0 - 180.0
    x = dl * np.cos(b0 * D2R); y = D['b'] - b0
    m = (np.abs(x) < box) & (np.abs(y) < box)
    idx = np.nonzero(m)[0]; x, y = x[idx], y[idx]
    c = V[-1] - V[0]; Lc = float(np.hypot(*c)); e1 = c / Lc; e2 = np.array([-e1[1], e1[0]])
    T = (V - V[0]) @ e1; N = (V - V[0]) @ e2
    assert np.all(np.diff(T) > 0), 'polyline must advance along its chord'
    t = (x - V[0, 0]) * e1[0] + (y - V[0, 1]) * e1[1]
    n = (x - V[0, 0]) * e2[0] + (y - V[0, 1]) * e2[1]
    d = n - np.interp(t, T, N)                                  # np.interp holds the end values constant: parallel continuation
    keep = (t >= -ext) & (t <= Lc + ext) & (np.abs(d) > DMIN) & (np.abs(d) <= DMAX)
    return t[keep], d[keep], idx[keep]


# ---------------------------------------------------------------- fitting
def design(d, s, zones=True, quad=False, step=None):
    cols = [np.ones_like(d), d, s]
    if quad:
        cols += [d * d, d * s, s * s]
    if zones:
        sg = np.sign(d); ad = np.abs(d)
        cols += [sg * (ad < ZIN), sg * ((ad >= ZIN) & (ad < ZMID)), sg * (ad >= ZMID)]
    if step is not None:
        cols.append(np.sign(d - step))
    return np.column_stack(cols)


def wls(X, y, e, sig2=None, passes=NPASS):
    """Weighted least squares with sig^2 re-estimated from the residuals; returns beta, cov, sig2, w, chi2."""
    n, p = X.shape
    s2 = 0.0 if sig2 is None else sig2
    for _ in range(passes if sig2 is None else 1):
        w = 1.0 / (e ** 2 + s2)
        XtW = X.T * w
        beta = np.linalg.solve(XtW @ X, XtW @ y)
        r = y - X @ beta
        if sig2 is None:
            s2 = max(0.0, float((r ** 2).sum() / max(n - p, 1) - np.mean(e ** 2)))
    w = 1.0 / (e ** 2 + s2)
    XtW = X.T * w
    cov = np.linalg.inv(XtW @ X)
    beta = cov @ (XtW @ y); r = y - X @ beta
    return beta, cov, s2, w, float((w * r ** 2).sum())


def chi2_fixed(X, y, w):
    XtW = X.T * w
    beta = np.linalg.solve(XtW @ X, XtW @ y)
    r = y - X @ beta
    return float((w * r ** 2).sum()), beta


def clip5(rm):
    med = np.median(rm); mad = np.median(np.abs(rm - med))
    return np.abs(rm - med) < 5 * mad


def window(V, l0, b0, D):
    s, d, idx = project(V, l0, b0, D)
    rm, e = D['rm'][idx], D['e'][idx]
    k = clip5(rm)
    return dict(s=s[k], d=d[k], rm=rm[k], e=e[k], idx=idx[k], L=float(np.hypot(*(V[-1] - V[0]))))


def zone_fit(Wn, quad=False):
    X = design(Wn['d'], Wn['s'], quad=quad)
    beta, cov, s2, w, chi2 = wls(X, Wn['rm'], Wn['e'])
    j = 6 if quad else 3
    A = beta[j:j + 3]; eA = np.sqrt(np.diag(cov)[j:j + 3])
    sstar = 1.0 if (A[0] + A[1]) >= 0 else -1.0
    Xr = design(Wn['d'], Wn['s'], zones=False, quad=quad)
    chi2_red, _ = chi2_fixed(Xr, Wn['rm'], w)
    ad = np.abs(Wn['d'])
    return dict(N=int(len(Wn['d'])), sig=float(np.sqrt(s2)), A=A, eA=eA, sstar=sstar, a=A * sstar, improvement=chi2_red - chi2,
                n_zone=[int((ad < ZIN).sum()), int(((ad >= ZIN) & (ad < ZMID)).sum()), int((ad >= ZMID).sum())], beta=beta, w=w)


def step_scan(Wn, sel=None):
    """Best step position d0 in -6..6, sigma fixed from the no-step fit; returns best d0, signed amplitude there, delta chi2 profile."""
    d, s, y, e = Wn['d'], Wn['s'], Wn['rm'], Wn['e']
    if sel is not None:
        d, s, y, e = d[sel], s[sel], y[sel], e[sel]
    X0 = design(d, s, zones=False)
    b0_, cov0, s2, w, chi0 = wls(X0, y, e)
    prof = np.empty(len(D0GRID)); amp = np.empty(len(D0GRID))
    for i, d0 in enumerate(D0GRID):
        st = np.sign(d - d0)
        if (st > 0).sum() < 5 or (st < 0).sum() < 5:           # a step with (almost) nothing on one side cannot be fitted
            prof[i] = -np.inf; amp[i] = 0.0; continue
        X1 = np.column_stack([X0, st])
        c1, beta = chi2_fixed(X1, y, w)
        prof[i] = chi0 - c1; amp[i] = beta[-1]
    if not np.isfinite(prof).any():
        return None, 0.0, prof
    j = int(np.argmax(prof))
    return float(D0GRID[j]), float(amp[j]), prof


def pinned_count(Wn):
    """Four equal-arc-length segments; returns (k = most common sign count among pinned segments, per-segment details)."""
    L = Wn['L']; edges = np.linspace(0, L, 5); det = []
    for q in range(4):
        lo, hi = edges[q], edges[q + 1]
        sel = (Wn['s'] >= lo if q > 0 else Wn['s'] >= -1e9) & (Wn['s'] < hi if q < 3 else Wn['s'] <= 1e9)
        if sel.sum() < 40:
            det.append(dict(n=int(sel.sum()), d0=None, A=None, pinned=False, sign=0)); continue
        d0, A, prof = step_scan(Wn, sel)
        if d0 is None:
            det.append(dict(n=int(sel.sum()), d0=None, A=None, pinned=False, sign=0)); continue
        pin = abs(d0) <= 1.0 + 1e-9
        det.append(dict(n=int(sel.sum()), d0=d0, A=A, dchi2=float(prof.max()), pinned=bool(pin), sign=int(np.sign(A)) if pin else 0))
    npos = sum(1 for t in det if t['pinned'] and t['sign'] > 0); nneg = sum(1 for t in det if t['pinned'] and t['sign'] < 0)
    return max(npos, nneg), det


def ang_sep(l1, b1, l2, b2):
    a1, a2 = np.array([l1, b1]) * D2R, np.array([l2, b2]) * D2R
    c = np.sin(a1[1]) * np.sin(a2[1]) + np.cos(a1[1]) * np.cos(a2[1]) * np.cos(a1[0] - a2[0])
    return float(np.degrees(np.arccos(np.clip(c, -1, 1))))


def rank_up(null, v):
    null = np.asarray(null); return float((1 + (null >= v).sum()) / (1 + len(null)))


def rank_lo(null, v):
    null = np.asarray(null); return float((1 + (null <= v).sum()) / (1 + len(null)))


def hist(v, lo, hi, w):
    c, _ = np.histogram(v, bins=np.arange(lo, hi + w, w)); return dict(lo=lo, width=w, counts=[int(x) for x in c])


def r4(x):
    return None if x is None else float(np.round(x, 4))


# ---------------------------------------------------------------- one cloud
def run_cloud(key, D, mids):
    C = CLOUDS[key]; t0 = time.time()
    l0, b0, V = make_frame(C['poly'])
    mids[key] = (l0, b0)
    other = [k for k in CLOUDS if k != key]
    real = window(V, l0, b0, D)
    out = dict(name=C['name'], dist_pc=C['dist_pc'], pc_per_deg=C['pc_per_deg'], midpoint_lb=[r4(l0), r4(b0)], chord_length_deg=r4(real['L']), arc_length_deg=r4(polyline_length(V)),
               vertices_frame=[[r4(a), r4(b)] for a, b in V])
    if len(real['d']) < MIN_SOURCES:
        out['decision'] = 'not tested: fewer than %d sources in the window' % MIN_SOURCES; out['N'] = int(len(real['d'])); return out, {}
    zf = zone_fit(real); zq = zone_fit(real, quad=True)
    k_real, det = pinned_count(real)
    # chord geometry (checks a and b)
    chord = np.array([V[0], V[-1]])
    chordW = window(chord, l0, b0, D); chordZ = zone_fit(chordW)
    # ------ null positions
    thetas = list(range(-90, 91, 15)); i0 = thetas.index(0)
    nulls = []
    dls = [s * v for v in np.arange(10, 60.0001, 2.5) for s in (1, -1)]
    for dl in dls:
        for db in (-4.0, 0.0, 4.0):
            ln, bn = (l0 + dl) % 360.0, b0 + db
            if abs(bn) > 89 or any(ang_sep(ln, bn, *mids_other) < 15.0 for mids_other in OTHER_MIDS[key]):
                continue
            Wn = window(V, ln, bn, D)
            if len(Wn['d']) < 0.5 * len(real['d']) or len(Wn['d']) < 100:
                continue
            z = zone_fit(Wn); zq_ = zone_fit(Wn, quad=True)
            k, _ = pinned_count(Wn)
            Cw = window(chord, ln, bn, D)
            zc = zone_fit(Cw) if len(Cw['d']) >= 100 else None
            imp_n = np.array([zone_fit(window(rotate(V, th), ln, bn, D))['improvement'] for th in thetas])      # EXPLORATORY, not in the frozen rule
            nulls.append(dict(dl=dl, db=db, N=z['N'], a=z['a'], aq=zq_['a'], k=k, curved_minus_straight=(z['improvement'] - zc['improvement']) if zc else None,
                              orient_rank=int(1 + (imp_n > imp_n[i0]).sum()), orient_excess=float(imp_n[i0] - np.delete(imp_n, i0).max())))
    nn = len(nulls)
    na = np.array([n['a'] for n in nulls]) if nn else np.zeros((0, 3))
    naq = np.array([n['aq'] for n in nulls]) if nn else np.zeros((0, 3))
    nk = np.array([n['k'] for n in nulls])
    out['N'] = zf['N']; out['n_zone'] = zf['n_zone']; out['sig_int'] = r4(zf['sig'])
    out['A'] = [r4(v) for v in zf['A']]; out['eA_formal'] = [r4(v) for v in zf['eA']]; out['sstar'] = int(zf['sstar'])
    out['a'] = [r4(v) for v in zf['a']]; out['a_quadratic_background'] = [r4(v) for v in zq['a']]
    out['null'] = dict(n_positions=nn, a_in=dict(mean=r4(na[:, 0].mean()) if nn else None, sd=r4(na[:, 0].std(ddof=1)) if nn > 1 else None),
                       a_mid=dict(mean=r4(na[:, 1].mean()) if nn else None, sd=r4(na[:, 1].std(ddof=1)) if nn > 1 else None),
                       a_out=dict(mean=r4(na[:, 2].mean()) if nn else None, sd=r4(na[:, 2].std(ddof=1)) if nn > 1 else None,
                                  values=[r4(v) for v in na[:, 2]], q=[r4(v) for v in np.quantile(na[:, 2], [0.025, 0.16, 0.5, 0.84, 0.975])] if nn else None),
                       k_counts=[int((nk == i).sum()) for i in range(5)],
                       window_counts=[n['N'] for n in nulls])
    if nn >= 40:
        out['p'] = dict(a_in_up=r4(rank_up(na[:, 0], zf['a'][0])), a_mid_up=r4(rank_up(na[:, 1], zf['a'][1])), a_out_up=r4(rank_up(na[:, 2], zf['a'][2])),
                        a_out_lo=r4(rank_lo(na[:, 2], zf['a'][2])), p_pin=r4(rank_up(nk, k_real)),
                        a_out_up_quadratic=r4(rank_up(naq[:, 2], zq['a'][2])))
    out['pinned'] = dict(k=int(k_real), segments=[{kk: (r4(v) if isinstance(v, float) else v) for kk, v in t.items()} for t in det])
    # ------ decision (as frozen)
    if nn < 40:
        out['decision'] = 'not tested: fewer than 40 null positions'
    else:
        p = out['p']; pinned = (k_real >= 3) and p['p_pin'] <= 0.05
        a_out = zf['a'][2]
        if not pinned:
            dec = 'not pinned'
        elif a_out > 0 and p['a_out_up'] <= 0.05:
            dec = 'current-like, pinned out to tens of parsecs'
        elif (0.20 < p['a_out_up'] < 0.80) or (a_out < 0 and p['a_out_lo'] <= 0.05):
            dec = 'bent-field-like'
        else:
            dec = 'unresolved'
        out['decision'] = dec
        out['quadratic_background_a_out_positive'] = bool(zq['a'][2] > 0)
    # ------ checks
    chk = {}
    # (a) replication of Entries 161-162, on the straight chord, two zones
    ad = np.abs(chordW['d']); sg = np.sign(chordW['d'])
    X2 = np.column_stack([np.ones_like(ad), chordW['d'], chordW['s'], sg * (ad < ZIN), sg * (ad >= ZIN)])
    be, cv, s2, w, _ = wls(X2, chordW['rm'], chordW['e'])
    d0b, ampb, profb = step_scan(chordW)
    chk['replication'] = dict(N=int(len(ad)), A_in=r4(be[3]), eA_in=r4(np.sqrt(cv[3, 3])), A_out=r4(be[4]), eA_out=r4(np.sqrt(cv[4, 4])),
                              best_d0=r4(d0b), best_amp=r4(ampb), nvss=C['nvss'],
                              same_sign_both=bool(be[3] > 0 and be[4] > 0), d0_within_1p5_of_nvss=bool(d0b is not None and abs(d0b - C['nvss']['d0']) <= 1.5))
    # (b) curved against straight
    dif = zf['improvement'] - chordZ['improvement']
    ndif = np.array([n['curved_minus_straight'] for n in nulls if n['curved_minus_straight'] is not None])
    chk['curved_vs_straight'] = dict(improvement_curved=r4(zf['improvement']), improvement_straight=r4(chordZ['improvement']), difference=r4(dif),
                                     p_up=r4(rank_up(ndif, dif)) if len(ndif) >= 40 else None, n_null=int(len(ndif)))
    # (c) orientation
    imps = []
    for th in thetas:
        Wt = window(rotate(V, th), l0, b0, D); imps.append(zone_fit(Wt)['improvement'])
    imps = np.array(imps); i0 = thetas.index(0)
    rk = int(1 + (imps > imps[i0]).sum()); ex = float(imps[i0] - np.delete(imps, i0).max()); nr = np.array([n['orient_rank'] for n in nulls]); ne = np.array([n['orient_excess'] for n in nulls])
    chk['orientation'] = dict(theta=thetas, improvement=[r4(v) for v in imps], rank_of_true=rk, of=len(thetas), excess_over_best_rotated=r4(ex))
    if nn >= 40:                                                # EXPLORATORY calibration, added after the results were seen; not part of the frozen rule
        chk['orientation']['exploratory_null'] = dict(n_positions=nn, rank_counts=[int((nr == i).sum()) for i in range(1, 14)],
                                                      p_rank_le_real=r4((1 + (nr <= rk).sum()) / (1 + nn)), p_excess_ge_real=r4(rank_up(ne, ex)),
                                                      null_excess_quantiles=[r4(v) for v in np.quantile(ne, [0.16, 0.5, 0.84, 0.95])])
    out['checks'] = chk
    # ------ profile for the figure: residual RM (after c0 + c2 s of the zone fit) in 1-degree bins of d
    be = zf['beta']; res = real['rm'] - be[0] - be[2] * real['s']
    edges = np.arange(-DMAX, DMAX + 0.001, 1.0); prof = []
    for a, b_ in zip(edges[:-1], edges[1:]):
        m = (real['d'] >= a) & (real['d'] < b_)
        if m.sum() < 3:
            prof.append(None); continue
        ww = 1.0 / (real['e'][m] ** 2 + zf['sig'] ** 2)
        prof.append([r4((a + b_) / 2), int(m.sum()), r4((ww * res[m]).sum() / ww.sum()), r4(1 / np.sqrt(ww.sum()))])
    out['profile_binned'] = dict(columns=['d_centre_deg', 'n', 'residual_rm_weighted_mean', 'formal_error'], rows=prof)
    out['seconds'] = round(time.time() - t0)
    return out, dict(a=zf['a'])


OTHER_MIDS = {}


def main():
    t0 = time.time()
    D, cnt = load()
    print('data', cnt, flush=True)
    mids = {}
    for k in CLOUDS:
        mids[k] = make_frame(CLOUDS[k]['poly'])[:2]
    for k in CLOUDS:
        OTHER_MIDS[k] = [mids[o] for o in CLOUDS if o != k] + [mids[k]]      # a null must stay 15 deg from both real midpoints
    res = dict(generated=time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime()), rule='Workbench section 166: frozen in 1955eae, Amendment 1 in c069c9f, Amendment 2 in 9a58d95', counts=cnt, clouds={})
    for k in CLOUDS:
        res['clouds'][k], _ = run_cloud(k, D, mids)
        c = res['clouds'][k]
        print(k, c.get('decision'), 'N', c.get('N'), 'a', c.get('a'), 'p', c.get('p'), 'pinned k', c.get('pinned', {}).get('k'), flush=True)
    res['seconds'] = round(time.time() - t0)
    json.dump(res, open(os.path.join(HERE, 'results.json'), 'w'), indent=1)
    print('written results.json', res['seconds'], 's')


if __name__ == '__main__':
    main()
