"""Entry 183: does the Faraday growth with quasar redshift survive at fixed radio luminosity?

Entries 171-174 found the RM scatter of quasar sight lines growing with redshift (union of 10,449 sight lines: t = +3.54).
Entry 181 showed intervening haloes cannot make it.  In a flux-limited sample the more distant quasars are more luminous, so
the growth could come from the sources: more luminous radio quasars in denser, more magnetised surroundings.  This test
separates the two.  Entry 174's union statistic is recomputed with the redshifts shuffled within cells of Galactic latitude
band x RM-error tercile x radio-luminosity quintile (60 cells) instead of the first two alone (12 cells).  A shuffle inside a
luminosity quintile keeps each redshift with sources of similar luminosity, so the null slope takes up whatever growth
luminosity makes; t then measures the growth beyond it.

Luminosity: L = 4 pi D_L^2 S (1+z)^(-1-alpha) at 888 MHz in the rest frame, alpha = -0.7, S the sum of total_I_flux over the
source's usable SPICE-RACS rows, Planck 2018 distances (H0 67.66, Om 0.3111, flat).
Data: as Entry 174 (header of quaia_growth.py, copied unchanged): SPICE-RACS DR2, Milliquas v8, DESI DR1 quasars, Quaia.
Run from the repository root:
  .venv/bin/python -I entries/183-faraday-luminosity/lum_growth.py synthetic   # simulated skies only; reads no data
  .venv/bin/python -I entries/183-faraday-luminosity/lum_growth.py run         # the test; writes results.json
"""
import gzip, json, os, sys, time
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.spatial import cKDTree
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import quaia_growth as q

ALPHA = -0.7
NQ = 5
H0, OM = 67.66, 0.3111
_zg = np.linspace(0, 6, 6001)
_dc = cumulative_trapezoid(299792.458 / H0 / np.sqrt(OM * (1 + _zg) ** 3 + 1 - OM), _zg, initial=0)     # Mpc


def lum(S_jy, z):
    """log10 L (W/Hz) at 888 MHz rest frame."""
    dl = np.interp(z, _zg, _dc) * (1 + z) * 3.0857e22                                                  # m
    return np.log10(4 * np.pi * dl ** 2 * S_jy * 1e-26 * (1 + z) ** (-1 - ALPHA))


def quintile(x, k=NQ):
    edges = np.quantile(x, np.linspace(0, 1, k + 1)[1:-1])
    return np.searchsorted(edges, x, side='right')


def analyse_cells(Sx, cell, nshuf, label, seed=q.CFG['seed'], log=print):
    """Entry 174's explore_union.analyse_cells: the slope b, and its null from redshifts shuffled within the given cells."""
    z, d = Sx['z'], Sx['d']
    real = q.bins_fit(z, d, np.random.default_rng(seed), Sx['e'])
    groups = [np.nonzero(cell == c)[0] for c in np.unique(cell)]
    R = np.random.default_rng(seed + 7)
    nulls = np.empty(nshuf)
    for s in range(nshuf):
        z2 = z.copy()
        for g in groups:
            z2[g] = z[g][R.permutation(len(g))]
        nulls[s] = q.bins_fit(z2, d, np.random.default_rng(seed + 1000 + s))['b']
    mu, sd = float(nulls.mean()), float(nulls.std(ddof=1))
    out = dict(label=label, n=int(len(z)), cells=len(groups), b=real['b'], null_mean=mu, null_sd=sd, t=float((real['b'] - mu) / sd),
               net=float(real['b'] - mu), p_upper=float((1 + (nulls >= real['b']).sum()) / (1 + nshuf)),
               bins=dict(zmed=real['zmed'], sigma=real['sigma'], err=real['err'], e2mean=real['e2mean'], n=real['n']))
    log('%s: n %d, cells %d, b %+.3f, null %+.3f +- %.3f, t %+.2f' % (label, out['n'], out['cells'], out['b'], mu, sd, out['t']))
    return out


def decide_183(t12, tL, Rz, RL):
    """Section 183, item 5.  Rz, RL: tL of the calibration skies with growth from redshift only and from luminosity only."""
    if t12 < 3:
        return 'not tested: the union growth is below t = 3 with this code'
    beyond = tL > np.percentile(RL, 97.5)          # luminosity alone rarely leaves this much
    lumi = tL < np.percentile(Rz, 2.5)             # redshift alone rarely leaves this little
    if beyond and lumi:
        return 'both: more growth than luminosity alone leaves at fixed luminosity, less than redshift alone leaves'
    if beyond:
        return 'the growth is not a luminosity effect: at fixed radio luminosity it stays above what luminosity alone leaves'
    if lumi:
        return 'the growth is a luminosity effect: at fixed radio luminosity it falls below what a redshift growth leaves'
    return 'unresolved: at fixed radio luminosity the growth is consistent with both'


GROWTH = 19.0          # noise-removed variance growth between median z 0.52 and 2.65 on the union (Entry 174)
S0SQ = 8.2 ** 2 - 4.0


def extra_variance(model, z, L):
    x = np.log10(1 + z)
    xr = np.log10(1.52), np.log10(3.65)
    ez = GROWTH * np.clip((x - xr[0]) / (xr[1] - xr[0]), 0, 1)
    if model == 'redshift':
        return ez
    lo, hi = np.quantile(L, [0.1, 0.9])
    el = np.clip((L - lo) / (hi - lo), 0, 1)
    qz = np.quantile(z, [0.25, 0.75])
    scale = (ez[z >= qz[1]].mean() - ez[z <= qz[0]].mean()) / (el[z >= qz[1]].mean() - el[z <= qz[0]].mean())
    return scale * el                               # matched: same rise between the lowest and highest redshift quarters


def _sim_one(args):
    model, k, z, L, e, b, nshuf, seed = args
    g = np.random.default_rng(seed)
    d = g.normal(0, 1, len(z)) * np.sqrt(S0SQ + extra_variance(model, z, L)) + g.normal(0, 1, len(z)) * e
    Sx = dict(z=z, d=d, e=e, b=b)
    base = q.cells_for(Sx)
    quiet = lambda s: None
    r12 = analyse_cells(Sx, base, nshuf, '', seed=seed, log=quiet)
    rL = analyse_cells(Sx, base * NQ + quintile(L), nshuf, '', seed=seed, log=quiet)
    return model, k, r12['t'], rL['t'], r12['net'], rL['net']


def calibrate(z, L, e, b, nsims=40, nshuf=100, seed=18300):
    """Section 183, item 4: tL (and t12) of skies with the sample's own z, L, e, b and simulated RM, growth from redshift or
    luminosity only, the same rise in variance between the lowest and highest redshift quarters."""
    from multiprocessing import Pool
    jobs = [(m, k, z, L, e, b, nshuf, seed + 1000 * j + k) for j, m in enumerate(('redshift', 'luminosity')) for k in range(nsims)]
    with Pool(8) as pool:
        res = pool.map(_sim_one, jobs)
    out = {}
    for m in ('redshift', 'luminosity'):
        r = [x for x in res if x[0] == m]
        out[m] = dict(t12=[x[2] for x in r], tL=[x[3] for x in r], net12=[x[4] for x in r], netL=[x[5] for x in r])
        tl = np.array(out[m]['tL'])
        print('calibration, growth from %s only: t12 %.2f +- %.2f, tL %.2f +- %.2f (2.5-97.5%%: %.2f to %.2f), share %.2f' % (
            m, np.mean(out[m]['t12']), np.std(out[m]['t12']), tl.mean(), tl.std(), np.percentile(tl, 2.5), np.percentile(tl, 97.5),
            np.mean(out[m]['netL']) / np.mean(out[m]['net12'])), flush=True)
    return out


# ---------------------------------------------------------------- simulated skies (no data)
def synthetic(nshuf=200):
    """Toy skies (growth 19 from luminosity, from redshift, or none), each judged against its own calibration; t12 is floored at 3
    so that the decision step is exercised."""
    out = {}
    for model in ('luminosity', 'redshift', 'none'):
        g = np.random.default_rng(183 + len(model))
        n = 10000
        z = np.clip(g.lognormal(np.log(1.3), 0.5, 4 * n), 0.1, 4.0)
        S = 0.008 * (1 - g.random(4 * n)) ** (-1 / 1.5)                    # N(>S) ~ S^-1.5 above 8 mJy
        L = lum(S, z)
        keep = np.arange(n)
        z, S, L = z[keep], S[keep], L[keep]
        b = g.uniform(20, 90, n) * g.choice([-1, 1], n)
        e = g.lognormal(np.log(2.0), 0.4, n)
        x = np.log10(1 + z)
        extra = extra_variance(model, z, L) if model != 'none' else np.zeros(n)
        d = g.normal(0, 1, n) * np.sqrt(S0SQ + extra) + g.normal(0, 1, n) * e
        Sx = dict(z=z, d=d, e=e, b=b)
        base = q.cells_for(Sx)
        r12 = analyse_cells(Sx, base, nshuf, 'synthetic %s, 12 cells' % model)
        rL = analyse_cells(Sx, base * NQ + quintile(L), nshuf, 'synthetic %s, 60 cells (luminosity)' % model)
        cal = calibrate(z, L, e, b, nsims=30, nshuf=80, seed=18300 + 100 * len(model))
        out[model] = dict(corr_logL_logz=float(np.corrcoef(L, x)[0, 1]), t12=r12['t'], tL=rL['t'], net12=r12['net'], netL=rL['net'],
                          calibration=cal, decision=decide_183(max(r12['t'], 3.0), rL['t'], cal['redshift']['tL'], cal['luminosity']['tL']))
        print('  -> %s' % out[model]['decision'], flush=True)
    json.dump(out, open(os.path.join(HERE, 'synthetic.json'), 'w'), indent=1)


# ---------------------------------------------------------------- the data
def extract_flux():
    """Stream SPICE-RACS once with Entry 171's row selection; save total_I_flux in the same row order as spice_rm.npz."""
    t0 = time.time()
    fl, ra = [], []
    with gzip.open(q.SPICE, 'rb') as f:
        q.read_header(f)
        h = q.read_header(f)
        dt, rowbytes = q.table_dtype(h)
        nrows = int(q.card(h, 'NAXIS2'))
        done = 0
        while done < nrows:
            n = min(100_000, nrows - done)
            a = np.frombuffer(f.read(n * rowbytes), dtype=dt)
            clear = np.ones(n, bool)
            for k in q.FLAGS:
                clear &= a[k] == b'F'
            rm, e = a['rm'].astype('<f8'), a['rm_err'].astype('<f8')
            ok = clear & np.isfinite(rm) & np.isfinite(e) & (e > 0)
            fl.append(a['total_I_flux'][ok].astype('<f8')); ra.append(a['ra'][ok].astype('<f8'))
            done += n
    os.makedirs(os.path.join(HERE, 'work'), exist_ok=True)
    np.savez(os.path.join(HERE, 'work', 'flux.npz'), flux=np.concatenate(fl), ra=np.concatenate(ra))
    print('flux extracted, %d rows, %.0f s' % (sum(len(x) for x in fl), time.time() - t0), flush=True)


def union():
    """Entry 174's union (explore_union.py), with each sight line's source index kept."""
    rows, src, _ = q.load_spice()
    quaia, _ = q.load_quaia(); mq = q.load_milliquas(); desi = q.load_desi()
    tree = cKDTree(q.unit(src['ra'], src['dec']))
    used = np.zeros(len(src['rm']), bool)
    m_mq = q.match(mq, rows); used[m_mq['src']] = True
    m_de = q.match(desi, rows); used[m_de['src']] = True
    S, _ = q.build(q.match(quaia, rows), src, quaia, q.CFG['fg_radius'], q.CFG['min_neigh'], tree, exclude=used)
    isQ = mq['cls'] == 'Q'
    mq_Q = {k: v[isQ] for k, v in mq.items()}
    M, _ = q.build(q.match(mq_Q, rows), src, mq_Q, q.CFG['fg_radius'], q.CFG['min_neigh'], tree)
    Dd, _ = q.build(m_de, src, desi, q.CFG['fg_radius'], q.CFG['min_neigh'], tree)
    donly = ~np.isin(Dd['src'], M['src'])
    keys = ('z', 'd', 'e', 'b', 'l', 'src')
    U = {k: np.r_[M[k], Dd[k][donly], S[k]] for k in keys}
    U['sample'] = np.r_[np.zeros(len(M['z']), int), np.ones(int(donly.sum()), int), np.full(len(S['z']), 2)]
    return U, rows


def run():
    t0 = time.time()
    U, rows = union()
    fp = os.path.join(HERE, 'work', 'flux.npz')
    if not os.path.exists(fp):
        extract_flux()
    F = np.load(fp)
    assert len(F['flux']) == len(rows['row_src']) and np.allclose(F['ra'], rows['row_ra'])
    Ssrc = np.bincount(rows['row_src'], weights=np.where(np.isfinite(F['flux']) & (F['flux'] > 0), F['flux'], 0.0))
    U['S'] = Ssrc[U['src']]
    ok = U['S'] > 0
    U = {k: v[ok] for k, v in U.items()}
    U['L'] = lum(U['S'], U['z'])
    base = q.cells_for(U)
    res = dict(n=int(len(U['z'])), dropped_no_flux=int((~ok).sum()), corr_logL_logz=float(np.corrcoef(U['L'], np.log10(1 + U['z']))[0, 1]),
               L_quintile_edges=[float(v) for v in np.quantile(U['L'], [0.2, 0.4, 0.6, 0.8])], median_S_mJy=float(1e3 * np.median(U['S'])))
    r12 = analyse_cells(U, base, 1000, 'union, 12 cells (Entry 174)')
    rL = analyse_cells(U, base * NQ + quintile(U['L']), 1000, 'union, 60 cells (radio luminosity quintiles)')
    cal = calibrate(U['z'], U['L'], U['e'], U['b'])
    res.update(t12=r12, tL=rL, calibration=cal, decision=decide_183(r12['t'], rL['t'], cal['redshift']['tL'], cal['luminosity']['tL']),
               share_beyond_luminosity=rL['net'] / r12['net'])
    print('decision:', res['decision'], '| net slope %.2f -> %.2f (share %.2f)' % (r12['net'], rL['net'], res['share_beyond_luminosity']), flush=True)
    # reported only (section 183, item 6)
    rep = {}
    rep['flux_quintiles'] = analyse_cells(U, base * NQ + quintile(np.log10(U['S'])), 300, 'union, 60 cells (flux-density quintiles)')
    rep['L_deciles'] = analyse_cells(U, base * 10 + quintile(U['L'], 10), 300, 'union, 120 cells (luminosity deciles)')
    for k in range(NQ):
        sel = quintile(U['L']) == k
        Sx = {kk: v[sel] for kk, v in U.items()}
        rep['L_quintile_%d' % k] = analyse_cells(Sx, q.cells_for(Sx), 300, 'luminosity quintile %d alone (z median %.2f)' % (k, np.median(Sx['z'])))
    # scatter against luminosity at fixed redshift: luminosity shuffled within redshift quintiles x the 12 cells
    Lsw = dict(U); Lsw['z'] = U['L'] - U['L'].min() + 0.01                 # the same binned statistic with log L in place of z
    rep['L_at_fixed_z'] = analyse_cells(Lsw, base * NQ + quintile(U['z']), 300, 'scatter against log L, shuffled within redshift quintiles')
    res['reported'] = rep
    res['seconds'] = round(time.time() - t0)
    json.dump(res, open(os.path.join(HERE, 'results.json'), 'w'), indent=1, default=float)


if __name__ == '__main__':
    import multiprocessing as mp
    mp.freeze_support()
    {'synthetic': synthetic, 'run': run}[sys.argv[1]]()
