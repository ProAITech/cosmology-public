"""Entry 178: has the magnetic field of galaxy clusters changed since redshift 0.4, relative to the thermal pressure of their gas?

Entry 176 found the Faraday-rotation excess of eRASS1 clusters at 0.5-2 R500 and read it as a field of about 2 uG at 0.75 R500,
with thermal pressure about 11 times the magnetic pressure in every mass quarter.  This entry asks whether that ratio is the same
for the nearer and the farther half of the clusters.  The rule was frozen in the Workbench (section 178) before this script read
any rotation measure.

Data: as Entry 176 (cluster_rm.py in this folder is Entry 176's, copied unchanged):
  data/spice-racs-dr2/spice-racs.dr2.fits.gz (read through Entry 171's cache, entries/171-spice-faraday-redshift/work/spice_rm.npz)
  data/erass1-clusters/erass1cl_primary_v3.2.fits   columns RA, DEC, BEST_Z, PCONT, R500, M500, MGAS500.

Method (section 178 has the frozen text):
  Sources   Entry 176's 92,169: SPICE-RACS DR2 at |b| >= 20 deg, dRM = RM - median of the other sources within 3 deg.
  Clusters  Entry 176's 9,050.  Halves: BEST_Z below and at-or-above the median of the 9,050 (0.2953).
  Distance  x = smallest theta / theta500 over all clusters, as in Entry 176; a source belongs to the half of that nearest cluster.
  Excess    In each half and in each of the ranges 0.5-1, 1-1.5 and 1.5-2 R500: V = (1.4826 MAD dRM)^2 - <e^2>, minus its mean over
            500 rotations of the whole cluster list in Galactic longitude (a rotated source takes the half of its rotated nearest
            cluster); the error is the standard deviation over the rotations.
  Model     Each cluster: n_e(r) = n_0 / (1 + r^2/r_c^2), r_c = 0.2 R500, normalised to MGAS500 (Entry 176's model.py); isothermal
            kT = mu m_p G M500 / (2 R500); a random field with reversal length 10 kpc whose energy density is a fixed share of the
            thermal pressure, B^2 / 8 pi = 1.92 n_e kT / beta_pl.  A source behind the cluster at projected radius R has
                sigma_RM^2 = 812^2 Lambda / 3 * Int n_e^2 B^2 dl * (1 + z)^-4,
            members half of that, sources in front nothing (front share F(z) = 1 - exp(-(z/1.085)^1.2), members from the source
            counts, both as in model.py).  Predicted excesses come from the field sources' dRM with the cluster term added, with the
            same random draws for every value of 1/beta_pl, so each half's chi^2 curve is smooth.
  Ratio     rho = B_far / B_near at equal thermal pressure = sqrt(beta_near / beta_far); its 68% and 95% ranges from the profile of
            chi^2_near + chi^2_far (Delta chi^2 = 1 and 3.84).  Because Lambda cancels only if it is the same in both halves, rho
            measures B sqrt(Lambda) at fixed pressure.
  Decision  decide_178() below; frozen text in section 178.

Run from the repository root:
  .venv/bin/python -I entries/178-cluster-field-evolution/field_z.py synthetic   # synthetic sky; reads no data
  .venv/bin/python -I entries/178-cluster-field-evolution/field_z.py run         # the test; writes results.json
"""
import json, os, sys, time
import numpy as np
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cluster_rm as c

KPC, MSUN, MP, MUE, MU, G = 3.0857e21, 1.989e33, 1.6726e-24, 1.17, 0.59, 6.674e-8
CFG = dict(rc_over_r500=0.2, lambda_kpc=10.0, ranges=[[0.5, 1.0], [1.0, 1.5], [1.5, 2.0]], inner=[0.0, 0.5], foot=[4.0, 8.0],
           nrot=500, seed=20261178, nreal=200, zsplit='median of the 9,050 clusters', F=(1.085, 1.2), band=(0.67, 1.5),
           q_grid='1/beta_pl from 0.001 to 2 (B^2 proportional to it)', rho_grid=(0.1, 4.0))
QGRID = np.unique(np.round(np.concatenate([np.arange(0.0, 0.02, 0.001), np.arange(0.02, 0.2, 0.0025), np.arange(0.2, 0.6, 0.01),
                                           np.arange(0.6, 2.0001, 0.05)]), 5))
ALL = [CFG['inner']] + CFG['ranges'] + [CFG['foot']]          # cell order: inner, three fitted ranges, footprint range


# ---------------------------------------------------------------- clusters with gas
def load_clusters():
    from astropy.io import fits
    with fits.open(c.CLUSTERS, memmap=True) as h:
        d = h[1].data
        ra, dec, z = np.asarray(d['RA'], float), np.asarray(d['DEC'], float), np.asarray(d['BEST_Z'], float)
        pc, r500, m500, mgas = (np.asarray(d[k], float) for k in ('PCONT', 'R500', 'M500', 'MGAS500'))
    ok = np.isfinite(ra) & np.isfinite(dec) & np.isfinite(z) & np.isfinite(r500) & (pc < 0.5) & (z >= 0.05) & (r500 > 0)
    l, b = c.galactic(ra[ok], dec[ok])
    cl = dict(l=l, b=b, z=z[ok], r500=r500[ok], m500=m500[ok], mgas=mgas[ok])
    return gas(cl)


def gas(cl):
    """n_0 from the gas mass (median gas fraction where none is listed), theta500, isothermal kT; as Entry 176's model.py."""
    cl = dict(cl)
    cl['th'] = cl['r500'] / 1000.0 / c.ang_diam_Mpc(cl['z']) / c.D2R
    mg = cl['mgas']
    good = np.isfinite(mg) & (mg > 0)
    fgas = np.nanmedian(mg[good] * 1e11 / (cl['m500'][good] * 1e13))
    cl['mgas'] = np.where(good, mg, fgas * cl['m500'] * 100.0)
    rc = CFG['rc_over_r500'] * cl['r500']
    X = 1.0 / CFG['rc_over_r500']
    cl['rc'] = rc
    cl['n0'] = cl['mgas'] * 1e11 * MSUN / (4 * np.pi * MUE * MP * (rc * KPC) ** 3 * (X - np.arctan(X)))
    cl['kT'] = MU * MP * G * cl['m500'] * 1e13 * MSUN / (2 * cl['r500'] * KPC)          # erg
    return cl


def los_n3(cl, ks, R, nstep=400):
    """Int n_e^3 dl (cm^-9 kpc) along the line of sight through cluster ks at projected radius R (kpc), out to 3 R500."""
    rmax = 3.0 * cl['r500'][ks]
    L = np.sqrt(np.maximum(rmax ** 2 - R ** 2, 0.0))
    t = (np.arange(nstep) + 0.5) / nstep
    ll = L[:, None] * t[None, :]
    ne = cl['n0'][ks][:, None] / (1.0 + (R[:, None] ** 2 + ll ** 2) / cl['rc'][ks][:, None] ** 2)
    return 2.0 * (ne ** 3).sum(1) * (L / nstep)


def sigma2_per_q(cl, ks, R):
    """sigma_RM^2 of a background source per unit 1/beta_pl (rad^2 m^-4)."""
    B2_per_ne = 8 * np.pi * 1.92 * cl['kT'][ks] * 1e12          # microgauss^2 per (cm^-3) at 1/beta = 1
    return 812.0 ** 2 * CFG['lambda_kpc'] / 3.0 * B2_per_ne * los_n3(cl, ks, R) * (1.0 + cl['z'][ks]) ** -4


# ---------------------------------------------------------------- cells
def cells(x, which, half_of, d, e):
    """V and count in each (half, range) cell, ranges ALL; half_of maps cluster index -> 0 (near) / 1 (far)."""
    V = np.full((2, len(ALL)), np.nan); n = np.zeros((2, len(ALL)), int)
    h = np.where(which >= 0, half_of[np.maximum(which, 0)], -1)
    for i in range(2):
        for j, (lo, hi) in enumerate(ALL):
            m = (h == i) & (x >= lo) & (x < hi)
            n[i, j] = int(m.sum())
            if m.sum() >= 15:
                V[i, j] = c.V_of(d[m], e[m])
    return V, n


def measure(S, cl, half_of, nrot, seed, log=print):
    tree = cKDTree(c.gal_unit(S['l'], S['b']))
    ns = len(S['d'])
    x, which = c.scaled_distance(tree, ns, cl)
    V, n = cells(x, which, half_of, S['d'], S['e'])
    g = np.random.default_rng(seed)
    Vr = np.empty((nrot, 2, len(ALL))); nr = np.empty((nrot, 2, len(ALL)))
    for i in range(nrot):
        xr, wr = c.scaled_distance(tree, ns, cl, shift=g.uniform(10.0, 350.0))
        Vr[i], nr[i] = cells(xr, wr, half_of, S['d'], S['e'])
        if (i + 1) % 100 == 0:
            log('  rotations %d' % (i + 1))
    mu, sd = np.nanmean(Vr, 0), np.nanstd(Vr, 0, ddof=1)
    return dict(x=x, which=which, V=V, n=n, null_mean=mu, null_sd=sd, null_n=nr.mean(0), excess=V - mu, t=(V - mu) / sd)


# ---------------------------------------------------------------- model fit
def fit_half(S, cl, M, i, field, g, F=CFG['F']):
    """chi^2 (over the three fitted ranges) against 1/beta for half i, with common random numbers."""
    dF, eF, Vf = field
    x, which = M['x'], M['which']
    h = np.where(which >= 0, M['half_of'][np.maximum(which, 0)], -1)
    foot = M['n'][:, -1].sum() / M['null_n'][:, -1].sum()          # share of rotated clusters still on the survey
    rows = []
    pred = np.zeros((len(QGRID), 3)); psd = np.zeros((len(QGRID), 3))
    for j, (lo, hi) in enumerate(CFG['ranges']):
        m = np.nonzero((h == i) & (x >= lo) & (x < hi))[0]
        ks = which[m]
        R = x[m] * cl['r500'][ks]
        s2 = sigma2_per_q(cl, ks, R)
        pf = 1.0 - np.exp(-(cl['z'][ks] / F[0]) ** F[1])
        expect = M['null_n'][i, j + 1] * foot
        fm = max(0.0, 1.0 - expect / max(len(m), 1))
        k = g.integers(0, len(dF), (CFG['nreal'], len(m)))
        u = g.uniform(size=(CFG['nreal'], len(m)))
        z0 = g.normal(size=(CFG['nreal'], len(m)))
        factor = np.where(u < fm, 0.5, np.where(u < fm + (1 - fm) * pf[None, :], 0.0, 1.0))
        for a, q in enumerate(QGRID):
            dsim = dF[k] + z0 * np.sqrt(s2[None, :] * factor * q)
            vs = np.array([c.V_of(dsim[r], eF[k[r]]) for r in range(CFG['nreal'])])
            pred[a, j] = vs.mean() - Vf; psd[a, j] = vs.std(ddof=1)
        rows.append(dict(range=[lo, hi], sources=int(len(m)), member_share=float(fm), front_share_mean=float(pf.mean()),
                         excess=float(M['excess'][i, j + 1]), error=float(M['null_sd'][i, j + 1]), t=float(M['t'][i, j + 1]),
                         median_sigma2_per_q=float(np.median(s2))))
    obs = M['excess'][i, 1:4]; nsd = M['null_sd'][i, 1:4]
    # The rotations give the error with no cluster field; a field adds its own sample variance, which the draws measure.
    err2 = nsd[None, :] ** 2 + np.maximum(psd ** 2 - psd[:1] ** 2, 0.0)
    chi = (((obs[None, :] - pred) ** 2 / err2) + np.log(err2 / nsd[None, :] ** 2)).sum(1)      # -2 ln L, Gaussian
    return chi, pred, rows


def summarise_half(cl, sel, chi, label):
    a = int(np.argmin(chi))
    in1 = QGRID[chi <= chi[a] + 1.0]
    in2 = QGRID[chi <= chi[a] + 3.84]
    ne075 = float(np.median(cl['n0'][sel])) / (1 + (0.75 / CFG['rc_over_r500']) ** 2)
    kT = float(np.median(cl['kT'][sel]))
    Bq = lambda q: float(np.sqrt(8 * np.pi * 1.92 * ne075 * kT * 1e12 * q))
    out = dict(label=label, clusters=int(sel.sum()), z_median=float(np.median(cl['z'][sel])), m500_median_1e13=float(np.median(cl['m500'][sel])),
               kT_median_keV=kT / 1.602e-9, n_e_0p75_median=ne075, inv_beta=float(QGRID[a]), inv_beta_68=[float(in1.min()), float(in1.max())],
               inv_beta_95=[float(in2.min()), float(in2.max())], chi2=float(chi[a]), dof=2,
               B_0p75=Bq(QGRID[a]), B_0p75_68=[Bq(in1.min()), Bq(in1.max())], beta_pl=(1 / QGRID[a]) if QGRID[a] > 0 else None,
               chi2_zero_field=float(chi[0]))
    return out


def ratio_profile(chi_n, chi_f):
    rho = np.exp(np.linspace(np.log(CFG['rho_grid'][0]), np.log(CFG['rho_grid'][1]), 400))
    qf = lambda q: np.interp(q, QGRID, chi_f, right=np.inf)
    prof = np.array([np.min(chi_n + qf(r * r * QGRID)) for r in rho])
    i0 = int(np.argmin(prof)); m = prof[i0]
    in1, in2 = rho[prof <= m + 1.0], rho[prof <= m + 3.84]
    edge = lambda arr: [float(arr.min()), float(arr.max())]
    return dict(rho=float(rho[i0]), rho_68=edge(in1), rho_95=edge(in2), chi2_min=float(m),
                chi2_at_rho_1=float(np.interp(0.0, np.log(rho), prof)), touches_grid_edge=bool(in2.min() <= rho[0] or in2.max() >= rho[-1]))


def decide_178(r):
    lo, hi = r['rho_95']
    if r['touches_grid_edge']:
        return 'unresolved: the 95% range of rho reaches the edge of the grid'
    if CFG['band'][0] <= lo and hi <= CFG['band'][1]:
        return 'no evolution beyond the band: at equal pressure the field of the farther half is within 0.67-1.5 times the nearer (95%)'
    if hi < 1.0:
        return 'weaker at higher redshift: at equal pressure the farther half has a weaker field (95% range below 1)'
    if lo > 1.0:
        return 'stronger at higher redshift: at equal pressure the farther half has a stronger field (95% range above 1)'
    return 'unresolved: rho = 1 is allowed but the 95% range reaches beyond 0.67-1.5'


def analyse(S, cl, nrot, seed, log=print, variants=True):
    zcut = float(np.median(cl['z']))
    half_of = (cl['z'] >= zcut).astype(int)
    M = measure(S, cl, half_of, nrot, seed, log=log)
    M['half_of'] = half_of
    field = ~np.isfinite(M['x'])
    dF, eF = S['d'][field], S['e'][field]
    fld = (dF, eF, c.V_of(dF, eF))
    out = dict(z_split=zcut, cells=dict(ranges=ALL, n=M['n'].tolist(), V=M['V'].tolist(), null_mean=M['null_mean'].tolist(),
                                        null_sd=M['null_sd'].tolist(), excess=M['excess'].tolist(), t=M['t'].tolist(),
                                        null_n=M['null_n'].tolist()))
    g = np.random.default_rng(seed + 7)
    chis, halves = [], []
    for i, name in ((0, 'near'), (1, 'far')):
        chi, pred, rows = fit_half(S, cl, M, i, fld, g)
        s = summarise_half(cl, half_of == i, chi, name)
        s['ranges'] = rows
        s['predicted_at_best'] = pred[int(np.argmin(chi))].tolist()
        chis.append(chi); halves.append(s)
        log('%s half: %d clusters, z %.3f, M500 %.1fe13, kT %.1f keV; excesses %s +- %s; 1/beta %.4f (68%% %.4f-%.4f), B(0.75 R500) %.2f uG, chi2 %.1f for 2 (zero field %.1f)' % (
            name, s['clusters'], s['z_median'], s['m500_median_1e13'], s['kT_median_keV'], [round(r['excess'], 1) for r in rows],
            [round(r['error'], 1) for r in rows], s['inv_beta'], *s['inv_beta_68'], s['B_0p75'], s['chi2'], s['chi2_zero_field']))
    out['halves'] = halves
    R = ratio_profile(chis[0], chis[1])
    R['decision'] = decide_178(R)
    out['ratio'] = R
    log('rho = B_far / B_near at equal pressure: %.2f (68%% %.2f-%.2f; 95%% %.2f-%.2f); chi2 at rho = 1 minus minimum %.2f' % (
        R['rho'], *R['rho_68'], *R['rho_95'], R['chi2_at_rho_1'] - R['chi2_min']))
    log('decision: %s' % R['decision'])
    if variants:
        rep = {}
        for lab, F in (('front_median_0p6', (0.814, 1.2)), ('front_median_1p0', (1.357, 1.2))):
            cc = [fit_half(S, cl, M, i, fld, np.random.default_rng(seed + 7), F=F)[0] for i in (0, 1)]
            rr = ratio_profile(*cc)
            rep[lab] = dict(rho=rr['rho'], rho_95=rr['rho_95'])
            log('  variant %s: rho %.2f (95%% %.2f-%.2f)' % (lab, rr['rho'], *rr['rho_95']))
        out['reported_only'] = rep
    return out


# ---------------------------------------------------------------- synthetic sky (no data read)
def synthetic():
    res = {}
    for name, q_near, q_far, seed in (('no_cluster_field', 0.0, 0.0, 41), ('same_beta', 0.09, 0.09, 42), ('same_beta_b', 0.09, 0.09, 44),
                                      ('same_beta_c', 0.09, 0.09, 45), ('same_beta_d', 0.09, 0.09, 46), ('far_weaker_rho_0p5', 0.09, 0.0225, 43),
                                      ('far_stronger_rho_1p5', 0.09, 0.2025, 47)):
        g = np.random.default_rng(seed)
        n_src = 120000
        dec = np.arcsin(g.uniform(np.sin(-90 * c.D2R), np.sin(49 * c.D2R), n_src)) / c.D2R; ra = g.uniform(0, 360, n_src)
        l, b = c.galactic(ra, dec)
        fg = 25 / (np.abs(np.sin(b * c.D2R)) + 0.3) * np.sin(l * c.D2R + 0.7) + 8 * np.cos(3 * l * c.D2R) * np.sin(2 * b * c.D2R)
        e = g.lognormal(np.log(1.7), 0.5, n_src)
        rm = fg + g.normal(0, 8, n_src) + g.normal(0, 1, n_src) * e
        ncl = 12000
        cl_l = g.uniform(180, 360, ncl); cl_b = np.arcsin(g.uniform(-1, 1, ncl)) / c.D2R
        keep = np.abs(cl_b) > 18; nk = int(keep.sum())
        z = g.uniform(0.05, 0.6, nk) ** 1.0
        m500 = g.lognormal(np.log(20.0), 0.6, nk) * (1 + z) ** 1.5                     # 1e13 Msun, rising with z as in a flux-limited sample
        Ez2 = 0.3089 * (1 + z) ** 3 + 0.6911
        rho_c = 2.775e11 * 0.6774 ** 2 * Ez2                                            # Msun / Mpc^3
        r500 = (3 * m500 * 1e13 / (4 * np.pi * 500 * rho_c)) ** (1 / 3) * 1000.0      # kpc
        cl = gas(dict(l=cl_l[keep], b=cl_b[keep], z=z, r500=r500, m500=m500, mgas=0.11 * m500 * 100.0))
        tree = cKDTree(c.gal_unit(l, b))
        x, which = c.scaled_distance(tree, n_src, cl)
        sel = np.nonzero(x < 3.0)[0]
        if q_near > 0:
            ks = which[sel]
            s2 = sigma2_per_q(cl, ks, x[sel] * cl['r500'][ks])
            q = np.where(cl['z'][ks] >= np.median(cl['z']), q_far, q_near)
            pf = 1.0 - np.exp(-(cl['z'][ks] / 1.085) ** 1.2)
            behind = g.uniform(size=len(sel)) >= pf
            rm[sel] = rm[sel] + g.normal(0, 1, len(sel)) * np.sqrt(s2 * q * behind)
        src = dict(ra=ra, dec=dec, l=l, b=b, rm=rm, e=e, simple=np.ones(n_src, bool))
        S, cnt = c.sample(src, c.CFG['fg_radius'], c.CFG['min_neigh'])
        out = analyse(S, cl, 60, seed, log=lambda s: print(name, s, flush=True), variants=False)
        res[name] = dict(true_rho=(float(np.sqrt(q_far / q_near)) if q_near > 0 else None), true_inv_beta=[q_near, q_far],
                         fitted_inv_beta=[h['inv_beta'] for h in out['halves']], fitted_inv_beta_68=[h['inv_beta_68'] for h in out['halves']],
                         rho=out['ratio']['rho'], rho_95=out['ratio']['rho_95'], decision=out['ratio']['decision'],
                         t_by_half_range=np.round(np.array(out['cells']['t'])[:, 1:4], 2).tolist(),
                         n_by_half_range=np.array(out['cells']['n'])[:, 1:4].tolist())
    print(json.dumps(res, indent=1))
    json.dump(res, open(os.path.join(HERE, 'synthetic.json'), 'w'), indent=1)


def run():
    t0 = time.time()
    log = lambda s: print(s, flush=True)
    rows, src, sinfo = c.load_spice()
    S, cnt = c.sample(src, c.CFG['fg_radius'], c.CFG['min_neigh'])
    cl = load_clusters()
    log('sources %s; clusters %d' % (cnt, len(cl['z'])))
    res = dict(generated=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), config=CFG, spice=sinfo, sources=cnt, clusters=int(len(cl['z'])))
    res.update(analyse(S, cl, CFG['nrot'], CFG['seed'], log=log))
    res['seconds'] = round(time.time() - t0)
    json.dump(res, open(os.path.join(HERE, 'results.json'), 'w'), indent=1, default=float)
    log('written results.json %d s' % res['seconds'])


if __name__ == '__main__':
    {'synthetic': synthetic, 'run': run}[sys.argv[1]]()
