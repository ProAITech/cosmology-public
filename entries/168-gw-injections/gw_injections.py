"""Entry 168: are gravitational waves stretched? Entry 79's likelihood with the LVK detection probabilities (Workbench section 168, rule frozen in 9747176).

Events: code/gw_rows.txt (Entry 79; GWOSC cumulative GWTC event API, queried 11 November 2026), restricted to O1-O4a (run codes 1-5): 155 binary black holes.
Injections: data/lvk-gwtc4-injections/mixture-semi_o1_o2-real_o3_o4a-cartesian_spins_20250503134659UTC.hdf (LVK GWTC-4.0 sensitivity estimates,
  Zenodo record 16740128, https://zenodo.org/records/16740128; 624,783,664 bytes; sha256 fdaa4743a169fef3a3ca4071b2fb939af3433df9de8fa9973305fad1186090f1),
  table 'events', 1,697,750 accepted injections of 1,130,738,485 generated.
Likelihood: imported unchanged from code/gw_redshift.py (Entry 79): nll, profile, summarise, logf, the 21 knots and the smoothness term; only Pdet changes.
Selection (rule items 2-3): run groups by time (O1+O2 < 1.2e9 GPS s, O3 < 1.3e9, O4a after). Found: O1+O2 semi-analytic network SNR >= 10; O3 and O4a
  false-alarm rate < 1 per year in any of that run's searches. Each found injection weighted by p_ref / p_draw x the file's weight field, with
  p_ref = uniform in (ln Mc_det, ln D) x q^1.1 on [0.1, 1] x isotropic spins with magnitudes uniform to the group's maximum x isotropic inclination,
  Jacobian D'(z) / (D m1^2); histogrammed on a 0.02 grid in ln D and ln Mc_det, smoothed with a Gaussian of 0.08 (rule) in both; scaled to a maximum of 1
  per group (the constant cancels in the likelihood).
Run from the repository root:  .venv/bin/python -I entries/168-gw-injections/gw_injections.py   (a few minutes; seeds fixed)
"""
import json, os, sys, time
import numpy as np
import h5py
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
WORK = os.path.join(HERE, 'work'); os.makedirs(WORK, exist_ok=True)
INJ = os.path.join(ROOT, 'data', 'lvk-gwtc4-injections', 'mixture-semi_o1_o2-real_o3_o4a-cartesian_spins_20250503134659UTC.hdf')
NMOCK = 50
SIG_RULE = 0.08
STEP = 0.02

# ---------------------------------------------------------------- events: Entry 79's table, O1-O4a only, loaded through Entry 79's own module
src = os.path.join(ROOT, 'code', 'gw_rows.txt')
hdr = [l for l in open(src) if l.startswith('#')]
rows = [r for l in open(src) if not l.startswith('#') for r in l.strip().split(';') if r]
keep = [r for r in rows if int(float(r.split(',')[0])) <= 5]
ROWS = os.path.join(WORK, 'gw_rows_o1_o4a.txt')
with open(ROWS, 'w') as f:
    f.writelines(hdr); f.write('# Entry 168: rows of runs 1-5 (O1-O4a) only.\n'); f.write(';'.join(keep) + '\n')
sys.argv = ['gw_redshift.py', ROWS]
sys.path.insert(0, os.path.join(ROOT, 'code'))
import gw_redshift as gw                                   # noqa: E402  (module-level code reads ROWS)
N = gw.N
GRP = np.where(gw.RUN <= 2, 0, np.where(gw.RUN <= 4, 1, 2))
GNAMES = ['O1+O2', 'O3', 'O4a']
ZZ = gw.z_of_D(gw.D)
ALL = np.ones(N, bool)


# ---------------------------------------------------------------- injections
def load_injections():
    f = h5py.File(INJ, 'r')
    e = f['events'][:]
    t = e['time_geocenter']
    g = np.where(t < 1.2e9, 0, np.where(t < 1.3e9, 1, 2))
    found = np.zeros(len(e), bool)
    found[g == 0] = e['semianalytic_observed_phase_maximized_snr_net'][g == 0] >= 10
    for k, pre in ((1, 'o3_'), (2, 'o4a_')):
        cols = [n for n in e.dtype.names if n.startswith(pre) and n.endswith('_far')]
        far = np.vstack([e[c] for c in cols]).T
        far = np.where(np.isfinite(far), far, np.inf).min(axis=1)
        found[g == k] = far[g == k] < 1.0
    m1, m2, z, D = e['mass1_source'], e['mass2_source'], e['redshift'], e['luminosity_distance']
    mc = (m1 * m2) ** 0.6 / (m1 + m2) ** 0.2 * (1 + z)
    q = m2 / m1
    lnh = np.where((q >= 0.1) & (q <= 1.0), 1.1 * np.log(np.clip(q, 1e-9, None)) + np.log(2.1 / (1 - 0.1 ** 2.1)), -np.inf)
    lns = np.zeros(len(e))
    amax = {}
    for k in range(3):
        a1 = np.sqrt(e['spin1x'] ** 2 + e['spin1y'] ** 2 + e['spin1z'] ** 2)
        amax[k] = float(max(a1[g == k].max(), np.sqrt(e['spin2x'] ** 2 + e['spin2y'] ** 2 + e['spin2z'] ** 2)[g == k].max()))
    am = np.array([amax[k] for k in range(3)])[g]
    for s in ('1', '2'):
        a = np.sqrt(e['spin%sx' % s] ** 2 + e['spin%sy' % s] ** 2 + e['spin%sz' % s] ** 2)
        lns += np.where(a <= am, -np.log(4 * np.pi * np.maximum(a, 1e-12) ** 2 * am), -np.inf)
    lninc = np.log(np.sin(e['inclination']) / 2)
    lnref = lnh + np.log(e['dluminosity_distance_dredshift']) - np.log(D) - 2 * np.log(m1) + lns + lninc
    lnw = lnref - e['lnpdraw_mass1_source_mass2_source_redshift_spin1x_spin1y_spin1z_spin2x_spin2y_spin2z'] - e['lnpdraw_inclination'] + np.log(e['weights'])
    info = dict(rows=int(len(e)), found_by_group=[int((found & (g == k)).sum()) for k in range(3)], rows_by_group=[int((g == k).sum()) for k in range(3)],
                spin_max_by_group=[round(amax[k], 4) for k in range(3)])
    return dict(g=g, found=found, lnD=np.log(D), lnM=np.log(mc), lnw=lnw), info


DE = np.arange(np.log(10.0), np.log(30000.0) + STEP, STEP)
ME = np.arange(np.log(0.5), np.log(1000.0) + STEP, STEP)
DC, MC = 0.5 * (DE[1:] + DE[:-1]), 0.5 * (ME[1:] + ME[:-1])


def build_pdet(I, sigma=SIG_RULE, subset=None):
    """Per group: smoothed weighted histogram of found injections -> interpolator of Pdet(ln D, ln Mc_det), max 1. Also N_eff maps."""
    out = []
    for k in range(3):
        m = I['found'] & (I['g'] == k) & np.isfinite(I['lnw'])
        if subset is not None:
            m &= subset
        w = np.exp(I['lnw'][m] - I['lnw'][m].max())
        H = np.histogram2d(I['lnD'][m], I['lnM'][m], bins=[DE, ME], weights=w)[0]
        H2 = np.histogram2d(I['lnD'][m], I['lnM'][m], bins=[DE, ME], weights=w * w)[0]
        sb = sigma / STEP
        S1 = gaussian_filter(H, sb, mode='constant', truncate=4.0)
        S2 = gaussian_filter(H2, sb / np.sqrt(2), mode='constant', truncate=4.0) / (4 * np.pi * sb * sb)
        # N_eff of the kernel sum at a point: (sum w K)^2 / sum (w K)^2, with K the unnormalised Gaussian (peak 1)
        k1 = S1 * (2 * np.pi * sb * sb)                      # sum w K   (gaussian_filter is normalised to unit sum)
        k2 = S2 * (4 * np.pi * sb * sb) * (np.pi * sb * sb)  # sum w^2 K^2: K^2 is a Gaussian of width sb/sqrt2 with peak 1 -> normalised sum x (pi sb^2)
        neff = np.where(k2 > 0, k1 ** 2 / np.maximum(k2, 1e-300), 0.0)
        P = S1 / S1.max()
        out.append(dict(P=RegularGridInterpolator((DC, MC), P, bounds_error=False, fill_value=0.0),
                        NE=RegularGridInterpolator((DC, MC), neff, bounds_error=False, fill_value=0.0), n_found=int(m.sum()),
                        n_eff_total=float(w.sum() ** 2 / (w * w).sum())))
    return out


def pdet(PD, Dv, M, grp, key='P'):
    Dv, M, grp = np.broadcast_arrays(np.asarray(Dv, float), np.asarray(M, float), np.asarray(grp))
    out = np.zeros(Dv.shape)
    for k in range(3):
        s = grp == k
        if s.any():
            out[s] = PD[k][key](np.column_stack([np.log(Dv[s]), np.log(M[s])]))
    return out


def grids(PD, Dv=None):
    Dv = gw.D if Dv is None else Dv
    PG = pdet(PD, Dv[:, None], gw.MG[None, :], GRP[:, None])
    PE = np.maximum(pdet(PD, gw.D, gw.MDET, GRP), 1e-6)
    return PG, PE


def fit(PG, PE, sel=ALL, alphas=gw.AL):
    return gw.profile(alphas, ZZ, PG, PE, sel)


def T_of(prof, alphas=gw.AL):
    return float(np.interp(0.0, alphas, prof) - np.interp(1.0, alphas, prof))


def main():
    t0 = time.time()
    I, iinfo = load_injections()
    print('injections', iinfo, round(time.time() - t0), 's', flush=True)
    PD = build_pdet(I)
    PG, PE = grids(PD)
    heavy, light = gw.MDET > 40, gw.MDET < 60
    res = dict(generated=time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime()), rule='Workbench section 168, frozen in 9747176', injections=iinfo,
               events=dict(n=int(N), by_group={GNAMES[k]: int((GRP == k).sum()) for k in range(3)}, heavy=int(heavy.sum()), light=int(light.sum())),
               pdet_found_used=[PD[k]['n_found'] for k in range(3)], pdet_n_eff_total=[round(PD[k]['n_eff_total'], 1) for k in range(3)])
    # N_eff of the selection estimate at each event's own point
    ne = pdet(PD, gw.D, gw.MDET, GRP, key='NE')
    res['n_eff_at_events'] = dict(min=round(float(ne.min()), 1), median=round(float(np.median(ne)), 1), heavy_min=round(float(ne[heavy].min()), 1),
                                  below_20=int((ne < 20).sum()))
    prof_all = fit(PG, PE); prof_h = fit(PG, PE, heavy); prof_l = fit(PG, PE, light)
    T = T_of(prof_all); h = float(np.interp(1.0, gw.AL, prof_h))
    res['all'] = gw.summarise(prof_all); res['all']['T'] = round(T, 3)
    res['heavy'] = gw.summarise(prof_h); res['heavy']['h'] = round(h, 3)
    res['light'] = gw.summarise(prof_l)
    res['profiles'] = dict(alpha=[float(a) for a in gw.AL], all=[round(float(v), 3) for v in prof_all], heavy=[round(float(v), 3) for v in prof_h],
                           light=[round(float(v), 3) for v in prof_l])
    print('real: T', round(T, 2), 'h', round(h, 2), res['all'], res['heavy'], flush=True)
    # ---- convergence: halves of the found injections by row parity
    idx = np.arange(len(I['g']))
    conv = {}
    for nm, sub in (('even', idx % 2 == 0), ('odd', idx % 2 == 1)):
        PDh = build_pdet(I, subset=sub); PGh, PEh = grids(PDh)
        conv[nm] = dict(T=round(T_of(fit(PGh, PEh)), 3), h=round(float(np.interp(1.0, gw.AL, fit(PGh, PEh, heavy))), 3))
    conv['dT'] = round(abs(conv['even']['T'] - conv['odd']['T']), 3); conv['dh'] = round(abs(conv['even']['h'] - conv['odd']['h']), 3)
    conv['converged'] = bool(conv['dT'] <= 1.0 and conv['dh'] <= 1.0)
    res['convergence'] = conv
    print('convergence', conv, flush=True)
    # ---- mocks (rule item 6): true distance from the quoted interval; mass from the fit under the hypothesis, through the selection at the true distance
    sig = 0.5 * (np.log(gw.D + gw.DHI) - np.log(np.maximum(gw.D - gw.DLO, 1))) / 1.645
    mocks = {}
    MDET0 = gw.MDET.copy()
    for a in (0.0, 1.0):
        th = minimize(gw.nll, np.zeros(len(gw.KNOTS)), args=(a, ZZ, PG, PE, ALL), method='L-BFGS-B', jac=True).x
        Ts, Hs, Ab, Nh = [], [], [], []
        for k in range(NMOCK):
            rg = np.random.default_rng(168000 + 1000 * int(a) + k)
            Dt = gw.D * np.exp(rg.normal(0, 1, N) * sig)
            PGt = pdet(PD, Dt[:, None], gw.MG[None, :], GRP[:, None])
            W = np.exp(gw.logf(np.log(gw.MG)[None, :] - a * np.log1p(gw.z_of_D(Dt))[:, None], th)) * PGt
            cdf = np.cumsum(W, 1); cdf /= cdf[:, -1:]
            Mm = np.array([np.interp(rg.random(), cdf[i], gw.MG) for i in range(N)]) * np.exp(rg.normal(0, 0.05, N))
            gw.MDET[:] = Mm
            PEm = np.maximum(pdet(PD, gw.D, Mm, GRP), 1e-6)
            pm = fit(PG, PEm)
            Ts.append(round(T_of(pm), 3)); Ab.append(float(gw.AL[int(np.argmin(pm))]))
            if a == 1.0:
                hm = Mm > 40; Nh.append(int(hm.sum()))
                Hs.append(round(float(np.interp(1.0, gw.AL, fit(PG, PEm, hm))), 3))
            gw.MDET[:] = MDET0
            if k % 10 == 9:
                print('mocks alpha=%d: %d done, %.0f s' % (a, k + 1, time.time() - t0), flush=True)
        mocks['alpha%d' % int(a)] = dict(T=Ts, alpha_best=Ab, **({'h_heavy': Hs, 'n_heavy': Nh} if a == 1.0 else {}))
    res['mocks'] = mocks
    T0, T1, H1 = np.array(mocks['alpha0']['T']), np.array(mocks['alpha1']['T']), np.array(mocks['alpha1']['h_heavy'])
    dec = dict(T=round(T, 3), unstretched_mocks=dict(mean=round(float(T0.mean()), 3), sd=round(float(T0.std(ddof=1)), 3), max=round(float(T0.max()), 3)),
               stretched_mocks=dict(mean=round(float(T1.mean()), 3), sd=round(float(T1.std(ddof=1)), 3), min=round(float(T1.min()), 3)),
               h=round(h, 3), heavy_stretched_mocks=dict(median=round(float(np.median(H1)), 3), q84=round(float(np.quantile(H1, 0.84)), 3), max=round(float(H1.max()), 3)))
    z0 = (T - T0.mean()) / T0.std(ddof=1); z1 = (T - T1.mean()) / T1.std(ddof=1)
    dec['z_vs_unstretched'] = round(float(z0), 2); dec['z_vs_stretched'] = round(float(z1), 2)
    p_h = (1 + int((H1 >= h).sum())) / (1 + len(H1)); dec['p_h'] = round(p_h, 4)
    if not conv['converged']:
        dec['stretch'] = dec['heavy'] = 'selection not converged'
    else:
        dec['stretch'] = ('stretch confirmed' if (T > T0.max() and z0 >= 3) else 'unstretched favoured' if (T < T1.min() and z1 <= -3) else 'not decided')
        dec['heavy'] = ('the heavy residual survives' if p_h <= 0.02 else 'the heavy residual is gone' if p_h > 0.05 else 'weakened, not decided')
    res['decision'] = dec
    print('decision', dec, flush=True)
    # ---- checks
    chk = {}
    old_all, _ = gw.run_emp('Entry 79 selection, O1-O4a events'); old_h, _ = gw.run_emp('Entry 79 selection, heavy', sel=gw.MDET > 40)
    chk['a_entry79_selection_same_events'] = dict(all=old_all, heavy=old_h, T=round(old_all['dlnL_alpha0'] - old_all['dlnL_alpha1'], 3), h=old_h['dlnL_alpha1'])
    PD2 = build_pdet(I, sigma=0.16); PG2, PE2 = grids(PD2)
    p2 = fit(PG2, PE2); p2h = fit(PG2, PE2, heavy)
    chk['b_smoothing_0p16'] = dict(T=round(T_of(p2), 3), h=round(float(np.interp(1.0, gw.AL, p2h)), 3), all=gw.summarise(p2), heavy=gw.summarise(p2h))
    chk['c_by_group'] = {GNAMES[k]: dict(n=int((GRP == k).sum()), **gw.summarise(fit(PG, PE, GRP == k))) for k in range(3) if (GRP == k).sum() >= 5}
    res['checks'] = chk
    res['seconds'] = round(time.time() - t0)
    json.dump(res, open(os.path.join(HERE, 'results.json'), 'w'), indent=1)
    print('written results.json', res['seconds'], 's')


if __name__ == '__main__':
    main()
