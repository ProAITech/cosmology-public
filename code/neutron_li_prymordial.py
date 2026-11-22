"""Workbench Entry 170, independent-code check: the neutron-injection scan rerun inside PRyMordial.

neutron_li.py, neutron_li_primat.py and neutron_li_shapes.py all solve Entry 53's network (bbn_net.py): its thermodynamics, its
n <-> p weak rates and its solver, with Kawano, PRIMAT or NACRE II nuclear rates.  Here the same question is put to a different code,
PRyMordial (Burns, Tait & Valli 2024, Eur. Phys. J. C 84, 86; arXiv:2307.07061).  That code has its own background (incomplete
neutrino decoupling, QED plasma corrections), its own weak rates (radiative, finite-nucleon-mass and thermal corrections), a
63-reaction network and its own stiff solver.  Only one thing is added: a neutron source term in dYn/dt during the mid- and
low-temperature eras.  It is the same Gaussian in ln T9 (width 0.15) as neutron_li.py, normalised so that A neutrons per baryon are
injected in all, and there is optionally a decaying-particle source, A exp(-t/tau)/tau, as in neutron_li_shapes.py.  Nothing else in
PRyMordial is changed.

Data and code: data/prymordial/PRyMordial-main.zip, https://github.com/vallima/PRyMordial (archive of the main branch),
    1,653,114 bytes, sha256 e7d41823af285241d8bb28c19f40f772132176781b723d99e583bebce81d8074, fetched 21 November 2026, unpacked
    beside the zip.  The script copies the unpacked package to data/prymordial/run/ (git-ignored) and runs it there, so the
    background and weak rates it caches never touch the download.  numba and Julia are not installed: plain-Python paths are used.
Settings: PRyMordial defaults (PRIMAT key rates, the 63-reaction network, tau_n = 878.4 s, incomplete decoupling), except
    omega_b h^2 = 0.02237 as in neutron_li.py.  Before the scan, the default run (omega_b h^2 = 0.02230) is compared with the output
    stored in the package's own demo notebook (PRyMdemoSM.ipynb), as a check of the installation.
Question (as in neutron_li.py): how many neutrons per baryon cut A = 7 (7Li + 7Be) to a third of the standard value, and by how much
    does D/H then rise?  Entry 170's claim is that no injection epoch does it with less than about +37% in D/H.

Run from code/:  ../.venv/bin/python -I neutron_li_prymordial.py          (about 30 minutes)
"""
import os, sys, json, time, shutil, types
HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.abspath(os.path.join(HERE, '..', 'data', 'prymordial', 'PRyMordial-main'))
RUN = os.path.abspath(os.path.join(HERE, '..', 'data', 'prymordial', 'run', 'PRyMordial-main'))
OUT = os.path.join(HERE, 'neutron_li_prymordial.json')
if not os.path.isdir(RUN):
    shutil.copytree(PKG, RUN)
os.chdir(RUN)                      # PRyM_init reads its tables from os.getcwd() at import
sys.path.insert(0, RUN)
import numpy as np
import PRyM.PRyM_init as PRyMini

PRyMini.numba_flag = False
PRyMini.julia_flag = False
PRyMini.smallnet_flag = False
PRyMini.verbose_flag = False

# ---- PRyM_main with a neutron source.  INJ is set from outside before each run.
src = open(os.path.join(RUN, 'PRyM', 'PRyM_main.py')).read()
def swap(a, b):
    global src
    assert src.count(a) == 1, (src.count(a), a)
    src = src.replace(a, b)
swap("from scipy.interpolate import interp1d\n", "from scipy.interpolate import interp1d, CubicSpline\n")
swap("        T_of_t_vec = np.vectorize(T_of_t)\n",
     "        T_of_t_vec = np.vectorize(T_of_t)\n"
     "        _lnT = CubicSpline(np.log(t_vec), np.log(Tg_vec)); _dlnT = _lnT.derivative()\n"
     "        def SRC(t):\n"
     "            A = INJ['A']\n"
     "            if A == 0.0:\n"
     "                return 0.0\n"
     "            if INJ['tau'] > 0.0:\n"
     "                return A * np.exp(-t / INJ['tau']) / INJ['tau']\n"
     "            lt = np.log(t)\n"
     "            T9 = np.exp(_lnT(lt)) * PRyMini.MeV_to_Kelvin / 1.e9\n"
     "            x = np.log(T9 / INJ['T9c']) / INJ['wid']\n"
     "            return A * np.exp(-0.5 * x * x) / (np.sqrt(2. * np.pi) * INJ['wid']) * (-_dlnT(lt)) / t\n"
     "        INJ['src'] = SRC; INJ['t_vec'] = np.array(t_vec)\n")
swap("dY = PRyMnucl.dYndtMT(Y,T_t,rhoBBN,nTOp_frwrd,nTOp_bkwrd),", "dY = PRyMnucl.dYndtMT(Y,T_t,rhoBBN,nTOp_frwrd,nTOp_bkwrd)+SRC(t),")
swap("dY = PRyMnucl.dYndtLT(Y,T_t,rhoBBN,nTOp_frwrd,nTOp_bkwrd),", "dY = PRyMnucl.dYndtLT(Y,T_t,rhoBBN,nTOp_frwrd,nTOp_bkwrd)+SRC(t),")
swap("method='BDF',jac=Jacobian_MT,rtol=1.e-6,atol=1.e-9)", "method='BDF',jac=Jacobian_MT,rtol=INJ['rtol_mt'],atol=1.e-9,**INJ['kw_mt'])")
swap("method='BDF',jac=Jacobian_LT,atol=1.e-15)", "method='BDF',jac=Jacobian_LT,atol=1.e-15,**INJ['kw_lt'])")
M = types.ModuleType('PRyM_main_inj')
M.__dict__['INJ'] = dict(A=0.0, T9c=0.5, wid=0.15, tau=0.0, rtol_mt=1e-6, kw_mt={}, kw_lt={})
exec(compile(src, 'PRyM_main_inj', 'exec'), M.__dict__)
INJ = M.INJ

D_ERR = 0.030e-5                   # Cooke et al. 2018, as in neutron_li.py
DEMO = [3.04438852, 5.69940717, 93.03360798, 0.24555799, 0.24688388, 2.45571913, 1.04109377, 5.4396027]   # PRyMdemoSM.ipynb
NAMES = ['Neff', 'Omeganu_rel_x1e6', 'Summnu_over_Omeganu', 'YP_CMB', 'YP_BBN', 'DoH_x1e5', 'He3oH_x1e5', 'Li7oH_x1e10']


def run(A=0.0, T9c=0.5, wid=0.15, tau=0.0, rtol_mt=1e-6, kw_mt=None, kw_lt=None):
    INJ.update(A=A, T9c=T9c, wid=wid, tau=tau, rtol_mt=rtol_mt, kw_mt=kw_mt or {}, kw_lt=kw_lt or {})
    r = M.PRyMclass().PRyMresults()
    return dict(Yp=float(r[4]), D=float(r[5]) * 1e-5, He3=float(r[6]) * 1e-5, Li7=float(r[7]) * 1e-10)


def injected(A, **kw):
    """Integral of the source over the run, on PRyMordial's own time grid (a check of the normalisation)."""
    run(A=A, **kw)
    t = np.logspace(np.log10(INJ['t_vec'][0]) + 0.01, np.log10(INJ['t_vec'][-1]) - 0.01, 40001)
    s = np.array([INJ['src'](x) for x in t])
    return float(np.trapezoid(s * t, np.log(t)))


def solve(base, target=1 / 3., **kw):
    lo, hi, A, hist = None, None, 1e-5, []
    for it in range(14):
        r = run(A=A, **kw)
        ratio = r['Li7'] / base['Li7']
        hist.append((A, ratio, r))
        if abs(np.log(ratio / target)) < 0.01:
            return A, r, hist
        if ratio > target:
            lo = (A, ratio)
        else:
            hi = (A, ratio)
        if lo and hi:
            la, lb = np.log(lo[0]), np.log(hi[0]); ra, rb = np.log(lo[1]), np.log(hi[1])
            A = float(np.exp(la + (np.log(target) - ra) * (lb - la) / (rb - ra)))
        elif lo:
            A *= 10.0 if ratio > 0.9 else 3.0
        else:
            A /= 3.0
        if A > 0.02:
            return None, min(hist, key=lambda h: h[1])[2], hist
    return A, r, hist


def row(label, A, r, hist, base, **kw):
    x = dict(kw, label=label, solved=A is not None, A_per_baryon=A, Li_ratio=r['Li7'] / base['Li7'], D=r['D'],
             D_rise=r['D'] / base['D'] - 1.0, D_pull=(r['D'] - base['D']) / D_ERR, He3_rise=r['He3'] / base['He3'] - 1.0,
             dYp=r['Yp'] - base['Yp'], runs=len(hist))
    print('%-8s %s: A %.3g, Li x%.3f, D %+.1f%% (%.0f sigma), 3He %+.1f%%, Yp %+.5f [%d runs, %d s]' % (
        label, kw, A or float('nan'), x['Li_ratio'], 100 * x['D_rise'], x['D_pull'], 100 * x['He3_rise'], x['dYp'], len(hist),
        time.time() - T0), flush=True)
    return x


def main():
    global T0
    T0 = time.time()
    if len(sys.argv) > 1 and sys.argv[1] == 'resume':          # steps 1-3 already in OUT; background and weak rates cached
        return resume()
    out = dict(code='PRyMordial main branch (see the header), run unmodified except for the neutron source term',
               omega_b_h2=0.02237, numba=False, julia=False, network='63 reactions')
    # 1. installation check against the demo notebook (defaults: omega_b h^2 = 0.02230); compute and cache background and weak rates
    PRyMini.compute_bckg_flag = True; PRyMini.save_bckg_flag = True
    PRyMini.compute_nTOp_flag = True; PRyMini.save_nTOp_flag = True
    PRyMini.Omegabh2 = 0.02230; PRyMini.eta0b = PRyMini.Omegabh2_to_eta0b * 0.02230
    INJ.update(A=0.0, kw_mt={}, kw_lt={})
    r = M.PRyMclass().PRyMresults()
    out['demo_check'] = dict(ours=dict(zip(NAMES, map(float, r))), notebook=dict(zip(NAMES, DEMO)),
                             max_rel_diff=float(np.max(np.abs(np.array(r) / np.array(DEMO) - 1.0))))
    print('demo check: max relative difference %.2e  (D/H %.5f vs %.5f, Li %.4f vs %.4f)  [%d s]' % (
        out['demo_check']['max_rel_diff'], r[5], DEMO[5], r[7], DEMO[7], time.time() - T0), flush=True)
    PRyMini.compute_bckg_flag = False; PRyMini.save_bckg_flag = False
    PRyMini.compute_nTOp_flag = False; PRyMini.save_nTOp_flag = False
    PRyMini.Omegabh2 = 0.02237; PRyMini.eta0b = PRyMini.Omegabh2_to_eta0b * 0.02237
    # 2. normalisation of the source
    out['normalisation'] = {str(T9c): injected(1e-5, T9c=T9c) / 1e-5 for T9c in (1.0, 0.5, 0.3)}
    out['normalisation']['decay_300s'] = injected(1e-5, tau=300.0) / 1e-5
    print('injected / requested:', {k: round(v, 4) for k, v in out['normalisation'].items()}, flush=True)
    # 3. the scan, PRIMAT key rates (PRyMordial's default)
    base = run()
    out['standard_primat'] = base
    print('standard (PRIMAT):', base, flush=True)
    rows = []
    for T9c in (1.0, 0.7, 0.6, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25):
        A, r, hist = solve(base, T9c=T9c)
        rows.append(row('primat', A, r, hist, base, T9c=T9c, keV=round(T9c * 86.17, 2), wid=0.15))
    out['primat_gaussian'] = rows
    json.dump(out, open(OUT, 'w'), indent=1)
    json.dump(out, open(OUT, 'w'), indent=1)
    rest(out)


def resume():
    out = json.load(open(OUT))
    PRyMini.compute_bckg_flag = False; PRyMini.save_bckg_flag = False
    PRyMini.compute_nTOp_flag = False; PRyMini.save_nTOp_flag = False
    PRyMini.Omegabh2 = 0.02237; PRyMini.eta0b = PRyMini.Omegabh2_to_eta0b * 0.02237
    out['resumed'] = 'steps 4-6 rerun on ' + time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()) + ' after a solver check with a 50 s step cap proved too slow'
    rest(out)


def rest(out):
    base = run()
    assert abs(base['D'] / out['standard_primat']['D'] - 1) < 1e-6, (base, out['standard_primat'])
    rows = out['primat_gaussian']
    best = min((x for x in rows if x['solved']), key=lambda x: x['D_rise'])
    # 4. solver check at the minimum: tighter tolerances and a step cap
    A_b = best['A_per_baryon']
    tight = run(A=A_b, T9c=best['T9c'], rtol_mt=1e-8, kw_lt=dict(rtol=1e-6))
    loose = run(A=A_b, T9c=best['T9c'])
    out['solver_check'] = dict(T9c=best['T9c'], A=A_b, default=loose, tight=tight,
                               D_rel_diff=tight['D'] / loose['D'] - 1.0, Li_rel_diff=tight['Li7'] / loose['Li7'] - 1.0)
    base_t = run(rtol_mt=1e-8, kw_lt=dict(rtol=1e-6))
    out['solver_check']['D_rise_tight'] = tight['D'] / base_t['D'] - 1.0
    out['solver_check']['Li_ratio_tight'] = tight['Li7'] / base_t['Li7']
    print('solver check:', {k: v for k, v in out['solver_check'].items() if k not in ('default', 'tight')}, flush=True)
    # 5. wider bursts and a decaying particle, PRIMAT rates
    rows = []
    for wid, T9s in ((0.3, (0.7, 0.55, 0.45)), (0.05, (0.5, 0.45, 0.4))):
        for T9c in T9s:
            A, r, hist = solve(base, T9c=T9c, wid=wid)
            rows.append(row('shape', A, r, hist, base, T9c=T9c, keV=round(T9c * 86.17, 2), wid=wid))
    for tau in (150.0, 300.0, 600.0):
        A, r, hist = solve(base, tau=tau)
        rows.append(row('decay', A, r, hist, base, tau_s=tau))
    out['shapes'] = rows
    json.dump(out, open(OUT, 'w'), indent=1)
    # 6. NACRE II key rates near the minimum
    PRyMini.nacreii_flag = True
    PRyMini.ReloadKeyRates()
    base_n = run()
    out['standard_nacreii'] = base_n
    print('standard (NACRE II):', base_n, flush=True)
    rows = []
    for T9c in (0.6, 0.5, 0.45, 0.4):
        A, r, hist = solve(base_n, T9c=T9c)
        rows.append(row('nacreii', A, r, hist, base_n, T9c=T9c, keV=round(T9c * 86.17, 2), wid=0.15))
    out['nacreii_gaussian'] = rows
    PRyMini.nacreii_flag = False
    PRyMini.ReloadKeyRates()
    allrows = [x for x in out['primat_gaussian'] + out['shapes'] + out['nacreii_gaussian'] if x['solved']]
    least = min(allrows, key=lambda x: x['D_rise'])
    out['least_D_rise'] = {k: least[k] for k in least}
    out['seconds'] = round(time.time() - T0)
    json.dump(out, open(OUT, 'w'), indent=1)
    print('least D/H rise over all runs: %+.1f%% (%s)' % (100 * least['D_rise'], {k: least.get(k) for k in ('label', 'T9c', 'wid', 'tau_s')}))
    print('written', OUT, out['seconds'], 's')


if __name__ == '__main__':
    main()
