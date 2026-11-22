"""Entry 184: Entry 113's helium neutron clock checked with PRyMordial.

Entry 113 read the neutron lifetime from primordial helium with a slope dYp/dtau = 2.039e-4 per second from Entry 53's small network,
anchored to a full-code prediction of Yp = 0.2467 at 878.4 s (eta10 = 6.124, three neutrino species), against the measured
Yp = 0.24572 +- 0.00126 (Entry 102), and found tau = 873.6 +- 6.3 s.  Here both the anchor and the slope come from PRyMordial
(Burns, Tait & Valli 2024, arXiv:2307.07061), run unmodified: PRIMAT key rates, 63-reaction network, incomplete neutrino decoupling,
omega_b h^2 = 0.02237, tau_n set to each value in turn (the code normalises its weak rates to 1 / (F_n tau_n)).
Observers' Yp = 4y/(1+4y) with y = n_He/n_H is PRyMordial's YPBBN (= 4 Y_4He); YPCMB, the true mass fraction, is listed too.
Code: data/prymordial/PRyMordial-main.zip (github.com/vallima/PRyMordial, 1,653,114 bytes, sha256 e7d41823...8074), run from the
copy in data/prymordial/run/ as in neutron_li_prymordial.py.  Run from code/:  ../.venv/bin/python -I neutron_clock_prymordial.py
"""
import os, sys, json, shutil
HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.abspath(os.path.join(HERE, '..', 'data', 'prymordial', 'PRyMordial-main'))
RUN = os.path.abspath(os.path.join(HERE, '..', 'data', 'prymordial', 'run', 'PRyMordial-main'))
if not os.path.isdir(RUN):
    shutil.copytree(PKG, RUN)
os.chdir(RUN)
sys.path.insert(0, RUN)
import numpy as np
import PRyM.PRyM_init as PRyMini
PRyMini.numba_flag = False; PRyMini.julia_flag = False; PRyMini.smallnet_flag = False; PRyMini.verbose_flag = False
import PRyM.PRyM_main as M

MEAS, MEAS_ERR, THEORY = 0.24572, 0.00126, 0.0002          # Entry 102's helium, and Entry 113's theory allowance
TAUS = [868.4, 877.83, 878.4, 887.7, 888.4]


def run(tau, first=False):
    PRyMini.compute_bckg_flag = first; PRyMini.save_bckg_flag = first
    PRyMini.compute_nTOp_flag = first; PRyMini.save_nTOp_flag = first
    PRyMini.Omegabh2 = 0.02237; PRyMini.eta0b = PRyMini.Omegabh2_to_eta0b * 0.02237
    PRyMini.tau_n = tau * PRyMini.second
    c = M.PRyMclass()
    c.PRyMresults()
    return dict(tau=tau, YPBBN=float(c.YPBBN_f), YPCMB=float(c.YPCMB_f), DoH=float(c.DoHx1e5_f), Li7oH=float(c.Li7oHx1e10_f))


def main():
    have_cache = os.path.exists(os.path.join(RUN, 'PRyMrates', 'thermo')) or True
    rows = []
    for k, tau in enumerate(TAUS):
        r = run(tau, first=False)
        rows.append(r)
        print('tau %.2f s: YPBBN %.5f, YPCMB %.5f, D/H %.4f e-5, Li/H %.3f e-10' % (tau, r['YPBBN'], r['YPCMB'], r['DoH'], r['Li7oH']), flush=True)
    t = np.array([r['tau'] for r in rows]); y = np.array([r['YPBBN'] for r in rows]); yc = np.array([r['YPCMB'] for r in rows])
    slope, icpt = np.polyfit(t, y, 1)
    resid = float(np.max(np.abs(y - (icpt + slope * t))))
    y878 = float(icpt + slope * 878.4)
    tau_hat = (MEAS - icpt) / slope
    err = np.hypot(MEAS_ERR, THEORY) / slope
    out = dict(rows=rows, slope_YPBBN_per_s=float(slope), slope_YPCMB_per_s=float(np.polyfit(t, yc, 1)[0]), max_residual_linear=resid,
               YPBBN_at_878p4=y878, entry113=dict(slope=2.039e-4, anchor=0.2467, tau=873.6, tau_err=6.3),
               tau_from_helium=dict(tau=float(tau_hat), err=float(err), measured=MEAS, measured_err=MEAS_ERR, theory=THEORY),
               sigma_from=dict(bottle_877p83=float((MEAS - (icpt + slope * 877.83)) / np.hypot(MEAS_ERR, THEORY)),
                               beam_887p7=float((MEAS - (icpt + slope * 887.7)) / np.hypot(MEAS_ERR, THEORY))))
    print('slope YPBBN %.4e per s (Entry 113: 2.039e-4); YPBBN at 878.4 s %.5f (Entry 113 anchor 0.2467); tau from helium %.1f +- %.1f s (Entry 113: 873.6 +- 6.3); '
          'measured vs bottle %+.2f sigma, vs beam %+.2f sigma' % (slope, y878, tau_hat, err, out['sigma_from']['bottle_877p83'], out['sigma_from']['beam_887p7']))
    json.dump(out, open(os.path.join(HERE, 'neutron_clock_prymordial.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
