"""Entry 187, after the frozen result (not pre-registered): fit the beam-smoothed shapes of a cooling to the whole radial profile.

Reads radial_<sample>_smica_l5.json (radial.py).  The model is excess_j = delta * f_j, where f_j is the bin mean of the
beam-smoothed profile of unit depth ('top': inside 0.2 Mpc; 'fade': falling linearly to 0 at 1.5 Mpc), averaged over the
sample's redshifts.  Generalised least squares with the bin covariance of the 1,000 rotations (corr x sd x sd), so the
large-scale pattern the bins share is weighted out.  If radial.py saved each rotation's profile (work/), the same weights
are applied to every rotation and their scatter is reported beside the Gaussian error.  delta < 0 is a cooling.
Writes radial_fit.json.  Run from the repository root:
  .venv/bin/python -I entries/187-cmb-spirals/radial_fit.py"""
import json
import numpy as np
from pathlib import Path
HERE = Path(__file__).resolve().parent
out = {}
for s in ('P', 'P_dense', 'A_2MRS'):
    r = json.load(open(HERE / f'radial_{s}_smica_l5.json'))
    d, sd, corr = np.array(r['excess']), np.array(r['null_sd']), np.array(r['corr'])
    Ci = np.linalg.inv(corr * np.outer(sd, sd))
    row = dict(n_used=r['n_used'], chi2_zero=float(d @ Ci @ d), dof=len(d), corr_bin1_bin2=float(corr[0, 1]), corr_bin1_bin6=float(corr[0, 5]))
    for k in ('top', 'fade'):
        f = np.array(r['predicted_minus15'][k]) / -15.0
        a, sa = (f @ Ci @ d) / (f @ Ci @ f), (f @ Ci @ f) ** -0.5
        m = d + 15.0 * f
        rot = HERE / 'work' / f'radial_rot_{s}_smica_l5.npy'
        if rot.exists():                                    # the same fit on every rotation: its empirical scatter
            Pr = np.load(rot); w = Ci @ f / (f @ Ci @ f); ar = (Pr - Pr.mean(0)) @ w
            # split halves: weights from the covariance of rotations 1-500, scatter measured on rotations 501-1000
            h = len(Pr) // 2; C1 = np.cov(Pr[:h].T); w1 = np.linalg.solve(C1, f) / (f @ np.linalg.solve(C1, f))
            a2 = (Pr[h:] - Pr[h:].mean(0)) @ w1
            emp = dict(p_le_observed=float((ar <= a).mean()), q05=float(np.quantile(ar, 0.05)), q95=float(np.quantile(ar, 0.95)),
                       sd_split_half=float(a2.std(ddof=1)), n=len(ar))
        else:
            emp = None
        row[k] = dict(delta=float(a), sd=float(sa), upper95_cooling=float(-(a - 1.645 * sa)), sigma_from_minus15=float((a + 15) / sa),
                      sigma_from_minus27=float((a + 27) / sa), chi2_minus15=float(m @ Ci @ m), rotations=emp)
    out[s] = row
    print(s, r['n_used'], f"chi2(0) {row['chi2_zero']:.2f}/6", ' | '.join(
        f"{k}: {row[k]['delta']:+.2f} +- {row[k]['sd']:.2f} (95% {row[k]['upper95_cooling']:.1f}; -15 at {row[k]['sigma_from_minus15']:.1f} sigma; chi2(-15) {row[k]['chi2_minus15']:.0f}; rotations {row[k]['rotations']})"
        for k in ('top', 'fade')))
json.dump(out, open(HERE / 'radial_fit.json', 'w'), indent=1)
