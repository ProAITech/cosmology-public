"""Entry 176, analysis after the frozen result: what the weak depolarization at 0.5-1 R500 says about the reversal length.

Between 0.5 and 1 R500 the rotation measure differs from source to source by sigma = sqrt(excess variance), yet the sources there
are almost as polarized as field sources and almost none are missing.  If the rotation varied across each source as much as it
does between sources, the sources would be depolarized by exp(-2 sigma^2 lambda^4).  They are not, so the screen is nearly uniform
across a source: its reversal length exceeds the projected size of the polarized emission.

Inputs (all from this folder's JSON files): the excess variance at 0.5-1 R500 (results.json), the flux-matched fractional
polarization there against rotated clusters (variants.json), the source count against expectation (model.json).
Screen: Gaussian correlation, xi(r) = sigma^2 exp(-r^2 / 2 Lambda^2); source: a uniform disc of diameter D at the cluster, for
which the variance inside the disc is sigma^2 [1 - <exp(-r^2 / 2 Lambda^2)>] with r the separation of two random points.
Run from the repository root:  .venv/bin/python -I entries/176-cluster-faraday-profile/coherence.py
"""
import json, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
res = json.load(open(os.path.join(HERE, 'results.json'))); var = json.load(open(os.path.join(HERE, 'variants.json'))); mod = json.load(open(os.path.join(HERE, 'model.json')))
p = res['reported_only']['profile'][1]                       # 0.5-1 R500
sigma2 = p['excess_variance']; sigma2_err = p['null_sd']
fr = var['fracpol']['ranges'][1]
ratio = fr['flux_matched_ratio'] / fr['ratio_null_mean']; ratio_err = fr['ratio_null_sd'] / fr['ratio_null_mean']
mem = mod['members']['ranges'][1]
lost_max = max(0.0, -(mem['surplus'] - 2 * np.sqrt(mem['expected_without_members'])) / mem['expected_without_members'])   # 2 sigma
# band 799.5-1087.5 MHz (SPICE-RACS DR2, from RACS-low3; Thomson et al. 2026, arXiv:2605.16917): mean of lambda^4 over the band,
# uniform in frequency.  Until 22 November 2026 this line used 744-1032 MHz, the band of the first data release (RACS-low1).
nu = np.linspace(799.5e6, 1087.5e6, 2001); lam4 = float(np.mean((299792458.0 / nu) ** 4))
dep_all = ratio * (1 - lost_max)                              # mean polarization kept, counting lost sources as fully depolarized
dep_low = (ratio - 2 * ratio_err) * (1 - lost_max)            # the same with the ratio at its 2 sigma lower edge
out = dict(sigma_between_sources=float(np.sqrt(sigma2)), excess_variance=sigma2, lambda4_mean=lam4, depolarization_if_uncorrelated=float(np.exp(-2 * sigma2 * lam4)),
           polarization_ratio_survivors=float(ratio), ratio_error=float(ratio_err), share_lost_at_most=float(lost_max), kept=float(dep_all), kept_lower=float(dep_low))
g = np.random.default_rng(1)
a, b = g.uniform(-0.5, 0.5, (2, 400000, 2))
ok = (np.hypot(a[:, 0], a[:, 1]) <= 0.5) & (np.hypot(b[:, 0], b[:, 1]) <= 0.5)
r2 = ((a[ok] - b[ok]) ** 2).sum(1)                           # squared separations in a disc of unit diameter
def inside_fraction(L_over_D):
    return float(1 - np.mean(np.exp(-r2 / (2 * L_over_D ** 2))))
grid = np.round(np.arange(0.2, 6.01, 0.05), 2)
frac = np.array([inside_fraction(v) for v in grid])
for name, kept in (('central', dep_all), ('conservative', dep_low)):
    s_in2 = -np.log(kept) / (2 * lam4)                         # variance allowed inside a source
    f_max = s_in2 / sigma2
    L = float(grid[np.argmax(frac <= f_max)]) if (frac <= f_max).any() else None
    out[name] = dict(polarization_kept=float(kept), variance_inside_source_at_most=float(s_in2), share_of_between_source_variance=float(f_max), lambda_over_diameter_at_least=L)
kpc_per_arcsec = 4.5                                          # at z = 0.30, the clusters' median redshift
out['in_kpc'] = {('%d arcsec' % d): {k: (round(out[k]['lambda_over_diameter_at_least'] * d * kpc_per_arcsec, 1) if out[k]['lambda_over_diameter_at_least'] else None)
                                      for k in ('central', 'conservative')} for d in (1, 2, 5, 10)}
print(json.dumps(out, indent=1))
json.dump(out, open(os.path.join(HERE, 'coherence.json'), 'w'), indent=1)
