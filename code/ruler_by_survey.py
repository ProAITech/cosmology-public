"""ruler_by_survey.py - Workbench Entry 58, check added 11 November 2026
Does the clustering ruler exclude the logarithmic law (steady expansion) in each survey on its own?
SDSS (BOSS and eBOSS, final 2020 values), DESI first release (2024), DESI second release (2025).
Law: across = K ln(1+z), along = K / (1+z), combined D_V = (z across^2 along)^(1/3); one free K.
Standard model: flat, matter share and scale free."""
import math, json
import numpy as np
from scipy import optimize, integrate

DR1 = [('V',0.295,7.93,0.15),('MH',0.510,13.62,0.25,20.98,0.61,-0.445),('MH',0.706,16.85,0.32,20.08,0.60,-0.420),('MH',0.930,21.71,0.28,17.88,0.35,-0.389),
       ('MH',1.317,27.79,0.69,13.82,0.42,-0.444),('V',1.491,26.07,0.67),('MH',2.330,39.71,0.94,8.52,0.17,-0.477)]
DR2 = [('V',0.295,7.942,0.075),('MH',0.510,13.588,0.167,21.863,0.425,-0.459),('MH',0.706,17.351,0.177,19.455,0.330,-0.404),('MH',0.934,21.576,0.152,17.641,0.193,-0.416),
       ('MH',1.321,27.601,0.318,14.176,0.221,-0.434),('MH',1.484,30.512,0.760,12.817,0.516,-0.500),('MH',2.330,38.988,0.531,8.632,0.101,-0.431)]
def sdss(rho):   # the two z = 2.33 rows are the forest alone and the forest against quasars; correlations are not in the summary table
    return [('V',0.15,4.47,0.17),('MH',0.38,10.23,0.17,25.00,0.76,rho),('MH',0.51,13.36,0.21,22.33,0.58,rho),('MH',0.698,17.86,0.33,19.33,0.53,rho),('V',0.85,18.33,0.60),
            ('MH',1.48,30.69,0.80,13.26,0.55,rho),('MH',2.33,37.6,1.9,8.93,0.28,rho),('MH',2.33,37.3,1.7,9.08,0.34,rho)]
zg = np.linspace(0, 2.6, 5201)
def lcdm(om):
    dc = integrate.cumulative_trapezoid(1 / np.sqrt(om * (1 + zg) ** 3 + 1 - om), zg, initial=0)
    return (lambda z: float(np.interp(z, zg, dc))), (lambda z: 1 / math.sqrt(om * (1 + z) ** 3 + 1 - om))
LOG = ((lambda z: math.log(1 + z)), (lambda z: 1 / (1 + z)))
def chi(data, model, K, skip_hi=False):
    dm, dh = model; c2 = 0.0
    for d in data:
        if skip_hi and d[1] > 2: continue
        if d[0] == 'V':
            c2 += ((d[2] - K * (d[1] * dm(d[1]) ** 2 * dh(d[1])) ** (1 / 3)) / d[3]) ** 2
        else:
            a = (d[2] - K * dm(d[1])) / d[3]; b = (d[4] - K * dh(d[1])) / d[5]; r = d[6]; c2 += (a * a + b * b - 2 * r * a * b) / (1 - r * r)
    return c2
def best(data, model, skip_hi=False):
    r = optimize.minimize_scalar(lambda K: chi(data, model, K, skip_hi), bounds=(5, 80), method='bounded', options={'xatol': 1e-8}); return float(r.fun), float(r.x)
def best_lcdm(data, skip_hi=False):
    r = optimize.minimize_scalar(lambda om: best(data, lcdm(om), skip_hi)[0], bounds=(0.1, 0.7), method='bounded'); return float(r.fun), float(r.x)
out = {}
for name, data in (("SDSS, correlations -0.4", sdss(-0.4)), ("SDSS, correlations 0", sdss(0.0)), ("DESI first release", DR1), ("DESI second release", DR2)):
    n = sum(1 if d[0] == 'V' else 2 for d in data); cl, K = best(data, LOG); cs, om = best_lcdm(data)
    n_lo = sum((1 if d[0] == 'V' else 2) for d in data if d[1] < 2); cl_lo, K_lo = best(data, LOG, True); cs_lo, om_lo = best_lcdm(data, True)
    hi = [d for d in data if d[1] > 2]
    print("%-26s %2d numbers: logarithmic law %.1f, standard %.1f (matter share %.3f) | below redshift 2, %d numbers: %.1f and %.1f" % (name, n, cl, cs, om, n_lo, cl_lo, cs_lo))
    row = dict(n=n, log=cl, lcdm=cs, om=om, n_lo=n_lo, log_lo=cl_lo, lcdm_lo=cs_lo, hi=[])
    for d in hi:
        pm, ph = K_lo * math.log(1 + d[1]), K_lo / (1 + d[1]); a, b = (d[2] - pm) / d[3], (d[4] - ph) / d[5]; r = d[6]
        F = d[2] / d[4]; eF = F * math.sqrt((d[3] / d[2]) ** 2 + (d[5] / d[4]) ** 2 - 2 * r * d[3] * d[5] / (d[2] * d[4]))
        print("      at %.2f, scale from the lower tracers: across %.2f +- %.2f against %.2f (%+.1f sigma); along %.3f +- %.3f against %.3f (%+.1f sigma); the pair: chi2 %.1f"
              % (d[1], d[2], d[3], pm, a, d[4], d[5], ph, b, (a * a + b * b - 2 * r * a * b) / (1 - r * r)))
        print("      across over along, no scale needed: %.2f +- %.2f; the law %.3f (%+.1f sigma); standard at matter share 0.30: %.3f (%+.1f sigma)"
              % (F, eF, (1 + d[1]) * math.log(1 + d[1]), (F - (1 + d[1]) * math.log(1 + d[1])) / eF, lcdm(0.30)[0](d[1]) / lcdm(0.30)[1](d[1]), (F - lcdm(0.30)[0](d[1]) / lcdm(0.30)[1](d[1])) / eF))
        row["hi"].append(dict(z=d[1], dm=d[2], dm_err=d[3], dm_pred=pm, dh=d[4], dh_err=d[5], dh_pred=ph, pair_chi=(a * a + b * b - 2 * r * a * b) / (1 - r * r), ratio=F, ratio_err=eF))
    out[name] = row
json.dump(out, open("ruler_by_survey.json", "w"), indent=1)
