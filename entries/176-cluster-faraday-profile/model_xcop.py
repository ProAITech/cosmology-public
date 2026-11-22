"""Entry 176, analysis after the frozen result: the field at 0.75 R500 with the X-COP universal density profile in place of the beta model.

Ghirardini et al. (2019, A&A 621, A41; arXiv:1805.00042, Table 2) fit the scaled density of the 12 X-COP clusters with piecewise power
laws, n_e E(z)^-2 = A x^B (x = r/R500, A in 1e-4 cm^-3), in eight radial ranges from 0.01 to 2 R500.  Two variants:
  (a) X-COP shape, normalised to each cluster's eRASS1 gas mass inside R500 (as model.py normalises its beta model);
  (b) X-COP absolute density, A x^B E(z)^2, ignoring the gas masses (X-COP clusters are massive, 3-12 x 10^14 Msun, with higher gas
      fractions than the median eRASS1 group, so this overstates the density of most clusters here).
Everything else as model_variants.py's standard case (eta = 0.5, Lambda = 10 kpc, sources in front, line of sight to 3 R500).
The text below is model_variants.py's header, kept for the shared setup.

model.py assumed a beta model with beta = 2/3 and core 0.2 R500, a share of sources in front of each cluster from a redshift
distribution of median 0.8, and B proportional to n_e^0.5.  Here each is changed in turn and the fit to 0.5-2 R500 repeated:
  gas profile     core radius 0.1, 0.2, 0.3 R500; beta 0.6, 2/3, 0.75 (each normalised to the cluster's gas mass inside R500)
  sources in front  none; the standard share; twice the standard share (capped at 0.9); a median source redshift of 0.5 and 1.2
  outer truncation  the line of sight cut at 2, 3 and 5 R500
The reversal length Lambda enters only as B proportional to Lambda^-0.5 and is not varied here.
Run from the repository root:  .venv/bin/python -I entries/176-cluster-faraday-profile/model_xcop.py     (about 2 minutes)
"""
import json, os, sys
import numpy as np
from scipy.spatial import cKDTree
from scipy.integrate import quad
from astropy.io import fits
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cluster_rm as c

log = lambda s: print(s, flush=True)
KPC, MSUN, MP, MUE = 3.0857e21, 1.989e33, 1.6726e-24, 1.17
LAMBDA, NREF, ETA = 10.0, 1e-3, 0.5
res = json.load(open(os.path.join(HERE, 'results.json')))
prof = res['reported_only']['profile']
obs = np.array([p['excess_variance'] for p in prof[:4]]); err = np.array([p['null_sd'] for p in prof[:4]])
USE = [1, 2, 3]
rows, src, _ = c.load_spice()
S, _ = c.sample(src, c.CFG['fg_radius'], c.CFG['min_neigh'])
with fits.open(c.CLUSTERS, memmap=True) as h:
    d = h[1].data
    ra, dec, z = np.asarray(d['RA'], float), np.asarray(d['DEC'], float), np.asarray(d['BEST_Z'], float)
    pc, r500, m500, mgas = np.asarray(d['PCONT'], float), np.asarray(d['R500'], float), np.asarray(d['M500'], float), np.asarray(d['MGAS500'], float)
ok = np.isfinite(ra) & np.isfinite(dec) & np.isfinite(z) & np.isfinite(r500) & (pc < 0.5) & (z >= 0.05) & (r500 > 0)
l, b = c.galactic(ra[ok], dec[ok])
cl = dict(l=l, b=b, th=r500[ok] / 1000.0 / c.ang_diam_Mpc(z[ok]) / c.D2R, z=z[ok], m500=m500[ok], r500=r500[ok], mgas=mgas[ok])
tree = cKDTree(c.gal_unit(S['l'], S['b']))
x, which = c.scaled_distance(tree, len(S['d']), cl)
field = ~np.isfinite(x)
dF, eF = S['d'][field], S['e'][field]
V_field = c.V_of(dF, eF)
sel = (x >= 0.5) & (x < 2.0)
xs, ks = x[sel], which[sel]
R = xs * cl['r500'][ks]
binidx = np.digitize(xs, [0.5, 1.0, 1.5])
g = np.random.default_rng(c.CFG['seed'] + 1776)
grid = np.round(np.concatenate([np.arange(0.4, 12.0, 0.2), np.arange(12.0, 30.01, 1.0)]), 2)



XCOP = [(0.01, 0.07, 13.00, -0.48), (0.07, 0.13, 3.44, -1.04), (0.13, 0.21, 2.60, -1.21), (0.21, 0.31, 2.85, -1.16),
        (0.31, 0.46, 1.84, -1.60), (0.46, 0.72, 1.63, -1.80), (0.72, 1.14, 1.42, -2.38), (1.15, 2.00, 1.53, -2.47)]


def xcop_shape(xr):
    """n_e E(z)^-2 in cm^-3 from Table 2 of Ghirardini et al. 2019; the first and last pieces extended below 0.01 and beyond 2 R500."""
    xr = np.maximum(np.asarray(xr, float), 1e-3)
    out = np.empty_like(xr)
    for i, (lo, hi, A, B) in enumerate(XCOP):
        m = (xr >= (0 if i == 0 else lo)) & (xr < (np.inf if i == len(XCOP) - 1 else XCOP[i + 1][0]))
        out[m] = A * 1e-4 * xr[m] ** B
    return out


E2 = c.CFG['Om'] * (1 + cl['z']) ** 3 + 1 - c.CFG['Om']
yy = np.linspace(1e-4, 1.0, 20001)
SHAPE_INT = float(np.trapezoid(xcop_shape(yy) * yy ** 2, yy))          # Int_0^1 f(x) x^2 dx


def fit_profile(mode, nreal=40):
    if mode == 'gas':
        scale = cl['mgas'] * 1e11 * MSUN / (4 * np.pi * MUE * MP * (cl['r500'] * KPC) ** 3 * SHAPE_INT)      # multiplies f(x)
    else:
        scale = E2.copy()
    L = np.sqrt(np.maximum((3.0 * cl['r500'][ks]) ** 2 - R ** 2, 0.0))
    t = (np.arange(300) + 0.5) / 300
    ll = L[:, None] * t[None, :]
    xr = np.sqrt(R[:, None] ** 2 + ll ** 2) / cl['r500'][ks][:, None]
    ne = scale[ks][:, None] * xcop_shape(xr)
    I = 2.0 * ((ne / NREF) ** (2 * ETA) * ne ** 2).sum(1) * (L / 300)
    z0 = 0.8 / np.log(2) ** (1 / 1.2)
    pf = np.minimum(0.9, 1.0 - np.exp(-(cl['z'][ks] / z0) ** 1.2))
    n_at = float(np.median(scale * xcop_shape(0.75)))
    chi = []
    for B0 in grid:
        sig2 = 812.0 ** 2 * LAMBDA * (B0 ** 2 / 3.0) * I * (1.0 + cl['z'][ks]) ** -4
        acc = np.zeros(4)
        for r_ in range(nreal):
            k = g.integers(0, len(dF), len(xs)); u = g.uniform(size=len(xs))
            dsim = dF[k] + g.normal(0, 1, len(xs)) * np.sqrt(sig2 * (u >= pf))
            for j in USE:
                m = binidx == j
                acc[j] += c.V_of(dsim[m], eF[k][m]) - V_field
        pr = acc / nreal
        chi.append(float((((obs[USE] - pr[USE]) / err[USE]) ** 2).sum()))
    chi = np.array(chi); i0 = int(np.argmin(chi)); inside = grid[chi <= chi[i0] + 1.0]
    conv = (n_at / NREF) ** ETA
    return dict(B_at_0p75_R500=float(grid[i0] * conv), lo=float(inside.min() * conv), hi=float(inside.max() * conv), chi2=float(chi[i0]),
                n_e_at_0p75_R500=n_at, B0_best=float(grid[i0]))


out = dict(source='Ghirardini et al. 2019, A&A 621, A41 (arXiv:1805.00042), Table 2, density', pieces=XCOP, shape_integral=SHAPE_INT)
beta_n = float(np.median(cl['mgas'] * 1e11 * MSUN / (4 * np.pi * MUE * MP * (0.2 * cl['r500'] * KPC) ** 3 * (5 - np.arctan(5))))) * (1 + (0.75 / 0.2) ** 2) ** -1
out['beta_model_n_e_at_0p75_R500'] = beta_n
for mode, name in (('gas', 'X-COP shape, normalised to each gas mass'), ('xcop', 'X-COP absolute density')):
    f = fit_profile(mode); f['case'] = name; out[mode] = f
    log('%-42s B at 0.75 R500 = %.2f (%.2f-%.2f) uG, chi2 %.1f; median n_e there %.2e (beta model %.2e)' % (
        name, f['B_at_0p75_R500'], f['lo'], f['hi'], f['chi2'], f['n_e_at_0p75_R500'], beta_n))
json.dump(out, open(os.path.join(HERE, 'model_xcop.json'), 'w'), indent=1, default=float)
