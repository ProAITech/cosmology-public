"""Entry 176, analysis after the frozen result: how much the field at 0.75 R500 depends on the model's assumptions.

model.py assumed a beta model with beta = 2/3 and core 0.2 R500, a share of sources in front of each cluster from a redshift
distribution of median 0.8, and B proportional to n_e^0.5.  Here each is changed in turn and the fit to 0.5-2 R500 repeated:
  gas profile     core radius 0.1, 0.2, 0.3 R500; beta 0.6, 2/3, 0.75 (each normalised to the cluster's gas mass inside R500)
  sources in front  none; the standard share; twice the standard share (capped at 0.9); a median source redshift of 0.5 and 1.2
  outer truncation  the line of sight cut at 2, 3 and 5 R500
The reversal length Lambda enters only as B proportional to Lambda^-0.5 and is not varied here.
Run from the repository root:  .venv/bin/python -I entries/176-cluster-faraday-profile/model_variants.py     (about 6 minutes)
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


def fit(rc_frac=0.2, beta=2 / 3, zmed=0.8, front_scale=1.0, rmax_frac=3.0, nreal=40):
    rc = rc_frac * cl['r500']
    X = 1.0 / rc_frac
    shape = quad(lambda y: y * y * (1 + y * y) ** (-1.5 * beta), 0, X)[0]          # Int_0^X y^2 (1+y^2)^(-3 beta/2) dy
    n0 = cl['mgas'] * 1e11 * MSUN / (4 * np.pi * MUE * MP * (rc * KPC) ** 3 * shape)
    L = np.sqrt(np.maximum((rmax_frac * cl['r500'][ks]) ** 2 - R ** 2, 0.0))
    t = (np.arange(300) + 0.5) / 300
    ll = L[:, None] * t[None, :]
    ne = n0[ks][:, None] * (1.0 + (R[:, None] ** 2 + ll ** 2) / rc[ks][:, None] ** 2) ** (-1.5 * beta)
    I = 2.0 * ((ne / NREF) ** (2 * ETA) * ne ** 2).sum(1) * (L / 300)
    z0 = zmed / np.log(2) ** (1 / 1.2)
    pf = np.minimum(0.9, front_scale * (1.0 - np.exp(-(cl['z'][ks] / z0) ** 1.2)))
    n_at = float(np.median(n0)) * (1 + (0.75 / rc_frac) ** 2) ** (-1.5 * beta)
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
                n_e_at_0p75_R500=n_at, n0_median=float(np.median(n0)), mean_share_in_front=float(pf.mean()))


out = dict(standard=None, variants=[])
cases = [('standard', {}), ('core 0.1 R500', dict(rc_frac=0.1)), ('core 0.3 R500', dict(rc_frac=0.3)), ('beta 0.6', dict(beta=0.6)), ('beta 0.75', dict(beta=0.75)),
         ('no sources in front', dict(front_scale=0.0)), ('twice the share in front', dict(front_scale=2.0)),
         ('median source redshift 0.5', dict(zmed=0.5)), ('median source redshift 1.2', dict(zmed=1.2)),
         ('line of sight to 2 R500', dict(rmax_frac=2.0)), ('line of sight to 5 R500', dict(rmax_frac=5.0))]
for name, kw in cases:
    f = fit(**kw); f['case'] = name
    out['variants'].append(f)
    log('%-28s B at 0.75 R500 = %.2f (%.2f-%.2f) uG, chi2 %.1f; n_e there %.2e; share in front %.2f' % (
        name, f['B_at_0p75_R500'], f['lo'], f['hi'], f['chi2'], f['n_e_at_0p75_R500'], f['mean_share_in_front']))
vals = np.array([f['B_at_0p75_R500'] for f in out['variants']])
out['standard'] = out['variants'][0]
out['range'] = dict(lowest=float(vals.min()), highest=float(vals.max()))
log('range over all variants: %.2f to %.2f uG' % (vals.min(), vals.max()))
json.dump(out, open(os.path.join(HERE, 'model_variants.json'), 'w'), indent=1, default=float)
