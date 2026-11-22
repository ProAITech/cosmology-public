"""Entry 176, correction to model.py: the field fit with the cluster term's own sample variance in the errors.

model.py fitted the field with chi^2 errors taken from the rotated clusters, the scatter of V where there is no cluster field.
Entry 178's synthetic test (entries/178-cluster-field-evolution/) showed that this understates the error where the cluster
term dominates. Path lengths and densities differ widely between sources, so V of a few hundred sources scatters far more with
the field than without it. The fit's intervals were therefore too narrow. Here the same model is refitted with:
    error^2 (range j, field B_0) = rotated-cluster variance + [variance of V over the draws at B_0 - that at B_0 = 0],
and -2 ln L = sum (obs - pred)^2 / error^2 + ln error^2 (Gaussian, errors depending on B_0). The same draws are used for every
B_0 (common random numbers), so the curves are smooth. Everything else is model.py's: gas from MGAS500, beta model with
r_c = 0.2 R500, Lambda = 10 kpc, front share F(z), member shares, ranges 0.5-1, 1-1.5, 1.5-2 R500 fitted, the four mass quarters
at 0.5-1 R500.

Run from the repository root:  .venv/bin/python -I entries/176-cluster-faraday-profile/model_errors.py     (about 10 minutes)
"""
import json, os, sys
import numpy as np
from scipy.spatial import cKDTree
from astropy.io import fits
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cluster_rm as c

log = lambda s: print(s, flush=True)
KPC, MSUN, MP, MUE, MU, G = 3.0857e21, 1.989e33, 1.6726e-24, 1.17, 0.59, 6.674e-8
RC_FRAC, LAMBDA, NREF, NREAL = 0.2, 10.0, 1e-3, 300
res = json.load(open(os.path.join(HERE, 'results.json')))
prof = res['reported_only']['profile']
model = json.load(open(os.path.join(HERE, 'model.json')))
expl = json.load(open(os.path.join(HERE, 'explore.json')))
obs = np.array([p['excess_variance'] for p in prof[:4]]); nsd = np.array([p['null_sd'] for p in prof[:4]])
f_mem = np.array([m['member_share'] for m in model['members']['ranges'][:4]])

rows, src, _ = c.load_spice()
S, _ = c.sample(src, c.CFG['fg_radius'], c.CFG['min_neigh'])
with fits.open(c.CLUSTERS, memmap=True) as h:
    d = h[1].data
    ra, dec, z = np.asarray(d['RA'], float), np.asarray(d['DEC'], float), np.asarray(d['BEST_Z'], float)
    pc, r500, m500, mgas = (np.asarray(d[k], float) for k in ('PCONT', 'R500', 'M500', 'MGAS500'))
ok = np.isfinite(ra) & np.isfinite(dec) & np.isfinite(z) & np.isfinite(r500) & (pc < 0.5) & (z >= 0.05) & (r500 > 0)
l, b = c.galactic(ra[ok], dec[ok])
th = r500[ok] / 1000.0 / c.ang_diam_Mpc(z[ok]) / c.D2R
mg = mgas[ok]; good = np.isfinite(mg) & (mg > 0)
fg_med = float(np.nanmedian(mg[good] * 1e11 / (m500[ok][good] * 1e13)))
mg = np.where(good, mg, fg_med * m500[ok] * 100.0)
cl = dict(l=l, b=b, th=th, z=z[ok], m500=m500[ok], r500=r500[ok], mgas=mg)
rc = RC_FRAC * cl['r500']; X = 1.0 / RC_FRAC
n0 = cl['mgas'] * 1e11 * MSUN / (4 * np.pi * MUE * MP * (rc * KPC) ** 3 * (X - np.arctan(X)))
tree = cKDTree(c.gal_unit(S['l'], S['b']))
x, which = c.scaled_distance(tree, len(S['d']), cl)
field = ~np.isfinite(x)
dF, eF = S['d'][field], S['e'][field]
V_field = c.V_of(dF, eF)
sel = x < 2.0
xs, ks = x[sel], which[sel]
R = xs * cl['r500'][ks]
binidx = np.digitize(xs, [0.5, 1.0, 1.5])
p_front = 1.0 - np.exp(-(cl['z'][ks] / 1.085) ** 1.2)


def los_integral(eta, nstep=400):
    rmax = 3.0 * cl['r500'][ks]
    L = np.sqrt(np.maximum(rmax ** 2 - R ** 2, 0.0))
    t = (np.arange(nstep) + 0.5) / nstep
    ne = n0[ks][:, None] / (1.0 + (R[:, None] ** 2 + (L[:, None] * t[None, :]) ** 2) / rc[ks][:, None] ** 2)
    return 2.0 * ((ne / NREF) ** (2 * eta) * ne ** 2).sum(1) * (L / nstep)


g = np.random.default_rng(c.CFG['seed'] + 1176)
K = g.integers(0, len(dF), (NREAL, len(xs))); U = g.uniform(size=(NREAL, len(xs))); Z0 = g.normal(size=(NREAL, len(xs)))
FACT = np.where(U < f_mem[binidx][None, :], 0.5, np.where(U < f_mem[binidx][None, :] + (1 - f_mem[binidx][None, :]) * p_front[None, :], 0.0, 1.0))


def sims(B0, I, mask=None, bins=range(4)):
    sig2 = 812.0 ** 2 * LAMBDA * (B0 ** 2 / 3.0) * I * (1.0 + cl['z'][ks]) ** -4
    ds = dF[K] + Z0 * np.sqrt(sig2[None, :] * FACT)
    out = np.empty((NREAL, len(bins)))
    for jj, j in enumerate(bins):
        m = (binidx == j) if mask is None else ((binidx == j) & mask)
        out[:, jj] = [c.V_of(ds[r][m], eF[K[r][m]]) - V_field for r in range(NREAL)]
    return out


grid = np.round(np.concatenate([np.arange(0.0, 8.0, 0.2), np.arange(8.0, 20.0, 0.5), np.arange(20.0, 40.01, 2.0)]), 2)
n_med = float(np.median(n0)); xq = 0.75
n_at = lambda xx: n_med / (1 + (xx / RC_FRAC) ** 2)
kT_med = MU * MP * G * float(np.median(cl['m500'])) * 1e13 * MSUN / (2 * float(np.median(cl['r500'])) * KPC)
USE = [1, 2, 3]


def fit(eta, I, label):
    m2l, m2l_old, preds = [], [], []
    for B0 in grid:
        v = sims(B0, I)
        mu, sd = v.mean(0), v.std(0, ddof=1)
        preds.append((mu, sd))
    sd0 = preds[0][1]
    for mu, sd in preds:
        e2 = nsd ** 2 + np.maximum(sd ** 2 - sd0 ** 2, 0.0)
        m2l.append(float(((obs[USE] - mu[USE]) ** 2 / e2[USE] + np.log(e2[USE] / nsd[USE] ** 2)).sum()))
        m2l_old.append(float((((obs[USE] - mu[USE]) / nsd[USE]) ** 2).sum()))
    m2l, m2l_old = np.array(m2l), np.array(m2l_old)
    conv = (n_at(xq) / NREF) ** eta
    out = {}
    for name, curve in (('corrected', m2l), ('rotation_errors_only', m2l_old)):
        i0 = int(np.argmin(curve)); in1 = grid[curve <= curve[i0] + 1.0]; in2 = grid[curve <= curve[i0] + 3.84]
        out[name] = dict(B0=float(grid[i0]), B0_68=[float(in1.min()), float(in1.max())], B0_95=[float(in2.min()), float(in2.max())],
                         B_0p75=float(grid[i0] * conv), B_0p75_68=[float(in1.min() * conv), float(in1.max() * conv)],
                         B_0p75_95=[float(in2.min() * conv), float(in2.max() * conv)], m2lnL=float(curve[i0]))
    i0 = int(np.argmin(m2l))
    mu, sd = preds[i0]
    out['errors_at_best'] = (np.sqrt(nsd ** 2 + np.maximum(sd ** 2 - sd0 ** 2, 0.0))).tolist()
    out['predicted_at_best'] = mu.tolist()
    Bg = out['corrected']['B_0p75'] * 1e-6
    out['plasma_beta_0p75'] = float(1.92 * n_at(xq) * kT_med / (Bg ** 2 / (8 * np.pi))) if Bg > 0 else None
    c_ = out['corrected']; o_ = out['rotation_errors_only']
    log('%s eta %.1f: B(0.75 R500) = %.2f uG, 68%% %.2f-%.2f, 95%% %.2f-%.2f (rotation errors only: %.2f, 68%% %.2f-%.2f); errors at best %s vs rotations %s' % (
        label, eta, c_['B_0p75'], *c_['B_0p75_68'], *c_['B_0p75_95'], o_['B_0p75'], *o_['B_0p75_68'],
        [round(v, 1) for v in out['errors_at_best']], [round(v, 1) for v in nsd]))
    return out


out = dict(note=__doc__.split('\n\n')[0], nreal=NREAL, fits={})
Is = {}
for eta in (0.0, 0.5, 1.0):
    Is[eta] = los_integral(eta)
    out['fits']['eta_%.1f' % eta] = fit(eta, Is[eta], 'all clusters')

# mass quarters, 0.5-1 R500, eta = 0.5
q = np.nanpercentile(cl['m500'], [0, 25, 50, 75, 100])
mq = []
for i in range(4):
    msel = (cl['m500'][ks] >= q[i]) & ((cl['m500'][ks] <= q[i + 1]) if i == 3 else (cl['m500'][ks] < q[i + 1]))
    meas = expl['mass'][i]['ranges'][1]
    curve, curve_old, pr = [], [], []
    for B0 in grid:
        v = sims(B0, Is[0.5], mask=msel, bins=[1])[:, 0]
        pr.append((v.mean(), v.std(ddof=1)))
    sd0 = pr[0][1]
    for mu, sd in pr:
        e2 = meas['null_sd'] ** 2 + max(sd ** 2 - sd0 ** 2, 0.0)
        curve.append((meas['excess_variance'] - mu) ** 2 / e2 + np.log(e2 / meas['null_sd'] ** 2))
        curve_old.append(((meas['excess_variance'] - mu) / meas['null_sd']) ** 2)
    curve = np.array(curve); i0 = int(np.argmin(curve)); in1 = grid[curve <= curve[i0] + 1.0]
    mcl = (cl['m500'] >= q[i]) & ((cl['m500'] <= q[i + 1]) if i == 3 else (cl['m500'] < q[i + 1]))
    nq = float(np.median(n0[mcl])) / (1 + (xq / RC_FRAC) ** 2)
    kTq = MU * MP * G * float(np.median(cl['m500'][mcl])) * 1e13 * MSUN / (2 * float(np.median(cl['r500'][mcl])) * KPC)
    conv = (nq / NREF) ** 0.5
    Bq = grid[i0] * conv
    pb = lambda B: float(1.92 * nq * kTq / ((B * 1e-6) ** 2 / (8 * np.pi))) if B > 0 else None
    mq.append(dict(quarter=i + 1, m500_median_1e13=expl['mass'][i]['m500_median_1e13'], measured=meas['excess_variance'], rotation_error=meas['null_sd'],
                   B_0p75=float(Bq), B_0p75_68=[float(in1.min() * conv), float(in1.max() * conv)], plasma_beta=pb(Bq),
                   plasma_beta_68=[pb(in1.max() * conv), pb(in1.min() * conv)]))
    log('mass quarter %d: B(0.75 R500) %.2f uG (68%% %.2f-%.2f); plasma beta %s (68%% %s)' % (
        i + 1, Bq, in1.min() * conv, in1.max() * conv, mq[-1]['plasma_beta'] and round(mq[-1]['plasma_beta'], 1),
        [round(v, 1) if v else None for v in mq[-1]['plasma_beta_68']]))
out['mass_quarters'] = mq
json.dump(out, open(os.path.join(HERE, 'model_errors.json'), 'w'), indent=1)
log('written model_errors.json')
