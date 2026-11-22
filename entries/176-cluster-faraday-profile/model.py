"""Entry 176, analysis after the frozen result: is the excess in the clusters' own radio galaxies, and what field does it imply?

Not a pre-registered test; the data were already open.  Two parts.

1. Members.  A cluster's own radio galaxies add to the source count near it; plasma in front of unrelated sources does not.
   The count in each range is set against the rotated clusters (corrected for the share of rotated clusters that leave the
   SPICE-RACS footprint, measured at 4-8 R500).  Where there is no surplus of sources, the excess scatter cannot be the members'.
   And a median-based scatter can rise only so far if a share f of the sources is replaced by sources of unlimited scatter:
   that ceiling is computed from the field sources.

2. Field strength.  Each cluster gets a gas density n_e(r) = n_0 / (1 + r^2/r_c^2) (a beta model, beta = 2/3, r_c = 0.2 R500),
   normalised to its own gas mass inside R500 from the eRASS1 catalogue (MGAS500).  The field is random with reversal length
   Lambda and strength B(r) = B_0 (n_e / 1e-3 cm^-3)^eta.  A source at projected radius R behind the cluster then has
       sigma_RM^2 = 812^2 Lambda (B_0^2 / 3) Int (n_e/1e-3)^(2 eta) n_e^2 dl (1 + z)^-4     [kpc, cm^-3, microgauss]
   (the standard random-cell formula).  Sources in front of the cluster get nothing; members get half the integral.  The
   chance that a source is in front follows an assumed redshift distribution of radio sources, F(z) = 1 - exp(-(z/1.085)^1.2),
   median 0.8; the share of members in each range comes from part 1.  For each B_0 and eta, synthetic dRM values are made by
   adding the cluster term to dRM values drawn from field sources, and the same statistic V is computed in the four ranges
   out to 2 R500.  chi^2 against the measured excesses (errors: the rotated-cluster standard deviations) gives B_0.
   The range inside 0.5 R500 is not fitted: 40% of its sources are members and the survivors there are the sources that did
   not depolarize.  The fit uses 0.5 to 2 R500 and is quoted as the field at 0.75 R500.

Run from the repository root:  .venv/bin/python -I entries/176-cluster-faraday-profile/model.py     (about 5 minutes)
"""
import json, os, sys
import numpy as np
from scipy.spatial import cKDTree
from astropy.io import fits
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cluster_rm as c

log = lambda s: print(s, flush=True)
KPC, MSUN, MP, MUE = 3.0857e21, 1.989e33, 1.6726e-24, 1.17
RC_FRAC, LAMBDA, NREF = 0.2, 10.0, 1e-3
res = json.load(open(os.path.join(HERE, 'results.json')))
prof = res['reported_only']['profile']
BINS = [[0.0, 0.5], [0.5, 1.0], [1.0, 1.5], [1.5, 2.0]]
obs = np.array([p['excess_variance'] for p in prof[:4]]); err = np.array([p['null_sd'] for p in prof[:4]])

# ---------------------------------------------------------------- data
rows, src, _ = c.load_spice()
S, _ = c.sample(src, c.CFG['fg_radius'], c.CFG['min_neigh'])
with fits.open(c.CLUSTERS, memmap=True) as h:
    d = h[1].data
    ra, dec, z = np.asarray(d['RA'], float), np.asarray(d['DEC'], float), np.asarray(d['BEST_Z'], float)
    pc, r500, m500, mgas = np.asarray(d['PCONT'], float), np.asarray(d['R500'], float), np.asarray(d['M500'], float), np.asarray(d['MGAS500'], float)
ok = np.isfinite(ra) & np.isfinite(dec) & np.isfinite(z) & np.isfinite(r500) & (pc < 0.5) & (z >= 0.05) & (r500 > 0)
l, b = c.galactic(ra[ok], dec[ok])
th = r500[ok] / 1000.0 / c.ang_diam_Mpc(z[ok]) / c.D2R
mg = mgas[ok]; fgas = mg * 1e11 / (m500[ok] * 1e13)
good = np.isfinite(mg) & (mg > 0)
fg_med = float(np.nanmedian(fgas[good]))
mg = np.where(good, mg, fg_med * m500[ok] * 100.0)          # no gas mass listed: the median gas fraction times M500
cl = dict(l=l, b=b, th=th, z=z[ok], m500=m500[ok], r500=r500[ok], mgas=mg)
log('clusters %d; gas mass listed for %d; median gas fraction %.3f; median n_0 below' % (len(th), good.sum(), fg_med))
rc = RC_FRAC * cl['r500']                                    # kpc
X = 1.0 / RC_FRAC
n0 = cl['mgas'] * 1e11 * MSUN / (4 * np.pi * MUE * MP * (rc * KPC) ** 3 * (X - np.arctan(X)))     # cm^-3
log('n_0: median %.2e cm^-3 (quartiles %.2e, %.2e); n_e at R500: median %.2e' % (np.median(n0), *np.percentile(n0, [25, 75]), np.median(n0) / (1 + X * X)))

tree = cKDTree(c.gal_unit(S['l'], S['b']))
n = len(S['d'])
x, which = c.scaled_distance(tree, n, cl)
field = ~np.isfinite(x)
dF, eF = S['d'][field], S['e'][field]
V_field = c.V_of(dF, eF)
out = dict(config=dict(rc_over_r500=RC_FRAC, lambda_kpc=LAMBDA, n_ref=NREF, beta=2 / 3, clusters=int(len(th)), gas_mass_listed=int(good.sum()),
                       median_gas_fraction=fg_med, n0_median=float(np.median(n0))))

# ---------------------------------------------------------------- 1. members
foot = sum(p['n'] for p in prof[6:8]) / sum(p['null_n_mean'] for p in prof[6:8])
absdev = np.sort(np.abs(dF - np.median(dF)))
mem = []
for p in prof[:6]:
    expect = p['null_n_mean'] * foot
    f = max(0.0, 1.0 - expect / p['n'])
    surplus_sigma = (p['n'] - expect) / np.sqrt(expect)
    q = min(0.999999, 0.5 / (1.0 - f)) if f < 0.5 else 0.999999
    ceiling = (1.4826 * absdev[int(q * (len(absdev) - 1))]) ** 2 - float(np.mean(eF ** 2)) - V_field
    mem.append(dict(lo=p['lo'], hi=p['hi'], n=p['n'], expected_without_members=float(expect), surplus=float(p['n'] - expect), surplus_sigma=float(surplus_sigma),
                    member_share=float(f), most_members_could_add=float(ceiling) if f < 0.5 else None, excess_measured=p['excess_variance']))
    log('members x %.1f-%.1f: %d sources, %.0f expected (surplus %+.0f, %+.1f sigma), member share %.2f; members could add at most %s; measured %+.1f' % (
        p['lo'], p['hi'], p['n'], expect, p['n'] - expect, surplus_sigma, f, ('%.1f' % ceiling) if f < 0.5 else 'no limit', p['excess_variance']))
out['members'] = dict(footprint_factor=float(foot), ranges=mem)
f_mem = np.array([m['member_share'] for m in mem[:4]])

# ---------------------------------------------------------------- 2. field strength
sel = x < 2.0
xs, ks = x[sel], which[sel]
R = xs * cl['r500'][ks]                                      # projected radius, kpc
def los_integral(eta, nstep=400):
    """Int (n_e/NREF)^(2 eta) n_e^2 dl in cm^-6 kpc, along the line of sight out to 3 R500 from the centre."""
    rmax = 3.0 * cl['r500'][ks]
    L = np.sqrt(np.maximum(rmax ** 2 - R ** 2, 0.0))
    t = (np.arange(nstep) + 0.5) / nstep
    ll = L[:, None] * t[None, :]
    ne = n0[ks][:, None] / (1.0 + (R[:, None] ** 2 + ll ** 2) / rc[ks][:, None] ** 2)
    return 2.0 * ((ne / NREF) ** (2 * eta) * ne ** 2).sum(1) * (L / nstep)
p_front = 1.0 - np.exp(-(cl['z'][ks] / 1.085) ** 1.2)
binidx = np.digitize(xs, [0.5, 1.0, 1.5])                    # 0..3
g = np.random.default_rng(c.CFG['seed'] + 177)
NREAL = 80
def predict(B0, eta, I, front=True):
    sig2 = 812.0 ** 2 * LAMBDA * (B0 ** 2 / 3.0) * I * (1.0 + cl['z'][ks]) ** -4
    acc = np.zeros(4)
    for r_ in range(NREAL):
        k = g.integers(0, len(dF), len(xs))
        u = g.uniform(size=len(xs))
        pf = p_front if front else np.zeros(len(xs))
        fm = f_mem[binidx]
        factor = np.where(u < fm, 0.5, np.where(u < fm + (1 - fm) * pf, 0.0, 1.0))
        dsim = dF[k] + g.normal(0, 1, len(xs)) * np.sqrt(sig2 * factor)
        for j in range(4):
            m = binidx == j
            acc[j] += c.V_of(dsim[m], eF[k][m]) - V_field
    return acc / NREAL
# The inner range (x < 0.5) is left out of the fit: 40% of its sources are cluster members, and its surviving sources are
# the ones that did not depolarize (their fractional polarization is 0.74 of the field's).  The fit uses 0.5-1, 1-1.5 and 1.5-2 R500,
# where there is no surplus of sources and little depolarization, and then says what it would predict for the inner range.
USE = [1, 2, 3]
grid = np.round(np.concatenate([np.arange(0.2, 8.0, 0.2), np.arange(8.0, 20.0, 0.5), np.arange(20.0, 40.01, 2.0)]), 2)
n_med = float(np.median(n0)); xq = 0.75
n_at = lambda xx: n_med / (1 + (xx / RC_FRAC) ** 2)
def run_fit(eta, front=True, label=''):
    I = los_integral(eta)
    chi, preds = [], []
    for B0 in grid:
        pr = predict(B0, eta, I, front=front)
        preds.append(pr); chi.append(float((((obs[USE] - pr[USE]) / err[USE]) ** 2).sum()))
    chi = np.array(chi); i0 = int(np.argmin(chi))
    inside = grid[chi <= chi[i0] + 1.0]
    conv = (n_at(xq) / NREF) ** eta                           # B at 0.75 R500 of the median cluster, per unit B_0
    f = dict(eta=eta, sources_in_front=front, B0=float(grid[i0]), B0_lo=float(inside.min()), B0_hi=float(inside.max()), chi2=float(chi[i0]), dof=len(USE) - 1,
             B_at_0p75_R500=float(grid[i0] * conv), B_at_0p75_lo=float(inside.min() * conv), B_at_0p75_hi=float(inside.max() * conv),
             predicted=[float(v) for v in preds[i0]])
    log('%seta %.1f: B_0 = %.2f (%.2f-%.2f) uG at n_e = 1e-3; B at 0.75 R500 = %.2f (%.2f-%.2f) uG; chi2 %.1f for %d; predicted %s, measured %s' % (
        label, eta, f['B0'], f['B0_lo'], f['B0_hi'], f['B_at_0p75_R500'], f['B_at_0p75_lo'], f['B_at_0p75_hi'], f['chi2'], f['dof'],
        [round(float(v), 1) for v in preds[i0]], [round(float(v), 1) for v in obs]))
    return f, I
fits_ = []
Is = {}
for eta in (0.0, 0.5, 1.0):
    f, Is[eta] = run_fit(eta)
    fits_.append(f)
nof, _ = run_fit(0.5, front=False, label='no sources in front: ')
out['field'] = dict(measured=obs.tolist(), errors=err.tolist(), bins=BINS, ranges_fitted=[BINS[j] for j in USE], fits=fits_, no_sources_in_front=nof,
                    n_e_at_0p75_R500=float(n_at(xq)), note='B scales as (Lambda / 10 kpc)^-0.5')
best = [f for f in fits_ if f['eta'] == 0.5][0]
# plasma beta at 0.75 R500 for the median cluster: thermal pressure over magnetic pressure, with kT = mu m_p G M500 / (2 R500)
G, MU = 6.674e-8, 0.59
M_med, R_med = float(np.median(cl['m500'])) * 1e13 * MSUN, float(np.median(cl['r500'])) * KPC
kT = MU * MP * G * M_med / (2 * R_med)                        # erg
P_th = 1.92 * n_at(xq) * kT
for f in fits_ + [nof]:
    Bg = f['B_at_0p75_R500'] * 1e-6
    f['plasma_beta_at_0p75_R500'] = float(P_th / (Bg ** 2 / (8 * np.pi)))
out['field']['kT500_keV_median'] = float(kT / 1.602e-9)
log('kT500 of the median cluster %.2f keV; n_e at 0.75 R500 %.2e; plasma beta: %s' % (kT / 1.602e-9, n_at(xq), [(f['eta'], round(f['plasma_beta_at_0p75_R500'], 1)) for f in fits_]))
# mass quarters, range 0.5-1 R500 only, eta = 0.5: the B_0 each quarter needs
e = json.load(open(os.path.join(HERE, 'explore.json')))
q = np.nanpercentile(cl['m500'], [0, 25, 50, 75, 100])
I05 = Is[0.5]
mq = []
for i in range(4):
    msel = (cl['m500'][ks] >= q[i]) & ((cl['m500'][ks] <= q[i + 1]) if i == 3 else (cl['m500'][ks] < q[i + 1])) & (xs >= 0.5) & (xs < 1.0)
    meas = e['mass'][i]['ranges'][1]
    chi = []
    for B0 in grid:
        sig2 = 812.0 ** 2 * LAMBDA * (B0 ** 2 / 3.0) * I05[msel] * (1.0 + cl['z'][ks][msel]) ** -4
        acc = 0.0
        for r_ in range(150):
            k = g.integers(0, len(dF), int(msel.sum())); u = g.uniform(size=int(msel.sum()))
            acc += c.V_of(dF[k] + g.normal(0, 1, len(k)) * np.sqrt(sig2 * (u >= p_front[msel])), eF[k]) - V_field
        chi.append(((meas['excess_variance'] - acc / 150) / meas['null_sd']) ** 2)
    chi = np.array(chi); i0 = int(np.argmin(chi)); inside = grid[chi <= chi[i0] + 1.0]
    mcl = (cl['m500'] >= q[i]) & ((cl['m500'] <= q[i + 1]) if i == 3 else (cl['m500'] < q[i + 1]))
    nq = float(np.median(n0[mcl])) / (1 + (xq / RC_FRAC) ** 2)
    kTq = MU * MP * G * float(np.median(cl['m500'][mcl])) * 1e13 * MSUN / (2 * float(np.median(cl['r500'][mcl])) * KPC)
    Bq = grid[i0] * (nq / NREF) ** 0.5
    mq.append(dict(quarter=i + 1, m500_median_1e13=e['mass'][i]['m500_median_1e13'], sources=int(msel.sum()), measured=meas['excess_variance'],
                   error=meas['null_sd'], B0=float(grid[i0]), B0_lo=float(inside.min()), B0_hi=float(inside.max()),
                   B_at_0p75_R500=float(Bq), B_at_0p75_lo=float(inside.min() * (nq / NREF) ** 0.5), B_at_0p75_hi=float(inside.max() * (nq / NREF) ** 0.5),
                   n_e_at_0p75_R500=nq, kT500_keV=float(kTq / 1.602e-9),
                   plasma_beta=float(1.92 * nq * kTq / ((Bq * 1e-6) ** 2 / (8 * np.pi))) if Bq > 0 else None))
    log('mass quarter %d (%.1fe13 Msun, kT500 %.1f keV): 0.5-1 R500 excess %+.0f +- %.0f from %d sources; B at 0.75 R500 = %.2f (%.2f-%.2f) uG; plasma beta %s' % (
        i + 1, mq[-1]['m500_median_1e13'], mq[-1]['kT500_keV'], meas['excess_variance'], meas['null_sd'], msel.sum(), Bq, mq[-1]['B_at_0p75_lo'], mq[-1]['B_at_0p75_hi'],
        ('%.0f' % mq[-1]['plasma_beta']) if mq[-1]['plasma_beta'] else 'undefined'))
out['mass_quarters'] = mq
# what the fit on the outer ranges predicts for the inner one, against what is measured there
out['inner_range'] = dict(measured=float(obs[0]), predicted_by_outer_fit={('eta_%.1f' % f['eta']): f['predicted'][0] for f in fits_},
                          note='the fit to 0.5-2 R500 predicts far more scatter inside 0.5 R500 than is measured: the missing sources are the depolarized ones')
log('inner range: measured %+.0f; the outer fit predicts %s' % (obs[0], out['inner_range']['predicted_by_outer_fit']))
json.dump(out, open(os.path.join(HERE, 'model.json'), 'w'), indent=1, default=float)
