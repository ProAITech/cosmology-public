"""Workbench Entry 129: electrons or layered orbits? A timing test for the wings of little red dots.

Entry 119 tested the Thomson-cocoon reading of little red dots with the one three-line spectrum then available. Since then:
Scholtz et al. (2026) fit 32 broad-line AGN and find exponential wings are not generally preferred, and that a broad-line region (BLR)
stratified in radius makes exponential-looking wings by itself; Sneppen et al. (2026) fit LRD spectra with a dense accreting cocoon;
Torralba et al. (2026) trace the wings of one bright LRD to electron scattering in a warm layer at about 7,000 K.
Both readings can make the same line shape. They predict opposite timing. In a stratified BLR the fastest gas is closest in, so the
wings respond first. In a scattering cocoon, wing photons have bounced off electrons and travelled farther, so the wings respond last.
This script (1) builds a stratified Keplerian BLR and checks how exponential its wings look; (2) runs a Monte Carlo of a central BLR
inside a Thomson cocoon, tracking each photon's path to get its delay; (3) gives the lag against velocity in both pictures; and
(4) recomputes the column a 7,000 K layer needs."""
import json, math, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(129)
KB = 1.380649e-23; ME = 9.1093837e-31; SIGT = 6.6524587e-25
def iso(n):
    v = rng.normal(size=(n, 3)); return v/np.linalg.norm(v, axis=1)[:, None]
def wingW(v, lo=50, hi=95):
    s = np.abs(v); a, b = np.percentile(s, [lo, hi]); m = (s >= a) & (s <= b)
    h, e = np.histogram(s[m], bins=30); c = 0.5*(e[1:] + e[:-1]); k = h > 0
    p = np.polyfit(c[k], np.log(h[k]), 1); res = np.log(h[k]) - np.polyval(p, c[k])
    pg = np.polyfit(c[k]**2, np.log(h[k]), 1); resg = np.log(h[k]) - np.polyval(pg, c[k]**2)
    return float(-1/p[0]), float(np.std(res)), float(np.std(resg))
out = {}
# 1. stratified Keplerian BLR: radii log-uniform over a factor X, line-of-sight speeds Gaussian with sigma = v1 r^-1/2
for X in (10, 30, 100):
    n = 400000; r = np.exp(rng.uniform(math.log(1/math.sqrt(X)), math.log(math.sqrt(X)), n)); sig = 1000*r**-0.5
    v = rng.normal(size=n)*sig; W, rexp, rgau = wingW(v)
    hist, e = np.histogram(np.abs(v), bins=np.linspace(0, 6000, 61)); out[f'strat_{X}'] = dict(W=W, rms_exp=rexp, rms_gauss=rgau)
    # lag against speed: delay = r (light-travel to the shell, geometric average)
    vb = np.linspace(0, 4000, 21); lag = [float(np.mean(r[(np.abs(v) >= a) & (np.abs(v) < b)])) if np.sum((np.abs(v) >= a) & (np.abs(v) < b)) > 200 else None for a, b in zip(vb[:-1], vb[1:])]
    core = float(np.mean(r[np.abs(v) < np.percentile(np.abs(v), 50)]))
    out[f'strat_{X}']['lag'] = [[float(0.5*(a+b)), (x/core if x else None)] for a, b, x in zip(vb[:-1], vb[1:], lag)]
    print(f"stratified BLR spanning x{X} in radius: exponential wing W = {W:.0f} km/s; log-residual of an exponential fit {rexp:.3f}, of a Gaussian {rgau:.3f}")
# 2. BLR shell (radius 1) inside a uniform Thomson cocoon of radius Rc and radial optical depth tau; photons tracked with path length
def cocoon(tau, Te, Rc=2.0, sig_core=300.0, n=200000, maxs=3000):
    s = math.sqrt(KB*Te/ME)/1e3; kap = tau/Rc                    # optical depth per unit length (lengths in BLR radii)
    pos = iso(n)*1.0; d = iso(n); path = np.zeros(n); shift = rng.normal(size=n)*sig_core
    t0 = np.ones(n)                                                # light-travel from the centre to the BLR shell
    ns = np.zeros(n, int); alive = np.ones(n, bool); exitp = np.zeros((n, 3)); exitd = np.zeros((n, 3))
    for it in range(maxs):
        idx = np.where(alive)[0]
        if not len(idx): break
        L = rng.exponential(size=len(idx))/kap
        p = pos[idx] + d[idx]*L[:, None]
        # distance to the cocoon surface along d
        b = np.sum(pos[idx]*d[idx], axis=1); cc = np.sum(pos[idx]**2, axis=1) - Rc**2; ds = -b + np.sqrt(b*b - cc)
        esc = L >= ds
        ie = idx[esc]; path[ie] += ds[esc]; exitp[ie] = pos[ie] + d[ie]*ds[esc][:, None]; exitd[ie] = d[ie]; alive[ie] = False
        isc = idx[~esc]; path[isc] += L[~esc]; pos[isc] = p[~esc]
        ve = rng.normal(size=(len(isc), 3))*s; dn = iso(len(isc)); shift[isc] += np.sum(ve*(dn - d[isc]), axis=1); d[isc] = dn; ns[isc] += 1
    delay = t0 + path - np.sum(exitp*exitd, axis=1)               # relative to light leaving the centre directly
    return shift, ns, delay
for tau in (3.0, 5.0):
    sh, ns, dl = cocoon(tau, 1e4)
    W, rexp, rgau = wingW(sh[ns > 0]) if np.any(ns > 0) else (None, None, None)
    core = float(np.mean(dl[np.abs(sh) < 300]))
    vb = np.linspace(0, 4000, 21); lag = []
    for a, b in zip(vb[:-1], vb[1:]):
        m = (np.abs(sh) >= a) & (np.abs(sh) < b); lag.append([float(0.5*(a+b)), float(np.mean(dl[m])/core) if m.sum() > 200 else None])
    out[f'cocoon_{tau}'] = dict(W=W, f_scat=float(np.mean(ns > 0)), core_lag=core, lag=lag)
    print(f"cocoon tau {tau}, 1e4 K, radius 2 BLR radii: scattered share {np.mean(ns>0):.2f}, wing W {W:.0f} km/s; core lag {core:.2f} BLR radii/c; lag at 1,500 and 3,000 km/s relative to the core: "
          f"{[x for v, x in lag if abs(v-1500) < 101][0]:.2f}, {[x for v, x in lag if abs(v-2900) < 101 or abs(v-3100) < 101][0]:.2f}")
for X in (30,):
    L = out[f'strat_{X}']['lag']; print(f"stratified BLR x{X}: lag at 1,500 and 3,000 km/s relative to the core: {[x for v, x in L if abs(v-1500) < 101][0]:.2f}, {[x for v, x in L if abs(v-2900) < 101 or abs(v-3100) < 101][0]:.2f}")
# 3. a warm layer at 7,000 K (Torralba et al.): W scales as sqrt(T); columns and thicknesses at n_H 1e9-1e10 cm^-3
W1e4 = {0.5: 487, 1: 580, 2: 800, 3: 1047, 5: 1559, 8: 2376}                # Entry 119 Monte Carlo, 1e4 K
ts = np.array(list(W1e4)); Ws = np.array(list(W1e4.values()))*math.sqrt(0.7)
out['warm'] = {}
for Wobs in (1000, 1500):
    t = float(np.interp(Wobs, Ws, ts)); N = t/SIGT
    out['warm'][Wobs] = dict(tau=t, N=N, L=[N/1e9, N/1e10])
    print(f"7,000 K layer: W = {Wobs} km/s needs tau {t:.1f}, an electron column of {N:.1e} cm^-2; at n_e 1e9-1e10 cm^-3 that is {N/1e10/1.496e13:.0f}-{N/1e9/1.496e13:.0f} au of fully ionized gas")
# 4. observed-frame lags for a typical LRD: BLR radius 10 light-days, z = 6
for R_ld, z in ((10, 6), (20, 5)):
    c3 = out['cocoon_3.0']; L = c3['lag']; wing = [x for v, x in L if abs(v-2900) < 101 or abs(v-3100) < 101][0]
    st = [x for v, x in out['strat_30']['lag'] if abs(v-2900) < 101 or abs(v-3100) < 101][0]
    out[f'obs_{R_ld}_{z}'] = (R_ld*(1+z), R_ld*(1+z)*c3['core_lag'], R_ld*(1+z)*c3['core_lag']*wing, R_ld*(1+z)*st)
    print(f"BLR {R_ld} light-days at z = {z}: observed core lag about {R_ld*(1+z)*c3['core_lag']:.0f} days; cocoon wings at 3,000 km/s about {R_ld*(1+z)*c3['core_lag']*wing:.0f} days; stratified wings about {R_ld*(1+z)*st:.0f} days")
json.dump(out, open(os.path.join(HERE, 'lrd_lags.json'), 'w'), indent=1)
