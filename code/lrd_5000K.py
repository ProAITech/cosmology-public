"""lrd_5000K.py - Workbench Entry 144
Entry 119's Thomson-scattering Monte Carlo (lrd_scatter.py, same photon transport) rerun at the electron temperatures that
Sneppen et al. (2026) find from Paschen and Brackett jumps in little red dots: about 5,000 K (5,400 +/- 300 K and 4,900 +/- 100 K
in the 'Rosetta Stone' dot). Question: what Thomson depth and electron column do exponential wings of 1,000-2,000 km/s then need,
and does that match the column their Balmer break implies (N_H ~ 1e25 cm^-2)?"""
import json, math, os
import numpy as np
rng = np.random.default_rng(144)
KB = 1.380649e-23; ME = 9.1093837e-31; SIGT = 6.6524587e-29
def iso(n):
    v = rng.normal(size=(n, 3)); return v/np.linalg.norm(v, axis=1)[:, None]
def run(tauR, Te, n=60000, maxs=4000):
    sig = math.sqrt(KB*Te/ME)/1e3
    pos = np.zeros((n, 3)); d = iso(n); shift = np.zeros(n); ns = np.zeros(n, int); alive = np.ones(n, bool)
    for it in range(maxs):
        idx = np.where(alive)[0]
        if not len(idx): break
        L = rng.exponential(size=len(idx)); p = pos[idx] + d[idx]*L[:, None]
        out = np.linalg.norm(p, axis=1) > tauR; alive[idx[out]] = False
        sc = idx[~out]
        if not len(sc): break
        pos[sc] = p[~out]; ve = rng.normal(size=(len(sc), 3))*sig; dn = iso(len(sc))
        shift[sc] += np.sum(ve*(dn - d[sc]), axis=1); d[sc] = dn; ns[sc] += 1
    return shift, ns
def wing(shift, ns):
    s = np.abs(shift[ns > 0]); lo, hi = np.percentile(s, [50, 95]); m = (s >= lo) & (s <= hi)
    h, e = np.histogram(s[m], bins=30); c = 0.5*(e[1:] + e[:-1]); k = h > 0
    return float(-1/np.polyfit(c[k], np.log(h[k]), 1)[0])
res = {}
for Te in (5000.0, 1e4):
    taus = [1, 2, 3, 4, 5, 6, 8, 10, 12]; Ws = []
    for t in taus:
        sh, ns = run(t, Te); Ws.append(wing(sh, ns))
    res[str(int(Te))] = dict(tau=taus, W=[round(x) for x in Ws])
    for Wobs in (1000, 1500, 2000):
        t = float(np.interp(Wobs, Ws, taus)); res[str(int(Te))][f'tau_W{Wobs}'] = round(t, 2); res[str(int(Te))][f'col_W{Wobs}'] = t/(SIGT*1e4)
    print(Te, list(zip(taus, [round(x) for x in Ws])), {k: v for k, v in res[str(int(Te))].items() if k.startswith(('tau_', 'col_'))})
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lrd_5000K.json'), 'w'), indent=1)
