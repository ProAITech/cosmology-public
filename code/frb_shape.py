"""frb_shape.py - Workbench Entry 96
Does the electron column to fast radio bursts grow with redshift as an expanding universe says, or as a static plasma does?
A burst's dispersion measure counts the free electrons on its path. In an expanding universe the mean density was (1+z)^3
higher in the past, so the intergalactic part grows as  DM_cos = K * Int_0^z (1+z')/E(z') dz'  (Macquart et al. 2020).
In a static universe filled with plasma of constant density the same electrons give DM = n_e Int dl / (1+z)^2: the local
delay goes as 1/frequency^2, and the burst was bluer by (1+z) where the delay was picked up. If time is stretched too, as the
supernova light curves require, one factor of (1+z) comes back. With dl from the distance law of each model:
   Lerner's straight line (d = cz/H0):     z/(1+z)          stretched: ln(1+z)
   tired light (z = exp(H0 d/c) - 1):     (1-(1+z)^-2)/2   stretched: z/(1+z)
The size of the column, K or n_e, is left free in every model: only the shape is tested.
Data: 92 localized bursts with host redshifts (Wang, Gao & Fan 2025, table 3) and FRB 20240304B at z = 2.148 (Caleb et al. 2025).
Likelihood as in Macquart et al. (2020): DM_ex = DM_obs - DM_MW,ISM (NE2001) - DM_halo (50) = DM_cos + DM_host,obs;
  p(DM_cos) with the published shape (alpha = beta = 3, sigma = F z^-1/2, mean = the model), DM_host log-normal in the rest frame,
  seen as DM_host / (1+z) where time is stretched and / (1+z)^2 where it is not. Free: K, F, host median and width.
  Bursts with DM_obs - DM_MW < 80 are left out, as in the source table's analysis.
"""
import json, math, os, sys
import numpy as np
from scipy import integrate, optimize
HERE = os.path.dirname(os.path.abspath(__file__))
OM = 0.3
rows = [l.split() for l in open(os.path.join(HERE, 'frb93.txt')) if not l.startswith('#') and l.strip()]
NAME = [r[0] for r in rows]; Z = np.array([float(r[1]) for r in rows]); DMO = np.array([float(r[2]) for r in rows])
MW = {'NE2001': np.array([float(r[5]) for r in rows]), 'YMW16': np.array([float(r[6]) for r in rows])}; TYP = np.array([int(r[7]) for r in rows])
zg = np.linspace(0, 3, 30001)
ex = integrate.cumulative_trapezoid((1+zg)/np.sqrt(OM*(1+zg)**3 + 1 - OM), zg, initial=0)
SHAPES = {   # name: (shape S(z), host factor exponent k: DM_host,obs = DM_host,rest / (1+z)^k)
 'expanding (flat, matter 0.3)':      (lambda z: np.interp(z, zg, ex), 1),
 "static, Lerner's straight line":    (lambda z: z/(1+z), 2),
 "static, Lerner's line, time stretched": (lambda z: np.log1p(z), 1),
 'static, tired light':               (lambda z: (1 - (1+z)**-2)/2, 2),
 'static, tired light, time stretched': (lambda z: z/(1+z), 1),
}
# Macquart p(Delta): A Delta^-3 exp(-(Delta^-3 - C0)^2 / (18 sigma^2)), with C0 set so that <Delta> = 1
_dl = np.linspace(1e-3, 8, 8000)
def _p(C0, sig):
    lp = -3*np.log(_dl) - (_dl**-3 - C0)**2/(18*sig*sig); return np.exp(lp - lp.max())
def pdelta_table(sig):
    C0 = optimize.brentq(lambda C0: np.trapezoid(_dl*_p(C0, sig), _dl)/np.trapezoid(_p(C0, sig), _dl) - 1, -60, 6)
    p = _p(C0, sig); return p/np.trapezoid(p, _dl)
SIGS = np.exp(np.linspace(math.log(0.05), math.log(6.0), 240)); PTAB = np.array([pdelta_table(sg) for sg in SIGS])   # p(Delta) for each sigma
DSTEP = _dl[1] - _dl[0]; U = np.linspace(0.0005, 0.9995, 700)
def pdelta_cache(sig): return PTAB[int(np.argmin(abs(SIGS - sig)))]
def loglike(params, S, k, dmex, z):
    lK, F, mu, sh = params
    if not (0.01 < F < 2.0 and 0.1 < sh < 2.5 and 2.0 < mu < 8.0) or np.any(dmex <= 1): return -1e30
    K = math.exp(lK); mean = K*S[:, None]
    sig = np.minimum(F/np.sqrt(z), 6.0); isg = np.clip(np.searchsorted(np.log(SIGS), np.log(sig)), 0, len(SIGS)-1)
    x = dmex[:, None]*U[None, :]
    t = (x/mean - _dl[0])/DSTEP; i0 = np.clip(t.astype(int), 0, len(_dl)-2); fr = np.clip(t - i0, 0, 1)
    rowsP = PTAB[isg]                                            # (N, len(_dl))
    p0 = np.take_along_axis(rowsP, i0, 1); p1 = np.take_along_axis(rowsP, i0+1, 1)
    pd = np.where(t < len(_dl)-1, (p0*(1-fr) + p1*fr), 0.0)/mean
    g = ((1 + z)**-k)[:, None]; h = (dmex[:, None] - x)/g
    ph = np.exp(-(np.log(h) - mu)**2/(2*sh*sh))/(h*sh*math.sqrt(2*math.pi))/g
    L = np.trapezoid(pd*ph, x, axis=1)
    return float(np.sum(np.log(np.maximum(L, 1e-300))))
def fit(name, dmex, z, Ffix=None):
    fS, k = SHAPES[name]; S = fS(z)
    if Ffix is not None:                                         # scatter held fixed
        best = None
        for lK0 in (math.log(max(np.median((dmex - 100/(1+z)**k)/S), 50.0)),):
            for dl in (0.0, 0.4, -0.4):
                rr = optimize.minimize(lambda p: -loglike([p[0], Ffix, p[1], p[2]], S, k, dmex, z), [lK0+dl, 4.5, 0.8], method='Nelder-Mead', options={'xatol': 1e-3, 'fatol': 1e-3, 'maxiter': 4000})
                if best is None or rr.fun < best.fun: best = rr
        return dict(lnL=round(-best.fun, 2), K=round(math.exp(best.x[0]), 1), F=Ffix, host_median=round(math.exp(best.x[1]), 1), host_width=round(best.x[2], 3))
    best = None
    lK1 = math.log(max(np.median((dmex - 100/(1+z)**k)/S), 50.0))
    for lK0 in (lK1, lK1 + 0.4, lK1 - 0.4):
        for F0 in (0.25, 1.0):
            r = optimize.minimize(lambda p: -loglike(p, S, k, dmex, z), [lK0, F0, 4.5, 0.8], method='Nelder-Mead', options={'xatol': 1e-3, 'fatol': 1e-3, 'maxiter': 4000})
            if best is None or r.fun < best.fun: best = r
    lK, F, mu, sh = best.x
    return dict(lnL=round(-best.fun, 2), K=round(math.exp(lK), 1), F=round(F, 3), host_median=round(math.exp(mu), 1), host_width=round(sh, 3))
def sample(mwkey='NE2001', halo=50.0, cut=80.0, drop=()):
    m = (DMO - MW[mwkey] > cut) & np.array([n not in drop for n in NAME])
    return DMO[m] - MW[mwkey][m] - halo, Z[m], [n for n, mm in zip(NAME, m) if mm]
out = {'n_table': len(Z)}
dmex, z, names = sample(); out['n_used'] = len(z); print('bursts used', len(z), 'of', len(Z))
out['main'] = {}
for nm in SHAPES:
    out['main'][nm] = fit(nm, dmex, z); print('%-40s' % nm, out['main'][nm])
ref = out['main']['expanding (flat, matter 0.3)']['lnL']
for nm in SHAPES: out['main'][nm]['dlnL_vs_expanding'] = round(out['main'][nm]['lnL'] - ref, 2)
# what the column's size means in the expanding reading: K = 3 c H0 Omega_b f_d / (8 pi G m_p), with f_d the share of baryons
# in diffuse ionized gas times electrons per baryon (0.875 for hydrogen and helium fully ionized)
c, G, mp, pc = 2.99792458e10, 6.674e-8, 1.6726e-24, 3.0857e18
for h, wb in ((0.6851, 0.0222), (0.6736, 0.02237)):
    H0 = 100*h*1e5/(pc*1e6); Ob = wb/h/h
    K1 = 3*H0*Ob*c/(8*math.pi*G*mp)/pc                      # n_b0 * c/H0 in pc cm^-3: one electron per baryon
    out.setdefault('K_full', {})['h=%.4f' % h] = round(K1*0.875, 1)
print('K if every baryon were diffuse and ionized:', out['K_full'])
# the expanding fit's column, profiled: log-likelihood at fixed K (F and host free), and the share of all baryons it implies
fS, k = SHAPES[EXN := 'expanding (flat, matter 0.3)']; Sx = fS(z); prof = []
for K in (600, 700, 760, 800, 840, 880, 955, 1050, 1150, 1300):
    bb = None
    for F0 in (0.25, 0.6, 1.0):
        rr = optimize.minimize(lambda p: -loglike([math.log(K), p[0], p[1], p[2]], Sx, k, dmex, z), [F0, 4.5, 0.9], method='Nelder-Mead', options={'xatol': 1e-3, 'fatol': 1e-3, 'maxiter': 3000})
        if bb is None or rr.fun < bb.fun: bb = rr
    prof.append([K, round(-bb.fun, 2)])
pk = np.array(prof); top = max(pk[:, 1].max(), out['main'][EXN]['lnL']); kk = np.linspace(600, 1300, 701); ll = np.interp(kk, pk[:, 0], pk[:, 1])
k2 = float(kk[np.argmax(ll >= top - 2)])
out['K_profile'] = dict(points=prof, K_at_2sigma_low=round(k2, 0), share_at_2sigma_low=round(k2/out['K_full']['h=0.6851'], 2), share_best=round(out['main'][EXN]['K']/out['K_full']['h=0.6851'], 2),
                        gain_above_all_baryons=round(top - float(np.interp(out['K_full']['h=0.6851'], pk[:, 0], pk[:, 1])), 2))
print('K profile', out['K_profile'])
# free exponent: S_q = Int (1+z)^q dz / E(z)^w  -- the static family has w = 0 and q = -2, -1, -3; test q free with w = 0 and host k = 1
def Sq(q): return lambda zz: np.interp(zz, zg, integrate.cumulative_trapezoid((1+zg)**q, zg, initial=0))
qs = np.round(np.arange(-3.0, 1.01, 0.25), 2); prof = []
for q in qs:
    SHAPES['_q'] = (Sq(q), 1); r = fit('_q', dmex, z); prof.append(r['lnL'])
prof = np.array(prof); ib = int(prof.argmax())
lo = [qq for qq, p in zip(qs, prof) if p >= prof.max() - 0.5]
cq = np.polyfit(qs[max(ib-2, 0):ib+3], prof[max(ib-2, 0):ib+3], 2); qb = -cq[1]/(2*cq[0]); qe = 1/math.sqrt(-2*cq[0])
out['free_q'] = dict(q=qs.tolist(), lnL=prof.round(2).tolist(), best=float(qs[ib]), range_1sigma=[float(min(lo)), float(max(lo))], quad_best=round(float(qb), 2), quad_err=round(float(qe), 2),
                     sigma_static={str(qq): round(math.sqrt(max(0.0, 2*(np.polyval(cq, qb) - float(prof[list(qs).index(qq)])))), 1) for qq in (-1.0, -2.0, -3.0)})
print('free exponent q:', out['free_q']['best'], out['free_q']['range_1sigma'])
del SHAPES['_q']
# robustness: other Galaxy model, halo 30 and 80, without the z = 2.15 burst, without the bursts with the largest host (repeaters in dwarfs), only z < 0.7
out['checks'] = {}
for lab, kw in (('YMW16 Galaxy', dict(mwkey='YMW16')), ('halo 30', dict(halo=30.0)), ('halo 80', dict(halo=80.0)), ('without FRB 20240304B', dict(drop=('FRB240304B',))),
                ('without repeaters in dwarfs', dict(drop=tuple(n for n, t in zip(NAME, TYP) if t == 1)))):
    d2, z2, n2 = sample(**kw); res = {}
    for nm in ('expanding (flat, matter 0.3)', "static, Lerner's line, time stretched", 'static, tired light, time stretched'):
        res[nm] = fit(nm, d2, z2)
    e = res['expanding (flat, matter 0.3)']['lnL']
    out['checks'][lab] = dict(n=len(z2), dlnL={nm: round(r['lnL'] - e, 2) for nm, r in res.items()})
    print(lab, out['checks'][lab])
res = {}
for nm in SHAPES: res[nm] = fit(nm, dmex, z, Ffix=out['main']['expanding (flat, matter 0.3)']['F'])
e = res['expanding (flat, matter 0.3)']['lnL']; out['checks']['scatter held at the expanding value'] = dict(n=len(z), dlnL={nm: round(r['lnL'] - e, 2) for nm, r in res.items()})
print('F fixed', out['checks']['scatter held at the expanding value'])
m = z < 0.7; res = {}
for nm in ('expanding (flat, matter 0.3)', "static, Lerner's line, time stretched", 'static, tired light, time stretched'): res[nm] = fit(nm, dmex[m], z[m])
e = res['expanding (flat, matter 0.3)']['lnL']; out['checks']['only z < 0.7'] = dict(n=int(m.sum()), dlnL={nm: round(r['lnL'] - e, 2) for nm, r in res.items()})
print('only z < 0.7', out['checks']['only z < 0.7'])
# binned medians for the figure
edges = [0.0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.9, 1.5, 2.5]; bins = []
for a, b in zip(edges[:-1], edges[1:]):
    s = (z >= a) & (z < b)
    if s.sum(): bins.append([round(float(np.median(z[s])), 3), round(float(np.median(dmex[s])), 1), int(s.sum()), round(float(np.percentile(dmex[s], 25)), 1), round(float(np.percentile(dmex[s], 75)), 1)])
out['bins'] = bins
out['points'] = [[round(float(a), 4), round(float(b), 1)] for a, b in zip(z, dmex)]
# model medians for the figure: median of cosmic + host at each z (Monte Carlo of the fitted distributions)
rng = np.random.default_rng(1); zz = np.round(np.concatenate([np.linspace(0.02, 0.3, 15), np.linspace(0.35, 2.3, 40)]), 3); curves = {}
for nm in SHAPES:
    p = out['main'][nm]; fS, k = SHAPES[nm]; cur = []
    for zi in zz:
        sig = min(p['F']/math.sqrt(zi), 6.0); tab = pdelta_cache(sig); cdf = np.cumsum(tab); cdf /= cdf[-1]
        dlt = np.interp(rng.random(40000), cdf, _dl); host = np.exp(rng.normal(math.log(p['host_median']), p['host_width'], 40000))/(1+zi)**k
        cur.append([float(zi), round(float(np.median(p['K']*float(fS(np.array([zi]))[0])*dlt + host)), 1)])
    curves[nm] = cur
out['curves'] = curves
json.dump(out, open(os.path.join(HERE, 'frb_shape.json'), 'w'), indent=1)
