"""Entry 87: is cluster gas held and heated by gravity, or by fields and currents?

beta_spec = mu m_p sigma^2 / kT compares the galaxies' line-of-sight velocity dispersion
(galaxies feel gravity only) with the gas temperature. One gravitational well holding both
gives beta near 1 (gravity_beta.py: 0.70-1.31 for sigma inside r200).
Data: SDSS DR18 dispersions measured here (sdss_sigma.js) for clusters with X-ray
temperatures from Migkas et al. 2020; and the Wu, Fang & Xu 1998 compilation (149 clusters).
"""
import json, math, numpy as np
from scipy import optimize, stats
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); os.chdir(os.path.dirname(os.path.abspath(__file__)))
import gravity_beta as gb
K = 0.6*1.67262e-24*1e10/1.60218e-9            # beta = K sigma[km/s]^2 / T[keV], mu = 0.6
rng = np.random.default_rng(87)

mig = {r.split(',')[0]: r.split(',') for r in open('migkas_cand.txt').read().split(';')}
rows = []
for r in open('sdss_sigma.txt').read().split(';'):
    f = r.split(',')
    name = f[0]; m = mig[name]; T, ET, eT = float(m[2]), float(m[3]), float(m[4])
    if f[3] == '' or f[5] == '': continue
    n, sig, lo, hi = int(f[3]), float(f[5]), float(f[6]), float(f[7])
    n1 = int(f[8]) if f[8] else 0; s1 = float(f[9]) if f[9] else None
    rows.append(dict(name=name, z=float(m[1]), T=T, dlnT=(ET+eT)/(2*T), n=n, sig=sig, dlns=(hi-lo)/(2*sig),
                     r200=float(f[2]), n1=n1, sig1=s1, sig600=float(f[13]) if f[13] else None))
sel = [r for r in rows if r['n'] >= 20]
for r in sel:
    r['beta'] = K*r['sig']**2/r['T']; r['dlnb'] = math.hypot(2*r['dlns'], r['dlnT'])
lnb = np.array([math.log(r['beta']) for r in sel]); e = np.array([r['dlnb'] for r in sel])
lnT = np.log([r['T'] for r in sel]); lns = np.log([r['sig'] for r in sel])
dlnT = np.array([r['dlnT'] for r in sel]); dlns = np.array([r['dlns'] for r in sel])

def robust(x):
    med = np.median(x); mad = 1.4826*np.median(abs(x-med)); return med, mad
med, mad = robust(lnb)
out_mask = abs(lnb-med) > 3*mad
outliers = [(sel[i]['name'], round(sel[i]['beta'], 2), sel[i]['n'], sel[i]['sig'], sel[i]['T']) for i in np.where(out_mask)[0]]
keep = ~out_mask

def boot_median(x, nb=4000):
    b = [np.median(rng.choice(x, len(x))) for _ in range(nb)]; return float(np.std(b))
def intrinsic(x, ex):
    # ML mean and intrinsic scatter of ln beta given measurement errors
    def nll(p):
        mu, s = p; v = ex**2+s**2; return 0.5*np.sum((x-mu)**2/v + np.log(v))
    p = optimize.minimize(nll, [np.mean(x), 0.2], method='Nelder-Mead').x
    return float(p[0]), float(abs(p[1]))
def fit_slope(x, ex, y, ey):
    # ln sigma = a + b ln(T/5 keV), effective variance with intrinsic scatter in ln sigma
    x = x-np.log(5)
    def nll(p):
        a, b, s = p; v = ey**2 + (b*ex)**2 + s**2; return 0.5*np.sum((y-a-b*x)**2/v + np.log(v))
    p = optimize.minimize(nll, [6.7, 0.5, 0.1], method='Nelder-Mead', options=dict(maxiter=4000, xatol=1e-6, fatol=1e-8)).x
    bs = []
    for _ in range(600):
        i = rng.integers(0, len(x), len(x))
        def nl2(q):
            a, b, s = q; v = ey[i]**2+(b*ex[i])**2+s**2; return 0.5*np.sum((y[i]-a-b*x[i])**2/v+np.log(v))
        bs.append(optimize.minimize(nl2, p, method='Nelder-Mead', options=dict(maxiter=3000)).x)
    bs = np.array(bs)
    return dict(a=float(p[0]), b=float(p[1]), scat=float(abs(p[2])), sig5=float(math.exp(p[0])),
                b_err=float(np.std(bs[:, 1])), sig5_err=float(math.exp(p[0])*np.std(bs[:, 0])))

res = {}
res['n_clusters'] = len(sel); res['n_fetched'] = len(rows)
res['beta_median_all'] = round(math.exp(med), 3); res['beta_median_all_err'] = round(math.exp(med)*boot_median(lnb), 3)
res['outliers'] = outliers
mu, s = intrinsic(lnb[keep], e[keep])
res['beta_mean_clean'] = round(math.exp(mu), 3); res['intrinsic_scatter_lnbeta'] = round(s, 3)
res['beta_median_clean'] = round(float(np.exp(np.median(lnb[keep]))), 3)
res['beta_median_clean_err'] = round(float(np.exp(np.median(lnb[keep]))*boot_median(lnb[keep])), 3)
res['fit_clean'] = fit_slope(lnT[keep], dlnT[keep], lns[keep], dlns[keep])
res['fit_all'] = fit_slope(lnT, dlnT, lns, dlns)
rho, p = stats.spearmanr(lnT, lns); res['spearman'] = [round(float(rho), 3), float('%.2g' % p)]
for lab, lo_, hi_ in (('T<3', 0, 3), ('3-6', 3, 6), ('T>6', 6, 99)):
    m_ = keep & (np.exp(lnT) >= lo_) & (np.exp(lnT) < hi_)
    res['beta_' + lab] = [int(m_.sum()), round(float(np.exp(np.median(lnb[m_]))), 3), round(float(np.exp(np.median(lnb[m_]))*boot_median(lnb[m_])), 3)]
# 1 Mpc aperture and 600 km/s gapper as checks
b1 = [K*r['sig1']**2/r['T'] for r in sel if r['sig1'] and r['n1'] >= 20]
res['beta_median_1Mpc'] = [len(b1), round(float(np.median(b1)), 3)]
b6 = [K*r['sig600']**2/r['T'] for r in sel if r['sig600']]
res['beta_median_gap600'] = [len(b6), round(float(np.median(b6)), 3)]
# gravity-only expectation
band = []
for c in (3, 4, 6):
    for cg in (c/2, c):
        for slope in (2.0, 2.2, 2.5):
            band.append(gb.model(c, cg, slope, 1.0)[0])
res['gravity_r200'] = dict(min=round(min(band), 2), max=round(max(band), 2), central=round(gb.model(4, 4, 2.2, 1.0)[0], 2))
# if galaxies and gas sat in unrelated wells: what share of clusters would land within the band by chance?
# Wu, Fang & Xu 1998 compilation
wu = []
for r in open('wu98_raw.txt').read().strip().split(';'):
    f = r.split(',')
    if f[0] == 'AC118': continue                     # same cluster as A2744
    T, s = float(f[2]), float(f[5]); wu.append((f[0], T, s, K*s*s/T))
wb = np.log([w[3] for w in wu]); wT = np.log([w[1] for w in wu]); ws = np.log([w[2] for w in wu])
wmed, wmad = robust(wb)
res['wu98'] = dict(n=len(wu), beta_median=round(float(np.exp(wmed)), 3), beta_median_err=round(float(np.exp(wmed)*boot_median(wb)), 3),
                   scatter_lnbeta_robust=round(float(wmad), 3),
                   slope_ols_bisector=None)
b1_ = np.polyfit(wT, ws, 1)[0]; b2_ = 1/np.polyfit(ws, wT, 1)[0]
res['wu98']['slope_ols_bisector'] = round(float((b1_*b2_-1+math.sqrt((1+b1_**2)*(1+b2_**2)))/(b1_+b2_)), 3)
res['wu98']['spearman'] = round(float(stats.spearmanr(wT, ws)[0]), 3)
# magnetic against thermal pressure in Coma (Bonafede et al. 2010; B ∝ n_e^0.5 so the ratio is ~constant)
n0, T_coma = 3.44e-3, 8.0
Pth = 1.929*n0*T_coma*1.60218e-9
res['coma'] = {f'B0={B}uG': round((B*1e-6)**2/(8*math.pi)/Pth, 4) for B in (3.9, 4.7, 5.4)}
res['coma']['B_equipartition_uG'] = round(math.sqrt(8*math.pi*Pth)*1e6, 1)
res['coma']['Pth_cgs'] = float('%.3g' % Pth)
# a list for the figure
res['points'] = [[r['name'], round(r['T'], 2), round(r['dlnT']*r['T'], 2), int(r['sig']), round(r['dlns']*r['sig']), int(r['n']), bool(o)] for r, o in zip(sel, out_mask)]
res['wu_points'] = [[w[0], w[1], int(w[2])] for w in wu]
json.dump(res, open('cluster_heat.json', 'w'), indent=1)
for k, v in res.items():
    if k not in ('points', 'wu_points'): print(k, v)
