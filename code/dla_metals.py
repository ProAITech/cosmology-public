"""Entry 89: does the gas get poorer in metals with distance?

Damped Lyman-alpha absorbers (DLAs) are gas clouds with N(HI) >= 2e20 cm^-2 seen against
quasars. Their metallicity [M/H] comes from ratios of column densities measured in the
absorber's own frame; no distance law enters. A universe that is the same at every distance
(a static universe with no cosmic evolution) predicts no trend of [M/H] with redshift.
Data: Rafelski et al. 2012 (VizieR J/ApJ/755/89; 242 rows, alpha-element or Fe+0.3 metallicity)
and De Cia et al. 2018 (J/A+A/611/A76; 235 dust-corrected [Fe/H]tot), largely the same absorbers.
"""
import json, math, os, numpy as np
from scipy import optimize, stats
os.chdir(os.path.dirname(os.path.abspath(__file__)))
rng = np.random.default_rng(89)
def load_r():
    rows = []; seen = set()
    for r in open('rafelski12.txt').read().split(';'):
        z, nhi, f, m, e = r.split(',')
        key = (z, nhi, m)
        if key in seen: continue                     # exact duplicate rows in the literature table
        seen.add(key); rows.append((float(z), float(nhi), int(f), float(m), max(float(e), 0.05)))
    return np.array(rows)
def load_d():
    return np.array([[float(x) for x in r.split(',')] for r in open('decia18.txt').read().split(';')])
R = load_r(); D = load_d()
def ml_line(z, y, e):
    def nll(p):
        a, b, s = p; v = e**2+s**2; return 0.5*np.sum((y-a-b*z)**2/v+np.log(v))
    p = optimize.minimize(nll, [-0.8, -0.2, 0.4], method='Nelder-Mead', options=dict(xatol=1e-7, fatol=1e-9, maxiter=5000)).x
    bs = []
    for _ in range(1000):
        i = rng.integers(0, len(z), len(z))
        def n2(q):
            a, b, s = q; v = e[i]**2+s**2; return 0.5*np.sum((y[i]-a-b*z[i])**2/v+np.log(v))
        bs.append(optimize.minimize(n2, p, method='Nelder-Mead').x)
    bs = np.array(bs)
    return dict(a=round(float(p[0]), 3), b=round(float(p[1]), 3), b_err=round(float(np.std(bs[:, 1])), 3), scatter=round(float(abs(p[2])), 3))
def weighted_bins(z, nhi, m, edges):
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (z >= lo) & (z < hi)
        if s.sum() < 5: continue
        def wm(idx):
            w = 10**nhi[idx]; return math.log10(np.sum(w*10**m[idx])/np.sum(w))
        idx = np.where(s)[0]; val = wm(idx)
        bs = [wm(rng.choice(idx, len(idx))) for _ in range(2000)]
        out.append(dict(zlo=lo, zhi=hi, n=int(s.sum()), zmed=round(float(np.median(z[s])), 3), mean=round(val, 3),
                        err=round(float(np.std(bs)), 3), median=round(float(np.median(m[s])), 3)))
    return out
res = {}
z, nhi, f, m, e = R.T
res['rafelski'] = dict(n=len(R), fit=ml_line(z, m, e), fit_measured_only=ml_line(z[f <= 2], m[f <= 2], e[f <= 2]),
                       bins=weighted_bins(z, nhi, m, [0, 1.5, 2.5, 3.0, 3.5, 4.0, 4.5, 5.3]),
                       spearman=[round(float(x), 4) for x in stats.spearmanr(z, m)])
zd, md, ed = D.T
res['decia'] = dict(n=len(D), fit=ml_line(zd, md, ed), spearman=[round(float(x), 4) for x in stats.spearmanr(zd, md)])
# medians nearby vs far
for lab, arr_z, arr_m in (('rafelski', z, m), ('decia', zd, md)):
    near = arr_m[arr_z < 1.5]; far = arr_m[arr_z >= 4.0]
    res[lab]['median_z_below_1p5'] = [int(len(near)), round(float(np.median(near)), 2)]
    res[lab]['median_z_above_4'] = [int(len(far)), round(float(np.median(far)), 2)]
    res[lab]['mannwhitney_p'] = float('%.2g' % stats.mannwhitneyu(near, far).pvalue)
# the most metal-poor star-forming galaxy known nearby: J0811+4730, 12+log(O/H) = 6.98 (Izotov et al. 2018); solar 8.69
local_floor = 6.98 - 8.69
res['local_floor_MH'] = round(local_floor, 2)
res['frac_far_below_local_floor'] = round(float(np.mean(m[z >= 4.0] < local_floor)), 3)
res['frac_near_below_local_floor'] = round(float(np.mean(m[z < 1.5] < local_floor)), 3)
# distances on the static laws for the redshift range
DH = 4.283
res['distances_Gpc'] = {str(zz): dict(linear=round(DH*zz, 1), tired=round(DH*math.log(1+zz), 1)) for zz in (0.5, 1.5, 4.0, 5.0)}
res['points'] = [[round(float(a), 3), round(float(b), 2), round(float(c), 2), int(d)] for a, b, c, d in zip(z, m, e, f)]
json.dump(res, open('dla_metals.json', 'w'), indent=1)
for k, v in res.items():
    if k != 'points': print(k, v)
