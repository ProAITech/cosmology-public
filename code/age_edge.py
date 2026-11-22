"""age_edge.py - Workbench Entry 140
Light older than the oldest stars? A static universe's light-travel paradox.
1. Where do local stellar ages stop? 133,970 subgiants with ages better than 8% (Xiang & Rix 2022, ages on a 0-20 Gyr isochrone grid;
   histograms by error class, built in a browser from VizieR J/other/Nat/603.599, saved in subgiant_hist.json).
   Model for the true age density near the old end: psi(t) = A_c (t/12)^beta for t <= T_max, zero beyond, plus a flat floor of
   outliers; observed ages are psi convolved with a Gaussian of width r_c * t for each error class c (r_c = class midpoint).
   Poisson likelihood on 0.25 Gyr bins from 9 to 19.5 Gyr (the grid ends at 20). Compared with the same model without an edge
   below 20 Gyr (what an unchanging universe with a long, steady history of star formation gives: beta ~ +0.4 for subgiants).
2. Light-travel time against redshift: Lerner's static law d = cz/H0, tired light d = (c/H0) ln(1+z), and the expanding LCDM
   lookback time (Planck 2018). Counts of DR16Q quasars (VizieR VII/289, read in the browser) beyond each static law's limit.
"""
import json, math, os
import numpy as np
from scipy.optimize import minimize
from scipy.integrate import quad
HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, 'subgiant_hist.json')))
edges = np.arange(0, 20.0001, 0.25); cen = 0.5*(edges[1:] + edges[:-1])
rc = [0.02, 0.04, 0.0575, 0.0725]
win = (cen > 9) & (cen < 19.5)
tt = np.linspace(4, 20, 1601); dt = tt[1] - tt[0]

def conv_matrix(r):
    # probability that true age tt lands in each observed bin
    s = r*tt
    from scipy.special import erf
    lo = (edges[:-1][:, None] - tt[None, :])/(np.sqrt(2)*s[None, :]); hi = (edges[1:][:, None] - tt[None, :])/(np.sqrt(2)*s[None, :])
    return 0.5*(erf(hi) - erf(lo))                                  # (bins, tt)
K = [conv_matrix(r) for r in rc]

def model_counts(params, H, edge=True):
    Tmax, beta, lf = params[0], params[1], params[2]
    psi = np.exp(np.clip(beta, -30, 30)*np.log(tt/12.0))
    if edge: psi = psi*(tt <= Tmax)
    out = []; nll = 0.0
    for c in range(4):
        shape = (K[c] @ psi)*dt
        sw = shape[win]; obs = H[c][win]
        floor = math.exp(lf)*np.ones_like(sw)
        # amplitude per class by maximum likelihood (analytic for the shape part with floor fixed -> solve numerically)
        a = max(obs.sum() - floor.sum(), 1.0)/max(sw.sum(), 1e-12)
        for _ in range(30):                                          # Poisson ML for amplitude a
            mu = a*sw + floor; g = (obs*sw/mu).sum() - sw.sum(); h_ = -(obs*sw*sw/mu**2).sum()
            a = max(a - g/h_, 1e-9)
        mu = a*sw + floor
        nll += (mu - obs*np.log(mu)).sum()
        out.append((a, mu))
    return nll, out

def fit(H, edge=True, beta_fixed=None):
    best = None
    for T0 in (12.5, 13.0, 13.5, 14.0, 15.0):
        for b0 in (-2.0, 0.0, 2.0):
            def f(p):
                q = [p[0] if edge else 25.0, beta_fixed if beta_fixed is not None else p[1], p[2]]
                return model_counts(q, H, edge)[0]
            res = minimize(f, [T0, b0, 0.0], method='Nelder-Mead', options=dict(xatol=1e-3, fatol=1e-4, maxiter=4000))
            if best is None or res.fun < best.fun: best = res
    p = best.x.copy()
    if not edge: p[0] = 25.0
    if beta_fixed is not None: p[1] = beta_fixed
    return best.fun, p

def profile_T(H, Ts):
    out = []
    for T in Ts:
        def f(p): return model_counts([T, p[0], p[1]], H, True)[0]
        best = min((minimize(f, [b0, 0.0], method='Nelder-Mead', options=dict(xatol=1e-3, fatol=1e-4)) for b0 in (-2.0, 0.0, 2.0)), key=lambda r: r.fun)
        out.append(best.fun)
    return np.array(out)

R = {}
for name in ('all', 'halo'):
    H = [np.array(h, float) for h in D[name]]
    nll_e, p_e = fit(H, True); nll_n, p_n = fit(H, False); nll_s, p_s = fit(H, False, beta_fixed=0.43)
    Ts = np.arange(12.0, 15.51, 0.05); prof = profile_T(H, Ts); ok = Ts[prof <= prof.min() + 2.0]   # 95% for one parameter
    # observed and predicted counts above 14 Gyr
    obs14 = sum(H[c][(cen > 14) & win].sum() for c in range(4))
    _, out_e = model_counts(p_e, H, True); _, out_s = model_counts(p_s, H, False)
    w14 = (cen[win] > 14)
    pred_e = sum(mu[w14].sum() for a, mu in out_e); pred_s = sum(mu[w14].sum() for a, mu in out_s)
    floor_frac = math.exp(p_e[2])*win.sum()*4/sum(H[c][win].sum() for c in range(4))
    R[name] = dict(n=int(sum(h.sum() for h in H)), Tmax=round(float(p_e[0]), 2), Tmax_95=[round(float(ok.min()), 2), round(float(ok.max()), 2)],
                   beta_edge=round(float(p_e[1]), 2), floor_share_window=round(floor_frac, 4),
                   dlnL_noedge_vs_edge=round(float(nll_n - nll_e), 1), beta_noedge=round(float(p_n[1]), 2),
                   dlnL_steady_vs_edge=round(float(nll_s - nll_e), 1),
                   obs_above14=int(obs14), pred_above14_edge=round(float(pred_e), 0), pred_above14_steady=round(float(pred_s), 0))
    # observed vs model histograms for the figure (summed over classes, window only)
    R[name]['fig'] = dict(t=[round(float(x), 3) for x in cen[win]], obs=[int(sum(H[c][win][i] for c in range(4))) for i in range(win.sum())],
                         edge=[round(float(sum(mu[i] for a, mu in out_e)), 1) for i in range(win.sum())],
                         steady=[round(float(sum(mu[i] for a, mu in out_s)), 1) for i in range(win.sum())])
    # model-free: counts in 9-14 and 14-19.5 Gyr, against a steady history (dN/dt ~ t^0.43 for subgiants, or flat)
    tot = sum(H[c] for c in range(4)); n914 = tot[(cen > 9) & (cen < 14)].sum(); n1419 = tot[(cen > 14) & (cen < 19.5)].sum()
    R[name].update(n_9_14=int(n914), n_14_19p5=int(n1419), ratio=round(float(n1419/n914), 4),
                   ratio_steady=round((19.5**1.43 - 14**1.43)/(14**1.43 - 9**1.43), 3), ratio_flat=round(5.5/5.0, 3),
                   per_bin_10_to_11=round(float(tot[(cen > 10) & (cen < 11)].mean()), 0), per_bin_16_to_19=round(float(tot[(cen > 16) & (cen < 19)].mean()), 0),
                   peak_bin=float(cen[(cen > 5)][np.argmax(tot[cen > 5])]))
    print(name, json.dumps({k: v for k, v in R[name].items() if k != 'fig'}))

# --- light-travel times ---
Gyr_per_H0 = lambda H0: 977.792/H0
h, Om = 0.674, 0.315; Or = 4.18e-5/h**2; OL = 1 - Om - Or
def Ea(a): return math.sqrt(Om/a**3 + Or/a**4 + OL)
def lookback(z): return quad(lambda a: 1/(a*Ea(a)), 1/(1 + z), 1)[0]*Gyr_per_H0(100*h)
age0 = quad(lambda a: 1/(a*Ea(a)), 0, 1)[0]*Gyr_per_H0(100*h)
T = R['all']['Tmax']; Thi = R['all']['Tmax_95'][1]
rows = []
for H0 in (67.4, 70.0, 73.0):
    tH = Gyr_per_H0(H0)
    rows.append(dict(H0=H0, z_lerner=round(T/tH, 3), z_tired=round(math.exp(T/tH) - 1, 3), z_lerner_16=round(16/tH, 3), z_tired_16=round(math.exp(16/tH) - 1, 3)))
objs = [('MoM-z14', 14.44), ('JADES-GS-z14-0', 14.32), ('quasar J0313-1806', 7.64), ('DR16Q median', None)]
tH70 = Gyr_per_H0(70.0)
lt = []
for name, z in objs[:3]:
    lt.append(dict(name=name, z=z, lerner_Gyr=round(z*tH70, 1), tired_Gyr=round(math.log(1 + z)*tH70, 1), lcdm_lookback_Gyr=round(lookback(z), 2),
                   lcdm_age_at_emission_Gyr=round(age0 - lookback(z), 2)))
curve = dict(z=[round(float(z), 3) for z in np.concatenate([np.linspace(0, 2, 41), np.linspace(2.1, 15, 130)])])
curve['lerner'] = [round(z*tH70, 2) for z in curve['z']]; curve['tired'] = [round(math.log(1 + z)*tH70, 2) for z in curve['z']]
curve['lcdm'] = [round(lookback(z), 3) for z in curve['z']]
R['light'] = dict(age_universe_LCDM=round(age0, 3), thresholds=rows, objects=lt, dr16q=D['dr16q'], curve=curve)
R['light']['lcdm_lookback_max_Gyr'] = round(age0, 3)
json.dump(R, open(os.path.join(HERE, 'age_edge.json'), 'w'), indent=1)
print(json.dumps(R['light']['thresholds'])); print(json.dumps(lt)); print('age0', age0)
