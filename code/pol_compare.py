"""Entry 92: is the CMB polarized the way sound in the primordial plasma requires?

Measured: positions and heights of the peaks and troughs of the TT, TE and EE spectra
(Planck 2018 results I, Table 5). Predicted: our compact Boltzmann code (Entry 12), now with
the E-mode source of Zaldarriaga & Seljak (1997), fed only non-CMB inputs (deuterium, DESI BAO
with BBN, atomic physics) and tilt 0.965. The code's polarization source uses Pi in the
Ma & Bertschinger normalisation, F_l = 4 Theta_l, so E is divided by 4 here.
Two tests: (1) model-free, from the data alone: are the EE peaks spaced like the TT peaks and
shifted by half a period, as the velocity of an oscillating fluid is? (2) does the plasma code,
with no CMB anisotropy input, put the EE and TE extrema where Planck finds them?
"""
import json, os, numpy as np
os.chdir(os.path.dirname(os.path.abspath(__file__)))
P = dict(
 TT_peaks=[(220.6,0.6),(538.1,1.3),(809.8,1.0),(1147.8,2.3),(1446.8,1.6),(1779,3)],
 TT_dips=[(416.3,1.1),(675.5,1.2),(1001.1,1.8),(1290.0,1.8),(1623.8,2.1)],
 TE_min=[(146.0,1.1),(471.1,0.7),(747.2,0.8),(1072.5,1.2),(1372.7,2.8),(1697.4,5.9)],
 TE_max=[(308.2,0.8),(595.8,0.9),(917.1,0.8),(1221.3,1.7),(1532.1,2.2)],
 EE_peaks=[(145,3),(398.3,1.0),(690.4,1.2),(993.1,1.8),(1296.4,4.3)],
 EE_dips=[(195.0,5.4),(522.0,1.1),(831.8,1.1),(1158.8,2.6)],
 EE_heights=[1.11,21.45,38.1,42.6,31.5], TT_heights=[5733,2586,2518,1227,799],
 TE_min_h=[-46.7,-74.1,-128.0,-79.1,-60.0,-27.2], TE_max_h=[117.1,27.8,59.0,3.5,8.9])
res = {'planck': P}
# (1) model-free phase test
tt = [p[0] for p in P['TT_peaks']]; ee = [p[0] for p in P['EE_peaks']]
phase = []
for e in ee[1:]:
    i = np.searchsorted(tt, e); lo, hi = tt[i-1], tt[i]
    phase.append(round((e-lo)/(hi-lo), 3))
res['EE_phase_between_TT_peaks'] = phase
res['TT_spacing'] = [round(b-a, 1) for a, b in zip(tt, tt[1:])]
res['EE_spacing'] = [round(b-a, 1) for a, b in zip(ee[1:], ee[2:])]
res['EE_vs_TT_dips'] = [(e, d[0], round(e/d[0]-1, 3)) for e, d in zip(ee[1:], P['TT_dips'])]
# (2) the plasma code with non-CMB inputs
M = json.load(open('cmbpol_noncmb.json')); run = M['runs']['0.965']
ells = np.array(M['ells']); TT = np.array(run['TT']); EE = np.array(run['EE'])/16; TE = np.array(run['TE'])/4
def ext(y, lmin):
    out = []
    for i in range(1, len(ells)-1):
        if ells[i] < lmin: continue
        for kind, cond in (('max', y[i] > y[i-1] and y[i] >= y[i+1]), ('min', y[i] < y[i-1] and y[i] <= y[i+1])):
            if cond:
                A = np.polyfit(ells[i-1:i+2], y[i-1:i+2], 2); lp = -A[1]/(2*A[0])
                out.append((kind, round(float(lp), 1), float(np.polyval(A, lp))))
    return out
pred = {'EE': ext(EE, 100), 'TE': ext(TE, 100), 'TT': ext(TT, 150)}
def match(meas, kind, series):
    rows = []
    cand = [p for p in pred[series] if p[0] == kind]
    for l, e in meas:
        best = min(cand, key=lambda p: abs(p[1]-l))
        rows.append(dict(planck=l, err=e, model=best[1], diff=round(best[1]-l, 1), pct=round(100*(best[1]/l-1), 2)))
    return rows
TEsign = 1 if np.interp(146, ells, TE) < 0 else -1          # fix the overall TE sign convention to Planck's (negative at l ~ 150)
if TEsign < 0: TE = -TE; pred['TE'] = [('max' if k == 'min' else 'min', l, -v) for k, l, v in pred['TE']]
res['TE_sign_flipped_by_convention'] = TEsign < 0
cmp = {'EE peaks': match(P['EE_peaks'], 'max', 'EE'), 'EE troughs': match(P['EE_dips'], 'min', 'EE'),
       'TE troughs': match(P['TE_min'], 'min', 'TE'), 'TE peaks': match(P['TE_max'], 'max', 'TE'),
       'TT peaks': match(P['TT_peaks'], 'max', 'TT'), 'TT troughs': match(P['TT_dips'], 'min', 'TT')}
res['model_vs_planck'] = cmp
allpct = [abs(r['pct']) for k, v in cmp.items() for r in v if r['planck'] > 300]
res['median_abs_pct_beyond_300'] = round(float(np.median(allpct)), 2); res['max_abs_pct_beyond_300'] = round(float(max(allpct)), 2)
# EE heights relative to EE peak 4 (normalisation-free) and the EE/TT ratio at the 4th EE peak vs the 1st TT peak
ee_model_peaks = [p for p in pred['EE'] if p[0] == 'max']
def hnear(l): return min(ee_model_peaks, key=lambda p: abs(p[1]-l))[2]
h4m = hnear(993.1); h4p = P['EE_heights'][3]
res['EE_rel_heights'] = [dict(l=l, planck=round(h/h4p, 3), model=round(hnear(l)/h4m, 3)) for (l, _), h in zip(P['EE_peaks'], P['EE_heights'])]
tt1m = max(TT[(ells > 150) & (ells < 300)])
res['EE4_over_TT1'] = dict(planck=round(h4p/P['TT_heights'][0], 5), model=round(h4m/tt1m, 5))
te_min_m = [p for p in pred['TE'] if p[0] == 'min']
res['TE_trough3_over_TT1'] = dict(planck=round(P['TE_min_h'][2]/P['TT_heights'][0], 4), model=round(min(te_min_m, key=lambda p: abs(p[1]-747))[2]/tt1m, 4))
# curves for the figure, scaled to Planck's first TT peak (5733 uK^2)
sc = 5733.0/tt1m
res['fig'] = dict(ells=ells.tolist(), EE=(EE*sc).round(3).tolist(), TE=(TE*sc).round(2).tolist())
json.dump(res, open('pol_compare.json', 'w'), indent=1)
for k, v in res.items():
    if k not in ('fig', 'planck'): print(k, v)
