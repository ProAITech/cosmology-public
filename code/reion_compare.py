"""Entry 94: the second plasma era. Large-angle E-mode polarization from reionization, computed with
the Entry 92 plasma code (cmbreion.py) for tau = 0, 0.054 and 0.09, against the measured optical
depths and the end of reionization seen in the Lyman-alpha forest."""
import json, os, math, numpy as np
os.chdir(os.path.dirname(os.path.abspath(__file__)))
runs = {k: json.load(open(f'cmbreion_{k}.json')) for k in ('tau0', 'tau054', 'tau09')}
ells = np.array(runs['tau054']['ells'])
out = {'z_re': {k: (round(v['z_re'], 2) if v['z_re'] else None) for k, v in runs.items()}}
# scale: E in Ma-Bertschinger normalisation -> divide EE by 16, TE by 4; physical units via TT first peak = 5733 uK^2
def spec(k):
    r = runs[k]['runs']['0.965']; TT = np.array(r['TT']); EE = np.array(r['EE'])/16; TE = np.array(r['TE'])/4
    return TT, EE, TE
TT0, EE0, _ = spec('tau0')
sc = 5733.0/max(TT0[(ells > 150) & (ells < 300)])                 # normalise all three runs on the tau = 0 first peak
low = ells <= 30
res = {}
for k in runs:
    TT, EE, TE = spec(k)
    e = EE*sc; t = TT*sc
    m = (ells >= 2) & (ells <= 10)
    res[k] = dict(EE_l2_10_mean=round(float(np.mean(e[m])), 4), EE_peak_low=round(float(max(e[low])), 4), l_peak_low=int(ells[low][np.argmax(e[low])]),
                  EE_l20=round(float(np.interp(20, ells, e)), 5), TT_peak1=round(float(max(t[(ells > 150) & (ells < 300)])), 1),
                  EE_low=[round(float(x), 5) for x in e[low]])
out['spectra'] = res
out['ells_low'] = [int(x) for x in ells[low]]
out['bump_ratio_09_to_054'] = round(res['tau09']['EE_l2_10_mean']/res['tau054']['EE_l2_10_mean'], 2)
out['tau_ratio_sq'] = round((0.09/0.054)**2, 2)
out['bump_over_no_reion'] = round(res['tau054']['EE_l2_10_mean']/max(res['tau0']['EE_l2_10_mean'], 1e-9), 1)
out['TT_peak_suppression_054'] = round(res['tau054']['TT_peak1']/res['tau0']['TT_peak1'], 4)
out['exp_minus_2tau'] = round(math.exp(-2*0.054), 4)
# ionization history for the figure (tau = 0.054 model), hydrogen + singly ionized helium, relative to hydrogen
Yp = 0.247; fHe = Yp/(4*(1-Yp)); zre = runs['tau054']['z_re']
def xe(z, zre):
    y = (1+z)**1.5; yre = (1+zre)**1.5; dy = 1.5*math.sqrt(1+zre)*0.5
    return (1+fHe)/2*(1+math.tanh((yre-y)/dy)) + fHe/2*(1+math.tanh((3.5-z)/0.5))
zs = np.linspace(0, 16, 161)
out['xe_curve'] = [[round(float(z), 2), round(xe(z, zre)/(1+2*fHe), 4)] for z in zs]          # as a fraction of full ionization
zre9 = runs['tau09']['z_re']
out['xe_curve_09'] = [[round(float(z), 2), round(xe(z, zre9)/(1+2*fHe), 4)] for z in zs]
# where the model is 95% ionized (hydrogen), i.e. the "end" in the Lyman-alpha-forest sense
zz = np.linspace(3, 15, 12001); xh = np.array([(1+math.tanh((((1+zre)**1.5)-(1+z)**1.5)/(1.5*math.sqrt(1+zre)*0.5)))/2 for z in zz])
out['z_95'] = round(float(zz[np.argmin(abs(xh-0.95))]), 2)
out['z_50'] = round(float(zz[np.argmin(abs(xh-0.5))]), 2)
# ---- optical depth of simple histories, in the cosmology of the runs (the non-CMB values of Entry 92) ----
from scipy.integrate import quad
h_, wb_, Om_ = runs['tau054']['h'], runs['tau054']['wb'], runs['tau054']['Om']
nH = (1-Yp)*wb_*1.87847e-29/1.6726e-24; sigT = 6.6524e-25; Mpc = 3.0857e24; c = 2.99792458e10; H0 = 100*h_*1e5/Mpc
Ez = lambda z: math.sqrt(Om_*(1+z)**3 + 1 - Om_)
def tau_of(xe, zmax=40): return quad(lambda z: xe(z)*nH*(1+z)**2/Ez(z)*sigT*c/H0, 0, zmax, limit=400, points=[3.5, 5.5, 6, 7.65, 11])[0]
he2 = lambda z: fHe if z < 3.5 else 0.0
def ramp(zmid, zend):                                    # ionized hydrogen fraction rises linearly from 0 to 1, complete at zend
    return lambda z: (1.0 if z <= zend else max(0.0, 1 - 0.5*(z - zend)/(zmid - zend)))*(1 + fHe) + he2(z)
H = {}
H['tanh_zre'] = round(tau_of(lambda z: xe(z, zre)), 4)
H['instant_6'] = round(tau_of(lambda z: (1 + fHe if z < 6 else 0) + he2(z)), 4)
H['instant_5.5'] = round(tau_of(lambda z: (1 + fHe if z < 5.5 else 0) + he2(z)), 4)
for zm in (7.0, 7.5, 8.0): H['ramp_end5.5_mid%.1f' % zm] = round(tau_of(ramp(zm, 5.5)), 4)
def mid_for(t, zend=5.5):                               # midpoint a linear ramp ending at zend needs to give tau = t
    lo, hi = zend + 0.05, 25.0
    for _ in range(60):
        m = (lo + hi)/2
        if tau_of(ramp(m, zend)) < t: lo = m
        else: hi = m
    return round((lo + hi)/2, 2)
H['ramp_end5.5_mid_for_0.054'] = mid_for(0.054); H['ramp_end5.5_mid_for_0.090'] = mid_for(0.090)
H['ramp_start_for_0.090'] = round(2*H['ramp_end5.5_mid_for_0.090'] - 5.5, 1)
ne0 = nH*(1 + 2*fHe)                                    # every electron free, today's density
H['static_ne0_cm3'] = float('%.3e' % ne0); H['static_tau_per_Gpc'] = float('%.3e' % (ne0*sigT*Mpc*1000))
H['static_path_Gpc_for_0.054'] = round(0.054/(ne0*sigT*Mpc*1000), 0)
out['histories'] = H
json.dump(out, open('reion_compare.json', 'w'), indent=1)
for k, v in out.items():
    if k not in ('xe_curve', 'xe_curve_09', 'spectra'): print(k, v)
for k, v in res.items(): print(k, {kk: vv for kk, vv in v.items() if kk != 'EE_low'})
