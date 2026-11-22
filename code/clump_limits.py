"""clump_limits.py - Workbench Entry 105: primordial-field clumping against the 2026 limits on electron-density fluctuations.

Lynch, Chluba, Battye & Vasil (2026) bound the variance of the electron density during recombination, <delta_e^2>, with
Planck + ACT DR6 + SPT-3G D1 ("SPA"): < 0.21 at 95%; < 0.30 with DESI DR2. This notebook's clumping model is Jedamzik & Pogosian's
three-zone M1 (f2 = 1/3, D1 = 0.1, D2 = 1, b = <Delta^2> - 1), run through the same hydrogen recombination as Entries 13, 92 and 94.
Here: the electron-density variance that M1 produces at each clumping b, through recombination, so that the two can be compared;
and the Hubble constant that each b gives in the notebook's own track (clump_track.json, Entry 13).
"""
import json, os, math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
HERE = os.path.dirname(os.path.abspath(__file__))
wb, Om, h = 0.02237, 0.315, 0.674
T0 = 2.7255; Yp = 0.247; Neff = 3.046
wg = 2.469e-5; wnu = 0.2271*Neff*wg; Or = (wg + wnu)/h**2; OL = 1 - Om - Or
H0 = h/2997.92458
def Ha(a): return H0*np.sqrt(Or/a**4 + Om/a**3 + OL)
kB=1.380649e-23; hP=6.62607015e-34; me=9.1093837e-31; c=2.99792458e8; sigT=6.6524587e-29
Eion=2.179872e-18; E2=Eion/4; Ealpha=0.75*Eion; lam_a=121.5668e-9; L2s=8.2245809; mH=1.6735575e-27; Mpc=3.0857e22
nH0 = (1-Yp)*wb*1.87847e-26/mH
def alphaB(T):
    t=T/1e4; return 1.14*1e-19*4.309*t**-0.6166/(1+0.6703*t**0.5300)
zs = np.linspace(1700, 600, 2201)
def solve_zone(D):
    def rhs(z, x):
        T=T0*(1+z); nH=D*nH0*(1+z)**3; xe=min(max(x[0],1e-12),1.0)
        al=alphaB(T); beta=al*(2*np.pi*me*kB*T/hP**2)**1.5*np.exp(-E2/(kB*T))
        Hs=Ha(1/(1+z))*c/Mpc; La=8*np.pi*Hs/(lam_a**3*nH*max(1-xe,1e-12)); Cr=(L2s+La)/(L2s+La+beta)
        return [Cr*(nH*xe*xe*al - beta*(1-xe)*np.exp(-Ealpha/(kB*T)))/(Hs*(1+z))]
    T=T0*1701; nH=D*nH0*1701**3; S=(2*np.pi*me*kB*T/hP**2)**1.5*np.exp(-Eion/(kB*T))/nH
    x0=(-S+np.sqrt(S*S+4*S))/2
    return solve_ivp(rhs,(1700,600),[x0],t_eval=zs,method='LSODA',rtol=1e-7,atol=1e-11).y[0]
def zones_M1(b):
    D1,D2,f2=0.1,1.0,1/3
    def resid(f3):
        f1=1-f2-f3; D3=(1-f1*D1-f2*D2)/f3; return f1*D1**2+f2*D2**2+f3*D3**2-(1+b)
    f3=brentq(resid,1e-4,2/3-1e-6); f1=1-f2-f3; D3=(1-f1*D1-f2*D2)/f3
    return [(f1,D1),(f2,D2),(f3,D3)]
def visibility(xe):
    a = 1/(1+zs); ne = xe*nH0*(1+zs)**3
    dtau_dz = ne*sigT*c/(Ha(a)*c/Mpc*(1+zs))/Mpc*Mpc   # per unit z: n sigma c dt/dz
    dtau_dz = ne*sigT*c/((Ha(a)*c/Mpc)*(1+zs))
    tau = np.concatenate([[0], np.cumsum(0.5*(dtau_dz[1:]+dtau_dz[:-1])*np.abs(np.diff(zs)))])
    tau = tau[-1] - tau                                      # optical depth from z down to 600 (integrated toward lower z)
    g = dtau_dz*np.exp(-tau); return g/np.trapezoid(g, -zs)
std = solve_zone(1.0); gstd = visibility(std); zpk = float(zs[np.argmax(gstd)])
rows = []
for b in [0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5]:
    Z = zones_M1(b); X = [solve_zone(D) for _, D in Z]
    ne = sum(f*D*x for (f, D), x in zip(Z, X)); ne2 = sum(f*(D*x)**2 for (f, D), x in zip(Z, X))
    var = ne2/ne**2 - 1                                      # electron-density variance through recombination
    g = visibility(ne)
    rows.append(dict(b=b, zones=[[round(f,3), round(D,3)] for f, D in Z], var_z1500=round(float(np.interp(-1500, -zs, var)), 3),
                     var_peak=round(float(np.interp(-zpk, -zs, var)), 3), var_vis=round(float(np.trapezoid(var*g, -zs)), 3),
                     zstar=round(float(zs[np.argmax(g)]), 1)))
out = dict(z_peak_standard=zpk, rows=rows)
T = json.load(open(os.path.join(HERE, 'clump_track.json')))
out['track'] = [dict(b=t['b'], H0=t['H0'], rd=t['rd']) for t in T]
# b at which the visibility-weighted variance reaches each bound, and the H0 the track gives there
bb = np.array([r['b'] for r in rows]); vv = np.array([r['var_vis'] for r in rows]); v15 = np.array([r['var_z1500'] for r in rows])
tb = np.array([t['b'] for t in T]); th = np.array([t['H0'] for t in T])
for lab, lim in (('SPA', 0.21), ('SPA_DESI', 0.30), ('SO_forecast', 0.10)):
    b_vis = float(np.interp(lim, vv, bb)); b_early = float(np.interp(lim, v15, bb))
    out['bound_'+lab] = dict(limit=lim, b_max_visibility=round(b_vis, 3), b_max_early=round(b_early, 3),
                             H0_at_b_max=round(float(np.interp(max(b_vis, b_early), tb, th)), 2))
out['entry36_b'] = [0.25, 0.14, 0.12]
out['H0_track_at_entry36_b'] = round(float(np.interp(0.25, tb, th)), 2)
json.dump(out, open(os.path.join(HERE, 'clump_limits.json'), 'w'), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != 'track'}, indent=0)[:3000]); print(out['track'])

# ---- the same variance, distributed log-normally (the form of the 2026 limits) against M1 ----
def zones_lognormal(V, n=16):
    s2 = math.log(1 + V); s = math.sqrt(s2)
    x, wts = np.polynomial.hermite_e.hermegauss(n); wts = wts/wts.sum()
    D = np.exp(s*x - s2/2); D = D/np.sum(wts*D)                     # mean density exactly 1
    return list(zip(wts, D))
def zstar_of(Z):
    X = [solve_zone(D) for _, D in Z]; ne = sum(f*D*x for (f, D), x in zip(Z, X)); g = visibility(ne)
    return float(zs[np.argmax(g)]), ne
# finer z grid near the peak for the shift
def zpk_fine(ne):
    g = visibility(ne); i = int(np.argmax(g)); y0, y1, y2 = g[i-1], g[i], g[i+1]
    dz = zs[1] - zs[0]; return float(zs[i] + 0.5*dz*(y0 - y2)/(y0 - 2*y1 + y2))
z0 = zpk_fine(std)
ln_rows = []
for V in [0.05, 0.1, 0.21, 0.3, 0.5]:
    _, ne = zstar_of(zones_lognormal(V)); ln_rows.append(dict(V=V, dzstar=round(zpk_fine(ne) - z0, 2)))
m1_rows = []
for b in [0.05, 0.1, 0.21, 0.25, 0.3, 0.5]:
    _, ne = zstar_of(zones_M1(b)); m1_rows.append(dict(b=b, dzstar=round(zpk_fine(ne) - z0, 2)))
out['shift_lognormal'] = ln_rows; out['shift_M1'] = m1_rows
mb = np.array([r['b'] for r in m1_rows]); md = np.array([r['dzstar'] for r in m1_rows])
for r in ln_rows: r['M1_b_same_shift'] = round(float(np.interp(r['dzstar'], md, mb)), 3)
out['b_equiv_SPA'] = [r for r in ln_rows if r['V'] == 0.21][0]['M1_b_same_shift']
out['b_equiv_SPA_DESI'] = [r for r in ln_rows if r['V'] == 0.3][0]['M1_b_same_shift']
out['H0_track_at_b_equiv_SPA'] = round(float(np.interp(out['b_equiv_SPA'], tb, th)), 2)
out['H0_track_at_b_equiv_SPA_DESI'] = round(float(np.interp(out['b_equiv_SPA_DESI'], tb, th)), 2)
json.dump(out, open(os.path.join(HERE, 'clump_limits.json'), 'w'), indent=1)
print('LN', ln_rows); print('M1', m1_rows); print(out['b_equiv_SPA'], out['b_equiv_SPA_DESI'], out['H0_track_at_b_equiv_SPA'], out['H0_track_at_b_equiv_SPA_DESI'])
