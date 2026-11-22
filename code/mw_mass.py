"""mw_mass.py - Workbench Entry 103: is the Milky Way light? Gaia's falling rotation curve against the halo's tracers.

Data
  ou24.txt: Ou, Eilers, Necib & Frebel (2024) Table 1, circular speed from Gaia DR3 + APOGEE red giants, 6.3-27.3 kpc
  Halo tracers (enclosed mass, all modelling the LMC):
    Brooks et al. 2026, outer halo stars (H3+SEGUE+MagE): M(<50 kpc) = 3.36 +- 0.15 e11
    Koposov et al. 2023, Orphan-Chenab stream:            M(<50 kpc) = 3.77 +0.26 -0.19 e11
    Vasiliev et al. 2021, Sagittarius stream:              M(<100 kpc) = 5.6 +- 0.4 e11
    Eadie & Juric 2019, globular clusters:                 M(<100 kpc) = 5.3 +2.1 -1.2 e11
  Baryons: McMillan 2017 main model (means): bulge 9.13e9; thin disc Sigma0 886.7 Msun/pc^2, Rd 2.53 kpc;
    thick disc 156.7, 3.38; HI Sigma0 53.1, Rd 7, hole 4 kpc; H2 Sigma0 2180, Rd 1.5, hole 12 kpc.
    Hot halo gas inside 50 kpc: 3.8e9 (Miller & Bregman 2015).
Models: Newtonian baryons + NFW halo; MOND (radial-acceleration-relation form) with and without an external field (1D approximation).
"""
import json, os, math
import numpy as np
from scipy.special import ellipk, i0, i1, k0, k1
from scipy import optimize
HERE = os.path.dirname(os.path.abspath(__file__))
G = 4.30091e-6                       # kpc (km/s)^2 / Msun
A0 = 1.2e-10*3.0857e19/1e6           # a0 in (km/s)^2/kpc
RHOC = 277.5*0.674**2                # Msun/kpc^3
D = np.loadtxt(os.path.join(HERE, 'ou24.txt'))
R, V, SP, SM, N = D.T

# ---- baryons: thin rings with a small softening for disc thickness ----
DISCS = {'thin': (886.7e6, 2.53, 0.0), 'thick': (156.7e6, 3.38, 0.0), 'HI': (53.1e6, 7.0, 4.0), 'H2': (2180e6, 1.5, 12.0)}
def sigma(a, S0, Rd, Rm): return S0*np.exp(-(Rm/np.maximum(a, 1e-6) if Rm else 0.0) - a/Rd)
ra = np.concatenate([np.arange(0.005, 60, 0.01), np.arange(60.05, 250, 0.1)]); da = np.concatenate([np.full(6000, 0.01), np.full(len(ra)-6000, 0.1)])
H = 0.15
def phi(Rr, m):
    s = (ra[None, :] + Rr[:, None])**2 + H**2
    return -(2*G/np.pi)*np.sum(m[None, :]*ellipk(4*ra[None, :]*Rr[:, None]/s)/np.sqrt(s), axis=1)
def gdisc(Rr, key):
    S0, Rd, Rm = DISCS[key]; m = 2*np.pi*ra*sigma(ra, S0, Rd, Rm)*da; e = 0.005
    return (phi(Rr + e, m) - phi(Rr - e, m))/(2*e)                       # inward acceleration, (km/s)^2/kpc
MB = 9.13e9
RG = np.unique(np.concatenate([R, np.linspace(5, 110, 421)]))
GD = {k: gdisc(RG, k) for k in DISCS}
MASS = {k: float(np.sum(2*np.pi*ra*sigma(ra, *DISCS[k])*da)) for k in DISCS}
# check against Freeman's exact thin exponential disc
y = 10/(2*2.53); vf2 = 4*np.pi*G*886.7e6*2.53*y**2*(i0(y)*k0(y) - i1(y)*k1(y)); vr2 = 10*np.interp(10, RG, GD['thin']) if 10 in RG else None
gt = gdisc(np.array([10.0]), 'thin')[0]*10
check = dict(freeman_v_thin_10kpc=round(math.sqrt(vf2), 2), ring_v_thin_10kpc=round(math.sqrt(gt), 2))
def gbar(Rr, ups=1.0):
    gs = ups*(np.interp(Rr, RG, GD['thin']) + np.interp(Rr, RG, GD['thick']) + G*MB/Rr**2)
    return gs + np.interp(Rr, RG, GD['HI']) + np.interp(Rr, RG, GD['H2'])
def vnfw2(Rr, M200, c):
    r200 = (3*M200/(4*np.pi*200*RHOC))**(1/3); rs = r200/c; x = Rr/rs
    return G*M200*(np.log(1 + x) - x/(1 + x))/(Rr*(np.log(1 + c) - c/(1 + c)))
def menc_nfw(r, M200, c): return vnfw2(r, M200, c)*r/G

# ---- tracers as circular speeds ----
TR = [dict(name='Outer halo stars (Brooks et al. 2026)', r=50, M=3.36e11, ep=0.15e11, em=0.15e11),
      dict(name='Orphan-Chenab stream (Koposov et al. 2023)', r=50, M=3.77e11, ep=0.26e11, em=0.19e11),
      dict(name='Sagittarius stream (Vasiliev et al. 2021)', r=100, M=5.6e11, ep=0.4e11, em=0.4e11),
      dict(name='Globular clusters (Eadie & Juric 2019)', r=100, M=5.3e11, ep=2.1e11, em=1.2e11)]
for t in TR:
    t['V'] = math.sqrt(G*t['M']/t['r']); t['Vep'] = math.sqrt(G*(t['M'] + t['ep'])/t['r']) - t['V']; t['Vem'] = t['V'] - math.sqrt(G*(t['M'] - t['em'])/t['r'])
BASE = [TR[0], TR[2]]                                                      # one tracer per radius: the two most precise

def chi_rc(vm, floor=3.0, syst_out=0.0):
    s = np.where(vm > V, SP, SM); s = np.sqrt(s**2 + floor**2 + (syst_out*V*(R > 22))**2)
    return float(np.sum(((vm - V)/s)**2))
def chi_tr(f, trs):
    out = 0.0
    for t in trs:
        v = f(t['r']); s = t['Vep'] if v > t['V'] else t['Vem']; out += ((v - t['V'])/s)**2
    return out
out = {'check': check, 'disc_masses': {k: float('%.3g' % v) for k, v in MASS.items()}, 'bulge': MB}
out['baryon_total'] = float('%.3g' % (MB + sum(MASS.values())))

# ---- 1. the spherical-equivalent mass the curve itself implies ----
def msph(v, r): return v**2*r/G
k = {}
for Rr in (21.72, 25.02, 27.31):
    i = int(np.argmin(abs(R - Rr))); k[str(Rr)] = dict(M=float('%.3g' % msph(V[i], R[i])), M_lo=float('%.3g' % msph(V[i] - SM[i], R[i])), M_hi=float('%.3g' % msph(V[i] + SP[i], R[i])), nstar=int(N[i]))
out['spherical_mass_from_curve'] = k
Mk = msph(V[-1], R[-1])
out['kepler_from_last_point'] = dict(V50=round(V[-1]*math.sqrt(R[-1]/50), 1), V100=round(V[-1]*math.sqrt(R[-1]/100), 1), M=float('%.3g' % Mk),
                                     shortfall_sigma_Brooks=round((TR[0]['M'] - Mk)/TR[0]['em'], 1))
out['tracers'] = [dict(name=t['name'], r=t['r'], M=t['M'], V=round(t['V'], 1), Vep=round(t['Vep'], 1), Vem=round(t['Vem'], 1)) for t in TR]
out['stars_beyond_22kpc'] = int(N[R > 22].sum())

# ---- 2. Newtonian baryons + NFW ----
def fit_nfw(trs, floor=3.0, syst_out=0.0, ups_free=False):
    def f(p):
        lm, c = p[0], p[1]; ups = p[2] if ups_free else 1.0
        if not (10 < lm < 13.5 and 1 < c < 60 and 0.5 < ups < 1.5): return 1e9
        vm = np.sqrt(gbar(R, ups)*R + vnfw2(R, 10**lm, c))
        return chi_rc(vm, floor, syst_out) + chi_tr(lambda r: math.sqrt(gbar(np.array([float(r)]), ups)[0]*r + vnfw2(r, 10**lm, c)), trs)
    best = None
    for lm0 in (11.3, 11.8, 12.2):
        for c0 in (5, 12, 25):
            r = optimize.minimize(f, [lm0, c0] + ([1.0] if ups_free else []), method='Nelder-Mead', options=dict(xatol=1e-6, fatol=1e-6, maxiter=8000))
            if best is None or r.fun < best.fun: best = r
    lm, c = best.x[0], best.x[1]; ups = best.x[2] if ups_free else 1.0
    vm = np.sqrt(gbar(R, ups)*R + vnfw2(R, 10**lm, c))
    res = dict(M200=float('%.3g' % 10**lm), c=round(c, 1), ups=round(ups, 3), chi2=round(best.fun, 1), chi2_rc=round(chi_rc(vm, floor, syst_out), 1), n_rc=len(R),
               V27=round(float(vm[-1]), 1), V20=round(float(np.interp(20.22, R, vm)), 1))
    for r in (50, 100):
        v = math.sqrt(gbar(np.array([float(r)]), ups)[0]*r + vnfw2(r, 10**lm, c)); res[f'V{r}'] = round(v, 1); res[f'M{r}'] = float('%.3g' % (v*v*r/G))
    res['tracer_sigma'] = {t['name'].split(' (')[0]: round((res['V%d' % t['r']] - t['V'])/(t['Vep'] if res['V%d' % t['r']] > t['V'] else t['Vem']), 1) for t in TR}
    return res
out['nfw'] = {'curve_only': fit_nfw([]), 'curve_plus_tracers': fit_nfw(BASE),
              'curve_only_syst15': fit_nfw([], syst_out=0.15), 'curve_plus_tracers_syst15': fit_nfw(BASE, syst_out=0.15),
              'curve_plus_tracers_stat_only': fit_nfw(BASE, floor=0.0), 'curve_only_stat_only': fit_nfw([], floor=0.0)}
a, b = out['nfw']['curve_only'], out['nfw']['curve_plus_tracers']
out['nfw']['dchi2_curve_from_tracers'] = round(b['chi2_rc'] - a['chi2_rc'], 1)
a, b = out['nfw']['curve_only_syst15'], out['nfw']['curve_plus_tracers_syst15']
out['nfw']['dchi2_curve_from_tracers_syst15'] = round(b['chi2_rc'] - a['chi2_rc'], 1)

# ---- 3. MOND, no dark matter ----
nu = lambda y: 1/(1 - np.exp(-np.sqrt(y)))
def gmond(gn, e=0.0):
    if e == 0: return nu(gn/A0)*gn
    ge = e*A0; return nu((gn + ge)/A0)*(gn + ge) - nu(ge/A0)*ge
def fit_mond(e=0.0, floor=3.0, a0fac=1.0):
    global A0
    A0s = A0; A0 = A0s*a0fac
    def f(p):
        if not (0.2 < p[0] < 2.5): return 1e9
        return chi_rc(np.sqrt(gmond(gbar(R, p[0]), e)*R), floor)
    r = optimize.minimize_scalar(lambda u: f([u]), bounds=(0.2, 2.5), method='bounded'); ups = r.x
    res = dict(e=e, a0fac=a0fac, ups=round(ups, 3), Mstar=float('%.3g' % (ups*(MB + MASS['thin'] + MASS['thick']))), chi2_rc=round(r.fun, 1),
               V27=round(float(np.sqrt(gmond(gbar(np.array([27.31]), ups), e)*27.31)[0]), 1))
    for rr in (50, 100):
        v = float(np.sqrt(gmond(gbar(np.array([float(rr)]), ups), e)*rr)[0]); res[f'V{rr}'] = round(v, 1); res[f'M{rr}'] = float('%.3g' % (v*v*rr/G))
    res['tracer_sigma'] = {t['name'].split(' (')[0]: round((res['V%d' % t['r']] - t['V'])/(t['Vep'] if res['V%d' % t['r']] > t['V'] else t['Vem']), 1) for t in TR}
    # at McMillan's own stellar mass
    v1 = {rr: float(np.sqrt(gmond(gbar(np.array([float(rr)]), 1.0), e)*rr)[0]) for rr in (27.31, 50, 100)}
    res['at_McMillan'] = {str(k_): round(v_, 1) for k_, v_ in v1.items()}
    res['at_McMillan_chi2_rc'] = round(chi_rc(np.sqrt(gmond(gbar(R, 1.0), e)*R), floor), 1)
    A0 = A0s
    return res
out['mond'] = [fit_mond(0.0), fit_mond(0.01), fit_mond(0.03)]

# ---- 4. what gravitates inside 50 kpc, against the baryons that are seen ----
bary50 = MB + sum(float(np.sum((2*np.pi*ra*sigma(ra, *DISCS[k_])*da)[ra < 50])) for k_ in DISCS) + 3.8e9
out['inside_50kpc'] = dict(baryons=float('%.3g' % bary50), gravitating=3.36e11, ratio=round(3.36e11/bary50, 1), missing=float('%.3g' % (3.36e11 - bary50)),
                           hot_gas=3.8e9, hot_gas_factor_needed=round((3.36e11 - bary50 + 3.8e9)/3.8e9))
# curves for the figure
rr = np.linspace(6, 105, 199)
fc = out['nfw']['curve_only']; ft = out['nfw']['curve_plus_tracers']; m0 = out['mond'][0]; m2 = out['mond'][2]
out['fig'] = dict(r=[round(x, 2) for x in rr],
                  nfw_curve=[round(float(np.sqrt(gbar(np.array([x]))[0]*x + vnfw2(x, fc['M200'], fc['c']))), 1) for x in rr],
                  nfw_all=[round(float(np.sqrt(gbar(np.array([x]))[0]*x + vnfw2(x, ft['M200'], ft['c']))), 1) for x in rr],
                  mond=[round(float(np.sqrt(gmond(gbar(np.array([x]), m0['ups']), 0.0)*x)[0]), 1) for x in rr],
                  mond_efe=[round(float(np.sqrt(gmond(gbar(np.array([x]), m2['ups']), 0.03)*x)[0]), 1) for x in rr],
                  bary=[round(float(np.sqrt(gbar(np.array([x]))[0]*x)), 1) for x in rr],
                  kepler=[round(float(V[-1]*math.sqrt(R[-1]/x)), 1) if x >= R[-1] else None for x in rr])
json.dump(out, open(os.path.join(HERE, 'mw_mass.json'), 'w'), indent=1)
for k_ in out:
    if k_ != 'fig': print(k_, json.dumps(out[k_]))

# ---- 5. MOND with the external-field strength free, fitted to the curve and the tracers together ----
def mond_joint(e, trs=BASE, floor=3.0):
    def f(u):
        vm = np.sqrt(gmond(gbar(R, u), e)*R)
        return chi_rc(vm, floor) + chi_tr(lambda r: float(np.sqrt(gmond(gbar(np.array([float(r)]), u), e)*r)[0]), trs)
    r = optimize.minimize_scalar(f, bounds=(0.2, 2.5), method='bounded'); u = r.x
    vm = np.sqrt(gmond(gbar(R, u), e)*R)
    v50 = float(np.sqrt(gmond(gbar(np.array([50.0]), u), e)*50)[0]); v100 = float(np.sqrt(gmond(gbar(np.array([100.0]), u), e)*100)[0])
    return dict(e=round(e, 4), ups=round(u, 3), Mstar=float('%.3g' % (u*(MB + MASS['thin'] + MASS['thick']))), chi2=round(r.fun, 1), chi2_rc=round(chi_rc(vm, floor), 1), V50=round(v50, 1), V100=round(v100, 1))
scan = [mond_joint(e) for e in np.linspace(0, 0.03, 61)]
best = min(scan, key=lambda d: d['chi2'])
ok = [d['e'] for d in scan if d['chi2'] <= best['chi2'] + 1]
out['mond_efe_scan'] = dict(best=best, e_1sigma=[min(ok), max(ok)], nfw_joint_chi2=out['nfw']['curve_plus_tracers']['chi2'], no_efe=scan[0])
# same with the stream value at 50 kpc instead of the halo stars
scanK = [mond_joint(e, [TR[1], TR[2]]) for e in np.linspace(0, 0.03, 61)]; bK = min(scanK, key=lambda d: d['chi2'])
out['mond_efe_scan_koposov'] = bK
json.dump(out, open(os.path.join(HERE, 'mw_mass.json'), 'w'), indent=1)
print('scan', json.dumps(out['mond_efe_scan'])); print('scanK', json.dumps(bK))

# ---- figure curves: MOND with the best external field; concentration and stellar fraction of the NFW fit ----
bE = out['mond_efe_scan']['best']
out['fig']['mond_efe'] = [round(float(np.sqrt(gmond(gbar(np.array([x]), bE['ups']), bE['e'])*x)[0]), 1) for x in rr]
out['fig']['mond_efe_e'] = bE['e']
fc = out['nfw']['curve_only']
c_dm = 10**(0.905 - 0.101*math.log10(fc['M200']*0.674/1e12))          # Dutton & Maccio 2014, NFW c200 at z = 0, scatter 0.11 dex
out['nfw_context'] = dict(c_expected=round(c_dm, 1), c_fit=fc['c'], dex_above=round(math.log10(fc['c']/c_dm), 2), sigma=round(math.log10(fc['c']/c_dm)/0.11, 1),
                          mstar_over_m200=round(54.3e9/fc['M200'], 3), baryon_over_m200=round(out['baryon_total']/fc['M200'], 3), cosmic=round(0.0493/0.315, 3))
json.dump(out, open(os.path.join(HERE, 'mw_mass.json'), 'w'), indent=1)
print('ctx', json.dumps(out['nfw_context']))
bE['V27'] = round(float(np.sqrt(gmond(gbar(np.array([27.31]), bE['ups']), bE['e'])*27.31)[0]), 1)
bE['tracer_sigma'] = {t['name'].split(' (')[0]: round((bE['V%d' % t['r']] - t['V'])/(t['Vep'] if bE['V%d' % t['r']] > t['V'] else t['Vem']), 1) for t in TR}
bE['chi2_tracers'] = round(bE['chi2'] - bE['chi2_rc'], 1)
out['mond_efe_scan']['best'] = bE
json.dump(out, open(os.path.join(HERE, 'mw_mass.json'), 'w'), indent=1)
print('best', json.dumps(bE)); print('nfwjoint', json.dumps(out['nfw']['curve_plus_tracers']))
