"""tau_from_xhi.py - Workbench Entry 95
The optical depth to reionization from the neutral fraction of the intergalactic gas, measured in the spectra of quasars
and galaxies, with no CMB data. Then: how high can tau go and still fit, and could an early burst of ionization beyond the
reach of these spectra (z = 15-30) supply the tau = 0.09 that DESI's ruler with the small-scale CMB asks for?

Data: the 19 neutral-fraction constraints compiled by Elbers (2025, arXiv:2508.21069, table 1), from
  Jin et al. 2023 (dark pixels in the Lyman-alpha and beta forests), Greig et al. 2024 (forest), Umeda et al. 2025 (damping
  wings of JWST galaxies), Mason et al. 2018, 2019 and Hoag et al. 2019 (Lyman-alpha emitters; these three measure the
  change since z = 6), Greig et al. 2022 (quasar damping wings, averaged), Bruton et al. 2023, Curtis-Lake et al. 2023.
Likelihood: two-sided errors as split normals; upper and lower limits as one-sided normals (a bare limit gets a soft edge of
  0.10); a point with a redshift range is compared with the model averaged over that range.
Histories (ionized fraction of hydrogen, helium's first electron with it, its second at z = 3.5 as in Planck):
  A  tanh in (1+z)^1.5, midpoint and width free (Planck's form has the width fixed at 0.5)
  B  power law in z from z_end (fully ionized) to z = 20 (neutral), end and exponent free (Planck's 'asymmetric' form)
Cosmology: the non-CMB values of Entries 92 and 94 (Omega_b h^2 = 0.0222, Omega_m = 0.2975, h = 0.6851, Y_p = 0.247).
"""
import json, math, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
WB, OM, H, YP = 0.0222, 0.2975, 0.6851, 0.247
fHe = YP/(4*(1-YP)); nH0 = (1-YP)*WB*1.87847e-29/1.6726e-24
sigT, Mpc, c = 6.6524e-25, 3.0857e24, 2.99792458e10; H0 = 100*H*1e5/Mpc
ZG = np.linspace(0, 40, 16001); DZ = ZG[1]-ZG[0]
KER = nH0*(1+ZG)**2/np.sqrt(OM*(1+ZG)**3 + 1 - OM)*sigT*c/H0          # d tau / dz for one free electron per hydrogen atom
HE2 = fHe*0.5*(1+np.tanh((3.5-ZG)/0.5))
def tau_of(Q):                                                       # Q: ionized fraction of hydrogen on ZG
    return float(np.trapezoid(KER*(Q*(1+fHe) + HE2), ZG))
def tau_range(z1, z2):                                               # fully ionized (H + He I) between z1 and z2
    m = (ZG >= z1) & (ZG <= z2); return float(np.trapezoid(KER[m]*(1+fHe), ZG[m]))

# ---- data: (z, z_lo, z_hi, kind, value, err_minus, err_plus, ref) ; kind: 'v' value, 'ul' upper limit, 'll' lower limit, 'd6' change since z = 6
D = [(5.5, 0.2, 0.2, 'ul', 0.09, 0.08, 0.08, 'Jin 2023'), (5.7, 0.2, 0.2, 'ul', 0.16, 0.14, 0.14, 'Jin 2023'), (5.8, 0.4, 0.5, 'v', 0.25, 0.20, 0.10, 'Umeda 2025'),
     (5.9, 0.2, 0.2, 'ul', 0.28, 0.08, 0.08, 'Jin 2023'), (5.90, 0.05, 0.05, 'ul', 0.21, None, None, 'Greig 2024'), (5.95, 0.05, 0.05, 'ul', 0.20, None, None, 'Greig 2024'),
     (6.05, 0.05, 0.05, 'ul', 0.21, None, None, 'Greig 2024'), (6.15, 0.05, 0.05, 'v', 0.20, 0.12, 0.14, 'Greig 2024'), (6.35, 0.05, 0.05, 'v', 0.29, 0.13, 0.14, 'Greig 2024'),
     (6.55, 0.05, 0.05, 'ul', 0.18, None, None, 'Greig 2024'), (6.9, 0.5, 0.5, 'd6', 0.59, 0.15, 0.11, 'Mason 2018'), (7.0, 0.4, 0.5, 'v', 0.65, 0.35, 0.27, 'Umeda 2025'),
     (7.29, 0.27, 0.27, 'v', 0.49, 0.14, 0.13, 'Greig 2022'), (7.6, 0.6, 0.6, 'd6', 0.88, 0.10, 0.05, 'Hoag 2019'), (7.9, 0.6, 0.6, 'll', 0.76, None, None, 'Mason 2019'),
     (8.6, 0.9, 1.0, 'll', 1.00, 0.20, 0.0, 'Umeda 2025'), (10.4, 0.8, 2.0, 'll', 1.00, 0.40, 0.0, 'Umeda 2025'), (10.6, 0.0, 0.0, 'ul', 0.88, None, None, 'Bruton 2023'),
     (11.49, 0.05, 0.04, 'v', 0.84, 0.27, 0.14, 'Curtis-Lake 2023')]
SOFT = 0.10
def xavg(XH, z, lo, hi):
    if lo + hi < 1e-6: return float(np.interp(z, ZG, XH))
    m = (ZG >= z-lo) & (ZG <= z+hi); return float(XH[m].mean())
def chi2(Q):
    XH = 1 - Q; x6 = float(np.interp(6.0, ZG, XH)); c2 = 0.0
    for z, lo, hi, kind, v, em, ep, ref in D:
        x = xavg(XH, z, lo, hi)
        if kind == 'd6': x = x - x6
        if kind in ('v', 'd6'): s = ep if x > v else em; c2 += ((x - v)/s)**2
        elif kind == 'ul' and x > v: c2 += ((x - v)/(ep if ep else SOFT))**2
        elif kind == 'll' and x < v: c2 += ((v - x)/(em if em else SOFT))**2
    return c2
def QA(zmid, dz):
    y = (1+ZG)**1.5; ym = (1+zmid)**1.5; dy = 1.5*math.sqrt(1+zmid)*dz
    return 0.5*(1+np.tanh((ym - y)/dy))
def QB(zend, alpha, zearly=20.0):
    q = np.clip((zearly - ZG)/(zearly - zend), 0, 1)**alpha; q[ZG <= zend] = 1.0; return q
out = {'data': [dict(z=d[0], kind=d[3], value=d[4], ref=d[7]) for d in D]}
fams = {'A': (QA, np.linspace(5.5, 10.5, 201), np.linspace(0.1, 4.0, 157), ('zmid', 'width')),
        'B': (QB, np.linspace(4.6, 7.0, 121), np.exp(np.linspace(math.log(0.3), math.log(15.0), 161)), ('zend', 'alpha'))}
for name, (fn, p1, p2, labels) in fams.items():
    C = np.zeros((len(p1), len(p2))); T = np.zeros_like(C); ZM = np.zeros_like(C)
    for i, a in enumerate(p1):
        for j, b in enumerate(p2):
            Q = fn(a, b); C[i, j] = chi2(Q); T[i, j] = tau_of(Q); ZM[i, j] = float(np.interp(0.5, Q[::-1], ZG[::-1]))
    W = np.exp(-0.5*(C - C.min())); W /= W.sum()
    ts = T.ravel(); o = np.argsort(ts); cw = np.cumsum(W.ravel()[o])
    q = lambda p: float(ts[o][np.searchsorted(cw, p)])
    zs = ZM.ravel(); o2 = np.argsort(zs); cw2 = np.cumsum(W.ravel()[o2]); qz = lambda p: float(zs[o2][np.searchsorted(cw2, p)])
    ib = np.unravel_index(C.argmin(), C.shape)
    hi = T >= 0.09; c09 = float(C[hi].min()) if hi.any() else None
    i09 = np.unravel_index(np.where(hi, C, np.inf).argmin(), C.shape)
    res = dict(chi2_min=round(float(C.min()), 2), n=len(D), best={labels[0]: round(float(p1[ib[0]]), 3), labels[1]: round(float(p2[ib[1]]), 3)},
               tau_best=round(float(T[ib]), 4), tau_median=round(q(0.5), 4), tau_16=round(q(0.16), 4), tau_84=round(q(0.84), 4), tau_2p5=round(q(0.025), 4), tau_97p5=round(q(0.975), 4),
               zmid_median=round(qz(0.5), 2), zmid_16=round(qz(0.16), 2), zmid_84=round(qz(0.84), 2),
               chi2_at_tau_ge_0p09=round(c09, 1) if c09 is not None else None, best_at_0p09={labels[0]: round(float(p1[i09[0]]), 3), labels[1]: round(float(p2[i09[1]]), 3), 'tau': round(float(T[i09]), 4)},
               chi2_at_tau_ge_0p054=round(float(C[T >= 0.054].min()), 2))
    # the history curves for the figure: best fit and the best history with tau >= 0.09
    res['curve_best'] = [[round(float(z), 2), round(float(1 - x), 4)] for z, x in zip(ZG[::50], fn(p1[ib[0]], p2[ib[1]])[::50]) if z <= 16]
    res['curve_09'] = [[round(float(z), 2), round(float(1 - x), 4)] for z, x in zip(ZG[::50], fn(p1[i09[0]], p2[i09[1]])[::50]) if z <= 16]
    out[name] = res
    print(name, {k: v for k, v in res.items() if not k.startswith('curve')})
# leave-one-method-out for family A
groups = {'forest (Jin, Greig 2024)': ('Jin 2023', 'Greig 2024'), 'emitters (Mason, Hoag)': ('Mason 2018', 'Mason 2019', 'Hoag 2019'),
          'JWST galaxy damping wings (Umeda, Bruton, Curtis-Lake)': ('Umeda 2025', 'Bruton 2023', 'Curtis-Lake 2023'), 'quasar damping wings (Greig 2022)': ('Greig 2022',)}
full = list(D); out['drop'] = {}
for g, refs in groups.items():
    D[:] = [d for d in full if d[7] not in refs]
    fn, p1, p2, labels = fams['A']; C = np.zeros((len(p1), len(p2))); T = np.zeros_like(C)
    for i, a in enumerate(p1):
        for j, b in enumerate(p2):
            Q = fn(a, b); C[i, j] = chi2(Q); T[i, j] = tau_of(Q)
    W = np.exp(-0.5*(C - C.min())); W /= W.sum(); ts = T.ravel(); o = np.argsort(ts); cw = np.cumsum(W.ravel()[o])
    q = lambda p: float(ts[o][np.searchsorted(cw, p)])
    out['drop'][g] = dict(n=len(D), tau_median=round(q(0.5), 4), tau_16=round(q(0.16), 4), tau_84=round(q(0.84), 4), tau_97p5=round(q(0.975), 4))
    print('without', g, out['drop'][g])
D[:] = full
# ---- floors and the early burst
out['floor'] = dict(instant_5p3=round(tau_of((ZG <= 5.3).astype(float)), 4), instant_6=round(tau_of((ZG <= 6.0).astype(float)), 4))
full1530 = tau_range(15, 30)
need = {k: round(0.090 - out[k]['tau_median'], 4) for k in ('A', 'B')}
out['flash'] = dict(tau_full_15_30=round(full1530, 4), need_from_A=need['A'], need_from_B=need['B'], fraction_needed_A=round(need['A']/full1530, 3),
                    planck_tau_15_30_95=0.020, ratio_to_limit_A=round(need['A']/0.020, 2))
print('floors', out['floor']); print('flash', out['flash'])
json.dump(out, open(os.path.join(HERE, 'tau_from_xhi.json'), 'w'), indent=1)
