"""Workbench Entry 122: how deep are cluster potential wells? Gravitational redshift against a cluster of plasma alone.

Light from a galaxy deep in a cluster climbs out of the cluster's gravity and arrives redder than light from the same galaxy outside it.
Galaxies in the outskirts sit higher in the well than the central (brightest) galaxy, so on average they look blueshifted relative to it,
by about 10 km/s for clusters of 10^14 Msun. Electric and magnetic forces do not do this: only the depth of the gravitational potential
does, and it does not depend on how galaxies move. So it weighs the well directly. If clusters hold only their plasma and stars, as
plasma cosmology has it, the well is about 6 times shallower than with dark matter.
The complication (Zhao, Peacock & Li 2013; Kaiser 2013): motions add their own shifts of the same size. Transverse Doppler and the
light-cone effect add +2.5 sigma_los^2/c; the selection of galaxies by apparent brightness (surface-brightness modulation) subtracts about
K sigma_los^2/c with K ~ 10 for SDSS; the brightest galaxy's own terms add about +2.2 km/s. These depend on the measured velocity
dispersion, which is the same in both readings, so they shift both predictions alike.
Model: NFW halo (M200 = 2e14 or 4e14 Msun, c = 5, z = 0.24, Planck 2018), galaxies tracing mass, isotropic orbits; shifts averaged
along the line of sight; bins at 0.6, 1.6, 3.3 and 5.2 Mpc weighted 1:1:3:3 as in Wojtak et al. 2011. Plasma-only clusters: the same
profile scaled by the baryon fraction 0.157, generous because real gas is more extended than the dark matter.
"""
import json, math, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
G = 4.30091e-9; C = 299792.458; h = 0.674; Om = 0.315; FB = 0.157
z = 0.24; E2 = Om*(1+z)**3 + 1 - Om; rho_c = 2.775e11*h*h*E2
BINS = np.array([0.6, 1.6, 3.3, 5.2]); WTS = np.array([1, 1, 3, 3.])
def nfw(M200, c=5.0, rmax=6.0):
    R200 = (3*M200/(4*math.pi*200*rho_c))**(1/3); rs = R200/c; m = math.log(1+c) - c/(1+c)
    A = G*M200/m                                     # G 4 pi rho_s rs^3
    r = np.geomspace(1e-3, rmax, 4000)
    x = r/rs; rho = 1/(x*(1+x)**2)                   # shape only
    Mr = A/G*(np.log(1+x) - x/(1+x))
    Phi = -A*np.log(1+x)/r; Phi0 = -A/rs
    # isotropic Jeans: rho sigma^2 = int_r^rmax rho G M / r'^2 dr'
    integrand = rho*G*Mr/r**2
    cum = np.concatenate([[0], np.cumsum(0.5*(integrand[1:]+integrand[:-1])*np.diff(r))])
    sig2 = (cum[-1] - cum)/rho
    return dict(R200=R200, rs=rs, r=r, rho=rho, Phi=Phi, Phi0=Phi0, sig2=sig2)
def project(H, R):
    zz = np.linspace(0, math.sqrt(max(6.0**2 - R**2, 1e-6)), 1500); rr = np.sqrt(R**2 + zz**2)
    rho = np.interp(rr, H['r'], H['rho']); Phi = np.interp(rr, H['r'], H['Phi']); s2 = np.interp(rr, H['r'], H['sig2'])
    W = np.trapezoid(rho, zz)
    return np.trapezoid(rho*Phi, zz)/W, np.trapezoid(rho*s2, zz)/W
out = {'bins': BINS.tolist(), 'models': {}}
Rg = np.concatenate([np.linspace(0.1, 1, 10), np.linspace(1.25, 6, 20)])
for M in (2e14, 4e14):
    H = nfw(M)
    gr = []; s2l = []
    for R in Rg:
        pb, s2 = project(H, R); gr.append(-(pb - H['Phi0'])/C); s2l.append(s2)
    gr = np.array(gr); s2l = np.array(s2l)
    sig_06 = math.sqrt(np.interp(0.6, Rg, s2l))
    print(f"M200 = {M:.0e}: R200 = {H['R200']:.2f} Mpc; gravitational shift at 0.6, 1.6, 3.3, 5.2 Mpc: " + ", ".join(f"{np.interp(b, Rg, gr):+.1f}" for b in BINS) + f" km/s; sigma_los(0.6 Mpc) = {sig_06:.0f} km/s")
    res = {}
    for K in (5.0, 10.0, 13.5):
        kin = (2.5 - K)*s2l/C + 2.2
        for lab, f in (('dark matter', 1.0), ('plasma only', FB)):
            tot = f*gr + kin
            binv = np.array([np.interp(b, Rg, tot) for b in BINS]); integ = float(np.sum(WTS*binv)/WTS.sum())
            res[f'{lab}|{K}'] = dict(profile=tot.round(3).tolist(), bins=binv.round(2).tolist(), integ=integ)
    out['models'][f'{M:.0e}'] = dict(R=Rg.round(3).tolist(), gr=gr.round(3).tolist(), sig_06=sig_06, res=res,
                                     gr_bins=[float(np.interp(b, Rg, gr)) for b in BINS])
    for K in (5.0, 10.0, 13.5):
        d, p = res[f'dark matter|{K}'], res[f'plasma only|{K}']
        print(f"   K = {K:4.1f}: dark matter, bins {d['bins']} -> {d['integ']:+.1f}; plasma only, bins {p['bins']} -> {p['integ']:+.1f} km/s")
# Kaiser's own numbers at 0.6 Mpc for comparison: GR -5.7, TD+LC +1.9, SB -10.3, BCG +2.2; observed -2.6 ± 6
meas = [('Wojtak et al. 2011, 7,800 SDSS clusters', -7.7, 3.0, 3.0), ('Sadeh et al. 2015, SDSS + BOSS', -11.0, 5.0, 7.0),
        ('Mpetha et al. 2021, SPIDERS X-ray clusters', -13.5, 4.7, 4.7), ('Rosselli et al. 2023, 3,058 SDSS clusters', -11.4, 3.3, 3.3)]
out['meas'] = meas
print("measured, integrated:")
for n, v, lo, hi in meas: print(f"   {n}: {v:+.1f} (−{lo}, +{hi})")
# compare; for each model the error used is the one on the side of the prediction
def tension(pred, v, lo, hi):
    s = hi if pred > v else lo
    return (pred - v)/s
out['tension'] = {}
for Mk in ('2e+14', '4e+14'):
    for K in (5.0, 10.0, 13.5):
        for lab in ('dark matter', 'plasma only'):
            p = out['models'][Mk]['res'][f'{lab}|{K}']['integ']
            t = [tension(p, v, lo, hi) for n, v, lo, hi in meas]
            out['tension'][f'{Mk}|{K}|{lab}'] = t
            print(f"   M {Mk}, K {K:4.1f}, {lab:11s}: prediction {p:+.1f}; offsets " + ", ".join(f"{x:+.1f}σ" for x in t))
# the samples overlap (three are SDSS), so we do not multiply them; the least-overlapping pair: Rosselli (SDSS optical) and Mpetha (X-ray chosen)
for Mk in ('2e+14', '4e+14'):
    for K in (5.0, 10.0, 13.5):
        for lab in ('dark matter', 'plasma only'):
            t = out['tension'][f'{Mk}|{K}|{lab}']; chi = t[2]**2 + t[3]**2
            out['tension'][f'{Mk}|{K}|{lab}|pair'] = chi
json.dump(out, open(os.path.join(HERE, 'grz.json'), 'w'), indent=1)
# depth factor f: the gravitational term scaled by f (1 = dark-matter halo, 0.157 = plasma and stars only), kinematic terms fixed
from math import erf, sqrt
def p2sig(chi, dof=2):
    p = math.exp(-chi/2) if dof == 2 else None
    lo, hi = 0.0, 10.0
    for _ in range(80):
        m = (lo+hi)/2; pm = 1 - erf(m/sqrt(2))
        lo, hi = (m, hi) if pm > p else (lo, m)
    return p, (lo+hi)/2
out['depth'] = {}
for Mk in ('2e+14', '4e+14'):
    M = out['models'][Mk]
    g_int = float(np.sum(WTS*np.array(M['gr_bins']))/WTS.sum())
    for K in (5.0, 10.0, 13.5):
        kin = M['res'][f'dark matter|{K}']['integ'] - g_int
        fs = []
        for n, v, lo, hi in meas:
            s = (lo+hi)/2; fs.append(((v - kin)/g_int, s/abs(g_int)))
        # the two least-overlapping samples: X-ray chosen (SPIDERS) and optical (Rosselli)
        w = np.array([1/fs[2][1]**2, 1/fs[3][1]**2]); fp = float((w[0]*fs[2][0] + w[1]*fs[3][0])/w.sum()); ep = float(1/math.sqrt(w.sum()))
        wa = np.array([1/f[1]**2 for f in fs]); fa = float(np.sum(wa*np.array([f[0] for f in fs]))/wa.sum()); ea = float(1/math.sqrt(wa.sum()))
        out['depth'][f'{Mk}|{K}'] = dict(g_int=g_int, kin=kin, each=fs, pair=(fp, ep), all=(fa, ea), z_plasma_pair=(fp - FB)/ep, z_plasma_all=(fa - FB)/ea)
        print(f"   M {Mk}, K {K:4.1f}: gravity term {g_int:+.1f}, motions {kin:+.1f} km/s; depth factor per survey " + ", ".join(f"{a:.2f}±{b:.2f}" for a, b in fs) +
              f"; SPIDERS+Rosselli {fp:.2f} ± {ep:.2f} (plasma only, 0.157, at {(fp-FB)/ep:.1f}σ); all four as if independent {fa:.2f} ± {ea:.2f} ({(fa-FB)/ea:.1f}σ)")
json.dump(out, open(os.path.join(HERE, 'grz.json'), 'w'), indent=1)
