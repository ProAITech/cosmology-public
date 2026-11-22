"""progenitors.py - Workbench Entry 101
How early can an expanding universe make a galaxy's stars? The ceiling from the merger history of its halo.
For a galaxy seen at z0 inside a halo of mass M0, extended Press-Schechter theory (Bond et al. 1991; Lacey & Cole 1993) gives the share of
M0 that sat in collapsed progenitors above a minimum mass M_min at an earlier redshift z1:
   f = erfc[ (delta_c(z1) - delta_c(z0)) / sqrt(2 (S(M_min) - S(M0))) ],   delta_c(z) = 1.686 / D(z),  S = sigma^2(M)
Stars formed by z1, in the galaxy's own halo or in any that later merged with it, cannot exceed f_b x f x M0. The ratio of the stars a galaxy is
said to have formed by then to that ceiling is the efficiency the expanding universe needs; above 1 is impossible.
Linear power: Eisenstein & Hu (1998) no-wiggle transfer, Planck 2018 (h 0.6736, Om 0.3153, Ob 0.0493, n_s 0.9649, sigma_8 0.811).
M_min: the atomic-cooling halo, virial temperature 1e4 K, ~1e8 Msun ((1+z)/10)^-1.5; and 1e6 Msun for the most generous molecular-cooling case.
M0: 1e13 and 3e13 Msun, both generous for galaxies of 2e10-2e11 Msun in stars at z = 2.7-7.3 (a 3e13 halo at z = 2.7 is a protocluster core).
"""
import json, math, os
import numpy as np
from scipy import integrate, special
HERE = os.path.dirname(os.path.abspath(__file__))
h, Om, Ob, ns, s8 = 0.6736, 0.3153, 0.0493, 0.9649, 0.811
fb = Ob/Om; rho_m = 2.775e11*h*h*Om                      # Msun / Mpc^3
def T_EH(k):                                             # k in 1/Mpc; Eisenstein & Hu 1998 no-wiggle
    th = 2.7255/2.7; omh2 = Om*h*h; obh2 = Ob*h*h; fbm = Ob/Om
    s = 44.5*math.log(9.83/omh2)/math.sqrt(1 + 10*obh2**0.75)
    aG = 1 - 0.328*math.log(431*omh2)*fbm + 0.38*math.log(22.3*omh2)*fbm**2
    G = Om*h*(aG + (1 - aG)/(1 + (0.43*k*s)**4))
    q = (k/h)*th**2/G                                      # EH98 eq. 28: k in h/Mpc over Gamma_eff = Om h (...)
    L0 = np.log(2*math.e + 1.8*q); C0 = 14.2 + 731/(1 + 62.5*q)
    return L0/(L0 + C0*q*q)
kk = np.logspace(-5, 3.5, 6000)
Pk = kk**ns*T_EH(kk)**2
def sigma(R, P=None):
    P = Pk if P is None else P
    x = kk*R; W = 3*(np.sin(x) - x*np.cos(x))/x**3
    return math.sqrt(np.trapezoid(P*W*W*kk**2, kk)/(2*math.pi**2))
norm = s8/sigma(8/h); Pk = Pk*norm**2
def sigM(M): return sigma((3*M/(4*math.pi*rho_m))**(1/3))
def D(z):
    a = 1/(1+z); E = lambda aa: math.sqrt(Om/aa**3 + 1 - Om)
    g = lambda aa: 2.5*Om*E(aa)*integrate.quad(lambda x: 1/(x*E(x))**3, 1e-8, aa)[0]
    return g(a)/g(1.0)
def age(z): tH = 977.79/(100*h); OL = 1 - Om; return 2/(3*math.sqrt(OL))*tH*math.asinh(math.sqrt(OL/Om)*(1+z)**-1.5)
from scipy.optimize import brentq
def z_at(t): return brentq(lambda z: age(z) - t, 0, 300)
dc = lambda z: 1.686/D(z)
def fcoll(z1, z0, M0, Mmin):
    S1 = sigM(Mmin)**2; S0 = sigM(M0)**2
    return float(special.erfc((dc(z1) - dc(z0))/math.sqrt(2*(S1 - S0))))
def Matom(z): return 1e8*((1+z)/10)**-1.5
# galaxies: name, z_obs, stellar mass now, list of (label, cosmic time in Gyr, fraction of the stellar mass formed by then)
G = [('Eridu, z = 2.675 (McConachie et al. 2025)', 2.675, 10**10.96, [('half its stars by 108 Myr (fiducial fit)', 0.108, 0.5), ('half by 242 Myr (continuity fit)', 0.242, 0.5),
                                                                       ('90% by 353 Myr (fiducial fit)', 0.353, 0.9), ('half by 496 Myr (double power law)', 0.496, 0.5)]),
     ('ZF-UDS-7329, z = 3.205 (Glazebrook et al. 2024)', 3.205, 10**11.26, [('half its stars ~1.5 Gyr before (z ~ 11)', None, 0.5)]),
     ('GS-9209, z = 4.658 (Carnall et al. 2023)', 4.658, 3.8e10, [('half by ~0.70 Gyr (z ~ 7.2)', 0.70, 0.5)]),
     ('RUBIES-UDS-QG-z7, z = 7.29 (Weibel et al. 2025)', 7.29, 10**10.23, [('half by z ~ 8.5', None, 0.5)])]
out = {'params': dict(h=h, Om=Om, Ob=Ob, ns=ns, s8=s8), 'sigma_checks': {('%.0e' % M): round(sigM(M), 3) for M in (1e6, 1e8, 1e10, 1e12, 1e13)},
       'D_checks': {str(z): round(D(z), 4) for z in (0, 2.675, 7, 12, 28)}, 'rows': []}
print('sigma(M):', out['sigma_checks'], ' D(z):', out['D_checks'])
for name, z0, Ms, events in G:
    for lab, t, frac in events:
        if t is None: t = age(z0) - 1.5 if 'ZF' in name else age(8.5)
        z1 = z_at(t); Mst = frac*Ms
        row = dict(galaxy=name, event=lab, t_Gyr=round(t, 3), z1=round(z1, 1), Mstar_formed=float('%.3g' % Mst))
        for M0 in (1e13, 3e13):
            for mlab, Mmin in (('atomic', Matom(z1)), ('molecular', 1e6)):
                f = fcoll(z1, z0, M0, Mmin); ceil = fb*f*M0
                row['eff_M0=%.0e_%s' % (M0, mlab)] = float('%.3g' % (Mst/ceil)); row['fcoll_M0=%.0e_%s' % (M0, mlab)] = float('%.3g' % f)
        out['rows'].append(row)
        print('%-46s %-42s z1 %5.1f  M* %.2e  eff needed (1e13, atomic) %.3g  (3e13, molecular) %.3g' % (name, lab, z1, Mst, row['eff_M0=1e+13_atomic'], row['eff_M0=3e+13_molecular']))
# the ceiling curve for the figure: for M0 = 1e13 at z0 = 2.675 and 4.658, the most stellar mass that can exist (efficiency 1) against cosmic time
curves = {}
for z0 in (2.675, 4.658):
    pts = []
    for t in np.round(np.concatenate([np.linspace(0.06, 0.2, 15), np.linspace(0.22, 1.2, 50)]), 3):
        if t >= age(z0): break
        z1 = z_at(t); pts.append([float(t), float('%.4g' % (fb*fcoll(z1, z0, 1e13, Matom(z1))*1e13)), float('%.4g' % (fb*fcoll(z1, z0, 3e13, 1e6)*3e13))])
    curves[str(z0)] = pts
out['curves'] = curves
json.dump(out, open(os.path.join(HERE, 'progenitors.json'), 'w'), indent=1)
