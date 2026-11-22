"""efe.py - Workbench Entry 104: does the Milky Way sit in the external field MOND needs?

Entry 103 (mw_mass.py) found that MOND fits Gaia's rotation curve and the halo tracers only with a Newtonian external
field e = g_Ne/a0 of about 0.003. Here the field is measured two other ways:
  1. From what pulls on the Milky Way: the Local Group's motion against the CMB, 627 +- 22 km/s toward (l, b) = (276, 30)
     (Kogut et al. 1993), read through linear theory with baryons tracing the matter, plus Andromeda's own pull
     (stellar mass 1.0-1.5e11, Tamm et al. 2012; HI ~5e9; distance 765 kpc) toward (l, b) = (121.2, -21.6).
  2. From the escape speed at the Sun: 523.7 +12.8 -13.5 km/s (Wu et al. 2025, LAMOST K giants) and 528 +24 -25 (Deason et al. 2019).
     In MOND an isolated galaxy has no escape speed; the external field sets it.
And the dark-halo fit of Entry 103 is checked against the same escape speed and against a published escape-speed halo mass.
"""
import json, os, math, sys, io, contextlib
import numpy as np
from scipy import optimize, integrate
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
with contextlib.redirect_stdout(io.StringIO()):
    import mw_mass as M
G, A0 = M.G, M.A0
out = {}

# ---- 1. the field from what pulls on the Milky Way ----
def unit(l, b):
    l, b = math.radians(l), math.radians(b); return np.array([math.cos(b)*math.cos(l), math.cos(b)*math.sin(l), math.sin(b)])
H0, Om, Ob = 67.4, 0.315, 0.0493
def lss(v=627.0, om=Om):
    f = om**0.55
    gtot = 3*H0*om*v/(2*f)/1000.0                     # (km/s)^2/kpc; H0 in km/s/Mpc, so /1000 for per kpc
    return gtot/A0, gtot*(Ob/om)/A0                   # total Newtonian field and the baryons' share, in units of a0
gt, gb = lss()
gb_lo = lss(605.0, 0.33)[1]; gb_hi = lss(649.0, 0.30)[1]
def m31(Mb, d=765.0): return G*Mb/d**2/A0
e31 = m31(1.2e11); e31_lo, e31_hi = m31(1.05e11), m31(1.55e11)
u_lss, u_31 = unit(276, 30), unit(121.17, -21.57)
cosang = float(u_lss @ u_31)
def total(a, b): return float(np.linalg.norm(a*u_lss + b*u_31))
e_pull = total(gb, e31)
out['pull'] = dict(g_total_LSS=round(gt, 4), e_LSS_baryons=round(gb, 5), e_LSS_range=[round(gb_lo, 5), round(gb_hi, 5)], e_M31=round(e31, 6), e_M31_range=[round(e31_lo, 6), round(e31_hi, 6)],
                   angle_deg=round(math.degrees(math.acos(cosang)), 1), e_total=round(e_pull, 5),
                   e_total_range=[round(total(gb_lo, e31_hi), 5), round(total(gb_hi, e31_lo), 5)],
                   mond_field=round(float(M.nu(e_pull)*e_pull), 4), years_to_build_627=round(627e3/(float(M.nu(e_pull)*e_pull)*1.2e-10)/3.156e7/1e9, 1))

# ---- 2. MOND fits and escape speed as functions of the external field ----
R0 = 8.2
RF = np.concatenate([np.linspace(R0, 110, 600), np.geomspace(110.5, 1e6, 600)])
Mtot = lambda ups: ups*(M.MB + M.MASS['thin'] + M.MASS['thick']) + M.MASS['HI'] + M.MASS['H2']
def g_along(ups, e):
    gn = np.where(RF <= 105, M.gbar(np.minimum(RF, 105), ups), G*Mtot(ups)/RF**2)
    return M.gmond(gn, e)
def vesc(ups, e):
    g = g_along(ups, e); phi = integrate.trapezoid(g, RF)          # potential difference from R0 to 'infinity' (1e6 kpc)
    tail = g[-1]*RF[-1]                                             # remaining integral of a 1/r^2 tail
    return math.sqrt(2*(phi + tail))
def fit_curve(e, trs=M.BASE):
    f = lambda u: M.chi_rc(np.sqrt(M.gmond(M.gbar(M.R, u), e)*M.R)) + M.chi_tr(lambda r: float(np.sqrt(M.gmond(M.gbar(np.array([float(r)]), u), e)*r)[0]), trs)
    r = optimize.minimize_scalar(f, bounds=(0.2, 2.5), method='bounded'); return r.x, r.fun
grid = []
for e in [float('%.3g' % x) for x in np.geomspace(3e-4, 1e-2, 31)]:
    u, c2 = fit_curve(e)
    v50 = float(np.sqrt(M.gmond(M.gbar(np.array([50.0]), u), e)*50)[0]); v100 = float(np.sqrt(M.gmond(M.gbar(np.array([100.0]), u), e)*100)[0])
    grid.append(dict(e=e, ups=round(u, 3), Mstar=float('%.3g' % (u*(M.MB + M.MASS['thin'] + M.MASS['thick']))), chi2=round(c2, 1), V50=round(v50, 1), V100=round(v100, 1),
                     vesc=round(vesc(u, e), 1), mond_field=round(float(M.nu(e)*e), 4)))
out['grid'] = grid
VE, VEE = 523.7, 13.2
ee = np.array([d['e'] for d in grid]); vv = np.array([d['vesc'] for d in grid]); cc = np.array([d['chi2'] for d in grid])
le = np.log(ee)
e_esc = float(np.exp(np.interp(-VE, -vv, le))); e_esc_lo = float(np.exp(np.interp(-(VE + VEE), -vv, le))); e_esc_hi = float(np.exp(np.interp(-(VE - VEE), -vv, le)))
out['escape'] = dict(measured=[VE, 12.8, 13.5], deason=[528, 24, 25], e_needed=round(e_esc, 5), e_range=[round(e_esc_lo, 5), round(e_esc_hi, 5)],
                     vesc_at_curve_e=round(float(np.interp(math.log(0.003), le, vv)), 1), vesc_at_pull_e=round(float(np.interp(math.log(e_pull), le, vv)), 1))
cmin = float(cc.min())
out['curve'] = dict(e_best=0.003, chi2_min=cmin, chi2_at_pull=round(float(np.interp(math.log(e_pull), le, cc)), 1), chi2_at_escape=round(float(np.interp(math.log(e_esc), le, cc)), 1))
out['curve']['dchi2_pull'] = round(out['curve']['chi2_at_pull'] - cmin, 1); out['curve']['dchi2_escape'] = round(out['curve']['chi2_at_escape'] - cmin, 1)
# the fit at the pull's field: which tracers miss and by how much
u, c2 = fit_curve(e_pull)
v50 = float(np.sqrt(M.gmond(M.gbar(np.array([50.0]), u), e_pull)*50)[0]); v100 = float(np.sqrt(M.gmond(M.gbar(np.array([100.0]), u), e_pull)*100)[0])
out['at_pull'] = dict(e=round(e_pull, 5), ups=round(u, 3), V50=round(v50, 1), V100=round(v100, 1), sig50=round((v50 - 170.0)/3.8, 1), sig100=round((v100 - 155.2)/(5.4 if v100 > 155.2 else 5.6), 1), vesc=round(vesc(u, e_pull), 1))
u, c2 = fit_curve(e_esc)
v50 = float(np.sqrt(M.gmond(M.gbar(np.array([50.0]), u), e_esc)*50)[0]); v100 = float(np.sqrt(M.gmond(M.gbar(np.array([100.0]), u), e_esc)*100)[0])
out['at_escape'] = dict(e=round(e_esc, 5), ups=round(u, 3), V50=round(v50, 1), V100=round(v100, 1), sig50=round((v50 - 170.0)/3.8, 1), sig100=round((v100 - 155.2)/(5.4 if v100 > 155.2 else 5.6), 1))
# one field for everything: curve + tracers + escape speed
def joint(e):
    u, c2 = fit_curve(e); return c2 + ((vesc(u, e) - VE)/(12.8 if vesc(u, e) > VE else 13.5))**2
r = optimize.minimize_scalar(lambda x: joint(math.exp(x)), bounds=(math.log(3e-4), math.log(0.01)), method='bounded')
ej = math.exp(r.x); u, c2 = fit_curve(ej)
out['joint'] = dict(e=round(ej, 5), chi2=round(r.fun, 1), chi2_curve_tracers=round(c2, 1), vesc=round(vesc(u, ej), 1), ups=round(u, 3))

# ---- 3. the dark halo of Entry 103 against the same escape speed ----
nc = M.out['nfw']['curve_only']; M200, c = nc['M200'], nc['c']
r200 = (3*M200/(4*np.pi*200*M.RHOC))**(1/3); rs = r200/c; fc = math.log(1 + c) - c/(1 + c)
phi_nfw = lambda r: -G*M200/fc*math.log(1 + r/rs)/r
def phi_bar(r):                                            # discs + bulge, potential in the plane relative to infinity
    if r > 105: return -G*Mtot(1.0)/r
    tot = -G*M.MB/r
    for k_ in M.DISCS:
        S0, Rd, Rm = M.DISCS[k_]; m = 2*np.pi*M.ra*M.sigma(M.ra, S0, Rd, Rm)*M.da; tot += float(M.phi(np.array([r]), m)[0])
    return tot
vesc_nfw = math.sqrt(2*((phi_nfw(2*r200) + phi_bar(2*r200)) - (phi_nfw(R0) + phi_bar(R0))))
vesc_nfw_inf = math.sqrt(-2*(phi_nfw(R0) + phi_bar(R0)))
out['halo'] = dict(M200=M200, c=c, r200=round(r200, 1), vesc_2r200=round(vesc_nfw, 1), vesc_inf=round(vesc_nfw_inf, 1),
                   roche2024=dict(M200c=6.4e11, c=13.9), wu2025_M200=9.0e11)
json.dump(out, open(os.path.join(HERE, 'efe.json'), 'w'), indent=1)
for k_, v_ in out.items(): print(k_, json.dumps(v_))
