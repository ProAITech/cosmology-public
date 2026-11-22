"""pmf_halos.py - Workbench Entry 131
Can a primordial magnetic field build JWST's brightest early galaxies? Where its extra matter power lands in halo mass.
All field strengths below are the rms field sqrt<B^2>, the variable in which Fairbairn et al. (2026) quote their bounds; B_1Mpc follows from it.

PMF-induced matter power (Ralegankar, Garaldi & Viel 2025, JCAP 08, 011; used by Fairbairn et al. 2026):
  linear:  Delta_lin(k,a) = 0.918e-4 xi^2(a) (k/Mpc^-1)^(2 n_B + 10) (B_1Mpc/nG)^4 G_nB          (n_B < -1.5)
  n_B=-2:  Delta(k, z=100) = Delta_lin(k, 0.01) [1 + (k lJ)^p]^(-6/p),      lJ = [0.1 kJ B]^(2/(n_B+5)) Mpc, matter kJ=1.17, p=1.94
  n_B=+2:  Delta(k, z=100) = gamma (k lJ)^7 / [1 + (k lJ)^p]^(7/p),          matter kJ=0.593, p=2.36, gamma=0.244
  later:   x [xi_m(a)/xi_m(0.01)]^2
xi_b, xi_DM from the paper's linear equations (2.8)-(2.9) with zero initial conditions at recombination, thermal pressure neglected
(it matters only well above k ~ 100/Mpc). G_{-2} is the paper's integral (2.24), evaluated here numerically.
Check: the baryon plateau at z = 100 should come out near the paper's ~0.5.
Halo abundance: Sheth-Tormen with a sharp-k filter (c = 2.5), so power on small scales does not leak into large masses;
the boost is the ratio PMF+CDM over CDM with the same filter.
"""
import json, math, os, sys, io
import numpy as np
from scipy import integrate
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'early'))
_o = sys.stdout; sys.stdout = io.StringIO()
import progenitors as P
sys.stdout = _o
h, Om, rho_m, D = P.h, P.Om, P.rho_m, P.D

# --- growth of magnetically sourced perturbations (paper's parameters) ---
Ob_, Om_, OL_, aeq, arec = 0.0482, 0.308, 0.691, 2.94e-4, 1/1091
Odm_ = Om_ - Ob_
def Oa(a): return Om_*(1 + aeq/a) + OL_*a**3
def rhs(N, y):
    a = math.exp(N); xb, xbN, xd, xdN = y; O = Oa(a)
    xbNN = -0.5*xbN + 1.5*(Ob_/O)*xb + Om_/O + 1.5*(Odm_/O)*xd
    xdNN = -0.5*xdN + 1.5*(Odm_/O)*xd + 1.5*(Ob_/O)*xb
    return [xbN, xbNN, xdN, xdNN]
Ngrid = np.linspace(math.log(arec), math.log(1/6.0), 4000)
sol = integrate.solve_ivp(rhs, (Ngrid[0], Ngrid[-1]), [0, 0, 0, 0], t_eval=Ngrid, rtol=1e-9, atol=1e-12)
xb_, xd_ = sol.y[0], sol.y[2]; xm_ = (Ob_*xb_ + Odm_*xd_)/Om_
def xi(a, which='m'):
    arr = {'b': xb_, 'd': xd_, 'm': xm_}[which]
    return float(np.interp(math.log(a), Ngrid, arr))

# --- G_{n_B} for n_B = -2, eq. (2.24) ---
def G_nB(n=-2.0):
    g = math.gamma((n + 3)/2)**2
    def inner(x):
        f = lambda y: 0.5*x**(n + 2)*(1 + x*x - 2*x*y)**(n/2 - 1)*(1 + 2*x*x + 4*y**4*x*x - 4*y*y*x*x - 4*y**3*x + y*y)
        return integrate.quad(f, -1, 1, limit=400, points=[1.0] if abs(x - 1) < 0.3 else None)[0]
    pieces = [(0, 0.5), (0.5, 0.9), (0.9, 0.99), (0.99, 1.0), (1.0, 1.01), (1.01, 1.1), (1.1, 2), (2, 10), (10, 1e3)]
    tot = sum(integrate.quad(inner, a_, b_, limit=400)[0] for a_, b_ in pieces)
    tot += integrate.quad(lambda u: inner(1/u)/u**2, 1e-6, 1e-3, limit=200)[0]   # x > 1e3
    return tot/g
G2 = G_nB(-2.0)

KD = {-2: 0.8, 2: 0.9}
def B1(Brms, nB):                                     # B_1Mpc from the rms field sqrt<B^2> (the variable the bounds use)
    n = nB; return (Brms*(0.1*KD[n])**((n + 3)/(n + 5)))**((n + 5)/2)
def lamJ(B, nB, which='m'):
    kap = {(-2, 'm'): 1.17, (-2, 'b'): 1.57, (2, 'm'): 0.593}[(nB, which)]
    return (0.1*kap*B)**(2/(nB + 5))
def Delta_pmf(k, Brms, nB, z):
    a = 1/(1 + z); grow = (xi(a)/xi(0.01))**2
    B = B1(Brms, nB); lJ = lamJ(B, nB)
    if nB == -2:
        lin = 0.918e-4*xi(0.01)**2*k**6*B**4*G2
        return lin*(1 + (k*lJ)**1.94)**(-6/1.94)*grow
    return 0.244*(k*lJ)**7/(1 + (k*lJ)**2.36)**(7/2.36)*grow
plateau_b100 = 0.918e-4*xi(0.01, 'b')**2*G2/(0.1*1.57)**4      # baryon plateau at z = 100 (B-independent)
plateau_m100 = 0.918e-4*xi(0.01, 'm')**2*G2/(0.1*1.17)**4

# --- halo abundance ---
kk = P.kk; Pcdm0 = P.Pk                              # z = 0, Mpc^3
def Delta_cdm(k, z): return k**3*np.interp(k, kk, Pcdm0)/(2*math.pi**2)*D(z)**2
c_sk = 2.5
def sig2_sharpk(R, Dk):                              # Dk on kk grid: dimensionless power
    m = kk <= 1/R
    return float(np.trapezoid(Dk[m]/kk[m], kk[m]))
def dndlnM(Ms, Dk):
    R = (3*Ms/(4*math.pi*rho_m))**(1/3)/c_sk
    s2 = np.array([sig2_sharpk(r, Dk) for r in R]); s = np.sqrt(s2)
    # sharp-k: dln sigma^2 / dln R = -Delta(1/R)/sigma^2 ; dln sigma/dln M = (1/6) * that
    DkR = np.interp(1/R, kk, Dk)
    dlns_dlnM = -DkR/(6*s2)
    nu = 1.686/s; A, qq, aa = 0.3222, 0.3, 0.707
    f = A*math.sqrt(2*aa/math.pi)*(1 + (aa*nu*nu)**-qq)*nu*np.exp(-aa*nu*nu/2)
    return rho_m/Ms*f*np.abs(dlns_dlnM)
lnM = np.linspace(math.log(1e6), math.log(1e14), 400); Ms = np.exp(lnM)
def ncum(dn):                                         # n(>M)
    c = np.zeros_like(dn); c[:-1] = np.cumsum((dn[:-1][::-1] + dn[1:][::-1])/2*np.diff(lnM)[::-1])[::-1]; return c

out = dict(xi=dict(xb_001=round(xi(0.01, 'b'), 3), xm_001=round(xi(0.01, 'm'), 3), xb_rec10=round(xi(10*arec, 'b'), 3),
                   xm_z12_over_001=round(xi(1/13)/xi(0.01), 2), D12_over_D100=round(D(12)/D(100), 2)),
           G_minus2=round(G2, 4), plateau_b_z100=round(plateau_b100, 3), plateau_m_z100=round(plateau_m100, 4))
print(out)
Bs = {-2: [0.05, 0.1, 0.18, 0.27, 0.5, 0.87, 1.5, 3.0], 2: [0.03, 0.06, 0.1, 0.18, 0.5, 0.85, 1.5, 3.0]}
Zs = [8, 10, 12, 14, 17]
Mq = [1e8, 1e9, 1e10, 10**10.5, 1e11]
res = {}; curves = {}
for z in Zs:
    Dc = Delta_cdm(kk, z); dn0 = dndlnM(Ms, Dc); n0 = ncum(dn0)
    res[z] = dict(cdm={('%.1f' % math.log10(M)): float('%.3g' % np.interp(math.log(M), lnM, n0)) for M in Mq})
    for nB, BL in Bs.items():
        for B in BL:
            Dt = Dc + Delta_pmf(kk, B, nB, z); dn = dndlnM(Ms, Dt); n1 = ncum(dn)
            key = f'nB{nB:+d}_B{B}'
            res[z][key] = {('%.1f' % math.log10(M)): round(float(np.interp(math.log(M), lnM, n1)/np.interp(math.log(M), lnM, n0)), 3) for M in Mq}
            if z == 14 and key in ('nB-2_B0.27', 'nB-2_B0.87', 'nB-2_B1.5', 'nB+2_B0.18', 'nB+2_B0.85', 'nB+2_B1.5'):
                curves[key] = [[round(float(l/math.log(10)), 3), round(float(r), 4)] for l, r in zip(lnM[::8], (n1/np.maximum(n0, 1e-300))[::8]) if 6.5 <= l/math.log(10) <= 12]
# scales: where PMF power equals CDM power at z=12, and the Jeans mass
scales = {}
for nB, BL in Bs.items():
    for B in BL:
        lJ = lamJ(B1(B, nB), nB); kJ = 1/lJ; MJ = 4*math.pi/3*rho_m*(c_sk/kJ)**3
        Dc = Delta_cdm(kk, 12); Dp = Delta_pmf(kk, B, nB, 12); cross = kk[np.argmax(Dp > Dc)] if np.any(Dp > Dc) else None
        scales[f'nB{nB:+d}_B{B}'] = dict(B1Mpc=float('%.3g' % B1(B, nB)), lJ_kpc=round(1e3*lJ, 1), kJ=round(kJ, 2), MJ=float('%.3g' % MJ),
                                       k_equal_z12=None if cross is None else round(float(cross), 2),
                                       M_equal_z12=None if cross is None else float('%.3g' % (4*math.pi/3*rho_m*(c_sk/cross)**3)),
                                       plateau_z12=float('%.3g' % Delta_pmf(np.array([50*kJ]), B, nB, 12)[0]))
# halo mass whose cumulative abundance matches the bright z~12-14 galaxies (~1e-5 per Mpc^3)
abund = {}
for z in (12, 14):
    Dc = Delta_cdm(kk, z); n0 = ncum(dndlnM(Ms, Dc))
    for nt in (1e-5, 4e-5):
        abund[f'z{z}_n{nt:g}'] = float('%.3g' % math.exp(np.interp(-math.log(nt), -np.log(np.maximum(n0, 1e-300)), lnM)))
# B needed to raise n(>M) at z=14 by x10 and x100 for M = 1e10, 10^10.5
need = {}
for nB in (-2, 2):
    grid = np.logspace(-1.5, 1.5, 90)
    for M in (1e10, 10**10.5):
        Dc = Delta_cdm(kk, 14); n0 = np.interp(math.log(M), lnM, ncum(dndlnM(Ms, Dc)))
        bo = []
        for B in grid:
            n1 = np.interp(math.log(M), lnM, ncum(dndlnM(Ms, Dc + Delta_pmf(kk, B, nB, 14)))); bo.append(n1/n0)
        bo = np.array(bo)
        for f in (10, 100):
            i = np.argmax(bo >= f)
            need[f'nB{nB:+d}_M{math.log10(M):.1f}_x{f}'] = None if bo[i] < f else float('%.3g' % grid[i])
# power spectra for the figure (z = 12)
kfig = np.logspace(-1, 2.3, 120)
spec = dict(k=[round(float(x), 4) for x in kfig], cdm=[float('%.3g' % v) for v in Delta_cdm(kfig, 12)])
for key, (nB, B) in {'m2_0.27': (-2, 0.27), 'm2_0.87': (-2, 0.87), 'p2_0.18': (2, 0.18), 'p2_0.85': (2, 0.85)}.items():
    spec[key] = [float('%.3g' % v) for v in Delta_pmf(kfig, B, nB, 12)]
out.update(boost=res, scales=scales, abundance_match=abund, B_needed=need, curves=curves, spec=spec)
json.dump(out, open(os.path.join(HERE, 'pmf_halos.json'), 'w'), indent=1)
print(json.dumps(scales, indent=0)); print(abund); print(need)
for z in Zs: print(z, res[z]['cdm'], {k: v for k, v in res[z].items() if k in ('nB-2_B0.27', 'nB-2_B0.87', 'nB+2_B0.18', 'nB+2_B0.85')})

# --- cosmic dawn: collapsed fraction in minihalos (H2 cooling, >1e6) and atomic-cooling halos ---
def fcoll(Dk, Mmin):
    dn = dndlnM(Ms, Dk); m = Ms >= Mmin
    return float(np.trapezoid((Ms*dn)[m], lnM[m])/rho_m)
def Matom(z): return 1e8*((1 + z)/10)**-1.5
zz = np.arange(10, 40.1, 0.5)
dawn = {}
cases = {'cdm': None, 'nB-2_B0.1': (-2, 0.1), 'nB-2_B0.18': (-2, 0.18), 'nB-2_B0.27': (-2, 0.27), 'nB+2_B0.1': (2, 0.1), 'nB+2_B0.18': (2, 0.18)}
for key, c in cases.items():
    rows = []
    for z in zz:
        Dk = Delta_cdm(kk, z) + (0 if c is None else Delta_pmf(kk, c[1], c[0], z))
        rows.append([float(z), fcoll(Dk, 1e6), fcoll(Dk, Matom(z))])
    dawn[key] = rows
def z_at_f(rows, col, f):                             # highest z at which the collapsed fraction reaches f
    for z, *fs in rows[::-1]:
        if fs[col - 1] >= f: return z
    return None
ref = {}
for col, lab in ((1, 'mini'), (2, 'atom')):
    for zr in (15, 20):
        f0 = [r for r in dawn['cdm'] if r[0] == zr][0][col]
        ref[f'{lab}_zref{zr}'] = dict(f_cdm=float('%.3g' % f0), **{k: z_at_f(v, col, f0) for k, v in dawn.items()})
out['dawn'] = {k: [[r[0], float('%.3g' % r[1]), float('%.3g' % r[2])] for r in v] for k, v in dawn.items()}
out['dawn_shift'] = ref
json.dump(out, open(os.path.join(HERE, 'pmf_halos.json'), 'w'), indent=1)
print(json.dumps(ref, indent=0))
for k, v in dawn.items(): print(k, [(r[0], '%.2g' % r[1], '%.2g' % r[2]) for r in v if r[0] in (12, 15, 20, 25, 30)])
