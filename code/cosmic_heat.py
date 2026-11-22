"""Workbench Entry 121: what keeps the universe's plasma hot?

Chiang et al. (2020) measured the thermal pressure of all the ionized gas in the universe, slice by slice in redshift, by cross-correlating
Planck's Compton-y maps with galaxy surveys. The density-weighted mean electron temperature rises from about 7 x 10^5 K at z = 1 to
2 x 10^6 K today, and the thermal energy density reaches Omega_th ~ 1.5 x 10^-8. Plasma cosmology has cosmic currents and fields; the
standard model has gravitational collapse plus feedback. Which can supply this heat?
 1. The measured temperatures (their Table 2, z < 0.7, plus their z = 1 value) and a fit of the trend; a static universe without evolution
    predicts no trend.
 2. A halo-model prediction from gravity alone: gas in Sheth-Tormen halos at its virial temperature (Planck 2018, Entry 101's spectrum).
 3. The energy budget: heat built since z = 1, against gravity, black holes, supernovae, and electric currents. For currents, the field they
    make is capped by Faraday rotation (filaments 30 nG, Entry 13/15 and Carretti et al.), which caps the current density; Ohmic heating at
    Spitzer resistivity, and the electric field J.E would need, follow."""
import os, sys, json, math, contextlib, io
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
CODE = HERE if os.path.exists(os.path.join(HERE, 'progenitors.py')) else os.path.abspath(os.path.join(HERE, '..', 'artifact-files', '4379c2ac-1195-4a41-978e-eb44b945a265', 'code'))
sys.path.insert(0, CODE)
with contextlib.redirect_stdout(io.StringIO()): import progenitors as P
# 1. measurements: z, T_e (1e6 K), error, Omega_th (1e-8)
T2 = [(0.07, 2.51, 0.68, 1.93), (0.11, 1.47, 0.62, 1.13), (0.14, 1.52, 0.52, 1.17), (0.18, 1.93, 0.50, 1.49), (0.22, 2.07, 0.43, 1.59),
      (0.27, 0.69, 0.43, 0.53), (0.31, 1.40, 0.51, 1.08), (0.36, 0.79, 0.39, 0.60), (0.40, 1.92, 0.30, 1.47), (0.45, 1.33, 0.25, 1.02),
      (0.50, 1.15, 0.26, 0.88), (0.55, 1.56, 0.24, 1.20), (0.61, 0.77, 0.26, 0.59), (0.66, 1.16, 0.32, 0.89)]
z = np.array([r[0] for r in T2]); T = np.array([r[1] for r in T2]); e = np.array([r[2] for r in T2])
w = 1/e**2
Tc = np.sum(w*T)/np.sum(w); chi_c = float(np.sum(w*(T-Tc)**2))
# T = A (1+z)^-a, weighted least squares in linear space by grid
best = None
for a in np.linspace(-2, 6, 801):
    f = (1+z)**-a; A = np.sum(w*T*f)/np.sum(w*f*f); c = float(np.sum(w*(T-A*f)**2))
    if best is None or c < best[0]: best = (c, a, A)
chi_s, a_fit, A_fit = best
aa = np.linspace(-2, 6, 801); cc = []
for a in aa:
    f = (1+z)**-a; A = np.sum(w*T*f)/np.sum(w*f*f); cc.append(float(np.sum(w*(T-A*f)**2)))
cc = np.array(cc); ok = aa[cc <= chi_s+1]; a_lo, a_hi = ok.min(), ok.max()
print(f"measured mean electron temperature, 14 slices at z = 0.07-0.66: constant fit {Tc:.2f} x 10^6 K, chi2 {chi_c:.1f} for 13; "
      f"T ∝ (1+z)^-{a_fit:.2f} (+{a_hi-a_fit:.2f}, -{a_fit-a_lo:.2f}), A = {A_fit:.2f}, chi2 {chi_s:.1f} for 12; no-trend rejected by Δχ² {chi_c-chi_s:.1f} ({math.sqrt(max(chi_c-chi_s,0)):.1f}σ)")
print(f"   with the paper's z = 1 value of 0.7 x 10^6 K the fit predicts {A_fit*2**-a_fit:.2f} x 10^6 K at z = 1")
# 2. halo model from gravity alone
G = 4.30091e-9                       # Mpc (km/s)^2 / Msun
KB = 1.380649e-23; MP = 1.67262e-27; MU = 0.59
h = P.h; Om = P.Om; Ob = P.Ob; fb = P.fb
rho_c0 = 2.775e11*h*h               # Msun / Mpc^3
def Ez(zz): return math.sqrt(Om*(1+zz)**3 + 1 - Om)
def dndlnM(M, zz):
    s = np.array([P.sigM(m) for m in M])*P.D(zz)
    nu = 1.686/s; a, p, A = 0.707, 0.3, 0.3222; an = math.sqrt(a)*nu
    f = A*math.sqrt(2/math.pi)*an*(1 + an**(-2*p))*np.exp(-an**2/2)
    return P.rho_m/M*f*np.abs(np.gradient(np.log(s), np.log(M)))
def Tvir(M, zz, beta=1.1):           # K; calibrated: M200 = 5.9e14 Msun at z = 0 -> 4.9 keV (Vikhlinin et al. 2009)
    R = (3*M/(4*math.pi*200*rho_c0*Ez(zz)**2))**(1/3)
    v2 = G*M/(2*R)*1e6                # m^2/s^2
    return beta*MU*MP*v2/KB
def fgas(M, mode):
    if mode == 'all': return fb*np.ones_like(M)
    return fb*0.9/(1 + (10**13.5/M)**0.7)      # groups lose gas to feedback
Ms = np.geomspace(1e11, 3e15, 240)
print(f"   check: T at M200 = 5.9e14, z = 0: {Tvir(np.array([5.9e14]),0)[0]*KB/1.602e-16:.2f} keV")
zs = [0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]; HM = {}
for mode in ('all', 'depleted'):
    for mmin in (1e12, 1e13):
        row = []
        for zz in zs:
            M = Ms[Ms >= mmin]; dn = dndlnM(M, zz)
            Tm = np.trapezoid(dn*fgas(M, mode)*M*Tvir(M, zz), np.log(M))/(P.rho_m*fb)
            row.append(Tm/1e6)
        HM[f'{mode}|{mmin:.0e}'] = row
        print(f"   halo model, gas {mode:8s}, halos above {mmin:.0e}: mean T (10^6 K) at z = " + ", ".join(f"{q}: {v:.2f}" for q, v in zip(zs, row)) + f"; z=0/z=1 = {row[0]/row[4]:.1f}")
# 3. energy budget (per comoving m^3)
c = 2.99792458e8; MPC = 3.0857e22; MSUN = 1.98847e30; GYR = 3.15576e16
rho_c_SI = 3*(67.4e3/MPC)**2/(8*math.pi*6.674e-11)                  # kg/m^3
u_c = rho_c_SI*c*c
u_th = 1.5e-8*u_c
t01 = (P.age(0) - P.age(1))*GYR
heat = 0.7*u_th/t01
print(f"   thermal energy today: Omega_th 1.5e-8 = {u_th:.2e} J/m^3; 70% built over {t01/GYR:.1f} Gyr since z = 1: {heat:.2e} W/m^3 on average")
budget = {'thermal energy of the cosmic plasma, measured': u_th}
# gravity: the virial expectation, thermal = -(2/3) x gas share of halo binding energy, i.e. what our halo model gives
Tm0 = HM['depleted|1e+12'][0]*1e6
n_tot = (Ob*rho_c_SI/(1.14*MP))*(1.14/MU)
budget['gravity: virial heat of gas in halos (our halo model)'] = u_th*Tm0/2.0e6         # scaled by temperature, on the paper's own definition
u_grav = budget['gravity: virial heat of gas in halos (our halo model)']
# black holes: local mass density 4.6e5 Msun/Mpc^3 (Marconi et al. 2004), radiative efficiency 0.1, 5% coupled
u_bh = 4.6e5*MSUN*c*c*0.1/MPC**3
budget['black holes: all radiated energy'] = u_bh; budget['black holes: 5% coupled to gas'] = 0.05*u_bh
# supernovae: stellar mass today 5.8e8 Msun/Mpc^3, return fraction 0.27, 0.0068 core-collapse SNe per Msun formed, 1e44 J each
u_sn = 5.8e8/(1-0.27)*0.0068*1e44/MPC**3
budget['supernovae: all explosion energy'] = u_sn
# currents: field capped at 30 nG over a filament radius of 1 Mpc; filaments fill 6% of space
B = 30e-13*1e0; B = 3.0e-12; r = 1.0*MPC; fV = 0.06
J = 2*B/(4e-7*math.pi*r)
def eta_sp(TK, lnL=35): return 5.2e-5*lnL/((KB*TK/1.602e-19)**1.5)
eta = eta_sp(2e6)
p_ohm = eta*J*J*fV                                                 # W per m^3, averaged over space
budget['currents: Ohmic heat in 8 Gyr, field 30 nG, Spitzer'] = p_ohm*t01
u_B = B*B/(2*4e-7*math.pi)*fV
budget['magnetic energy in filaments, 30 nG'] = u_B
need_fil = heat/fV
J_need = math.sqrt(need_fil/eta); B_need = 4e-7*math.pi*J_need*r/2
E_need = need_fil/J; J_drive = E_need/eta
nE = 1e-4*1e6                                                       # electrons per m^3 in filaments
E_D = nE*(1.602e-19)**3*35/(4*math.pi*(8.854e-12)**2*KB*2e6)
print(f"   Spitzer resistivity at 2e6 K: {eta:.1e} Ohm m; 30 nG over 1 Mpc allows J = {J:.1e} A/m^2; Ohmic heat {p_ohm:.1e} W/m^3, {heat/p_ohm:.1e} times short")
print(f"   to heat filaments Ohmically at Spitzer resistivity needs J = {J_need:.1e} A/m^2, a field of {B_need*1e4:.2f} G (measured: 3e-8 G), {B_need/B:.1e} times too strong")
print(f"   with J capped by the field, J.E needs E = {E_need:.1e} V/m along the filaments ({E_need*10*MPC:.1e} V over 10 Mpc); that field would drive {J_drive:.1e} A/m^2, {J_drive/J:.1e} times the allowed current, and is {E_need/E_D:.1e} times the Dreicer field ({E_D:.1e} V/m), so electrons run away instead of heating")
# thin current sheets: J = B/(mu0 l) in sheets of thickness l filling l/r of the filament
l_sheet = eta*B*B/((4e-7*math.pi)**2*r*need_fil)
wpe = 56.4*math.sqrt(nE); skin = c/wpe
print(f"   thin current sheets would need thickness {l_sheet:.1e} m, {skin/l_sheet:.1e} times thinner than the electron skin depth ({skin/1e3:.0f} km)")
# the only generator named: a galaxy's rotation in its own field, 5.3e32 W for the Milky Way (Entry 107), 0.005-0.01 such galaxies per Mpc^3
for ng in (0.005, 0.01):
    pg = 5.3e32*ng/MPC**3
    print(f"   galactic generators at the Milky Way's power, {ng} per Mpc^3: {pg:.1e} W/m^3, {pg/heat*100:.1f}% of the heating")
out_extra = dict(l_sheet=l_sheet, skin=skin, gen=[5.3e32*ng/MPC**3/heat for ng in (0.005, 0.01)])
print(f"   resistivity needed with J capped: {need_fil/J**2:.1e} Ohm m, {need_fil/J**2/eta:.1e} times Spitzer")
for k, v in budget.items(): print(f"   {k:58s} {v:.2e} J/m^3  ({v/u_th:.1e} of the heat)")
slopes = {k: math.log(v[0]/v[4])/math.log(2) for k, v in HM.items()}
print("   halo-model slope a in T ∝ (1+z)^-a between z = 0 and 1: " + ", ".join(f"{k} {v:.2f}" for k, v in slopes.items()))
# cluster cores: 5 microgauss over 100 kpc, 5e7 K, against radiative cooling n^2 Lambda with n = 0.05 cm^-3, Lambda = 1e-23 erg cm^3/s
Bc = 5e-10; rc = 0.1*MPC; Jc = 2*Bc/(4e-7*math.pi*rc); etac = eta_sp(5e7)
pc = etac*Jc**2; cool = (0.05)**2*1e-23*0.1                     # erg/cm^3/s -> W/m^3 (x 0.1)
print(f"   cluster core: J = {Jc:.1e} A/m^2, Ohmic {pc:.1e} W/m^3 against cooling {cool:.1e} W/m^3: {cool/pc:.1e} times short")
out_extra.update(slopes=slopes, core_short=cool/pc)
out = dict(T2=T2, Tc=Tc, chi_c=chi_c, a_fit=float(a_fit), a_lo=float(a_lo), a_hi=float(a_hi), A_fit=float(A_fit), chi_s=chi_s, zs=zs, HM=HM,
           u_th=u_th, heat=heat, t01=t01/GYR, budget=budget, eta=eta, J=J, p_ohm=p_ohm, J_need=J_need, B_need=B_need, E_need=E_need, J_drive=J_drive, E_D=E_D,
           eta_need=need_fil/J**2, **out_extra)
json.dump(out, open(os.path.join(HERE, 'cosmic_heat.json'), 'w'), indent=1)
