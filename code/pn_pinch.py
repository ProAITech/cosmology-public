"""pn_pinch.py - Workbench Entry 136
Hourglass nebulae as z-pinch discharges, or shaped outflows? NGC 6302, the best-measured case.
Inputs: Meaburn et al. (2008) and Szyszka et al. (2011): proper motions grow linearly with distance from the star,
  0.445 mas/yr per arcsec (+0.58 mas/yr), distance 1.17 kpc, knots to >= 600 km/s, age 2,200-2,250 yr.
  Wright et al. (2011): lobes n_e ~ 2,000 cm^-3, T_e ~ 19,400 K; ionized mass 1.82 Msun; central star 220,000 K, ~14,300 Lsun;
  the star-powered photoionization model gives H-beta 852 against 776 observed (1e-12 erg/cm2/s).
  Greaves (2002), Sabin et al. (2007): toroidal fields in 3 of 4 nebulae; OH Zeeman fields of a few mG at 3e15-2e16 cm.
Bennett pinch: B^2/2mu0 at the column edge balances the column's gas pressure, I = 2 pi a B / mu0.
"""
import json, math
k, mp, mu0 = 1.380649e-23, 1.6726e-27, 4e-7*math.pi
pc, AU = 3.0857e16, 1.496e11
d_kpc = 1.17
def v_kms(r_arcsec): return 4.74*((0.445*r_arcsec + 0.58)/1000)*d_kpc*1000      # km/s from proper motion
ne, Te = 2000e6, 19400.0
P_th = 2*ne*k*Te                                     # electrons + ions (He neglected in the count)
B_pinch = math.sqrt(2*mu0*P_th)
a = 0.5*(1/60)*math.pi/180*d_kpc*1e3*pc               # half of a 1-arcmin lobe width
I_pinch = 2*math.pi*a*B_pinch/mu0
rho = ne*1.4*mp
cs = math.sqrt(2*k*Te/(0.6*mp))                       # isothermal sound speed of the ionized gas
vA = B_pinch/math.sqrt(mu0*rho)
res = dict(v_at={str(r): round(v_kms(r), 1) for r in (30, 60, 102.6, 150, 200, 240)},
           age_yr=round(1000/0.445),                                   # r / (dr/dt) for a linear flow
           P_th=P_th, B_pinch_mG=B_pinch*1e7, a_pc=a/pc, I_pinch=I_pinch, cs_kms=cs/1e3, vA_kms=vA/1e3,
           P_ram={str(v): rho*(v*1e3)**2 for v in (100, 263, 600)},
           ram_over_pinch={str(v): rho*(v*1e3)**2/P_th for v in (100, 263, 600)},
           KE_lobes_J=0.5*1.82*1.989e30*(300e3)**2, KE_over_age_W=0.5*1.82*1.989e30*(300e3)**2/(2250*3.156e7),
           Lstar_W=14300*3.828e26,
           I_near_star={f'{B}mG_at_{r:g}cm': 2*math.pi*(r/100)*(B*1e-7)/mu0 for B, r in ((1, 3e15), (3, 1e16), (3, 2e16))})
res['KE_power_over_Lstar'] = res['KE_over_age_W']/res['Lstar_W']
json.dump(res, open('pn_pinch.json', 'w'), indent=1)
print(json.dumps(res, indent=1))
