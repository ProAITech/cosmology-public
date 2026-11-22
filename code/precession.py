"""Precession: the standard torque account versus motion around a current filament or a Sirius companion (Workbench Entry 16)."""
import numpy as np
AS=np.degrees(1)*3600; yr=365.25*86400; G=6.674e-11; Msun=1.989e30; AU=1.496e11; ly=9.461e15; c=2.998e8
out={}
# ---------- 1. Earth's precession from satellite J2 and the seismological moment of inertia (no astronomy of the pole used)
J2=1.08263e-3            # from satellite orbits
C_MR2=0.3308             # PREM seismological Earth model (independent of precession)
H=J2/C_MR2
eps=np.radians(23.4393); w=2*np.pi/86164.0905
n_s=2*np.pi/(365.25636*86400); e_s=0.016709
n_m=2*np.pi/(27.321661*86400); e_m=0.0549; i_m=np.radians(5.145); mu=1/82.3006   # M_moon/(M_earth+M_moon)
sun=1.5*H*n_s**2/w*np.cos(eps)*(1-e_s**2)**-1.5
moon=1.5*H*n_m**2/w*np.cos(eps)*mu*(1-e_m**2)**-1.5*(1-1.5*np.sin(i_m)**2)
tot=(sun+moon)*yr*AS
print(f'Earth: H = J2/(C/MR2) = {H:.6f} (1/{1/H:.1f});  Sun {sun*yr*AS:.2f} + Moon {moon*yr*AS:.2f} = {tot:.2f} arcsec/yr;  period {360*3600/tot:,.0f} yr;  observed lunisolar 50.38, general 50.29')
out['earth']=dict(H=H,sun=sun*yr*AS,moon=moon*yr*AS,total=tot,period=360*3600/tot)
# torque needed, and what a magnetic torque could supply
Cabs=C_MR2*5.972e24*6.371e6**2; L=Cabs*w; Om=np.radians(50.38/3600)/yr
tau=L*Om*np.sin(eps); tau_mag=8.0e22*5e-9
print(f'Torque needed: {tau:.2e} N m; Earth dipole (8e22 A m^2) in a 5 nT field: {tau_mag:.1e} N m -> short by {tau/tau_mag:.1e}')
out['torque']=dict(needed=tau,magnetic=tau_mag,ratio=tau/tau_mag)
# ---------- 2. Mars: same formula, different planet
J2m=1.9566e-3; Cm=0.3639; Hm=J2m/Cm; epsm=np.radians(25.19); wm=2*np.pi/88642.66; nm=2*np.pi/(686.98*86400); em=0.0934
mars=1.5*Hm*nm**2/wm*np.cos(epsm)*(1-em**2)**-1.5*yr*AS
print(f'Mars: predicted {mars:.2f} arcsec/yr from solar torque (H = {Hm:.5f}); measured 7.608 +- 0.002; if precession were a motion of the whole solar system it would be 50.3 like Earth')
out['mars']=dict(pred=mars,meas=7.6076,H=Hm)
# ---------- 3. A Sun-Sirius orbit with a 25,772-year period
P=25772; d=8.60*ly; Mtot=1+2.06+1.02
a_kep=(Mtot*P**2)**(1/3); P_grav=np.sqrt((d/AU)**3/Mtot)
r=d/2; v=2*np.pi*r/(P*yr); acc=v**2/r
print(f'Gravity: a {P}-yr binary of {Mtot:.2f} Msun has separation {a_kep:,.0f} AU = {a_kep*AU/ly:.3f} ly; at Sirius\'s real distance the period is {P_grav/1e6:.0f} Myr')
print(f'Kinematics: opposite sides of a circle {d/ly:.1f} ly across in {P} yr -> speed {v/1e3:.0f} km/s each, acceleration {acc:.2e} m/s^2 ({acc/2.2e-10:,.0f}x the Sun\'s Galactic acceleration)')
pm=v/1e3/(4.74*10)
print(f'  nearby stars not sharing the motion would streak past at ~{pm:.1f} arcsec/yr at 10 pc; typical measured proper motions there are 0.1-1 arcsec/yr (speeds ~20-40 km/s)')
# Sirius's measured motion relative to the Sun
pm_s=np.hypot(0.54601,1.22307); vt=4.74*pm_s*2.637; vr=-5.5
print(f'  Sirius measured: proper motion {pm_s:.3f} arcsec/yr -> {vt:.1f} km/s across the sky, {vr} km/s toward us; total {np.hypot(vt,vr):.1f} km/s')
# Lorentz force that would be needed
F=Msun*acc; Bism=3e-10; q=F/(v*Bism); E_surf=8.988e9*q/6.957e8**2
print(f'Force on the Sun {F:.1e} N; with a 3 microgauss field this needs a net charge {q:.1e} C (surface field {E_surf:.1e} V/m); estimates of the Sun\'s real net charge are ~1e2 C')
out['sirius']=dict(a_kep_AU=a_kep,P_grav_Myr=P_grav/1e6,v_kms=v/1e3,acc=acc,pm10=pm,sirius_v=float(np.hypot(vt,vr)),q=q)
# ---------- 4. Does Sirius precess? ecliptic longitude now, and what 50.29"/yr implies for 138 AD
def ecl(ra,dec,e=np.radians(23.4393)):
    ra,dec=np.radians(ra),np.radians(dec); lam=np.arctan2(np.sin(ra)*np.cos(e)+np.tan(dec)*np.sin(e),np.cos(ra)); bet=np.arcsin(np.sin(dec)*np.cos(e)-np.cos(dec)*np.sin(e)*np.sin(ra)); return np.degrees(lam)%360,np.degrees(bet)
lam,bet=ecl(101.2872,-16.7161)
back=lam-(2000-138)*50.29/3600
print(f'Sirius ecliptic longitude J2000: {lam:.2f} deg, latitude {bet:.2f}; precessed back to 138 AD: {back:.2f} deg = Gemini {back-60:.1f}; Ptolemy lists Gemini 17 2/3 (77.67), latitude -39 1/6')
# B1950 -> J2000 change in right ascension expected from precession alone
m,n=3.075,1.336; ra=np.radians(101.0); dec=np.radians(-16.68)
dra=(m+n*np.sin(ra)*np.tan(dec))*50
print(f'  catalogue right ascension B1950 06h42m56.7s -> J2000 06h45m08.9s = +132.2 s; precession formula gives {dra:.1f} s, proper motion {-0.0379*50:.1f} s')
out['sirius_precess']=dict(lam2000=lam,beta=bet,lam138=back,dra=dra)
# ---------- 5. The Sun's apogee against the equinox (both inside the solar system)
h,b,now=65.5,82.2833,103.0   # Hipparchus c. -140; al-Battani c. 883; today
r1=(b-h)*3600/(883+140); r2=(now-b)*3600/(2026-883); r3=(now-h)*3600/(2026+140)
print(f'Solar apogee: Hipparchus 65.5 -> al-Battani 82.28 -> today 103.0 deg: {r1:.1f}, {r2:.1f}, overall {r3:.1f} arcsec/yr; standard = 50.3 (axis) + 11.6 (orbit) = 61.9; if the axis did not turn relative to the orbit: 11.6')
out['apogee']=dict(rates=[r1,r2,r3],std=61.9,alt=11.6)
# ---------- 6. Nutation: the same constant predicts the 18.6-yr nod tied to the Moon's node
nod=18.6*yr
print(f'Nutation: lunar node period 18.61 yr; observed nod of the pole 9.20 arcsec at exactly that period (the principal term of 1,365 in the standard model)')
import json; json.dump(out,open('precession.json','w'),indent=1,default=float)

# ====================== Part B: is electromagnetism the dominant force? ======================
print('\n--- force budget ---')
k=8.988e9; e=1.602e-19; mp=1.6726e-27; me=9.109e-31; mu0=4e-7*np.pi; kB=1.381e-23; kpc=3.086e19
print(f'electric/gravity for a proton and an electron: {k*e*e/(G*mp*me):.1e}; for two protons: {k*e*e/(G*mp*mp):.1e}')
qm=np.sqrt(G/k); frac=qm/(e/mp)
Qsun=Msun/mp*e
print(f'charge imbalance at which electric force equals gravity: {qm:.1e} C/kg = one unbalanced charge per {1/frac:.1e}')
print(f'Sun: {Qsun:.1e} C of each sign; estimated net charge ~77 C = one part in {Qsun/77:.0e}; that is {qm*Msun/77:.0e} times below the level that would rival gravity')
deb=lambda T,n: 69*np.sqrt(T/n)
print(f'Debye length (scale beyond which plasma is neutral): warm interstellar gas {deb(8000,3e4):.0f} m; cosmic-web filament {deb(1e6,10)/1e3:.0f} km; structures are 1e19-1e23 m')
rows=[('Milky Way disk gas',6e-10,1e6,220e3),('Milky Way halo gas',0.73e-10,1e3,200e3),('Cluster gas',5e-10,1e3,1000e3),('Cosmic-web filament',13e-13,10,300e3)]
out['energy']=[]
for name,B,n,v in rows:
    um=B*B/(2*mu0); uk=0.5*n*mp*v*v; print(f'  {name:22s} B={B*1e10:8.4f} uG  magnetic/kinetic energy = {um/uk:.1e}'); out['energy'].append((name,um/uk))
Fm=(6e-10)**2/(2*mu0)*np.pi*6.957e8**2; Fg=Msun*(220e3)**2/(8.2*kpc)
print(f'Sun in the Galaxy: magnetic pressure of 6 uG on its cross-section {Fm:.1e} N vs gravity holding its orbit {Fg:.1e} N -> {Fg/Fm:.0e} times larger')
print(f'Milky Way baryons: stars ~5.5e10 Msun, gas ~1.0e10 Msun -> {5.5/6.5:.0%} of the ordinary mass is in stars')
out['force']=dict(ratio_pe=k*e*e/(G*mp*me),imbalance=frac,sun_frac=77/Qsun,sun_force_ratio=Fg/Fm)
json.dump(out,open('precession.json','w'),indent=1,default=float)
