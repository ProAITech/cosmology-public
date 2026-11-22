"""Workbench Entry 21: a variable speed of light as the cause of redshift (Unzicker's Dirac-Sciama-Dicke cosmology,
arXiv:0708.3518). In absolute time t the paper has: c ~ t^-1/2, atomic wavelengths (rulers) ~ t^-1/4, atomic
frequencies ~ t^-1/4 (clock periods ~ t^1/4), matter at rest in flat space. What does an observer using atomic
rulers and clocks measure, and does it pass the timing test, the temperature test and DESI's shape test?"""
import json, numpy as np
from scipy.optimize import minimize_scalar
out={}
# ---------- 1. translate to measured units
# rulers shrink as t^-p, clock periods grow as t^q:  measured distances between fixed galaxies grow as a ~ t^p,
# measured time t' = integral dt / t^q ~ t^(1-q)  ->  a ~ t'^n with n = p/(1-q)
p,q=0.25,0.25; n=p/(1-q)
print(f"Rulers ~ t^-{p}, clock periods ~ t^{q}: measured scale factor a ~ t'^n with n = {n:.4f} (= 1/3)")
# redshift: photon wavelength fixed in flight, rulers shrink: 1+z = (t0/te)^p = a0/ae
# timing: c dt is the same for two pulses on the same path -> dt_obs/dt_emit = c_e/c_o = (t0/te)^(1/2) = (1+z)^2 in absolute time;
# in local clock units divide by the clock period ratio (t0/te)^q = (1+z): measured stretch = (1+z)^(2-1) = (1+z)
print("Measured time dilation: (1+z)^(1/(2p)) / (1+z)^(q/p) = (1+z)^%.1f" % (1/(2*p)-q/p))
print("Photon frequency against local atomic frequency: (1+z)^-2 / (1+z)^-1 = (1+z)^-1 -> a blackbody stays a blackbody with T ~ (1+z)")
out["n"]=n; out["time_dilation_exponent"]=1/(2*p)-q/p

# ---------- 2. the shape test (Alcock-Paczynski): D_M/D_H from DESI DR2, no ruler and no galaxy evolution in it
data=[(0.510,0.622,0.017),(0.706,0.892,0.021),(0.934,1.223,0.019),(1.321,1.948,0.045),(1.484,2.386,0.136),(2.330,4.518,0.097)]
def F_power(z,n):                      # a ~ t^n: H = H0 (1+z)^(1/n)
    if abs(n-1)<1e-9: return (1+z)*np.log(1+z)
    return ((1+z)**(1/n)-(1+z))/(1/n-1)
def F_lcdm(z,Om=0.315,N=4000):
    E=lambda x: np.sqrt(Om*(1+x)**3+1-Om); zz=(np.arange(N)+0.5)*z/N
    return np.sum(1/E(zz))*z/N*E(z)
chi=lambda f: sum(((o-f(z))/s)**2 for z,o,s in data)
c_vsl=chi(lambda z:F_power(z,n)); c_coast=chi(lambda z:F_power(z,1.0)); c_l=chi(F_lcdm)
best=minimize_scalar(lambda m: chi(lambda z:F_power(z,m)),bounds=(0.3,3),method='bounded')
# 1-sigma range of n
grid=np.linspace(0.6,1.4,1601); cc=np.array([chi(lambda z,m=m:F_power(z,m)) for m in grid]); ok=grid[cc<best.fun+1]
print("\nShape test, chi-squared for 6 DESI points:")
print(f"  LCDM (Om = 0.315): {c_l:.1f}")
print(f"  coasting, a ~ t (n = 1): {c_coast:.1f}")
print(f"  best power law: n = {best.x:.3f} (+{ok.max()-best.x:.3f}/-{best.x-ok.min():.3f}), chi2 = {best.fun:.1f}")
print(f"  variable speed of light as published (n = 1/3): {c_vsl:.0f}")
print("  predicted D_M/D_H at each redshift (measured in brackets):")
for z,o,s in data: print(f"    z = {z:.3f}: n=1/3 gives {F_power(z,n):6.2f}; LCDM {F_lcdm(z):.3f}  [{o} +- {s}]  -> {abs(F_power(z,n)-o)/s:.0f} sigma")
o=21.863/8.632; s=o*np.hypot(0.425/21.863,0.101/8.632)
r_vsl=(3.33/1.51)**(1/n); E=lambda z,Om=0.315: np.sqrt(Om*(1+z)**3+1-Om)
print(f"  ruler-free radial test D_H(0.51)/D_H(2.33): measured {o:.2f} +- {s:.2f}; n = 1/3 gives {r_vsl:.2f} ({(r_vsl-o)/s:.0f} sigma); LCDM {E(2.33)/E(0.51):.2f}")
out["shape"]={"chi2_lcdm":c_l,"chi2_coasting":c_coast,"best_n":best.x,"best_n_range":[ok.min(),ok.max()],"chi2_best":best.fun,"chi2_vsl":c_vsl,
              "pred_vsl":[F_power(z,n) for z,_,_ in data],"pred_lcdm":[F_lcdm(z) for z,_,_ in data],"DH_ratio":{"obs":o,"err":s,"vsl":r_vsl,"lcdm":E(2.33)/E(0.51)}}

# ---------- 3. other consequences of H = H0 (1+z)^3
GYR_PER_H0=977.8   # 1/H0 in Gyr for H0 = 1 km/s/Mpc
for H0 in (67.4,73.0): print(f"\nAge in atomic-clock time, t0 = n/H0: {n*GYR_PER_H0/H0:.1f} Gyr for H0 = {H0}")
dL=lambda z,m: (1+z)*(1-(1+z)**(1-1/m))/(1/m-1)
def dL_l(z,Om=0.315,N=4000):
    zz=(np.arange(N)+0.5)*z/N; return (1+z)*np.sum(1/np.sqrt(Om*(1+zz)**3+1-Om))*z/N
for z in (0.5,1.0): print(f"Supernova at z = {z}: {5*np.log10(dL(z,n)/dL_l(z)):+.2f} mag against LCDM (negative = would look brighter)")
print(f"Deceleration parameter q = 1/n - 1 = {1/n-1:.0f} (LCDM today: about -0.5)")
out["age_Gyr"]={"67.4":n*GYR_PER_H0/67.4,"73":n*GYR_PER_H0/73.0}; out["sn_offset_mag"]={"0.5":5*np.log10(dL(0.5,n)/dL_l(0.5)),"1.0":5*np.log10(dL(1.0,n)/dL_l(1.0))}

# ---------- 4. what unit laws would be needed instead
print(f"\nTo mimic the best-fitting power law (n = {best.x:.2f}) with Dicke-type units one needs p/(1-q) = {best.x:.2f}, e.g. rulers ~ t^-{best.x*(1-q):.2f} with clock periods ~ t^{q}.")
print("No single power law gives the change from slowing to speeding up that the data show; LCDM fits better by delta chi2 = %.1f" % (best.fun-c_l))
# ---------- 5. cross-check in the model's own absolute frame (no translation to "expansion" at all)
# flat static space, c(t) = c0 (t0/t)^(1/2), 1+z = (t0/t)^(1/4). A pattern of size r with no preferred direction at epoch t1:
#   across the sky it subtends theta = r / D with D = integral of c dt from t1 to t0 (numerical);
#   along the line of sight its far end is seen at a time earlier by dt = r / c(t1), a redshift step dz = |dz/dt| dt (numerical derivative).
# D_M/D_H = dz / theta, and r cancels: the shape test does not care how big the pattern is or whether its size changes with time.
def F_abs(z,N=200000):
    t1=(1+z)**-4.0                                   # in units of t0, with c0 = 1
    tt=t1+(np.arange(N)+0.5)*(1-t1)/N; D=np.sum(tt**-0.5)*(1-t1)/N
    h=1e-6*t1; dzdt=((t1-h)**-0.25-(t1+h)**-0.25)/(2*h)
    return dzdt/(t1**-0.5)*D
print("\nCross-check, light paths integrated in absolute time (should equal the n = 1/3 column):")
for z,o,s_ in data: print(f"    z = {z:.3f}: {F_abs(z):6.3f}  vs  {F_power(z,n):6.3f}")
out["abs_frame_check"]=[F_abs(z) for z,_,_ in data]
# time from the beginning to z = 1100 in atomic-clock time: t'/t0' = (1+z)^-3
for H0 in (67.4,73.0): print(f"Atomic-clock time at z = 1100: {n*GYR_PER_H0/H0*1e9*1101.0**-3:.1f} years (H0 = {H0}); standard model: about 380,000 years")
out["t_cmb_years"]=n*GYR_PER_H0/67.4*1e9*1101.0**-3
# curves for the figure
zs=np.linspace(0.05,2.6,52)
out["curve"]={"z":list(zs),"vsl":[F_power(z,n) for z in zs],"lcdm":[F_lcdm(z) for z in zs],"coast":[F_power(z,1.0) for z in zs],"best":[F_power(z,best.x) for z in zs]}
json.dump(out,open("vsl_results.json","w"),indent=1,default=float)
