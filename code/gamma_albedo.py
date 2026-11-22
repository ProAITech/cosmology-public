"""Workbench Entry 25, Part 3: do stars hide the gamma rays?
Making deuterium with cosmic rays also makes neutral pions, which decay to gamma rays. The model's own paper puts these at about 250 times
the observed gamma-ray background, and answers that the cosmic rays hit stellar surfaces, which swallow most of the gamma rays.
There is one star where this is measured: the Sun glows in gamma rays because cosmic rays hit its surface (Fermi-LAT, Abdo et al. 2011).
How much of the cosmic-ray power arriving at the Sun comes back out as gamma rays, and what does that efficiency give for the model?"""
import numpy as np, json
eV=1.602176634e-12; GeV=1.602176634e-3; c=2.99792458e10; AU=1.495979e13; Rsun=6.957e10
F=4.6e-7; Flo,Fhi=4.6e-7-0.8e-7,4.6e-7+1.0e-7           # solar disk photons /cm^2/s above 100 MeV
def meanE(index,E0=0.1,E1=10.0):                         # mean photon energy (GeV) of a power law between E0 and E1
    E=np.logspace(np.log10(E0),np.log10(E1),4000); w=E**-index*np.gradient(E); return float((w*E).sum()/w.sum())
out={}
print("Sun: gamma-ray luminosity above 100 MeV")
rows=[]
for idx in (2.1,2.3):
    L=4*np.pi*AU**2*F*meanE(idx)*GeV; rows.append(L); print(f"  photon index {idx}: mean energy {meanE(idx):.2f} GeV, L = {L:.2e} erg/s")
Lg=np.mean(rows)
print("cosmic-ray power arriving at the solar surface, if the intensity there were the same as near Earth: P = pi R^2 c u")
alb=[]
for u in (0.5,0.7,1.0):                                  # cosmic-ray energy density, eV/cm^3
    P=np.pi*Rsun**2*c*u*eV; alb.append(Lg/P); print(f"  u = {u} eV/cm^3: P = {P:.2e} erg/s -> fraction returned as gamma rays {Lg/P:.2e}")
a_lo,a_mid,a_hi=min(alb),alb[1],max(alb)
# the gamma-ray background
I=7.2e-6; E=np.logspace(-1,3,6000); w=E**-2.32*np.exp(-E/279)*np.gradient(E); mE=float((w*E).sum()/w.sum())
u_bg=4*np.pi*I*mE*1e9/c
print(f"\nisotropic gamma-ray background above 100 MeV: {I:.1e} /cm^2/s/sr, mean energy {mE:.2f} GeV -> {u_bg:.2e} eV/cm^3")
# the model: cosmic rays carry 0.5% (+-0.2%) of the thermonuclear energy; its paper quotes 2e-15 erg/cm^3
u_cr=2e-15/eV; print(f"model: cosmic-ray energy released per volume 2e-15 erg/cm^3 = {u_cr:.2e} eV/cm^3")
for fg,lab in ((1/3,"one third of the energy to gamma rays (the paper's implied figure)"),(0.10,"one tenth (GeV protons)")):
    print(f"  thin gas, {lab}: {fg*u_cr:.1e} eV/cm^3 = {fg*u_cr/u_bg:.0f} x the background; suppression available from a Sun-like surface: {fg/a_hi:.0f} to {fg/a_lo:.0f}")
print(f"  stopped in Sun-like surfaces: {a_lo*u_cr:.1e} to {a_hi*u_cr:.1e} eV/cm^3 = {a_lo*u_cr/u_bg:.1f} to {a_hi*u_cr/u_bg:.1f} x the background (central {a_mid*u_cr/u_bg:.1f})")
print(f"  with the model's upper cosmic-ray share (0.7% and the standard gas density): up to {a_hi*u_cr*(0.7/0.5)*(0.409/0.259)/u_bg:.1f} x")
out={"L_gamma":Lg,"albedo":[a_lo,a_mid,a_hi],"u_bg":u_bg,"u_cr":u_cr,"ratio":[a_lo*u_cr/u_bg,a_mid*u_cr/u_bg,a_hi*u_cr/u_bg],"ratio_max":a_hi*u_cr*(0.7/0.5)*(0.409/0.259)/u_bg,"thin":[0.10*u_cr/u_bg,u_cr/3/u_bg]}
# depth at which deuterium is made, against the depth a wind removes
lam=1/( (0.75/1.6726e-24)*30e-27 + (0.25/6.646e-24)*100e-27 ); print(f"\nmean depth of a proton's first nuclear collision in hydrogen-helium gas: {lam:.0f} g/cm^2")
for R,name,md in ((3,"a 5 Msun star, wind 1e-9 Msun/yr",1e-9),(10,"a 20 Msun star, wind 1e-7 Msun/yr",1e-7)):
    m=4*np.pi*(R*Rsun)**2*300/1.989e33; print(f"  {name}: the outer 300 g/cm^2 weighs {m:.1e} Msun and is blown off in {m/md:.0f} yr")
out["lambda"]=lam
json.dump(out,open("gamma_albedo.json","w"),indent=1)
