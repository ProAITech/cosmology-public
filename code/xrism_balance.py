"""Workbench Entry 118: is the hot plasma of galaxy clusters in balance? XRISM's measured gas motions against simulations and mass bias.

Cluster masses from X-rays assume the hot plasma is held up by its own thermal pressure (hydrostatic equilibrium). Turbulence, bulk flows,
magnetic fields and cosmic rays add pressure the X-ray temperature does not see, so hydrostatic masses come out low by about the non-thermal
share. XRISM/Resolve now measures gas motions directly from line widths and shifts. Ota et al. (2026) compile 45 measurements in 19 clusters
with the non-thermal pressure fraction alpha = P_nt/P_tot for each region; 44 are transcribed here (Perseus Reg3 was cut off in the source).
Compared with: the non-radiative simulations of Nelson et al. (2014), fit valid for 0.1 <= r/r200m <= 1.5 and used inward of that as a floor;
the X-COP measurement at R500 (Eckert et al. 2019); and the hydrostatic bias Planck's cluster counts would need, 1-b = 0.58 +- 0.04."""
import os, json
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
# cluster, region, class (CC/NCC), centre/offset, l_eff kpc, sigma_v km/s, v_bulk km/s, kT keV, alpha %
ROWS="""Virgo Reg1 CC C 3 262 39 1.63 21.3
Virgo Reg2 CC C 7 64 -7 2.08 1.2
Virgo Reg3 CC C 13 62 -33 2.26 1.2
Virgo E CC C 11 81 -26 2.00 2.1
Virgo SW CC C 11 91 -65 2.14 2.8
Centaurus Center CC C 27 117 -128 2.11 5.5
Perseus Reg1 CC C 21 172 152 3.62 6.2
Perseus Reg2 CC C 27 135 53 3.98 3.0
Perseus Reg4 CC C 56 186 -147 5.75 4.4
Perseus Reg5 CC O 75 115 -33 6.22 1.4
Perseus Reg6 CC O 103 188 197 7.82 3.8
Perseus E+NE CC O 163 349 282 5.47 14.7
Perseus N CC O 114 150 22 6.84 2.1
Perseus M3 CC O 54 171 -90 6.18 3.2
Perseus O3 CC O 85 275 95 6.32 7.3
Coma Center NCC C 155 208 -450 8.55 7.6
Coma South NCC O 171 202 -730 7.44 15.8
Coma North NCC O 172 167 -199 8.14 3.1
Ophiuchus Inner CC C 56 115 8 5.80 1.4
Ophiuchus Outer CC O 66 186 -104 8.40 2.8
A496 Center CC C 30 78 -69 3.23 1.5
A3571 R1 CC C 44 118 19 6.62 1.3
A3571 R2 CC O 75 103 -11 6.61 1.0
A3395S Center NCC C 71 124 264 4.96 4.7
A754 Main NCC C 107 279 220 8.10 6.9
A754 Sub NCC O 186 195 -436 9.48 6.4
Hydra-A Center CC C 103 160 -37 3.34 4.7
A2319 Full NCC C 87 232 -76 9.48 3.6
A3667 Center NCC C 154 180 -198 5.94 4.6
A3667 CF_IN NCC O 192 292 125 4.89 10.5
CygnusA Core CC C 28 261 120 5.53 7.7
A1795 R0-0.75 CC C 28 114 22 3.92 2.1
A1795 R0.75-1.5 CC C 56 115 7 6.28 1.3
A1795 R1.5-3 CC O 123 86 -32 5.53 0.9
A1795 R3-4.5 CC O 209 68 47 4.23 0.8
A2029 Center CC C 54 148 13 6.62 2.1
A2029 N1 CC O 129 58 -218 7.49 1.6
A2029 N2 CC O 259 94 -93 8.26 0.9
PKS0745-19 Center CC C 28 78 33 4.55 0.9
PKS0745-19 Outer CC O 68 110 68 7.85 1.1
A2034 Full NCC C 117 470 70 8.10 14.9
A1914 Blue NCC C 93 190 416 9.50 5.9
A1914 Red NCC C 93 245 641 9.80 11.3
H1821+643 Core CC C 44 281 51 6.00 7.8"""
R=[l.split() for l in ROWS.split("\n")]
name=[r[0]+" "+r[1] for r in R]; cl=np.array([r[2] for r in R]); leff=np.array([float(r[4]) for r in R]); sv=np.array([float(r[5]) for r in R])
vb=np.array([float(r[6]) for r in R]); kT=np.array([float(r[7]) for r in R]); al=np.array([float(r[8]) for r in R])/100
out={"n":len(R)}
# check: alpha from sigma_v alone (turbulent part), alpha_turb = 3 sigma^2 / (3 sigma^2 + 3 kT/(mu m_p)) with mu = 0.6
mp=1.6726e-27; keV=1.602e-16; cs2=kT*keV/(0.6*mp)/1e6                 # (km/s)^2
aturb=sv**2/(sv**2+cs2)
print(f"{len(R)} regions in {len(set(r[0] for r in R))} clusters. Check: turbulent share from line width alone, median {np.median(aturb)*100:.1f}%; the paper's total, median {np.median(al)*100:.1f}%")
for lab,m in (("all",np.ones(len(R),bool)),("cool-core",cl=="CC"),("non-cool-core",cl=="NCC")):
    q=np.percentile(al[m],[25,50,75]); out[lab]=dict(n=int(m.sum()),q25=float(q[0]),median=float(q[1]),q75=float(q[2]),mean=float(al[m].mean()),max=float(al[m].max()))
    print(f"   {lab:14s} n = {m.sum():2d}: alpha median {q[1]*100:.1f}% (quartiles {q[0]*100:.1f}-{q[2]*100:.1f}%), mean {al[m].mean()*100:.1f}%, max {al[m].max()*100:.1f}%")
# radius in units of r500 from the temperature (r500 = 1.10 Mpc (kT/5 keV)^0.5, a standard M-T scaling), and of r200m (r500c ~ 0.37 r200m)
Tcl={}
for r in R:
    Tcl.setdefault(r[0],[]).append(float(r[7]))
r500=np.array([1100*(max(Tcl[r[0]])/5.0)**0.5 for r in R])                  # kpc, from each cluster's hottest region
x500=leff/r500; x200m=0.37*x500
out["x500_range"]=[float(x500.min()),float(x500.max())]
print(f"   radii: {x500.min():.3f} to {x500.max():.2f} of r500 ({x200m.min():.4f} to {x200m.max():.3f} of r200m)")
nel=lambda x: 1-0.452*(1+np.exp(-(x/0.841)**1.628))
an=nel(np.maximum(x200m,0.0)); out["nelson_at_regions_median"]=float(np.median(an)); out["ratio_median"]=float(np.median(al/an))
print(f"   Nelson et al. (2014) non-radiative simulations at these radii (extrapolated inward): median {np.median(an)*100:.1f}%; measured/simulated median {np.median(al/an):.2f}")
print(f"   the simulations' own floor at the centre, 1 - 2A = {(1-2*0.452)*100:.1f}%; at r500 {nel(0.37)*100:.0f}%, at r200m {nel(1.0)*100:.0f}%")
out["nelson_r500"]=float(nel(0.37)); out["nelson_floor"]=1-2*0.452
# trend with radius: fit log alpha against log r/r500 (all, and cool cores)
for lab,m in (("all",np.ones(len(R),bool)),("CC",cl=="CC")):
    b=np.polyfit(np.log10(x500[m]),np.log10(al[m]),1); out[f"slope_{lab}"]=float(b[0]); out[f"alpha_at_r500_{lab}"]=float(10**np.polyval(b,0.0))
    print(f"   power-law trend ({lab}): alpha ∝ (r/r500)^{b[0]:.2f}; extrapolated to r500: {10**np.polyval(b,0.0)*100:.0f}% (an extrapolation over a decade)")
# hydrostatic bias implications
need=0.58; xcop=(0.059,0.031,0.94,0.04)
print(f"\nHydrostatic bias: constant alpha gives M_HSE/M = 1 - alpha. XRISM cores: about {1-np.median(al):.2f}; X-COP at R500: 0.94 ± 0.04 (alpha {xcop[0]*100:.1f}%);")
print(f"   Planck's cluster counts with Planck's CMB cosmology need 1-b = {need} ± 0.04, i.e. alpha ≈ {1-need:.2f} at R500; the X-COP authors reject it at 4.4σ")
out["xcop"]=dict(alpha=0.059,mhse=0.94,err=0.04); out["need"]=need
# magnetic pressure in Coma's centre as an example: B0 = 4.7 microG, n_e0 = 3.44e-3 cm^-3, kT = 8.25 keV
B0=4.7e-6; ne=3.44e-3; kTc=8.25
Pth=1.92*ne*kTc*1.602e-9; PB=B0**2/(8*np.pi); out["coma_magnetic_share"]=float(PB/(Pth+PB))
print(f"   magnetic pressure in Coma's centre (4.7 µG): {PB/(Pth+PB)*100:.1f}% of the total")
out["regions"]=[[float(a),float(b),c] for a,b,c in zip(x500,al,cl)]
xg=np.geomspace(0.02,1.5,80); out["nelson_curve"]=[[float(x),float(nel(0.37*x))] for x in xg]
json.dump(out,open(os.path.join(HERE,"xrism_balance.json"),"w"),indent=1)
