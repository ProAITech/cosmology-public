"""Workbench Entry 19: is there a universal direction for galactic currents?
Published counts for jets, our signed stack of background rotation measures around nearby inclined galaxies
(numbers from halo_sign_stack.js, run in a browser tab), and what a Milky-Way-like halo would show."""
import json, numpy as np
from scipy import stats
out={}
# ---------- 1. Jets, counted one by one (Contopoulos, Christodoulou, Kazanas & Gabuzda 2009)
p_core=stats.binom.sf(21,29,0.5); p_far=0.5**7
print(f"Jets: 22 of 29 gradients nearest the core are clockwise: chance of 22 or more = {100*p_core:.2f}%")
print(f"      the 7 measured further out are all counter-clockwise: chance = {100*p_far:.2f}% (one-sided)")
out["jets"]={"core":[22,29,p_core],"outer":[7,7,p_far]}

# ---------- 2. Our signed stack. Model: RRM = -A sin(2 theta); A > 0 for outward currents. [n sightlines, n galaxies, A, error]
stack={"NVSS":{"0-30 kpc":[26,22,-2.80,3.48],"30-60 kpc":[80,37,0.41,2.77],"60-100 kpc":[178,55,-1.94,2.13],"100-200 kpc":[1053,133,0.17,0.72],
               "all within 100 kpc":[284,75,-1.44,1.65],"all within 200 kpc":[1337,143,-0.18,0.63]},
       "LoTSS":{"30-60 kpc":[21,9,-0.50,3.55],"60-100 kpc":[36,15,-0.20,1.64],"100-200 kpc":[184,37,-0.28,0.68],
               "all within 100 kpc":[59,18,0.00,1.72],"all within 200 kpc":[243,42,-0.21,0.71]}}
mirror={"NVSS <100":[284,0.60,1.30],"NVSS <200":[1337,0.50,0.62],"LoTSS <100":[59,1.98,1.44],"LoTSS <200":[243,0.59,0.59]}
print("\nSigned stack (rad/m2); 95% upper limit on an outward pattern = A + 2 sigma")
for cat,rows in stack.items():
    for k,(n,g,A,e) in rows.items(): print(f"  {cat:5s} {k:20s} {n:5d} sightlines, {g:3d} galaxies: A = {A:+.2f} +- {e:.2f}   limit {A+2*e:.1f}")
print("  mirror control (position angles reflected):", {k:f"{v[1]:+.2f} +- {v[2]:.2f}" for k,v in mirror.items()})
print("  injection: adding A = 5 to the NVSS sightlines within 100 kpc returns 3.56 = 5 - 1.44")
# NVSS and LoTSS see many of the same sources, so they are not combined.
out["stack"]=stack; out["mirror"]=mirror
print("  sightlines closer than 10, 15, 20, 30 kpc: NVSS 3, 7, 13, 26; LoTSS 0, 1, 1, 2   (228 galaxies)")

# ---------- 3. What the Milky Way's rings would look like from outside, edge-on
# Ring field: 0.73 uG peak at 8 kpc, Gaussian width 5.31 kpc (Xu & Han 2024), taken as uniform for 1 < |z| < 5 kpc.
# Electrons: thick disk 0.0113 cm-3 with 1.67 kpc scale height (YMW16), out to 15 kpc.
def rm_edge(x,z,B0=0.73,R0=7.97,RT=5.31,n0=0.0113,H=1.673,Rmax=15.0):
    if not (1<abs(z)<5) or abs(x)>=Rmax: return 0.0
    s=np.linspace(-np.sqrt(Rmax**2-x**2),np.sqrt(Rmax**2-x**2),2001); R=np.hypot(x,s)
    Bpar=B0*np.exp(-(R-R0)**2/(2*RT**2))*abs(x)/R
    return 0.812*n0*np.exp(-abs(z)/H)*np.trapezoid(Bpar,s*1000)
print("\nMilky-Way-like rings seen edge-on from outside, rad/m2 at (distance along the disk, height):")
for z in (1.5,2,3,4):
    print(f"  height {z} kpc:", {x: round(rm_edge(x,z),1) for x in (3,6,9,12)})
out["mw_like"]={f"{x},{z}":rm_edge(x,z) for x in (3,6,9,12) for z in (1.5,2,3,4)}
# average over the patch 3<|x|<13 kpc, 1.5<|z|<4 kpc
xs=np.linspace(3,13,21); zs=np.linspace(1.5,4,11); patch=np.mean([[rm_edge(x,z) for x in xs] for z in zs])
print(f"  mean over 3-13 kpc along the disk and 1.5-4 kpc above it: {patch:.0f} rad/m2")
out["mw_patch_mean"]=patch

# ---------- 4. What it takes: sightlines needed through the inner halo
sig=10.0   # scatter of one background source after foreground removal, rad/m2 (NVSS-like 14-16; LoTSS-like 6-8)
for target in (5,3):
    n=(sig/ (patch/target))**2*2   # factor 2: sin^2(2 theta) averages 1/2
    print(f"To see a {patch:.0f} rad/m2 pattern at {target} sigma: about {n:.0f} sightlines through such patches")
area=4*(10*2.5)   # four patches of 10 kpc x 2.5 kpc
for D in (5,10,20):
    kpc_per_deg=D*1000*np.pi/180
    print(f"  galaxy at {D} Mpc: the four patches cover {area/kpc_per_deg**2:.4f} deg2 -> {40*area/kpc_per_deg**2:.2f} sources at 40 per deg2; {1.1*area/kpc_per_deg**2:.3f} at NVSS's 1.1")
out["forecast"]={"patch_area_kpc2":area}
json.dump(out,open("sign_results.json","w"),indent=1)
