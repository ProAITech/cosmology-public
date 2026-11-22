"""Workbench Entry 20: can a force that acts on gas but not on stars (magnetic support) hold up the Galaxy's
rotation curve? Measured numbers come from gas_vs_stars.js (run in a browser tab on the maser parallaxes of
Reid+ 2019 and the cluster catalogue of Hunt & Reffert 2023)."""
import json, numpy as np
from scipy.special import i0, i1, k0, k1
G=4.30091e-6          # kpc (km/s)^2 / Msun
KK=1.0227e-3          # (km/s/kpc) -> 1/Myr
out={}

# ---------- measured
gas={"n":129,"R_range":[5.0,13.5],"V_at_8.15":232.1,"slope":-1.16,"scatter":9.9,"vR":[-4.9,1.3],"vR_scatter":11.5}
ages=[[0,20],[20,50],[50,120],[120,300],[300,1000],[1000,20000]]
obs={"n":[154,134,308,505,360,102],"age_med":[9.8,30.7,91.6,173.7,500.4,1352.1],"R_med":[8.3,8.2,8.1,8.2,8.6,9.1],
     "vR":[-2.5,0.4,0.1,-2.3,-2.5,2.1],"vR_err":[0.7,0.7,0.9,0.8,1.4,2.6],"vR_scatter":[7.2,6.5,12.8,15.1,20.7,21.2],
     "lead":[1.9,3.4,4.0,3.5,-0.5,1.1],"lead_err":[0.7,0.9,0.6,0.5,0.8,1.8],"lead_scatter":[6.7,8.6,8.6,9.7,12.7,14.3]}
# predictions from the same clusters' ages and radii (ages smeared by 0.2 dex): [median vR, rms vR, median lead]
pred={"0.33":[[19.3,25.1,-2.3],[47.2,46.8,-20.3],[18.3,43.2,-63.8],[-10.7,42.0,-47.5],[-1.8,42.2,-41.5],[0.2,41.7,-42.2]],
      "0.50":[[26.8,36.2,-2.7],[72.4,72.7,-26.1],[56.9,72.5,-98.7],[-17.8,67.1,-85.5],[0.3,67.6,-66.3],[3.3,66.7,-66.0]],
      "0.05":[[3.2,4.0,-0.4],[7.0,6.8,-4.1],[0.5,5.9,-8.8],[-1.1,5.8,-5.9],[0.2,5.8,-5.6],[0.1,5.9,-5.7]],
      "0.02":[[1.3,1.6,-0.2],[2.8,2.7,-1.6],[0.2,2.3,-3.5],[-0.6,2.3,-2.4],[0.0,2.3,-2.3],[-0.1,2.3,-2.3]]}
radial=[[5,7,54,236.6,1.6,33,240.9,2.0],[7,8,14,231.9,1.6,147,234.8,0.9],[8,9,23,229.4,1.6,389,235.1,0.5],[9,10,20,230.4,2.9,26,234.5,1.8]]
out.update({"gas":gas,"obs":obs,"pred":pred,"radial":radial})

# ---------- 1. the physics: a star born at the gas speed when gravity alone supports less
def epicycle(f,R=8.2,Vg=232.0):
    Vc=Vg*np.sqrt(1-f); dv=Vg-Vc; kap=np.sqrt(2)*Vc/R*KK
    return Vc,dv,kap,2*np.pi/kap
for f in (0.33,0.05):
    Vc,dv,kap,P=epicycle(f)
    print(f"f = {f:.2f}: gravity supports {Vc:.0f} km/s of the gas's 232; excess {dv:.1f} km/s; radial swing up to {np.sqrt(2)*dv:.0f} km/s; epicycle period {P:.0f} Myr")
    for t in (10,30,90):
        print(f"    age {t:3d} Myr: outward {np.sqrt(2)*dv*np.sin(kap*t):+6.1f} km/s, rotation relative to gas {dv*(np.cos(kap*t)-1):+6.1f} km/s")

# ---------- 2. how much non-gravitational support the dark-matter-free picture needs at these radii
def vdisk2(R,M,Rd):
    y=R/(2*Rd); return 2*G*M/Rd*y*y*(i0(y)*k0(y)-i1(y)*k1(y))
def vbar(R,Mstar=4.5e10,Rd=2.6,Mbulge=1.5e10,Mgas=1.0e10,Rg=5.0):
    return np.sqrt(vdisk2(R,Mstar,Rd)+G*Mbulge/R+vdisk2(R,Mgas,Rg))
print("\nBaryons alone (assumed: stellar disk 4.5e10 Msun with 2.6 kpc scale length, bulge 1.5e10, gas 1.0e10):")
need={}
for R in (6,8.15,10,12):
    Vo=232.1-1.16*(R-8.15); vb=vbar(R); vb_hi=vbar(R,Mstar=5.5e10); f=1-(vb/Vo)**2; f_hi=1-(vb_hi/Vo)**2
    need[R]=[vb,f,vb_hi,f_hi]
    print(f"  R = {R:5.2f} kpc: baryons {vb:.0f} km/s of {Vo:.0f} observed -> missing share of the centripetal force {100*f:.0f}% ({100*f_hi:.0f}% with a 20% heavier disk)")
out["needed"]=need

# ---------- 3. limits on a gas-only force from the three signatures (linear in f for small f)
k_vR=pred["0.05"][1][0]/0.05          # outward streaming at 20-50 Myr per unit f
k_lag=-pred["0.05"][2][2]/0.05        # lag at 50-120 Myr per unit f
k_rms=pred["0.05"][1][1]/0.05         # scatter at 20-50 Myr per unit f
sys_vR, sys_lag = 3.0, 5.0            # allowances: solar motion and the gas's own radial drift; the masers' known lag behind circular rotation
f1=(obs["vR"][1]+2*np.hypot(obs["vR_err"][1],sys_vR))/k_vR
f2=(-obs["lead"][2]+2*np.hypot(obs["lead_err"][2],sys_lag))/k_lag
f3=obs["vR_scatter"][1]/k_rms
print(f"\nLimits on the share f of the gas's centripetal force that does not act on stars (95%):")
print(f"  outward streaming at 20-50 Myr: measured {obs['vR'][1]:+.1f} +- {obs['vR_err'][1]} km/s, {k_vR:.0f} km/s per unit f -> f < {100*f1:.1f}% (statistical only: {100*(obs['vR'][1]+2*obs['vR_err'][1])/k_vR:.1f}%)")
print(f"  lag behind the gas at 50-120 Myr: measured {obs['lead'][2]:+.1f} +- {obs['lead_err'][2]} km/s (stars lead), {k_lag:.0f} km/s per unit f -> f < {100*f2:.1f}%")
print(f"  scatter at 20-50 Myr: measured {obs['vR_scatter'][1]} km/s in total, {k_rms:.0f} km/s per unit f -> f < {100*f3:.1f}% even if all the scatter came from this")
out["limits"]={"streaming":f1,"lag":f2,"scatter":f3}
print(f"  predicted at f = 0.33: streaming {pred['0.33'][1][0]:+.0f} km/s ({(pred['0.33'][1][0]-obs['vR'][1])/np.hypot(obs['vR_err'][1],sys_vR):.0f} sigma from the measurement), lag {pred['0.33'][2][2]:.0f} km/s, scatter {pred['0.33'][1][1]:.0f} km/s")

# ---------- 4. Cepheids out to 20 kpc (Gaia Collaboration, Drimmel+ 2023, VizieR J/A+A/674/A37/table2; motions and
# radial velocities from Gaia DR3). 1,901 stars with radial velocities, ages 17-620 Myr (median 100).
# [R range, n, median R, median age, v_phi, err, v_R, err, scatter of v_R, lead over the gas line, err,
#  share needed, predicted median v_R, predicted rms v_R, predicted lead]
cep=[[5,7,270,6.5,72.4,236.1,1.2,-1.5,1.7,21.7,2.0,1.2,0.23,8.7,31.2,-46.6],[7,9,490,8.0,85.0,232.0,0.7,-5.2,0.8,14.6,-0.2,0.7,0.30,19.1,39.1,-58.4],
     [9,11,403,10.1,101.9,226.3,0.8,-3.1,1.0,16.0,-3.4,0.9,0.38,40.1,52.1,-74.4],[11,13,298,12.0,112.8,228.8,1.0,5.1,1.0,13.9,1.1,1.0,0.46,62.2,64.9,-84.4],
     [13,16,327,14.1,128.5,230.1,0.8,-1.8,0.9,12.8,5.2,0.8,0.54,81.9,79.6,-87.9],[16,20,89,16.8,135.5,223.6,1.7,-7.9,1.5,11.6,1.0,1.5,0.62,99.1,93.0,-81.6]]
out["cepheids"]=cep
print("\nCepheids: scatter of radial velocity against the prediction, and the cap on f if all scatter came from a gas-only force (rms ~ V f / 2)")
for c in cep:
    print(f"  {c[0]:2d}-{c[1]:2d} kpc: n={c[2]:3d}, rotation {c[5]:.0f}+-{c[6]}, streaming {c[7]:+.1f}+-{c[8]}, scatter {c[9]:.1f} (predicted {c[14]:.0f} for the {100*c[12]:.0f}% needed) -> f < {100*2*c[9]/c[5]:.0f}%")
json.dump(out,open("gas_vs_stars_results.json","w"),indent=1,default=float)
