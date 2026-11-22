"""Workbench Entry 22: is intergalactic space a radio fog?  Derived numbers.
Lerner's reading of the CMB needs the intergalactic medium to absorb or scatter radio waves between roughly 150 MHz and 200 GHz.
His stated evidence and numbers (Lerner 1990, 1993; restated in 'A Reassessment of the Galactic Origin of Light Elements Hypothesis'):
  at fixed infrared luminosity, radio luminosity falls with distance as D^-0.41+-0.06 (1990), D^-0.32+-0.04 (1993, 1.4 GHz, z < 0.07),
  D^-0.42+-0.087 (LOFAR 150 MHz, to 600 Mpc); a factor 10 out to 300 Mpc; 8+-4 % of radio from beyond 600 Mpc arrives;
  optical depth 7.6 by z = 2.2+-1.0 and 10 by z = 4.3+-2.
This script turns the catalogue measurements in fog_results.json (made by radio_fog.js) into comparisons with those numbers."""
import json, math
import numpy as np
R=json.load(open("fog_results.json")); out={}
ln10=math.log(10)
print("== 1. What the claimed law means")
for k in (0.32,0.41,0.42): print(f"  radio ~ D^-{k}: {100*(1-10**-k):.0f}% lost per tenfold distance; optical depth {k*ln10:.2f} per tenfold; from 1 to 300 Mpc a factor {300**k:.1f}")
print(f"  'only 8% from beyond 600 Mpc' corresponds to an exponent {-math.log10(0.08)/math.log10(600):.2f} from 1 Mpc")

print("\n== 2. Local test: slope of q (infrared/radio, log) with log distance at fixed infrared luminosity (fog predicts +k)")
def cmp(name,c,se,claims):
    s=f"  {name}: {c:+.3f} +- {se:.3f}; 2-sigma upper limit {c+2*se:+.3f} -> at most {100*(1-10**-(max(c+2*se,0))):.0f}% lost per tenfold distance"
    for k,ke in claims: s+=f"; vs {k}+-{ke}: {(k-c)/math.hypot(se,ke):.1f} sigma ({(k-c)/se:.1f} on our error alone)"
    print(s); return {"c":c,"se":se,"upper2":c+2*se,"max_loss_per_decade":1-10**-(max(c+2*se,0)),"vs":[[k,ke,(k-c)/math.hypot(se,ke),(k-c)/se] for k,ke in claims]}
n=R["nvss"]["f2"]["fixedL"]; l=R["lotss"]["f06"]["fixedL"]
out["nvss"]=cmp("NVSS 1.4 GHz, 1,387 galaxies",n["c"],n["se"],[(0.32,0.04),(0.41,0.06)])
out["lotss_all"]=cmp("LoTSS 144 MHz, all 5,619 galaxies from 29 Mpc (biased low, see section 6)",l["c"],l["se"],[(0.42,0.087),(0.32,0.04)])
lf=R["lotss"]["f06"]["fixedL_by_distance"][">= 100 Mpc"]
out["lotss"]=cmp("LoTSS 144 MHz, 3,665 galaxies beyond 100 Mpc",lf[0],lf[1],[(0.42,0.087),(0.32,0.04)])
w=R["wenss"]["fixedL"]; print(f"  WENSS 325 MHz, 329 galaxies (17% undetected): {w['c']:+.2f} +- {w['se']:.2f} -- cannot tell the two apart")
i=R["nvss"]["f2"]["inject"]; print(f"  injection check (NVSS): a D^-0.32 fog added to the data is recovered as {i['0.32'][0]:+.3f} (expected {n['c']+0.32:+.3f}); D^-0.41 as {i['0.41'][0]:+.3f}")
# raw trend: the far bins against the 100-155 Mpc bins (LoTSS), and what a fog would give between the same distances
b=R["lotss"]["f06"]["byD"]; base=np.mean([b[2][2],b[3][2]]); far=b[-1]
dq=far[2]-base; print(f"  LoTSS raw medians: q = {base:.3f} at 100-155 Mpc, {far[2]:.3f} at {far[0]:.0f} Mpc: +{dq:.2f} dex, i.e. {100*10**-dq:.0f}% of the radio if all of it were fog;")
for k in (0.32,0.42): print(f"     a D^-{k} fog predicts +{k*math.log10(far[0]/128):.2f} dex between those distances ({100*10**-(k*math.log10(far[0]/128)):.0f}%)")
s=R["lotss"]["f06"]["sameL_11_115_byD"]; print(f"  same luminosity (10^11.0-11.5 Lsun) at {s[0][0]} and {s[2][0]} Mpc: q = {s[0][2]:.2f}+-{s[0][3]:.2f} and {s[2][2]:.2f}+-{s[2][3]:.2f}; fog (0.42) predicts +{0.42*math.log10(s[2][0]/s[0][0]):.2f}")
out["raw_far"]={"dq":dq,"fraction":10**-dq}

print("\n== 3. Far test: optically selected quasars (2QZ) detected by NVSS above 2.5 mJy")
def tau(z,slow=False):                      # linear interpolation between Lerner's stated anchor points
    pts=[(0.14,-math.log(0.08)),(3.2 if slow else 2.2,7.6),(6.3 if slow else 4.3,10.0)]
    for (z0,t0),(z1,t1) in zip(pts,pts[1:]):
        if z<=z1: return t0+(t1-t0)*(z-z0)/(z1-z0)
    return 10.0
Q=R["quasars"]["twoqz"]; S=np.log10(Q["flux_cum_mJy"]); Nc=np.log10(Q["flux_cum_n"])
def G(k):                                   # fraction of detected quasars that would stay above 2.5 mJy if dimmed k-fold (measured flux distribution)
    x=math.log10(2.5*k)
    if x>=S[-1]: return 10**(Nc[-1]+(Nc[-1]-Nc[-2])/(S[-1]-S[-2])*(x-S[-1])-Nc[0])
    return 10**(np.interp(x,S,Nc)-Nc[0])
bins=Q["bins"]; z0,N0,d0,c0=bins[0]; f0=(d0-c0)/N0; rows=[]
print("   z     N    net detections   fraction        fog, even dimming   fog, seen through holes   (expected detections)")
tot={"obs":0,"A":0,"B":0,"As":0,"Bs":0}
for z,N,d,c in bins:
    f=(d-c)/N; e=math.sqrt(d)/N; dt=tau(z)-tau(z0); dts=tau(z,True)-tau(z0,True)
    A=f0*G(math.exp(dt)); B=f0*math.exp(-dt); As=f0*G(math.exp(dts)); Bs=f0*math.exp(-dts)
    rows.append([z,N,d-c,f,e,A,B,As,Bs,tau(z)])
    if z>1.5: tot["obs"]+=d-c; tot["A"]+=A*N; tot["B"]+=B*N; tot["As"]+=As*N; tot["Bs"]+=Bs*N
    print(f"  {z:.2f} {N:5d}   {d-c:6.1f}       {100*f:.2f}+-{100*e:.2f}%      {100*A:.2f}%              {100*B:.3f}%          ({A*N:.1f}, {B*N:.2f})   tau={tau(z):.1f}")
print(f"  z > 1.5: {tot['obs']:.0f} detections; fog predicts {tot['A']:.0f} (even dimming) or {tot['B']:.1f} (holes); with the slowest reading of his numbers {tot['As']:.0f} or {tot['Bs']:.1f}")
f_lo=sum(d-c for z,N,d,c in bins[:2])/sum(N for z,N,d,c in bins[:2]); f_hi=sum(d-c for z,N,d,c in bins[4:])/sum(N for z,N,d,c in bins[4:])
e_lo=math.sqrt(sum(d for z,N,d,c in bins[:2]))/sum(N for z,N,d,c in bins[:2]); e_hi=math.sqrt(sum(d for z,N,d,c in bins[4:]))/sum(N for z,N,d,c in bins[4:])
print(f"  detected fraction: z < 0.9 {100*f_lo:.2f}+-{100*e_lo:.2f}%, z > 1.5 {100*f_hi:.2f}+-{100*e_hi:.2f}%; ratio {f_hi/f_lo:.2f}+-{f_hi/f_lo*math.hypot(e_lo/f_lo,e_hi/f_hi):.2f}")
L=R["quasars"]["lbqs"]["bins"]; print("  LBQS (brighter quasars):", ", ".join(f"z~{z}: {100*(d-c)/N:.1f}+-{100*math.sqrt(d)/N:.1f}%" for z,N,d,c in L))
out["quasars"]={"rows":rows,"z_gt_1p5":tot,"f_lo":f_lo,"f_hi":f_hi,"ratio":f_hi/f_lo,"lbqs":[[z,(d-c)/N,math.sqrt(d)/N] for z,N,d,c in L]}

print("\n== 4. Other consequences of the stated optical depths")
for name,z,E in (("FRB 20220610A",1.016,2e42),("FRB 20240304B",2.148,1e41)):
    t=tau(z); print(f"  {name} at z = {z}: tau = {t:.1f}, flux cut {math.exp(t):.0f}-fold; quoted energy {E:.0e} erg would become {E*math.exp(t):.0e} erg")
for z in (0.52,1.0,1.91): print(f"  SZ shadow of a cluster at z = {z}: tau = {tau(z):.1f}; {100*math.exp(-tau(z)):.1f}% of the contrast survives scattering")
zq=2.2; print(f"  infrared/radio ratio at z = {zq}: fog predicts q higher by {(tau(zq)-tau(0.14))/ln10:.1f} dex than at 600 Mpc; Delhaize et al. (2017) measure q = 2.88(1+z)^-0.19, i.e. {2.88*(1+zq)**-0.19-2.88:+.2f}")
out["tau"]={str(z):tau(z) for z in (0.14,0.48,0.52,1.0,1.016,1.91,2.148,2.2,2.25,2.62)}
json.dump(out,open("fog_model_out.json","w"),indent=1,default=float)
print("\n== 5. Cluster shadows (SZ) against X-ray masses: fitted fog strength s (1 = the stated fog, 0 = none)")
Z=R["sz"]
def tau_local(z):                           # the stated local law (8% arrives from 600 Mpc) joined to the far anchors
    return 0.395*math.log(4282.7*z) if z<=0.14 else tau(z)
for name,key in (("ACT x MCXC","act_mcxc"),("Planck x MCXC","psz_mcxc")):
    d=Z[key]; print(f"  {name}: {d['n']} clusters, s = {d['s']:+.2f} +- {d['se']:.2f}; the stated fog is {(1-d['s'])/d['se']:.1f} sigma away; 2-sigma upper limit {d['s']+2*d['se']:.2f}")
    for v,(s_,e_,n_) in d["variants"].items(): print(f"      {v}: {s_:+.2f} +- {e_:.2f} ({n_})")
d=Z["act_xcs"]; print(f"  ACT x XCS temperatures: {d['n']} clusters, s = {d['s']:+.2f} +- {d['se']:.2f}")
print("  what the stated fog would do to the SZ mass (and to the shadow itself), relative to z = 0.07:")
for z in (0.17,0.28,0.44,0.55,0.83,1.1):
    dt=tau_local(z)-tau_local(0.07); print(f"    z = {z}: optical depth +{dt:.2f}; shadow x{math.exp(-dt):.2f}; SZ mass x{math.exp(-0.6*dt):.2f} ({-0.6*dt/ln10:+.2f} dex)")
a=Z["act_mcxc"]["bins"]; p=Z["psz_mcxc"]["bins"]
print(f"  measured change from z = {p[0][0]} to {p[-1][0]} (Planck): {p[-1][2]-p[0][2]:+.3f} +- {math.hypot(p[-1][3],p[0][3]):.3f} dex; fog: {-0.6*(tau_local(p[-1][0])-tau_local(p[0][0]))/ln10:+.2f}")
print(f"  measured change from z = {a[1][0]} to {a[-1][0]} (ACT):    {a[-1][2]-a[1][2]:+.3f} +- {math.hypot(a[-1][3],a[1][3]):.3f} dex; fog: {-0.6*(tau_local(a[-1][0])-tau_local(a[1][0]))/ln10:+.2f}")
c=Z["act_counts"]; dt=tau_local(1.1)-tau_local(0.14)
print(f"  ACT lists {c['>=1'][0]} clusters at z >= 1 with median mass {c['>=1'][1]}e14; behind the stated fog their shadows would be {math.exp(dt):.0f}x deeper at the source, and their masses {math.exp(0.6*dt):.1f}x larger ({c['>=1'][1]*math.exp(0.6*dt):.0f}e14)")
out["sz"]={"fog_curve":[[z,-0.6*(tau_local(z)-tau_local(0.07))/ln10] for z in [0.05,0.07,0.1,0.14,0.2,0.3,0.4,0.5,0.6]]}
json.dump(out,open("fog_model_out.json","w"),indent=1,default=float)
print("\n== 6. Spectra between 144 MHz and 1.4 GHz (same galaxies in both surveys)")
P=R["spectra"]; f=P["f12"]
aL=f["alpha_vs_L_at_fixed_D_L>=10.4"]; print(f"  at fixed distance the spectral index flattens with luminosity: {aL[0]:+.3f} +- {aL[1]:.3f} per tenfold ({aL[2]} galaxies above 10^10.4 Lsun)")
q1=f["q14_vs_L_at_fixed_D_L>=10.4"]; q2=f["q144_vs_L_at_fixed_D_L>=10.4"]; print(f"  so q rises with luminosity at fixed distance by {q2[0]:+.3f} +- {q2[1]:.3f} at 144 MHz and {q1[0]:+.3f} +- {q1[1]:.3f} at 1.4 GHz")
r=f["raw_slopes_L>=10.4"]; print(f"  ignoring luminosity, q rises with distance by {r['q144'][0]:+.2f} +- {r['q144'][1]:.2f} (144 MHz) and {r['q14'][0]:+.2f} +- {r['q14'][1]:.2f} (1.4 GHz) per tenfold; luminosity rises as D^{r['dlogL_dlogD']:.2f} in that sample")
aD=P["f06"]["alpha_vs_D_at_fixed_L"]; print(f"  survey systematic: at fixed luminosity the 144 MHz flux, relative to 1.4 GHz, rises with distance as D^{aD[0]*math.log10(1400/144):+.3f}+-{aD[1]*math.log10(1400/144):.3f} ({aD[2]} galaxies)")
shift=aD[0]*math.log10(1400/144); c=lf[0]+shift; e=math.hypot(lf[1],aD[1]*math.log10(1400/144))
print(f"  if all of that is 144 MHz flux lost in nearer galaxies, the corrected 144 MHz slope beyond 100 Mpc is {c:+.3f} +- {e:.3f}: {(0.42-c)/math.hypot(e,0.087):.1f} sigma below 0.42+-0.087, {(0.32-c)/math.hypot(e,0.04):.1f} sigma below 0.32+-0.04")
print(f"  Lerner's exponents (0.42 at 150 MHz, 0.32 at 1.4 GHz) would make the spectral index fall with distance by {-(0.42-0.32)/math.log10(1400/144):+.2f} per tenfold; measured {aD[0]:+.3f} +- {aD[1]:.3f}")
out["spectra"]={"corrected_lotss_slope":[c,e],"sigma_vs_042":(0.42-c)/math.hypot(e,0.087),"sigma_vs_032":(0.32-c)/math.hypot(e,0.04)}
json.dump(out,open("fog_model_out.json","w"),indent=1,default=float)


