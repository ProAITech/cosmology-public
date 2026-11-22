"""Workbench Entry 30: the nitrogen that comes with helium made in intermediate-mass stars.
The helium-from-stars model (Lerner's GOLE paper, section 4) makes the helium in stars of 4-12 solar masses, which return it to the gas in winds.
It allows for the carbon these stars release (one carbon atom per 6,000-30,000 helium nuclei, 'consistent with low-metallicity models') and does not
mention nitrogen. In stars of this mass the carbon brought up from the core is burned to nitrogen at the bottom of the envelope before it is blown off.
This script takes published yields of such stars and asks what gas with a quarter of its mass in helium made this way would contain.
Yields (net, solar masses, lost in the wind over the star's life):
  Karakas (2010, MNRAS 403, 1413; VizieR J/MNRAS/403/1413), metallicity Z = 0.0001: initial mass, final (white dwarf) mass, 4He, 12C, 13C, 14N, 16O
  Gil-Pons et al. (2021, A&A 645, A10; VizieR J/A+A/645/A10), Z = 0.00001, Vassiliadis-Wood winds: initial mass, hydrogen and helium ejected, 4He, 12C, 13C, 14N, 16O"""
import numpy as np, json
K10=[(4.0,0.874,2.03e-1,5.02e-3,4.28e-4,2.63e-2,6.49e-4),(4.5,0.897,2.65e-1,3.44e-3,3.51e-4,2.24e-2,4.77e-4),(5.0,0.929,3.44e-1,2.73e-3,4.05e-4,2.40e-2,4.14e-4),(5.5,0.966,4.21e-1,2.50e-3,4.05e-4,2.79e-2,3.72e-4),(6.0,1.008,5.18e-1,1.86e-3,3.87e-4,3.39e-2,6.38e-4)]
G21=[(4.0,1.674472+1.256955,4.924313e-1,5.910526e-3,1.207540e-3,1.273353e-1,3.238156e-3),(5.0,1.916038+1.933389,9.368315e-1,5.355807e-3,1.282432e-3,1.467588e-1,3.002503e-3),
     (6.0,2.519705+2.274205,1.044621e+0,4.082616e-3,1.182496e-3,1.490042e-1,2.590566e-3),(7.0,3.121133+2.628525,1.175831e+0,4.039098e-3,1.085353e-3,1.093838e-1,2.086780e-3)]
obs=dict(leoP=dict(HeH=0.0837,OH=10**(7.17-12),NO=10**-1.36),plateau_NO=10**-1.60,sun=dict(HeH=10**(10.93-12),NH=10**(7.83-12),CH=10**(8.43-12),OH=10**(8.69-12)))
obs["leoP"]["NH"]=obs["leoP"]["OH"]*obs["leoP"]["NO"]
out={"models":[]}
print("helium nuclei made per atom of nitrogen, carbon and oxygen released (by number), and new helium as a share of the star's mass:")
rows=[]
for name,tab,z in (("Karakas 2010, Z = 1e-4",K10,1e-4),("Gil-Pons et al. 2021, Z = 1e-5",G21,1e-5)):
    for r in tab:
        M=r[0]; he,c12,c13,n14,o16=r[2:]; nHe=he/4; nN=n14/14; nC=c12/12+c13/13; nO=o16/16
        rem=r[1] if name.startswith("Karakas") else M-(r[1]+n14+c12+c13+o16)
        rows.append(dict(set=name,Z=z,M=M,rem=rem,y=he/M,HeN=nHe/nN,HeC=nHe/nC,HeO=nHe/nO,he=he,n14=n14))
        print(f"  {name:32s} {M:.1f} Msun: He/N = {nHe/nN:5.1f}   He/C = {nHe/nC:5.0f}   He/O = {nHe/nO:6.0f}   new helium {100*he/M:4.1f}% of the star, remnant {rem:.2f} Msun ({100*rem/M:.0f}%)")
out["models"]=rows
HeN=np.array([r["HeN"] for r in rows]); HeC=np.array([r["HeC"] for r in rows]); HeO=np.array([r["HeO"] for r in rows])
print(f"\nrange over the nine models: one nitrogen atom per {HeN.min():.0f}-{HeN.max():.0f} helium nuclei; one carbon atom per {HeC.min():.0f}-{HeC.max():.0f}; one oxygen atom per {HeO.min():.0f}-{HeO.max():.0f}")
print("the model's allowance: one carbon atom per 6,000-30,000 helium nuclei; nitrogen not mentioned")
# Salpeter-weighted mean over each set
for name in ("Karakas 2010, Z = 1e-4","Gil-Pons et al. 2021, Z = 1e-5"):
    s=[r for r in rows if r["set"]==name]; w=np.array([r["M"]**-2.35 for r in s]); he=sum(wi*r["he"] for wi,r in zip(w,s)); n=sum(wi*r["n14"] for wi,r in zip(w,s)); print(f"  weighted by a Salpeter mass function, {name}: He/N = {(he/4)/(n/14):.1f}")
    out.setdefault("salpeter",{})[name]=(he/4)/(n/14)
HeH=obs["leoP"]["HeH"]; S=obs["sun"]; L=obs["leoP"]
NH=[HeH/HeN.max(),HeH/HeN.min()]; CH=[HeH/HeC.max(),HeH/HeC.min()]; OH=[HeH/HeO.max(),HeH/HeO.min()]
print(f"\ngas with He/H = {HeH} (a quarter helium by mass) made entirely this way would hold:")
print(f"  N/H = {NH[0]:.1e} to {NH[1]:.1e}   = {NH[0]/S['NH']:.0f} to {NH[1]/S['NH']:.0f} times the Sun's ({S['NH']:.2e}); {NH[0]/L['NH']:.0f} to {NH[1]/L['NH']:.0f} times Leo P's ({L['NH']:.2e})")
print(f"  C/H = {CH[0]:.1e} to {CH[1]:.1e}   = {CH[0]/S['CH']:.2f} to {CH[1]/S['CH']:.1f} times the Sun's")
print(f"  O/H = {OH[0]:.1e} to {OH[1]:.1e}   = {OH[0]/S['OH']:.2f} to {OH[1]/S['OH']:.2f} times the Sun's;  N/O = {NH[0]/OH[1]:.0f} to {NH[1]/OH[0]:.0f}; observed N/O: {obs['plateau_NO']:.3f} (metal-poor plateau), {L['NO']:.3f} (Leo P), {S['NH']/S['OH']:.2f} (Sun)")
out["predicted"]=dict(NH=NH,CH=CH,OH=OH); out["obs"]=obs
print("\nturned round: the helium that can have passed through such stars, if ALL the nitrogen seen came from them:")
for nm,o in (("Leo P",L),("the Sun",S)):
    f=[o["NH"]*HeN.min()/o["HeH"],o["NH"]*HeN.max()/o["HeH"]]; print(f"  {nm:8s} N/H = {o['NH']:.2e}, He/H = {o['HeH']:.4f}: at most {100*f[0]:.3f}% to {100*f[1]:.3f}% of its helium"); out.setdefault("max_share",{})[nm]=f
# the model's own account: first generation leaves He/H = 0.11 by mass (Y about 0.10); later generations, which do contain helium, add the rest
later=(0.248-0.10)/0.248; print(f"\nif the first, helium-free generation made no nitrogen at all, the later generations still make {100*later:.0f}% of the helium: N/H = {later*NH[0]:.1e} to {later*NH[1]:.1e}, {later*NH[0]/S['NH']:.0f} to {later*NH[1]/S['NH']:.0f} times solar")
out["later_share"]=later
# how much carbon against the model's allowance
print(f"carbon: the model allows one atom per 6,000-30,000 helium nuclei; these yields give one per {HeC.min():.0f}-{HeC.max():.0f}, i.e. {6000/HeC.max():.0f} to {30000/HeC.min():.0f} times more")
out["carbon_excess"]=[6000/HeC.max(),30000/HeC.min()]
# lock-up in white dwarfs
print("\nwhite dwarfs left behind while the gas is brought to a quarter helium (all gas cycled through stars of one mass):")
lock=[]
for r in rows:
    f_ret=1-r["rem"]/r["M"]; dY=r["y"]/f_ret; n=np.log(1-0.25)/np.log(1-dY); locked=1-f_ret**n; lock.append(locked)
    print(f"  {r['set'][:12]:12s} {r['M']:.1f} Msun: returned gas gains {dY:.3f} in helium fraction per pass; {n:.1f} passes; {100*locked:.0f}% of the matter ends in white dwarfs")
out["locked"]=[float(min(lock)),float(max(lock))]; print(f"  range: {100*min(lock):.0f}% to {100*max(lock):.0f}%")
# the top of the model's mass range: stars that explode. Limongi & Chieffi (2018, ApJS 237, 13; VizieR J/ApJS/237/13, table 8), no rotation, [Fe/H] = -3:
# total ejected mass (wind and explosion) of H, 4He, 12C, 14N, 16O for 13 and 15 solar masses. Initial helium fraction 0.240.
print("\nstars that explode (Limongi & Chieffi 2018, [Fe/H] = -3, no rotation):")
out["massive"]=[]
for M,H,He,C12,N14,O16 in ((13,6.230,4.448,1.652e-1,1.017e-4,2.905e-1),(15,6.979,4.629,2.724e-1,1.046e-4,8.096e-1)):
    newHe=He-0.240*(H+He+C12+O16)      # helium in the ejecta beyond what the gas started with
    nHe=newHe/4; hO=nHe/(O16/16); hC=nHe/(C12/12); hN=nHe/(N14/14)
    print(f"  {M} Msun: new helium {newHe:.2f} Msun ({100*newHe/M:.0f}% of the star); one oxygen atom per {hO:.0f} helium nuclei, one carbon per {hC:.0f}, one nitrogen per {hN:.0f}")
    print(f"           gas with a quarter helium made this way: O/H = {HeH/hO:.1e} ({HeH/hO/S['OH']:.1f} times solar, {HeH/hO/L['OH']:.0f} times Leo P), C/H = {HeH/hC:.1e} ({HeH/hC/S['CH']:.1f} times solar)")
    out["massive"].append(dict(M=M,newHe=newHe,HeO=hO,HeC=hC,OH=HeH/hO,OH_solar=HeH/hO/S["OH"],OH_leoP=HeH/hO/L["OH"],CH_solar=HeH/hC/S["CH"]))
json.dump(out,open("nitrogen_check.json","w"),indent=1,default=float)
