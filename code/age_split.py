"""Workbench Entry 51: can older stars hide the mass that static distances need?
The published stellar masses cap each galaxy's age at the age of the expanding universe at its redshift.  A static universe has no such cap.
But the cap can only bind for galaxies whose spectra look old.  Strong H-delta absorption means the light is dominated by stars a billion
or two years old, far below any cap.  So: split the LEGA-C sample of Entry 49 by H-delta and ask whether the static excess is confined
to the old-looking galaxies."""
import os, json
import numpy as np
from scipy import integrate, stats
HERE=os.path.dirname(os.path.abspath(__file__))
L=[l.split() for l in open(os.path.join(HERE,"legac_box.txt")) if not l.startswith("#")]
C=[l.split() for l in open(os.path.join(HERE,"cosmos_box.txt")) if not l.startswith("#")]
f=lambda s: np.nan if s=="NA" else float(s)
lid=np.array([int(r[0]) for r in L]); lra,lde,lz,ls,les,lhd,ln,lmv=np.array([[f(x) for x in r[1:]] for r in L]).T
cra,cde,czp,cms,csf,cks=np.array([[f(x) for x in r] for r in C]).T
G=[]
for i in sorted(set(lid)):
    q=lid==i; G.append([i,lra[q][0],lde[q][0],lz[q].mean(),ls[q].mean(),np.nanmean(lhd[q]) if np.any(~np.isnan(lhd[q])) else np.nan,ln[q][0],lmv[q].mean()])
M=[]
for g in G:
    d=np.hypot((cra-g[1])*np.cos(np.radians(g[2])),cde-g[2])*3600; j=np.argmin(d)
    if d[j]<1.0: M.append(list(g)+[cms[j],csf[j]])
M=np.array(M); z,sig,hd,n,lmv,lms,lsfr=M[:,3],M[:,4],M[:,5],M[:,6],M[:,7],M[:,8],M[:,9]
Om=0.3; DM=np.vectorize(lambda x: integrate.quad(lambda t:1/np.sqrt(Om*(1+t)**3+1-Om),0,x)[0]); FL=lambda x:(1+x)**2*DM(x)/x
LOCAL=-np.log10(0.59); zl=0.06; y=lmv-lms-LOCAL; ys=y+np.log10(FL(z))-np.log10(FL(np.array([zl])))[0]
ok=~np.isnan(hd); quiet=(lsfr-lms)<-11.0; out={}
print(f"{len(M)} matched galaxies, {ok.sum()} with an H-delta measurement; age of the expanding universe at z = 0.6, 0.8, 1.0: "+", ".join(f"{977.8/70*integrate.quad(lambda t:1/((1+t)*np.sqrt(0.3*(1+t)**3+0.7)),zz,np.inf)[0]:.1f}" for zz in (0.6,0.8,1.0))+" Gyr")
print("\nsplit by H-delta absorption (Angstrom): stronger means younger stars dominate the light"); out["bins"]=[]
for lo,hi,lab in ((-9,1,"below 1: old"),(1,3,"1 to 3"),(3,5,"3 to 5"),(5,99,"above 5: young")):
    q=ok&(hd>=lo)&(hd<hi)
    if q.sum()<3: continue
    m,e=y[q].mean(),y[q].std(ddof=1)/q.sum()**.5; ms=ys[q].mean()
    print(f"   H-delta {lab:16s} n={q.sum():2d} (quiet {int((q&quiet).sum()):2d}), median z {np.median(z[q]):.2f}, median dispersion {np.median(sig[q]):.0f}, log M* {np.median(lms[q]):.2f}: standard x{10**m:.2f} (+-{100*(10**e-1):.0f}%), Lerner x{10**ms:.2f}")
    out["bins"].append(dict(label=lab,n=int(q.sum()),hd=float(np.median(hd[q])),std=float(m),err=float(e),lerner=float(ms)))
r,p=stats.spearmanr(hd[ok],y[ok]); b=np.polyfit(hd[ok],y[ok],1); rng=np.random.default_rng(51); bs=[np.polyfit(hd[ok][j],y[ok][j],1)[0] for j in (rng.integers(0,ok.sum(),ok.sum()) for _ in range(3000))]
print(f"\ntrend of log(dynamical/stellar) with H-delta, all {ok.sum()}: slope {b[0]:+.3f} +- {np.std(bs):.3f} dex per Angstrom; rank correlation {r:+.2f} (chance {p:.2f})")
out["trend"]=[float(b[0]),float(np.std(bs)),float(r),float(p)]
young=ok&(hd>=3); old=ok&(hd<1)
d=ys[young].mean()-ys[old].mean(); ed=np.hypot(ys[young].std(ddof=1)/young.sum()**.5,ys[old].std(ddof=1)/old.sum()**.5)
print(f"young-looking (H-delta 3 or more, n={young.sum()}) against old-looking (below 1, n={old.sum()}): Lerner ratio x{10**ys[young].mean():.2f} and x{10**ys[old].mean():.2f}; difference {d:+.2f} +- {ed:.2f} dex")
print(f"if the cap hid a factor 3.3 in the old-looking galaxies only, the young-looking ones would sit at the nearby value: they sit {ys[young].mean()/ (ys[young].std(ddof=1)/young.sum()**.5):.0f} standard errors above it")
out["young_old"]=[float(ys[young].mean()),int(young.sum()),float(ys[old].mean()),int(old.sum()),float(d),float(ed)]
out["pts"]=[[round(float(a),2),round(float(b_),3),int(c)] for a,b_,c in zip(hd[ok],ys[ok],quiet[ok])]
out["bins_l"]=[[float(np.median(hd[q])),float(ys[q].mean()),float(ys[q].std(ddof=1)/q.sum()**.5),int(q.sum())] for q in (ok&(hd<1),ok&(hd>=1)&(hd<3),ok&(hd>=3))]
print("three bins under Lerner:",[(round(a_,2),round(10**b_,2),round(c_,3),n_) for a_,b_,c_,n_ in out["bins_l"]])
print("colour range g-r 0.5 to 0.9 moves the K-band mass-to-light ratio by a factor",round(10**(0.197*0.4),2),"and the r-band one by",round(10**(1.097*0.4),2))
# what colour allows: published colour - mass-to-light relation (Bell et al. 2003): log(M/L_r) = -0.306 + 1.097 (g - r)
for dm in (3.3,2.6,2.0): print(f"to raise the stellar mass {dm} times at fixed light, the relation needs g - r redder by {np.log10(dm)/1.097:.2f} mag")
json.dump(out,open(os.path.join(HERE,"age_split.json"),"w"),indent=1)
