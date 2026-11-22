"""Workbench Entry 49: the line-width test on a second, independent sample at redshift 0.6 to 1.
LEGA-C (van der Wel et al. 2021) gives a stellar velocity dispersion and a virial mass for each galaxy; COSMOS2015 (Laigle et al. 2016)
gives stellar masses from photometry for the same field.  One box of sky, RA 150.0-150.2, Dec 2.1-2.4, typed in from VizieR.
Matched by position.  Then the same question as Entry 45: dynamical over stellar mass under each distance law, against nearby galaxies."""
import os, json
import numpy as np
from scipy import integrate
HERE=os.path.dirname(os.path.abspath(__file__))
L=[l.split() for l in open(os.path.join(HERE,"legac_box.txt")) if not l.startswith("#")]
C=[l.split() for l in open(os.path.join(HERE,"cosmos_box.txt")) if not l.startswith("#")]
f=lambda s: np.nan if s=="NA" else float(s)
lid=np.array([int(r[0]) for r in L]); lra,lde,lz,ls,les,lhd,ln,lmv=np.array([[f(x) for x in r[1:]] for r in L]).T
cra,cde,czp,cms,csf,cks=np.array([[f(x) for x in r] for r in C]).T
print(f"{len(L)} spectra of {len(set(lid))} galaxies; {len(C)} photometric galaxies above 10^10.8 solar masses")
# one row per galaxy: spectra of the same galaxy averaged (dispersion and virial mass)
G=[]
for i in sorted(set(lid)):
    q=lid==i; G.append([i,lra[q][0],lde[q][0],lz[q].mean(),ls[q].mean(),np.nanmean(lhd[q]) if np.any(~np.isnan(lhd[q])) else np.nan,ln[q][0],lmv[q].mean(),q.sum(),np.ptp(lmv[q])])
G=np.array(G); rep=G[G[:,8]>1]; print(f"galaxies observed twice: {len(rep)}; their two virial masses differ by a median {np.median(rep[:,9]):.2f} dex (largest {rep[:,9].max():.2f})")
# match by position within 1 arcsec
M=[]
for g in G:
    d=np.hypot((cra-g[1])*np.cos(np.radians(g[2])),cde-g[2])*3600; j=np.argmin(d)
    if d[j]<1.0: M.append(list(g[:8])+[cms[j],csf[j],czp[j],d[j]])
M=np.array(M); print(f"matched within 1 arcsec: {len(M)}; median offset {np.median(M[:,11]):.2f} arcsec; photometric minus spectroscopic redshift: median {np.median(M[:,10]-M[:,3]):+.3f}, scatter {np.std(M[:,10]-M[:,3]):.3f}")
z,sig,hd,n,lmv,lms,lsfr=M[:,3],M[:,4],M[:,5],M[:,6],M[:,7],M[:,8],M[:,9]
ssfr=lsfr-lms; y=lmv-lms
Om=0.3; DM=np.vectorize(lambda x: integrate.quad(lambda t:1/np.sqrt(Om*(1+t)**3+1-Om),0,x)[0])
LAW={"expanding":lambda x:np.ones_like(x),"lerner":lambda x:(1+x)**2*DM(x)/x,"tired":lambda x:(1+x)**2*DM(x)/np.log(1+x),"tired_t":lambda x:(1+x)*DM(x)/np.log(1+x)}
LOCAL=-np.log10(0.59); zl=0.06; rng=np.random.default_rng(49); out={"n_matched":int(len(M))}
def stats(q,label):
    res={"n":int(q.sum()),"zmed":float(np.median(z[q])),"laws":{}}
    line=f"{label:52s} n={q.sum():3d}, median z {np.median(z[q]):.2f}, median dispersion {np.median(sig[q]):.0f} km/s, log M* {np.median(lms[q]):.2f}:"
    for k,F in LAW.items():
        yy=y[q]+np.log10(F(z[q])); loc=LOCAL+np.log10(F(np.array([zl])))[0]; m=yy.mean()-loc; e=yy.std(ddof=1)/q.sum()**0.5
        res["laws"][k]=[float(m),float(e),float(np.median(yy)-loc)]; line+=f" {k} x{10**m:.2f} ({m:+.2f}+-{e:.2f})"
    res["scatter"]=float(y[q].std(ddof=1)); res["over"]=int(np.sum(y[q]<0)); print(line+f" | scatter {y[q].std(ddof=1):.2f} dex; more stars than mass: {int(np.sum(y[q]<0))}"); return res
quiet=ssfr<-11.0
print(f"\nquiescent by photometry (specific star formation below 10^-11 per year): {quiet.sum()} of {len(M)}; by spectrum (H-delta below 2): {np.nansum(hd<2):.0f}")
out["samples"]={}
for lab,q in (("all matched",np.ones(len(M),bool)),("quiescent (photometry)",quiet),("quiescent, M* >= 10^11",quiet&(lms>=11.0)),("quiescent, Sersic n >= 2.5",quiet&(n>=2.5)),
              ("quiescent and old by spectrum (H-delta < 2)",quiet&(hd<2)),("H-delta < 2 alone",hd<2),("star-forming (photometry)",~quiet),
              ("quiescent, z < 0.8",quiet&(z<0.8)),("quiescent, z >= 0.8",quiet&(z>=0.8)),("quiescent, 0.80 < z < 0.87",quiet&(z>0.80)&(z<0.87))):
    out["samples"][lab]=stats(q,lab)
# trend with redshift inside the quiescent sample, per law
lx=np.log10(1+z)
print("\ntrend inside the quiescent sample, power of (1+z):"); out["trend"]={}
for k,F in LAW.items():
    q=quiet; yy=(y+np.log10(F(z)))[q]; X=np.c_[np.ones(q.sum()),lx[q]]; b=np.linalg.lstsq(X,yy,rcond=None)[0][1]; bs=[]
    for _ in range(3000):
        j=rng.integers(0,q.sum(),q.sum()); bs.append(np.linalg.lstsq(X[j],yy[j],rcond=None)[0][1])
    print(f"   {k:10s} {b:+.2f} +- {np.std(bs):.2f}"); out["trend"][k]=[float(b),float(np.std(bs))]
out["points"]=[[round(float(a),4),round(float(b-LOCAL),3),int(c)] for a,b,c in zip(z,y,quiet)]
json.dump(out,open(os.path.join(HERE,"legac_test.json"),"w"),indent=1)
