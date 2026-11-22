"""Entry 24, second version: how much of the trend is the selection of the sample?
The clusters of Sanders et al. were picked because their CMB shadow was strong (SPT), so at any temperature the ones with an unusually deep shadow
are the ones in the sample. Regressing shadow depth on temperature in such a sample is biased. Two answers:
 (a) turn the regression round: at a given shadow depth, the temperature is not affected by a cut on the shadow;
 (b) simulate the cut under each hypothesis and push the simulated samples through the same fit."""
import numpy as np, json
from scipy import integrate
from scipy.special import ndtr
cH=299792.458/70.; Ez=lambda z: np.sqrt(0.3*(1+z)**3+0.7)
DA=np.vectorize(lambda z: cH*integrate.quad(lambda x:1/Ez(x),0,z)[0]/(1+z)); d_static=lambda z: cH*np.log(1+z)
q=np.array([[float(x) for x in l.split()] for l in open("act_q_curve.txt") if not l.startswith("#")])
logQ=lambda th: np.interp(np.log10(np.clip(th,q[0,0],q[-1,0])),np.log10(q[:,0]),np.log10(q[:,1]))
am=180*60/np.pi
th_exp=lambda T,z,da: 1148*(T/5)**0.51/Ez(z)/(da*1000)*am
th_sta=lambda T,z: 1148*(T/5)**0.51/(d_static(z)*1000)*am
D=[l.split() for l in open("clusters_y_T.txt") if not l.startswith("#")]
z=np.array([float(r[0]) for r in D]); T=np.array([float(r[1]) for r in D]); src=np.array([r[2] for r in D]); lf=np.array([float(r[3]) for r in D]); da=DA(z)
lT=np.log10(T/5); lE=np.log10(Ez(z)); ys=lf-logQ(th_sta(T,z)); ye=lf-logQ(th_exp(T,z,da))
def ols(X,y):
    b=np.linalg.lstsq(X,y,rcond=None)[0]; r=y-X@b; s2=r@r/(len(y)-X.shape[1]); return b,np.sqrt(np.diag(s2*np.linalg.pinv(X.T@X))),np.sqrt(s2)
out={}
# local relation (clusters picked by their X-rays, shadow well above the ACT limit)
loc=(z<0.2)&(T>=4); b,e,rms=ols(np.c_[np.ones(loc.sum()),lT[loc]],ys[loc]); a0,al0=b; print(f"local relation (z<0.2, kT>=4, n={loc.sum()}): log signal = {a0:.3f} + {al0:.2f} log(T/5), scatter {rms:.3f} dex"); out["local"]=dict(n=int(loc.sum()),a=float(a0),alpha=float(al0),rms=float(rms))
S=src=="S"; n=S.sum(); zs=z[S]; das=da[S]
print(f"\nSanders et al. sample: {n} clusters, raw signal: lowest {lf[S].min():.2f}, 10th percentile {np.percentile(lf[S],10):.2f}, median {np.median(lf[S]):.2f} (log10 of 1e-4)")
def both(y,lTx,lEx):
    N=len(y); bf,ef,_=ols(np.c_[np.ones(N),lTx,lEx],y); bi,ei,_=ols(np.c_[np.ones(N),y,lEx],lTx); return bf[2],ef[2],bi[2],ei[2],bi[1],-bi[2]/bi[1]
for name,y in (("static sizes",ys),("expanding sizes",ye)):
    mf,sf,c,sc,bb,mi=both(y[S],lT[S],lE[S]); print(f"  {name:16s} forward (signal on T): m = {mf:+.2f}+-{sf:.2f} | reverse (T on signal): T ~ E^{c:+.2f}+-{sc:.2f} at fixed signal, signal power {bb:.2f}, implied m = {mi:+.2f}")
    out["sanders_"+name.split()[0]]=dict(m_fwd=float(mf),m_fwd_se=float(sf),c=float(c),c_se=float(sc),b=float(bb),m_inv=float(mi))
# (a) the temperature the static universe needs: take the local relation, read off the temperature that gives each cluster's size-corrected signal
Tneed=5*10**((ys[S]-a0)/al0); r=np.log10(Tneed/T[S]); bb,ee,_=ols(np.c_[np.ones(n),lE[S]],r)
print("\n(a) static universe: temperature needed to cast the observed shadow (local relation), against the temperature measured:")
for lo,hi in ((0.25,0.45),(0.45,0.65),(0.65,0.9),(0.9,1.3)):
    k=(zs>=lo)&(zs<hi); print(f"    z {lo}-{hi}: n={k.sum():2d}  needed {np.median(Tneed[k]):5.1f} keV, measured {np.median(T[S][k]):4.1f} keV, ratio {10**np.mean(r[k]):.2f} +- {10**np.mean(r[k])*np.log(10)*r[k].std(ddof=1)/np.sqrt(k.sum()):.2f}")
    out.setdefault("Tneed",[]).append([lo,hi,int(k.sum()),float(np.median(Tneed[k])),float(np.median(T[S][k])),float(10**np.mean(r[k])),float(10**np.mean(r[k])*np.log(10)*r[k].std(ddof=1)/np.sqrt(k.sum()))])
print(f"    trend of needed/measured with E(z): power {bb[1]:+.2f}+-{ee[1]:.2f}  (static universe: 0)")
out["Tneed_trend"]=[float(bb[1]),float(ee[1])]
# (b) simulation of the selection
rng=np.random.default_rng(2)
def sim(model,gamma,thr,sint=0.10,nrep=300,width=0.06):
    res=[]
    for _ in range(nrep):
        Ts=[];Ls=[]
        for zi,dai in zip(zs,das):
            while True:
                lTc=np.log10(2/5)+rng.exponential(1/(gamma*np.log(10)),400)        # dn/dlnT ~ T^-gamma above 2 keV
                Tc=5*10**lTc; Tc=np.minimum(Tc,25.0); lTc=np.log10(Tc/5)
                if model=="static": sig=a0+al0*lTc+logQ(th_sta(Tc,zi))
                else: sig=a0+al0*lTc+0.92*np.log10(Ez(zi))+logQ(th_exp(Tc,zi,dai))
                sig=sig+rng.normal(0,np.hypot(sint,0.05),len(sig)); ok=rng.random(len(sig))<ndtr((sig-thr)/width)
                if ok.any(): j=np.flatnonzero(ok)[0]; Ts.append(Tc[j]*10**rng.normal(0,0.05)); Ls.append(sig[j]); break
        Ts=np.array(Ts); Ls=np.array(Ls); lTs=np.log10(Ts/5); lEs=np.log10(Ez(zs)); y=Ls-(logQ(th_sta(Ts,zs)) if model=="static" else logQ(th_exp(Ts,zs,das)))
        mf,_,c,_,bb,mi=both(y,lTs,lEs); res.append((mf,c,np.median(Ts),np.median(Ts[zs>0.75])/np.median(Ts[zs<0.45])))
    return np.array(res)
print("\n(b) simulated samples cut on the raw signal like the real one, pushed through the same fits (mean, spread over 300 samples):")
obs=out["sanders_static"]; obe=out["sanders_expanding"]; hot=np.median(T[S][zs>0.75])/np.median(T[S][zs<0.45])
print(f"    observed: forward m = {obs['m_fwd']:+.2f} (static sizes), {obe['m_fwd']:+.2f} (expanding sizes); reverse T ~ E^{obs['c']:+.2f}; median kT {np.median(T[S]):.1f} keV; median kT at z>0.75 over z<0.45 = {hot:.2f}")
out["observed_Tratio"]=float(hot); out["sim"]={}
for model in ("static","expanding"):
    for gamma in (3.0,5.0):
        for thr in (-0.05,0.05):
            r=sim(model,gamma,thr); key=f"{model}_g{gamma:.0f}_t{thr:+.2f}"; out["sim"][key]=dict(m_fwd=[float(r[:,0].mean()),float(r[:,0].std())],c=[float(r[:,1].mean()),float(r[:,1].std())],Tmed=float(r[:,2].mean()),Tratio=[float(r[:,3].mean()),float(r[:,3].std())])
            print(f"    {model:9s} temperature function T^-{gamma:.0f}, cut at {thr:+.2f}: forward m = {r[:,0].mean():+.2f}+-{r[:,0].std():.2f} | reverse T ~ E^{r[:,1].mean():+.2f}+-{r[:,1].std():.2f} | median kT {r[:,2].mean():4.1f} keV | kT(z>0.75)/kT(z<0.45) = {r[:,3].mean():.2f}+-{r[:,3].std():.2f}")
print("\n    the same with more scatter between shadow and temperature (0.15 dex instead of 0.10):")
for model in ("static","expanding"):
    for gamma in (3.0,5.0):
        r=sim(model,gamma,0.0,sint=0.15); key=f"{model}_g{gamma:.0f}_s15"; out["sim"][key]=dict(m_fwd=[float(r[:,0].mean()),float(r[:,0].std())],c=[float(r[:,1].mean()),float(r[:,1].std())],Tmed=float(r[:,2].mean()),Tratio=[float(r[:,3].mean()),float(r[:,3].std())])
        print(f"    {model:9s} temperature function T^-{gamma:.0f}: forward m = {r[:,0].mean():+.2f}+-{r[:,0].std():.2f} | reverse T ~ E^{r[:,1].mean():+.2f}+-{r[:,1].std():.2f} | median kT {r[:,2].mean():4.1f} keV | kT(z>0.75)/kT(z<0.45) = {r[:,3].mean():.2f}+-{r[:,3].std():.2f}")
# distance of the observation from each simulated case, in units of the simulated spread
print("\n    distance of the observed values from each simulated case (sigma): forward m, reverse power, temperature ratio")
for key,v in out["sim"].items():
    d=[(obs["m_fwd"]-v["m_fwd"][0])/v["m_fwd"][1],(obs["c"]-v["c"][0])/v["c"][1],(hot-v["Tratio"][0])/v["Tratio"][1]]; v["dist"]=[float(x) for x in d]; print(f"      {key:24s} {d[0]:+.1f} {d[1]:+.1f} {d[2]:+.1f}")
json.dump(out,open("selection_sim.json","w"),indent=1)
