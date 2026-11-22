"""Workbench Entry 24, second version: were galaxy clusters denser in the past?
Two things are measured for each cluster:
  T  : the temperature of its gas, from the shape of its X-ray spectrum (no distance needed);
  y~ : the depth of the shadow it casts on the CMB, as ACT measures it through a fixed filter matched to clusters 2.4 arcmin in radius.
CORRECTION TO THE FIRST VERSION. The first version rebuilt the shadow depth from the ACT catalogue as E(z)^2 x (uncorrected mass)^1.08 and called it
"the measured signal". It is not. That quantity is the measured signal divided by a filter factor Q, which allows for the cluster being smaller or larger
on the sky than the filter, and ACT computes Q from the angular size the standard model gives a cluster of that mass. Q falls from 1.4 for nearby clusters
to 0.4-0.5 beyond z = 1, so the correction itself grows with redshift roughly as E(z)^1.5. The test therefore needs an angular size, and is not free of
assumptions the way the first version claimed. This version does the comparison openly: it takes the raw signal, and corrects it for size twice, once with
the size an expanding universe gives a cluster of that temperature and once with the size a static universe gives it.
Data: first_sample_act.txt (125 clusters of the first version), temps_matched.txt (new temperatures: S = Sanders et al. 2018, Chandra, SPT clusters;
L = Lovisari et al. 2020, XMM; R = Reichert et al. 2011 compilation; GC/GX = Migkas et al. 2020, Chandra / XMM), act_signals.txt, act_q_curve.txt."""
import numpy as np, json
from scipy import integrate
cH=299792.458/70.
Ez=lambda z: np.sqrt(0.3*(1+z)**3+0.7)
DA=np.vectorize(lambda z: cH*integrate.quad(lambda x:1/Ez(x),0,z)[0]/(1+z))
d_static=lambda z: cH*np.log(1+z)
q=np.array([[float(x) for x in l.split()] for l in open("act_q_curve.txt") if not l.startswith("#")])
def logQ(theta): return np.interp(np.log10(np.clip(theta,q[0,0],q[-1,0])),np.log10(q[:,0]),np.log10(q[:,1]))
i1,i2=np.searchsorted(q[:,0],[1.2,2.3]); print(f"filter factor Q: {q[0,1]:.2f} at {q[0,0]:.1f} arcmin, 1.0 at 2.5 arcmin, peak {q[:,1].max():.2f}; below 2.3 arcmin Q ~ size^{np.polyfit(np.log10(q[1:7,0]),np.log10(q[1:7,1]),1)[0]:.2f}")
# ---------- assemble one record per ACT cluster
C={}
for l in open("first_sample_act.txt"):
    if l.startswith("#"): continue
    z,T,ly,snr,src,i,lf,lc,Mu,Mc,es=l.split(); c=C.setdefault(int(i),dict(z=float(z),ly=float(ly),snr=float(snr),lf=float(lf),lc=float(lc),T={})); c["T"].setdefault(src,[]).append((float(T),np.nan))
sig={int(l.split()[0]):l.split() for l in open("act_signals.txt") if not l.startswith("#")}
for l in open("temps_matched.txt"):
    if l.startswith("#"): continue
    src,i,z,T,eT,ly,snr,T2=l.split(); i=int(i); s=sig[i]
    c=C.setdefault(i,dict(z=float(z),ly=float(ly),snr=float(snr),lf=float(s[1]),lc=float(s[2]),T={})); c["T"].setdefault(src,[]).append((float(T),float(eT)))
prio=["S","M","L","A","GC","XC","R","XS","GX","B"]
rec=[]
for i,c in C.items():
    for s in prio:
        if s in c["T"]:
            good=[t for t,e in c["T"][s] if not (e>0.30)]
            if good: rec.append((i,c["z"],10**np.mean(np.log10(good)),s,c["ly"],c["lf"],c["lc"],c["snr"])); break
rec.sort(key=lambda r:r[1])
idx=np.array([r[0] for r in rec]); z=np.array([r[1] for r in rec]); T=np.array([r[2] for r in rec]); src=np.array([r[3] for r in rec]); ly=np.array([r[4] for r in rec]); lf=np.array([r[5] for r in rec]); lc=np.array([r[6] for r in rec]); snr=np.array([r[7] for r in rec])
lE=np.log10(Ez(z)); lT=np.log10(T/5); n=len(z)
print(f"{n} clusters with a usable temperature; z {z.min():.3f}-{z.max():.3f}; {int((z>=0.45).sum())} at z >= 0.45, {int((z>=0.7).sum())} at z >= 0.7; sources",{s:int((src==s).sum()) for s in prio if (src==s).any()})
open("clusters_y_T.txt","w").write("# z, kT keV, temperature source, log10 raw fixed-filter signal (1e-4), log10 signal corrected by ACT for size = log10(E^2 Mu^1.08), log10 optimal-filter central y (1e-4), ACT S/N, ACT catalogue row\n"+"\n".join(f"{a:.3f} {b:.2f} {c} {d:.3f} {e:.3f} {f:.3f} {g:.1f} {h}" for a,b,c,d,e,f,g,h in zip(z,T,src,lf,ly,lc,snr,idx))+"\n")
# angular sizes from the temperature (R500 = 1148 kpc (T/5 keV)^0.51 today; Vikhlinin et al. 2009), in arcmin
R0=1148*(T/5)**0.51
th_exp=R0/Ez(z)/(DA(z)*1000)*180*60/np.pi      # expanding: clusters of a given temperature are smaller by 1/E(z), distance = standard D_A
th_sta=R0/(d_static(z)*1000)*180*60/np.pi      # static: same size at all times, distance = (c/H0) ln(1+z)
print(f"angular size at fixed temperature, static / expanding: min {np.min(th_sta/th_exp):.3f}, max {np.max(th_sta/th_exp):.3f}  (the two models predict nearly the same size on the sky)")
def ols(X,y):
    b=np.linalg.lstsq(X,y,rcond=None)[0]; r=y-X@b; s2=r@r/(len(y)-X.shape[1]); cov=s2*np.linalg.pinv(X.T@X); return b,np.sqrt(np.diag(cov)),np.sqrt(s2)
def design(mask,offsets=True):
    ss=[s for s in prio if ((src==s)&mask).sum()>0]
    if offsets and len(ss)>1: D=np.array([(src[mask]==s).astype(float) for s in ss]).T
    else: D=np.ones((mask.sum(),1))
    return np.c_[D,lT[mask],lE[mask]]
def fit(y,mask=None,offsets=True):
    mask=np.ones(n,bool) if mask is None else mask; X=design(mask,offsets); b,e,rms=ols(X,y[mask]); return dict(n=int(mask.sum()),alpha=float(b[-2]),alpha_se=float(e[-2]),m=float(b[-1]),m_se=float(e[-1]),rms=float(rms))
def show(label,r,extra=""): print(f"  {label:60s} n={r['n']:3d}  temperature power {r['alpha']:.2f}+-{r['alpha_se']:.2f}  m = {r['m']:+.2f}+-{r['m_se']:.2f} {extra}")
out={"n":n,"sources":{s:int((src==s).sum()) for s in prio if (src==s).any()},"n_hi":int((z>=0.45).sum()),"n_hi7":int((z>=0.7).sum()),"zmax":float(z.max()),"size_ratio":[float(np.min(th_sta/th_exp)),float(np.max(th_sta/th_exp))]}
y_exp=lf-logQ(th_exp); y_sta=lf-logQ(th_sta)
print("\nfit: log(signal) = constant per temperature catalogue + alpha log T + m log E(z)")
print("1. raw signal through the fixed filter, no size correction:")
out["raw"]=fit(lf); show("observed",out["raw"])
# what each model predicts for the raw signal, cluster by cluster, then fitted the same way
a_obs=out["raw"]["alpha"]; pe=fit(1.62*lT+0.92*lE+logQ(th_exp)); ps=fit(1.62*lT+logQ(th_sta)); out["raw_pred_exp"]=pe; out["raw_pred_sta"]=ps
print(f"     predicted for the raw signal: expanding, denser clusters m = {pe['m']:+.2f};  static, unchanged clusters m = {ps['m']:+.2f}")
print("2. signal corrected for angular size (expanding model expects m = 0.92; static model expects m = 0):")
out["exp"]=fit(y_exp); show("size from the temperature, expanding model",out["exp"])
out["sta"]=fit(y_sta); show("size from the temperature, static model",out["sta"])
out["act"]=fit(ly); show("size as ACT assumed it (first version's quantity)",out["act"])
out["opt"]=fit(lc); show("ACT's best-fitting filter size, no model for the size",out["opt"])
out["exp_single"]=fit(y_exp,offsets=False); show("expanding sizes, one constant for all catalogues",out["exp_single"])
out["sta_single"]=fit(y_sta,offsets=False); show("static sizes, one constant for all catalogues",out["sta_single"])
print("3. subsets (static-model sizes; the static model expects m = 0 in every row):")
sub={}
for key,label,mask in [("sanders","Sanders et al. only (one Chandra analysis, z 0.28-1.29)",src=="S"),("chandra","Chandra temperatures only (S, M, A, GC)",np.isin(src,["S","M","A","GC"])),("xmm","XMM temperatures only (L, XC, XS, GX)",np.isin(src,["L","XC","XS","GX"])),
                       ("noR","without the two literature compilations (R, B)",~np.isin(src,["R","B"])),("T5","kT >= 5 keV",T>=5),("snr8","ACT signal-to-noise >= 8",snr>=8),("z01","z >= 0.1",z>=0.1),("zlt1","z < 1",z<1),("zlt06","z < 0.6",z<0.6)]:
    sub[key]=dict(sta=fit(y_sta,mask),exp=fit(y_exp,mask)); show(label,sub[key]["sta"],f"(expanding sizes: {sub[key]['exp']['m']:+.2f}+-{sub[key]['exp']['m_se']:.2f})")
out["subsets"]=sub
m,se=out["sta"]["m"],out["sta"]["m_se"]; me,see=out["exp"]["m"],out["exp"]["m_se"]
print(f"\nstatic reading:    m = {m:+.2f}+-{se:.2f}, expected 0:    {abs(m)/se:.1f} sigma away")
print(f"expanding reading: m = {me:+.2f}+-{see:.2f}, expected 0.92: {abs(me-0.92)/see:.1f} sigma away")
rng=np.random.default_rng(11); bs=[]; X=design(np.ones(n,bool))
for _ in range(3000):
    j=rng.integers(0,n,n); bs.append(np.linalg.lstsq(X[j],y_sta[j],rcond=None)[0][-1])
bs=np.array(bs); print(f"bootstrap (static sizes): m = {bs.mean():+.2f}+-{bs.std():.2f}; share of resamples with m <= 0: {np.mean(bs<=0):.5f}"); out["boot"]=[float(bs.mean()),float(bs.std()),float(np.mean(bs<=0))]
# how wrong would the sizes have to be to bring the static reading to m = 0?
for f in (0.8,0.7,0.6,0.5):
    r=fit(lf-logQ(th_sta*f**((z-z.min())/(1.0-z.min())))); out.setdefault("shrink",{})[str(f)]=r["m"]
print("if static-model clusters at z = 1 were really smaller on the sky than assumed by a factor f (growing smoothly with z): m =",{k:round(v,2) for k,v in out["shrink"].items()})
# binned residuals at alpha = fitted, static sizes, for the figure
b,e,_=ols(design(np.ones(n,bool)),y_sta); Dm=design(np.ones(n,bool)); res=y_sta-Dm[:,:-1]@b[:-1]     # remove catalogue constants and temperature term only
edges=[0.03,0.12,0.22,0.32,0.45,0.6,0.8,1.3]; bins=[]
print("binned signal at fixed temperature (static sizes), relative to the lowest-redshift bin:")
for lo,hi in zip(edges,edges[1:]):
    k=(z>=lo)&(z<hi); bins.append([float(np.median(z[k])),int(k.sum()),float(res[k].mean()),float(res[k].std(ddof=1)/np.sqrt(k.sum())),float(lE[k].mean())])
for bn in bins: print(f"   z~{bn[0]:.2f} n={bn[1]:2d}  {bn[2]-bins[0][2]:+.3f} +- {bn[3]:.3f}   expanding prediction {0.92*(bn[4]-bins[0][4]):+.3f}, static prediction +0.000")
out["bins"]=bins
zz=np.linspace(0.03,1.3,52); out["curve"]={"z":zz.tolist(),"logE":np.log10(Ez(zz)).tolist()}
out["qcurve"]=q[:,:2].tolist()
json.dump(out,open("cluster_density_test.json","w"),indent=1,default=float)
