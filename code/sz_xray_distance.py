"""Workbench Entry 27: the shadow-plus-X-ray distance test.
A cluster's hot gas is seen twice. Its shadow on the CMB (ACT DR5, Hilton et al. 2021) measures the pressure summed along the line of sight
and does not care how far away the cluster is. Its X-ray glow (Chandra profiles of SPT clusters, Sanders et al. 2018) measures density squared
times depth, and turning that into a pressure needs the distance and the rule for how surface brightness fades with redshift.
Sanders et al. assumed the standard expanding model. If that model is right, shadow / X-ray prediction is the same at every redshift.
If the universe is static and light simply tires, the same X-ray data imply less gas at high redshift, and the ratio computed with the standard
model must fall as F(z) = sqrt(d_static/D_A) * (1+z)^-1.5.
Inputs (same folder, or a folder given as the first argument): sanders_act_profiles.txt, act_signals.txt, temps_matched.txt, written from the browser session by raw_chunks.py."""
import numpy as np, json, sys
from scipy import integrate, optimize
dd=sys.argv[1] if len(sys.argv)>1 else "."
H0=70.; cH=299792.458/H0
Ez=lambda z: np.sqrt(0.3*(1+z)**3+0.7)
def DA(z): return cH*integrate.quad(lambda x:1/Ez(x),0,z)[0]/(1+z)
DA=np.vectorize(DA)
d_log=lambda z: cH*np.log(1+z)          # static, tired light, distance = (c/H0) ln(1+z)
d_lin=lambda z: cH*z                    # static, Hubble law taken as exact at all z (Lerner et al. 2014)
# ---------- universal pressure profile (Arnaud et al. 2010) shape integrals
c5,ga,al,be=1.177,0.3081,1.0510,5.4905
p=lambda x: 1/((c5*x)**ga*(1+(c5*x)**al)**((be-ga)/al))
def los(b,xmax=5.0):
    if b>=xmax: return 0.0
    return 2*integrate.quad(lambda l:p(np.hypot(b,l)),0,np.sqrt(xmax**2-b**2),limit=200)[0]
I0=los(1e-6,5.0)
Ksph=lambda a: 4*np.pi*integrate.quad(lambda x:p(x)*x*x,0,a)[0]/I0
def Kcyl(a,xmax):
    return integrate.quad(lambda b:2*np.pi*b*los(b,xmax),0,a,limit=200)[0]/I0
print(f"profile shape: central line integral {I0:.3f} (x R500 x P500); sphere/(y0 R500^2) = {Ksph(1):.4f} inside R500, {Ksph(0.5):.4f} inside half of it")
# ---------- data
P=[l.split() for l in open(dd+"/sanders_act_profiles.txt") if not l.startswith("#")]
S={l.split()[0]:l.split() for l in open(dd+"/act_signals.txt") if not l.startswith("#")}
T={l.split()[1]:l.split() for l in open(dd+"/temps_matched.txt") if l.startswith("S ")}
k=1.3018e-27*3.0857e21   # sigma_T/(m_e c^2) in cm^2/keV times cm per kpc: (keV cm^-3 kpc^n) -> kpc^(n-1)
rows=[]
for r in P:
    i=r[0]; s=S[i]; t=T[i]; z=float(t[2]); E=Ez(z); Mu=float(s[3]); Mc=float(s[4])
    y0=4.95e-5*E*E*(Mu/3)**1.08; y0c=4.95e-5*E*E*(Mc/3)**1.08         # central Compton y of the template ACT fitted (uncorrected / corrected mass)
    R=float(r[1]); Rc=float(r[2]); am=float(r[13]); rmax=float(r[14]); th5=R/am
    X=dict(yc=k*float(r[3]),Y5=k*float(r[4]),Y1=k*float(r[5]),err=float(r[6]),cov=float(r[7]),Y1b=k*float(r[8]),Y1n=k*float(r[9]),ycn=k*float(r[10]),Y1c=k*float(r[11]),Y5c=k*float(r[12]),
           C1=k*float(r[15]),C2=k*float(r[16]),CR=k*float(r[17]))
    xm=rmax/R
    A=dict(Y1=y0*R*R*Ksph(1),Y5=y0*R*R*Ksph(0.5),Y1c=y0c*Rc*Rc*Ksph(1),C1=y0*R*R*Kcyl(am/R,xm),C2=y0*R*R*Kcyl(2*am/R,xm),CR=y0*R*R*Kcyl(1.0,xm),yc=y0)
    rows.append(dict(i=int(i),name=r[20],z=z,E=E,T=float(t[3]),eT=float(t[4]),snr=float(t[6]),Mu=Mu,Mc=Mc,th5=th5,R=R,Rh=(float(r[18]) if r[18]!="-" else np.nan),Mg=float(r[19]),esz=float(s[5]),X=X,A=A,cov=X["cov"],
                     lf=float(s[1]),lc=float(s[2])))
n=len(rows); z=np.array([r["z"] for r in rows]); lz=np.log10(1+z); Tk=np.array([r["T"] for r in rows]); snr=np.array([r["snr"] for r in rows])
Fs=np.sqrt(d_log(z)/DA(z))*(1+z)**-1.5; Fl=np.sqrt(d_lin(z)/DA(z))*(1+z)**-1.5
Ft=np.sqrt(d_log(z)/DA(z))*(1+z)**-1.0     # static space in which photons also arrive stretched in time (as supernova light curves require)
print(f"{n} clusters, z {z.min():.2f}-{z.max():.2f} (median {np.median(z):.2f}); kT median {np.median(Tk):.1f} keV; ACT S/N median {np.median(snr):.1f}")
print("static prediction for the ratio, F(z) at z = 0.3, 0.6, 1.0, 1.3:",[round(float(np.sqrt(d_log(x)/DA(x))*(1+x)**-1.5),3) for x in (0.3,0.6,1.0,1.3)],"(log distance);",[round(float(np.sqrt(d_lin(x)/DA(x))*(1+x)**-1.5),3) for x in (0.3,0.6,1.0,1.3)],"(linear)")
def ols(Xm,y,w=None):
    Xm=np.asarray(Xm,float); y=np.asarray(y,float)
    if w is None: w=np.ones(len(y))
    W=np.sqrt(w)[:,None]; b=np.linalg.lstsq(Xm*W,y*W[:,0],rcond=None)[0]; r=y-Xm@b; s2=(w*r*r).sum()/(len(y)-Xm.shape[1]); cov=s2*np.linalg.inv((Xm*W).T@(Xm*W)); return b,np.sqrt(np.diag(cov)),np.sqrt((r*r).mean())
def trend(lr,label,mask=None,quiet=False):
    m=np.isfinite(lr) if mask is None else (mask&np.isfinite(lr))
    b,e,rms=ols(np.c_[np.ones(m.sum()),lz[m]],lr[m]); bk,ek,_=ols(np.c_[np.ones(m.sum()),np.log10(Fs[m])],lr[m]); bl,el,_=ols(np.c_[np.ones(m.sum()),np.log10(Fl[m])],lr[m])
    med=10**np.median(lr[m])
    if not quiet: print(f"  {label:58s} n={m.sum():2d} median ratio {med:.2f} scatter {rms:.3f} dex | power of (1+z): {b[1]:+.2f}+-{e[1]:.2f} | k = {bk[1]:+.2f}+-{ek[1]:.2f} (0 expanding, 1 static) | linear-distance static k = {bl[1]:+.2f}+-{el[1]:.2f}")
    return dict(n=int(m.sum()),median=float(med),rms=float(rms),p=float(b[1]),p_se=float(e[1]),k=float(bk[1]),k_se=float(ek[1]),klin=float(bl[1]),klin_se=float(el[1]))
g=lambda key,side: np.array([r[side][key] for r in rows])
out={"n":n,"zrange":[float(z.min()),float(z.max())],"shape":{"I0":I0,"Ksph1":Ksph(1),"Ksph05":Ksph(0.5)}}
pz=ols(np.c_[np.ones(n),lz],np.log10(Fs))[0][1]; pzl=ols(np.c_[np.ones(n),lz],np.log10(Fl))[0][1]
print(f"over these redshifts the static prediction is close to a power law: ratio ~ (1+z)^{pz:.2f} (log distance), (1+z)^{pzl:.2f} (linear distance); expanding: (1+z)^0")
out["static_power"]=[float(pz),float(pzl)]
print("\nratio = shadow (ACT) / X-ray prediction (Sanders et al. profiles, standard distances):")
res={}
res["Y1"]=trend(np.log10(g("Y1","A")/g("Y1","X")),"pressure summed inside R500 (sphere)")
res["Y5"]=trend(np.log10(g("Y5","A")/g("Y5","X")),"pressure summed inside half of R500")
res["C2"]=trend(np.log10(g("C2","A")/g("C2","X")),"cylinder of radius 2 arcmin on the sky")
res["C1"]=trend(np.log10(g("C1","A")/g("C1","X")),"cylinder of radius 1 arcmin")
res["CR"]=trend(np.log10(g("CR","A")/g("CR","X")),"cylinder of radius R500")
res["yc"]=trend(np.log10(g("yc","A")/g("yc","X")),"central line of sight only")
res["Y1_mbeta"]=trend(np.log10(g("Y1","A")/g("Y1b","X")),"inside R500, second X-ray model (modified beta)")
res["Y1_free"]=trend(np.log10(g("Y1","A")/g("Y1n","X")),"inside R500, X-ray temperatures free in each shell")
res["Y1_corr"]=trend(np.log10(g("Y1c","A")/g("Y1c","X")),"inside R500, ACT signal corrected for selection boost")
base=np.log10(g("Y1","A")/g("Y1","X"))
print("cuts on the main version (pressure inside R500):")
res["snr8"]=trend(base,"ACT signal-to-noise >= 8",snr>=8); res["snr10"]=trend(base,"ACT signal-to-noise >= 10",snr>=10)
res["err25"]=trend(base,"X-ray pressure half-range below 25%",np.array([r["X"]["err"] for r in rows])<0.25)
res["T6"]=trend(base,"kT >= 6 keV",Tk>=6); res["noEG"]=trend(base,"without the five largest signals",snr<30)
res["zlt1"]=trend(base,"z < 1",z<1); res["zgt04"]=trend(base,"z > 0.4",z>0.4)
# temperature and mass as extra terms
b,e,_=ols(np.c_[np.ones(n),np.log10(Fs),np.log10(Tk/6)],base); print(f"with a temperature term: k = {b[1]:+.2f}+-{e[1]:.2f}, temperature power {b[2]:+.2f}+-{e[2]:.2f}"); res["withT"]=dict(k=float(b[1]),k_se=float(e[1]),a=float(b[2]),a_se=float(e[2]))
Mu=np.array([r["Mu"] for r in rows]); b,e,_=ols(np.c_[np.ones(n),np.log10(Fs),np.log10(Mu/4)],base); print(f"with an ACT-mass term:   k = {b[1]:+.2f}+-{e[1]:.2f}, mass power {b[2]:+.2f}+-{e[2]:.2f}"); res["withM"]=dict(k=float(b[1]),k_se=float(e[1]),a=float(b[2]),a_se=float(e[2]))
# weighted by measurement error (shadow error and X-ray half-range in quadrature, plus 0.06 dex intrinsic)
sig=np.sqrt((np.array([r["esz"] for r in rows])*1.08/np.log(10))**2+(np.array([r["X"]["err"] for r in rows])/np.log(10))**2+0.06**2); b,e,_=ols(np.c_[np.ones(n),np.log10(Fs)],base,1/sig**2); print(f"weighted by measurement error: k = {b[1]:+.2f}+-{e[1]:.2f}"); res["weighted"]=dict(k=float(b[1]),k_se=float(e[1]))
rng=np.random.default_rng(5); bs=[]
for _ in range(4000):
    j=rng.integers(0,n,n); bs.append(ols(np.c_[np.ones(n),np.log10(Fs[j])],base[j])[0][1])
bs=np.array(bs); print(f"bootstrap: k = {bs.mean():+.2f}+-{bs.std():.2f}; share of resamples with k >= 1: {np.mean(bs>=1):.5f}; with k >= 0.5: {np.mean(bs>=0.5):.5f}")
res["boot"]=[float(bs.mean()),float(bs.std()),float(np.mean(bs>=1)),float(np.mean(bs>=0.5))]
bt,et,_=ols(np.c_[np.ones(n),np.log10(Ft)],base); ptd=ols(np.c_[np.ones(n),lz],np.log10(Ft))[0][1]
print(f"static space with time-stretched photons predicts ratio ~ (1+z)^{ptd:.2f}; fitted k for that model = {bt[1]:+.2f}+-{et[1]:.2f}: {abs(1-bt[1])/et[1]:.1f} sigma from 1")
res["timedil"]=dict(power=float(ptd),k=float(bt[1]),k_se=float(et[1]))
allk=[res[v]["k"] for v in ("Y1","Y5","C2","C1","CR","yc","Y1_mbeta","Y1_free","Y1_corr")]; print(f"nine versions of the comparison: k from {min(allk):+.2f} to {max(allk):+.2f}; closest to the static value: {max(allk):+.2f}")
res["k_range"]=[float(min(allk)),float(max(allk))]
bl=np.array([ols(np.c_[np.ones(n),np.log10(Fl[j])],base[j])[0][1] for j in (rng.integers(0,n,n) for _ in range(4000))]); print(f"bootstrap, linear-distance static: k = {bl.mean():+.2f}+-{bl.std():.2f}; share >= 1: {np.mean(bl>=1):.5f}"); res["boot_lin"]=[float(bl.mean()),float(bl.std()),float(np.mean(bl>=1))]
# the radio fog of Entry 22 would also thin the shadows (X-rays pass through it): shadow x exp(-tau(z)), Lerner's stated optical depths
tau=lambda zz: np.where(zz<=0.14,0.395*np.log(4282.7*np.maximum(zz,1e-3)),2.526+(7.6-2.526)*(zz-0.14)/(2.2-0.14))
bf,ef,_=ols(np.c_[np.ones(n),-tau(z)/np.log(10)],base); bsf=np.array([ols(np.c_[np.ones(n),-tau(z[j])/np.log(10)],base[j])[0][1] for j in (rng.integers(0,n,n) for _ in range(2000))])
print(f"fog strength from these shadows (1 = as stated, 0 = none): {bf[1]:+.3f}+-{ef[1]:.3f} (bootstrap +-{bsf.std():.3f}); stated fog would thin the shadow {np.exp(tau(1.29)-tau(0.28)):.0f}-fold between z = 0.28 and 1.29")
res["fog"]=dict(s=float(bf[1]),s_se=float(ef[1]),s_boot=float(bsf.std()),cut=float(np.exp(tau(1.29)-tau(0.28))))
K=res["Y1"]; print(f"\nmain result: k = {K['k']:+.2f}+-{K['k_se']:.2f}: {abs(K['k'])/K['k_se']:.1f} sigma from the expanding value 0, {abs(1-K['k'])/K['k_se']:.1f} sigma from the static value 1")
out["ratios"]=res
# binned ratio for the figure
edges=[0.25,0.38,0.5,0.62,0.78,1.0,1.35]; bins=[]
print("binned (pressure inside R500):  z range, n, mean log ratio, error, static prediction relative to the lowest bin")
for lo,hi in zip(edges,edges[1:]):
    m=(z>=lo)&(z<hi); bins.append([float(np.median(z[m])),int(m.sum()),float(base[m].mean()),float(base[m].std(ddof=1)/np.sqrt(m.sum())),float(np.log10(Fs[m]).mean()),float(np.log10(Fl[m]).mean())])
for bn in bins: print(f"   z~{bn[0]:.2f} n={bn[1]:2d}  {bn[2]:+.3f} +- {bn[3]:.3f}   static: {bn[4]-bins[0][4]:+.3f} (log), {bn[5]-bins[0][5]:+.3f} (linear)")
out["bins"]=bins
# ---------- is the template size ACT assumed the size the X-ray gas actually has?  (pure angles, no distances)
sx=np.log10(g("Y5","X")/g("Y1","X")); sa=np.log10(Ksph(0.5)/Ksph(1)); b,e,rms=ols(np.c_[np.ones(n),lz],sx-sa)
print(f"\nshape check: X-ray share of the R500 pressure that lies inside half of R500, relative to the template: median {10**np.median(sx-sa):.2f}, trend with (1+z): power {b[1]:+.2f}+-{e[1]:.2f}")
cx=np.log10(g("C1","X")/g("C2","X")); ca=np.log10(g("C1","A")/g("C2","A")); b2,e2,_=ols(np.c_[np.ones(n),lz],cx-ca)
print(f"             1-arcmin / 2-arcmin cylinder, X-ray relative to template: median {10**np.median(cx-ca):.2f}, trend with (1+z): power {b2[1]:+.2f}+-{e2[1]:.2f}")
Rh=np.array([r["Rh"] for r in rows]); R5=np.array([r["R"] for r in rows]); ok=np.isfinite(Rh); b3,e3,_=ols(np.c_[np.ones(ok.sum()),lz[ok]],np.log10(Rh[ok]/R5[ok]))
print(f"             X-ray (hydrostatic) R500 / ACT R500: median {np.median(Rh[ok]/R5[ok]):.2f}, trend with (1+z): power {b3[1]:+.2f}+-{e3[1]:.2f}  ({ok.sum()} clusters)")
out["shape_check"]=dict(inner_share_ratio=float(10**np.median(sx-sa)),inner_trend=[float(b[1]),float(e[1])],cyl_ratio=float(10**np.median(cx-ca)),cyl_trend=[float(b2[1]),float(e2[1])],R500_ratio=float(np.median(Rh[ok]/R5[ok])),R500_trend=[float(b3[1]),float(e3[1])])
# ---------- what the static reading needs instead: gas clumping or a missing-gas factor that changes with redshift
print("\nif the universe were static, the X-ray analysis redone with static distances gives shadow/X-ray =",[round(float(10**(K0)),2) for K0 in (np.median(base[z<0.45])-np.log10(np.median(Fs[z<0.45])),np.median(base[z>0.75])-np.log10(np.median(Fs[z>0.75])))],"at z<0.45 and z>0.75: the shadow would be deeper than the gas seen in X-rays can cast")
out["static_ratio_lowz_highz"]=[float(10**(np.median(base[z<0.45])-np.log10(np.median(Fs[z<0.45])))),float(10**(np.median(base[z>0.75])-np.log10(np.median(Fs[z>0.75]))))]
out["points"]=[[r["z"],float(base[i]),r["snr"],r["name"]] for i,r in enumerate(rows)]
zz=np.linspace(0.25,1.32,44); out["curve"]={"z":zz.tolist(),"static_log":np.log10(np.sqrt(d_log(zz)/DA(zz))*(1+zz)**-1.5).tolist(),"static_lin":np.log10(np.sqrt(d_lin(zz)/DA(zz))*(1+zz)**-1.5).tolist()}
json.dump(out,open("sz_xray_distance.json","w"),indent=1,default=float)
