"""Workbench Entry 28: the shadow-against-glow distance test, repeated with two other pairs of instruments.
Entry 27 compared ACT shadows with Chandra profiles (Sanders et al. 2018) for 66 clusters. Its weak points were one shadow catalogue, whose
size correction leans on the standard model, and a spread between ways of adding up the outer gas. Two repeats:
 A. South Pole Telescope shadows (Bleem et al. 2015: Compton Y inside a fixed 0.75 arcmin circle, from SPT's own maps) against the same Chandra
    profiles, summed in the same 0.75 arcmin column. 83 clusters, z 0.28-1.22. No template size from the standard model enters, and the inner
    gas is used, where the X-ray profiles are best measured.
 B. Planck shadows (PSZ2) against XMM-Newton gas masses and temperatures (Lovisari et al. 2020). 117 clusters, z 0.06-0.55. Different shadow
    telescope, different X-ray telescope, different clusters.
In each case the X-ray prediction was made with standard distances. A static universe with tired light needs measured/predicted to fall as
F(z) = sqrt(d_static/D_A) (1+z)^-1.5; an expanding one needs it level.  k = fitted power of F: 0 expanding, 1 static."""
import numpy as np, json, sys
from scipy import integrate
dd=sys.argv[1] if len(sys.argv)>1 else "."
cH=299792.458/70.; Ez=lambda z: np.sqrt(0.3*(1+z)**3+0.7)
DA=np.vectorize(lambda z: cH*integrate.quad(lambda x:1/Ez(x),0,z)[0]/(1+z))
F_log=lambda z: np.sqrt(cH*np.log(1+z)/DA(z))*(1+z)**-1.5      # static, tired light, d = (c/H0) ln(1+z)
F_lin=lambda z: np.sqrt(cH*z/DA(z))*(1+z)**-1.5                # static, Lerner's d = cz/H0
F_td =lambda z: np.sqrt(cH*np.log(1+z)/DA(z))*(1+z)**-1.0      # static, photons also stretched in time
def ols(X,y,w=None):
    X=np.asarray(X,float); y=np.asarray(y,float); w=np.ones(len(y)) if w is None else w; W=np.sqrt(w)[:,None]
    b=np.linalg.lstsq(X*W,y*W[:,0],rcond=None)[0]; r=y-X@b; s2=(w*r*r).sum()/(len(y)-X.shape[1]); return b,np.sqrt(np.diag(s2*np.linalg.inv((X*W).T@(X*W)))),np.sqrt((r*r).mean())
rng=np.random.default_rng(7)
def kfit(z,lr,label,mask=None,quiet=False,nboot=4000):
    m=np.isfinite(lr) if mask is None else mask&np.isfinite(lr); zz=z[m]; y=lr[m]; n=m.sum(); o=dict(n=int(n),median=float(10**np.median(y)))
    b,e,rms=ols(np.c_[np.ones(n),np.log10(1+zz)],y); o.update(p=float(b[1]),p_se=float(e[1]),rms=float(rms))
    for name,F in (("k",F_log),("klin",F_lin),("ktd",F_td)):
        x=np.log10(F(zz)); b,e,_=ols(np.c_[np.ones(n),x],y); bs=[]
        if nboot:
            for _ in range(nboot):
                j=rng.integers(0,n,n); bs.append(ols(np.c_[np.ones(n),x[j]],y[j])[0][1])
            bs=np.array(bs); o[name+"_boot"]=float(bs.std()); o[name+"_ge1"]=float(np.mean(bs>=1))
        o[name]=float(b[1]); o[name+"_se"]=float(e[1])
    if not quiet: print(f"  {label:54s} n={n:3d} median {o['median']:.2f} scatter {o['rms']:.3f} dex | (1+z)^{o['p']:+.2f}+-{o['p_se']:.2f} | k = {o['k']:+.2f}+-{o['k_se']:.2f}"+(f" (boot {o['k_boot']:.2f}; {o['k_ge1']*100:.2f}% reach 1)" if nboot else "")+f" | Lerner law k = {o['klin']:+.2f}+-{o['klin_se']:.2f} | stretched time k = {o['ktd']:+.2f}+-{o['ktd_se']:.2f}")
    return o
out={}
# ---------- universal pressure profile, used only to allow for the X-ray profiles stopping at a finite radius
c5,ga,al,be=1.177,0.3081,1.0510,5.4905
p=lambda x: 1/((c5*x)**ga*(1+(c5*x)**al)**((be-ga)/al))
def los(b,xmax): return 0.0 if b>=xmax else 2*integrate.quad(lambda l:p(np.hypot(b,l)),0,np.sqrt(xmax**2-b**2),limit=200)[0]
Kcyl=lambda a,xmax: integrate.quad(lambda b:2*np.pi*b*los(b,xmax),0,a,limit=200)[0]
# ================= A. SPT shadows against Chandra profiles
A=[l.split() for l in open(dd+"/sanders_spt.txt") if not l.startswith("#")]
name=np.array([r[0] for r in A]); z=np.array([float(r[1]) for r in A]); xi=np.array([float(r[2]) for r in A]); thc=np.array([float(r[3]) for r in A])
Ysz=np.array([float(r[4]) for r in A])*1e-6; eY=np.array([float(r[5]) for r in A])*1e-6; zb=np.array([float(r[6]) for r in A]); am=np.array([float(r[7]) for r in A]); rmax=np.array([float(r[8]) for r in A])
kk=1.3018e-27*3.0857e21; Yx=kk*np.array([float(r[9]) for r in A])/am**2; ex=np.array([float(r[10]) for r in A]); Yxm=kk*np.array([float(r[11]) for r in A])/am**2; Yxb=kk*np.array([float(r[12]) for r in A])/am**2; T=np.array([float(r[13]) for r in A])
R500=1148*(T/5)**0.51/Ez(z); a=0.75*am/R500; xm=rmax/R500; corr=np.array([Kcyl(ai,8.0)/Kcyl(ai,xi_) for ai,xi_ in zip(a,xm)])
print(f"A. SPT against Chandra: {len(z)} clusters, z {z.min():.2f}-{z.max():.2f} (median {np.median(z):.2f}); 0.75 arcmin is {np.median(a):.2f} R500 (range {a.min():.2f}-{a.max():.2f}); allowance for the X-ray profiles ending at {np.median(xm):.1f} R500: +{100*(np.median(corr)-1):.1f}% (max +{100*(corr.max()-1):.1f}%)")
lr=np.log10(Ysz/(Yx*corr)); zok=np.abs(z-zb)<=0.03*(1+z)
print(f"   redshift revised by more than 3% of (1+z) since the X-ray analysis: {int((~zok).sum())} clusters ({', '.join(name[~zok])})")
resA={}
resA["all"]=kfit(z,lr,"all 83")
resA["zok"]=kfit(z,lr,"redshift not revised since",zok)
resA["nocorr"]=kfit(z,np.log10(Ysz/Yx),"no allowance for the profile ending",zok,nboot=0)
resA["mbeta"]=kfit(z,np.log10(Ysz/(Yxm*corr)),"second X-ray model",zok,nboot=0)
resA["free"]=kfit(z,np.log10(Ysz/(Yxb*corr)),"X-ray temperature free in each shell",zok,nboot=0)
resA["xi8"]=kfit(z,lr,"SPT significance of 8 or more",zok&(xi>=8),nboot=0)
resA["xi10"]=kfit(z,lr,"SPT significance of 10 or more",zok&(xi>=10),nboot=0)
resA["err25"]=kfit(z,lr,"X-ray half-range below 25%",zok&(ex<0.25),nboot=0)
resA["core05"]=kfit(z,lr,"SPT core radius 0.5 arcmin or less",zok&(thc<=0.5),nboot=0)
resA["core75"]=kfit(z,lr,"SPT core radius 0.75 arcmin or more",zok&(thc>=0.75),nboot=0)
resA["zgt04"]=kfit(z,lr,"z > 0.4",zok&(z>0.4),nboot=0)
m=zok; b,e,_=ols(np.c_[np.ones(m.sum()),np.log10(F_log(z[m])),np.log10(T[m]/6)],lr[m]); print(f"   with a temperature term: k = {b[1]:+.2f}+-{e[1]:.2f}; temperature power {b[2]:+.2f}+-{e[2]:.2f}"); resA["withT"]=dict(k=float(b[1]),k_se=float(e[1]),a=float(b[2]),a_se=float(e[2]))
b,e,_=ols(np.c_[np.ones(m.sum()),np.log10(F_log(z[m])),np.log10(xi[m]/9)],lr[m]); print(f"   with a significance term: k = {b[1]:+.2f}+-{e[1]:.2f}; significance power {b[2]:+.2f}+-{e[2]:.2f}"); resA["withXi"]=dict(k=float(b[1]),k_se=float(e[1]),a=float(b[2]),a_se=float(e[2]))
sig=np.sqrt((eY/Ysz/np.log(10))**2+(ex/np.log(10))**2+0.05**2); b,e,_=ols(np.c_[np.ones(m.sum()),np.log10(F_log(z[m]))],lr[m],1/sig[m]**2); print(f"   weighted by measurement error: k = {b[1]:+.2f}+-{e[1]:.2f}"); resA["weighted"]=dict(k=float(b[1]),k_se=float(e[1]))
edges=[0.25,0.38,0.5,0.62,0.78,1.0,1.3]; binsA=[]
for lo,hi in zip(edges,edges[1:]):
    q=m&(z>=lo)&(z<hi); binsA.append([float(np.median(z[q])),int(q.sum()),float(lr[q].mean()),float(lr[q].std(ddof=1)/np.sqrt(q.sum())),float(np.log10(F_log(z[q])).mean()),float(np.log10(F_lin(z[q])).mean())])
print("   binned:",", ".join(f"z~{b[0]:.2f} (n={b[1]}) {b[2]:+.3f}+-{b[3]:.3f}" for b in binsA),f"| static predicts a change of {binsA[-1][4]-binsA[0][4]:+.3f} from first to last; measured {binsA[-1][2]-binsA[0][2]:+.3f}+-{np.hypot(binsA[-1][3],binsA[0][3]):.3f}")
resA["bins"]=binsA; out["A"]=resA
ptsA=[[float(a_),float(b_),n_,float(c_)] for a_,b_,n_,c_ in zip(z[m],lr[m],name[m],xi[m])]
# ================= B. Planck shadows against XMM gas
B=np.array([[float(x) if x!="-" else np.nan for x in l.split()] for l in open(dd+"/planck_lovisari.txt") if not l.startswith("#")])
zB,Msz,Y5,eY5,snr,Mx,Mg,eMg,Tin,Tex,eTex,Lx=B.T
Tuse=np.where(np.isfinite(Tex),Tex,Tin); good=(eTex/Tuse<0.3)|~np.isfinite(eTex); good&=(eMg/Mg<0.3)
Ysz_B=1e-4*10**-0.19*(Msz/6)**1.79*Ez(zB)**(2/3)                 # Mpc^2, the Planck signal inside R500 as encoded in its mass proxy
Yx_B=1.4262e-5*(Mg/10)*Tuse                                      # Mpc^2, gas mass x temperature inside the XMM R500
Y5phys=(Y5*1e-3/1.79)*(np.pi/180/60)**2*DA(zB)**2                # Mpc^2, from the blind Y5R500 column (size not fixed by the X-rays)
lrB=np.log10(Ysz_B/Yx_B)
print(f"\nB. Planck against XMM: {len(zB)} clusters matched, {int(good.sum())} with gas mass and temperature good to 30%; z {zB[good].min():.3f}-{zB[good].max():.3f} (median {np.median(zB[good]):.2f}); {int((zB[good]>0.3).sum())} beyond z = 0.3")
resB={}
resB["all"]=kfit(zB,lrB,"all with good X-ray values",good)
resB["Tin"]=kfit(zB,np.log10(Ysz_B/(1.4262e-5*(Mg/10)*Tin)),"temperature with the core included",good,nboot=0)
resB["blind"]=kfit(zB,np.log10(Y5phys/Yx_B),"Planck's blind Y5R500 instead of its mass proxy",good,nboot=0)
resB["snr10"]=kfit(zB,lrB,"Planck signal-to-noise of 10 or more",good&(snr>=10),nboot=0)
resB["snr12"]=kfit(zB,lrB,"Planck signal-to-noise of 12 or more",good&(snr>=12),nboot=0)
resB["z01"]=kfit(zB,lrB,"z > 0.1",good&(zB>0.1),nboot=0)
resB["hot"]=kfit(zB,lrB,"kT of 6 keV or more",good&(Tuse>=6),nboot=0)
b,e,_=ols(np.c_[np.ones(good.sum()),np.log10(F_log(zB[good])),np.log10(Tuse[good]/6)],lrB[good]); print(f"   with a temperature term: k = {b[1]:+.2f}+-{e[1]:.2f}; temperature power {b[2]:+.2f}+-{e[2]:.2f}"); resB["withT"]=dict(k=float(b[1]),k_se=float(e[1]),a=float(b[2]),a_se=float(e[2]))
edges=[0.05,0.1,0.16,0.21,0.29,0.36,0.56]; binsB=[]
for lo,hi in zip(edges,edges[1:]):
    q=good&(zB>=lo)&(zB<hi); binsB.append([float(np.median(zB[q])),int(q.sum()),float(lrB[q].mean()),float(lrB[q].std(ddof=1)/np.sqrt(q.sum())),float(np.log10(F_log(zB[q])).mean()),float(np.log10(F_lin(zB[q])).mean())])
print("   binned:",", ".join(f"z~{b[0]:.2f} (n={b[1]}) {b[2]:+.3f}+-{b[3]:.3f}" for b in binsB),f"| static predicts a change of {binsB[-1][4]-binsB[0][4]:+.3f} from first to last; measured {binsB[-1][2]-binsB[0][2]:+.3f}+-{np.hypot(binsB[-1][3],binsB[0][3]):.3f}")
resB["bins"]=binsB; out["B"]=resB
ptsB=[[float(a_),float(b_),float(c_)] for a_,b_,c_ in zip(zB[good],lrB[good],snr[good])]
# ================= joint: one k, a separate level for each pair of instruments
zj=np.r_[z[m],zB[good]]; yj=np.r_[lr[m],lrB[good]]; s1=np.r_[np.ones(m.sum()),np.zeros(good.sum())]; s2=1-s1; nj=len(zj); resJ={}
print(f"\njoint fit, SPT+Chandra and Planck+XMM, {nj} clusters, z {zj.min():.2f}-{zj.max():.2f}, a separate level for each pair:")
for nm,F,lab in (("k",F_log,"static, tired light"),("klin",F_lin,"static, Lerner's distance law"),("ktd",F_td,"static, photons also stretched in time")):
    x=np.log10(F(zj)); b,e,_=ols(np.c_[s1,s2,x],yj); bs=[]
    for _ in range(4000):
        j=rng.integers(0,nj,nj); bs.append(ols(np.c_[s1[j],s2[j],x[j]],yj[j])[0][2])
    bs=np.array(bs); resJ[nm]=dict(k=float(b[2]),k_se=float(e[2]),boot=float(bs.std()),ge1=float(np.mean(bs>=1)),sigma=float((1-b[2])/max(e[2],bs.std())))
    print(f"  {lab:40s} k = {b[2]:+.2f}+-{e[2]:.2f} (boot {bs.std():.2f}); {100*np.mean(bs>=1):.3f}% of resamples reach 1; {(1-b[2])/max(e[2],bs.std()):.1f} sigma from the static value, {abs(b[2])/max(e[2],bs.std()):.1f} from expanding")
x=np.log10(1+zj); b,e,_=ols(np.c_[s1,s2,x],yj); print(f"  as a power of (1+z): {b[2]:+.2f}+-{e[2]:.2f}   (expanding 0; static tired light about -1.0; Lerner's law about -0.8; stretched time about -0.5)"); resJ["p"]=[float(b[2]),float(e[2])]
# least favourable cut found: the strongest detections only, in both pairs
hi1=m&(xi>=10); hi2=good&(snr>=12); zh=np.r_[z[hi1],zB[hi2]]; yh=np.r_[lr[hi1],lrB[hi2]]; t1=np.r_[np.ones(hi1.sum()),np.zeros(hi2.sum())]; xh=np.log10(F_log(zh)); b,e,_=ols(np.c_[t1,1-t1,xh],yh)
bs=np.array([ols(np.c_[t1[j],1-t1[j],xh[j]],yh[j])[0][2] for j in (rng.integers(0,len(zh),len(zh)) for _ in range(4000))])
print(f"  strongest detections only (SPT xi >= 10, Planck S/N >= 12), {len(zh)} clusters: k = {b[2]:+.2f}+-{e[2]:.2f} (boot {bs.std():.2f}); {100*np.mean(bs>=1):.2f}% reach 1, {100*np.mean(bs<=0):.2f}% reach 0")
resJ["strong"]=dict(n=int(len(zh)),k=float(b[2]),k_se=float(e[2]),boot=float(bs.std()),ge1=float(np.mean(bs>=1)),le0=float(np.mean(bs<=0)))
Th=np.r_[T[hi1],Tuse[hi2]]; b,e,_=ols(np.c_[t1,1-t1,xh,np.log10(Th/6)],yh); print(f"     the same with a temperature term: k = {b[2]:+.2f}+-{e[2]:.2f}; temperature power {b[3]:+.2f}+-{e[3]:.2f}"); resJ["strongT"]=dict(k=float(b[2]),k_se=float(e[2]),a=float(b[3]),a_se=float(e[3]))
Tj=np.r_[T[m],Tuse[good]]; b,e,_=ols(np.c_[s1,s2,np.log10(F_log(zj)),np.log10(Tj/6)],yj); print(f"  all {nj}, with a temperature term: k = {b[2]:+.2f}+-{e[2]:.2f}; temperature power {b[3]:+.2f}+-{e[3]:.2f}"); resJ["allT"]=dict(k=float(b[2]),k_se=float(e[2]),a=float(b[3]),a_se=float(e[3]))
out["joint"]=resJ
out["pointsA"]=ptsA; out["pointsB"]=ptsB
zz=np.linspace(0.05,1.3,51); out["curve"]={"z":zz.tolist(),"static_log":np.log10(F_log(zz)).tolist(),"static_lin":np.log10(F_lin(zz)).tolist(),"static_td":np.log10(F_td(zz)).tolist()}
json.dump(out,open("distance_test_repeat.json","w"),indent=1)
