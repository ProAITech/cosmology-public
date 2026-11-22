"""Workbench Entry 32: the shadow-against-glow distance test with a third X-ray data set, reaching redshift 1.5.
Entries 27 and 28 compared the Sunyaev-Zel'dovich shadows of galaxy clusters with the shadow predicted from their X-ray gas, using Chandra
profiles (to z = 1.2) and XMM-Newton gas masses (to z = 0.55). Open since then: the distant half rested on one X-ray telescope.
Here: South Pole Telescope shadows (Bleem et al. 2015: Compton Y inside a 0.75 arcmin circle) against XMM-Newton gas masses and temperatures
of 59 SPT clusters at z = 0.20 to 1.48 (Bulbul et al. 2019). 42 of the 59 are not in the Chandra sample of Entry 28.
The X-ray side gives gas mass x temperature inside the sphere of radius R500; a pressure-profile shape turns that into the shadow inside the
0.75 arcmin column. The shape is an assumption, so it is varied.
As before the X-ray prediction is made with standard distances; static tired light needs measured/predicted to fall as
F(z) = sqrt(d_static/D_A) (1+z)^-1.5, an expanding universe needs it level.  k = fitted power of F: 0 expanding, 1 static.
Run: python3 spt_xmm_distance.py data   (data/ holds bulbul_spt.txt and sanders_spt.txt; entry28_points.json, the per-cluster ratios of Entry 28, sits beside this script)."""
import numpy as np, json, sys, os
from scipy import integrate
dd=sys.argv[1] if len(sys.argv)>1 else "."
H0,Om=67.74,0.304                                             # the cosmology in which the X-ray quantities were derived
cH=299792.458/H0; Ez=lambda z: np.sqrt(Om*(1+z)**3+1-Om)
DA=np.vectorize(lambda z: cH*integrate.quad(lambda x:1/Ez(x),0,z)[0]/(1+z))
F_log=lambda z: np.sqrt(cH*np.log(1+z)/DA(z))*(1+z)**-1.5      # static, tired light, d = (c/H0) ln(1+z)
F_lin=lambda z: np.sqrt(cH*z/DA(z))*(1+z)**-1.5                # static, Lerner's d = cz/H0
F_td =lambda z: np.sqrt(cH*np.log(1+z)/DA(z))*(1+z)**-1.0      # static, photons also stretched in time
def ols(X,y,w=None):
    X=np.asarray(X,float); y=np.asarray(y,float); w=np.ones(len(y)) if w is None else w; W=np.sqrt(w)[:,None]
    b=np.linalg.lstsq(X*W,y*W[:,0],rcond=None)[0]; r=y-X@b; s2=(w*r*r).sum()/(len(y)-X.shape[1]); return b,np.sqrt(np.diag(s2*np.linalg.inv((X*W).T@(X*W)))),np.sqrt((r*r).mean())
rng=np.random.default_rng(11)
def kfit(z,lr,label,mask=None,quiet=False,nboot=4000):
    m=np.isfinite(lr) if mask is None else mask&np.isfinite(lr); zz=z[m]; y=lr[m]; n=int(m.sum()); o=dict(n=n,median=float(10**np.median(y)))
    b,e,rms=ols(np.c_[np.ones(n),np.log10(1+zz)],y); o.update(p=float(b[1]),p_se=float(e[1]),rms=float(rms))
    for name,F in (("k",F_log),("klin",F_lin),("ktd",F_td)):
        x=np.log10(F(zz)); b,e,_=ols(np.c_[np.ones(n),x],y)
        if nboot:
            bs=np.array([ols(np.c_[np.ones(n),x[j]],y[j])[0][1] for j in (rng.integers(0,n,n) for _ in range(nboot))]); o[name+"_boot"]=float(bs.std()); o[name+"_ge1"]=float(np.mean(bs>=1))
        o[name]=float(b[1]); o[name+"_se"]=float(e[1])
    if not quiet: print(f"  {label:58s} n={n:3d} median {o['median']:.2f} scatter {o['rms']:.3f} dex | (1+z)^{o['p']:+.2f}+-{o['p_se']:.2f} | k = {o['k']:+.2f}+-{o['k_se']:.2f}"+(f" (boot {o['k_boot']:.2f}; {o['k_ge1']*100:.2f}% reach 1)" if nboot else "")+f" | Lerner {o['klin']:+.2f}+-{o['klin_se']:.2f} | stretched {o['ktd']:+.2f}+-{o['ktd_se']:.2f}")
    return o
# ---- pressure-profile shapes (generalised NFW, Arnaud et al. 2010): universal, cool-core, disturbed
SHAPES={"universal":(1.177,0.3081,1.0510,5.4905),"cool-core":(1.128,0.7736,1.2223,5.4905),"disturbed":(1.083,0.3798,1.4063,5.4905)}
def make(shape):
    c5,ga,al,be=SHAPES[shape]; p=lambda x: 1/((c5*x)**ga*(1+(c5*x)**al)**((be-ga)/al))
    los=lambda b,xmax: 0.0 if b>=xmax else 2*integrate.quad(lambda l:p(np.hypot(b,l)),0,np.sqrt(xmax**2-b**2),limit=200)[0]
    Kcyl=lambda a,xmax=8.0: integrate.quad(lambda b:2*np.pi*b*los(b,xmax),0,a,limit=200)[0]
    Ksph=lambda x1: integrate.quad(lambda x:4*np.pi*x*x*p(x),0,x1,limit=200)[0]
    return Kcyl,Ksph
G={s:make(s) for s in SHAPES}
# ---- data
B=[l.split() for l in open(dd+"/bulbul_spt.txt") if not l.startswith("#")]
name=np.array([r[0] for r in B]); f=lambda i: np.array([float(r[i]) for r in B])
z=f(6); xi=f(2); thc=f(3); Ysz=f(4)*1e-6; eY=f(5)*1e-6; R500=f(7); Tin=f(8); Tin_lo=f(9); Tin_hi=f(10); Tex=f(11); Mg=f(14)/10; Mg_lo=f(15)/10; Mg_hi=f(16)/10; Yx_print=f(17); M500=f(19)
chk=Mg*10*Tin/10/Yx_print; dup=np.array([n=="0330-5228" for n in name]); ok=(np.abs(chk-1)<0.10)&~dup
print(f"{len(z)} clusters, z {z.min():.2f}-{z.max():.2f}; dropped: {', '.join(name[~ok])} (printed gas mass and Y_X disagree, or the row repeats another cluster's values)")
am=DA(z)*1e3*np.pi/180/60                                         # kpc per arcmin
a=0.75*am/R500                                                    # the SPT aperture in units of R500
arcmin2=(DA(z)*np.pi/180/60)**2                                   # Mpc^2 per arcmin^2
def predict(T,shape="universal",scale=1.0):
    """X-ray predicted Y inside the 0.75 arcmin column, arcmin^2. scale: R500(true)/R500(printed), for the aperture mismatch check."""
    Kcyl,Ksph=G[shape]; g=np.array([Kcyl(ai/scale)/Ksph(1/scale) for ai in a]); return 1.4262e-5*Mg*T*g/arcmin2, g
Yp,g=predict(Tin); lr=np.log10(Ysz/Yp)
print(f"0.75 arcmin is {np.median(a):.2f} R500 (range {a.min():.2f}-{a.max():.2f}); share of the R500 sphere's signal inside the column: median {np.median(g):.2f} (range {g.min():.2f}-{g.max():.2f})")
out={"n_all":int(len(z)),"dropped":name[~ok].tolist()}; res={}
res["main"]=kfit(z,lr,"all usable, temperature with core, universal shape",ok)
res["Tex"]=kfit(z,np.log10(Ysz/predict(Tex)[0]),"temperature with the core cut out",ok,nboot=0)
res["cc"]=kfit(z,np.log10(Ysz/predict(Tin,"cool-core")[0]),"cool-core shape for all",ok,nboot=0)
res["md"]=kfit(z,np.log10(Ysz/predict(Tin,"disturbed")[0]),"disturbed shape for all",ok,nboot=0)
# least favourable shape mix for expansion: nearby clusters cool-core, distant ones disturbed, and the reverse
lr_cc=np.log10(Ysz/predict(Tin,"cool-core")[0]); lr_md=np.log10(Ysz/predict(Tin,"disturbed")[0]); zmid=np.median(z[ok])
res["mixA"]=kfit(z,np.where(z<zmid,lr_cc,lr_md),"cool-core below the median redshift, disturbed above",ok,nboot=0)
res["mixB"]=kfit(z,np.where(z<zmid,lr_md,lr_cc),"disturbed below the median redshift, cool-core above",ok,nboot=0)
Rm=(3*M500*1e14/(4*np.pi*500*2.775e11*(H0/100)**2*Ez(z)**2))**(1/3)*1e3; sc=Rm/R500
print(f"   printed R500 is {np.median(1/sc):.3f} times the R500 that the printed mass implies (range {(1/sc).min():.3f}-{(1/sc).max():.3f})")
res["scale"]=kfit(z,np.log10(Ysz/np.array([1.4262e-5*Mg[i]*Tin[i]*G['universal'][0](a[i]/sc[i])/G['universal'][1](1/sc[i])/arcmin2[i] for i in range(len(z))])),"profile scaled to the R500 that the printed mass implies",ok,nboot=0)
res["xi6"]=kfit(z,lr,"SPT significance of 6 or more",ok&(xi>=6),nboot=0)
res["xi8"]=kfit(z,lr,"SPT significance of 8 or more",ok&(xi>=8),nboot=0)
res["xi10"]=kfit(z,lr,"SPT significance of 10 or more",ok&(xi>=10),nboot=0)
res["core05"]=kfit(z,lr,"SPT core radius 0.5 arcmin or less",ok&(thc<=0.5),nboot=0)
res["core75"]=kfit(z,lr,"SPT core radius 0.75 arcmin or more",ok&(thc>=0.75),nboot=0)
res["zgt03"]=kfit(z,lr,"z > 0.3",ok&(z>0.3),nboot=0)
res["zlt1"]=kfit(z,lr,"z < 1",ok&(z<1.0),nboot=0)
eT=(Tin_lo+Tin_hi)/2/Tin; eM=(Mg_lo+Mg_hi)/2/Mg; ex=np.hypot(eT,eM)
res["err30"]=kfit(z,lr,"X-ray error below 30%",ok&(ex<0.30),nboot=0)
m=ok; n=int(m.sum())
b,e,_=ols(np.c_[np.ones(n),np.log10(F_log(z[m])),np.log10(Tin[m]/6)],lr[m]); print(f"   with a temperature term: k = {b[1]:+.2f}+-{e[1]:.2f}; temperature power {b[2]:+.2f}+-{e[2]:.2f}"); res["withT"]=dict(k=float(b[1]),k_se=float(e[1]),a=float(b[2]),a_se=float(e[2]))
b,e,_=ols(np.c_[np.ones(n),np.log10(F_log(z[m])),np.log10(xi[m]/9)],lr[m]); print(f"   with a significance term: k = {b[1]:+.2f}+-{e[1]:.2f}; significance power {b[2]:+.2f}+-{e[2]:.2f}"); res["withXi"]=dict(k=float(b[1]),k_se=float(e[1]),a=float(b[2]),a_se=float(e[2]))
b,e,_=ols(np.c_[np.ones(n),np.log10(F_log(z[m])),np.log10(a[m]/0.25)],lr[m]); print(f"   with an aperture term (a shape error would show here): k = {b[1]:+.2f}+-{e[1]:.2f}; aperture power {b[2]:+.2f}+-{e[2]:.2f}"); res["withA"]=dict(k=float(b[1]),k_se=float(e[1]),a=float(b[2]),a_se=float(e[2]))
sig=np.sqrt((eY/Ysz/np.log(10))**2+(ex/np.log(10))**2+0.05**2); b,e,_=ols(np.c_[np.ones(n),np.log10(F_log(z[m]))],lr[m],1/sig[m]**2); print(f"   weighted by measurement error: k = {b[1]:+.2f}+-{e[1]:.2f}"); res["weighted"]=dict(k=float(b[1]),k_se=float(e[1]))
edges=[0.19,0.30,0.40,0.47,0.60,0.80,1.5]; bins=[]
for lo,hi in zip(edges,edges[1:]):
    q=m&(z>=lo)&(z<hi); bins.append([float(np.median(z[q])),int(q.sum()),float(lr[q].mean()),float(lr[q].std(ddof=1)/np.sqrt(q.sum())),float(np.log10(F_log(z[q])).mean()),float(np.log10(F_lin(z[q])).mean())])
print("   binned:",", ".join(f"z~{b[0]:.2f} (n={b[1]}) {b[2]:+.3f}+-{b[3]:.3f}" for b in bins),f"| static predicts a change of {bins[-1][4]-bins[0][4]:+.3f} from first to last; measured {bins[-1][2]-bins[0][2]:+.3f}+-{np.hypot(bins[-1][3],bins[0][3]):.3f}")
res["bins"]=bins
# the distant third on its own
hi=m&(z>=0.55); lo=m&(z<0.55); d=lr[hi].mean()-lr[lo].mean(); ed=np.hypot(lr[hi].std(ddof=1)/np.sqrt(hi.sum()),lr[lo].std(ddof=1)/np.sqrt(lo.sum())); pl=np.log10(F_log(z[hi])).mean()-np.log10(F_log(z[lo])).mean(); pn=np.log10(F_lin(z[hi])).mean()-np.log10(F_lin(z[lo])).mean(); pt=np.log10(F_td(z[hi])).mean()-np.log10(F_td(z[lo])).mean()
print(f"   {int(hi.sum())} clusters at z >= 0.55 (median {np.median(z[hi]):.2f}) against {int(lo.sum())} nearer (median {np.median(z[lo]):.2f}): measured change {d:+.3f}+-{ed:.3f} dex; static tired light predicts {pl:+.3f}, Lerner's law {pn:+.3f}, stretched time {pt:+.3f}")
res["split"]=dict(n_hi=int(hi.sum()),n_lo=int(lo.sum()),z_hi=float(np.median(z[hi])),z_lo=float(np.median(z[lo])),d=float(d),ed=float(ed),pred_log=float(pl),pred_lin=float(pn),pred_td=float(pt))
out["C"]=res
# ---- the 17 clusters also in the Chandra sample: do the two X-ray routes agree, and does their ratio drift with redshift?
S={l.split()[0]:l.split() for l in open(dd+"/sanders_spt.txt") if not l.startswith("#")}; kk=1.3018e-27*3.0857e21
cz=[];cr=[];cn=[]
cHs=299792.458/70.; Es=lambda t: np.sqrt(0.3*(1+t)**3+0.7)
for i,nm in enumerate(name):
    if nm in S and ok[i]:
        s=S[nm]; zs=float(s[1]); ams=float(s[7]); Ych=kk*float(s[9])/ams**2       # Chandra cylinder, arcmin^2, in the h = 0.7 cosmology of Entry 28
        if abs(zs-float(s[6]))>0.03*(1+zs): continue
        # put the Chandra value on this script's distance scale: predicted y scales as sqrt(D_A)
        das=cHs*integrate.quad(lambda x:1/Es(x),0,zs)[0]/(1+zs); Ych*=np.sqrt(DA(z[i])/das)
        cz.append(z[i]); cr.append(np.log10(Yp[i]/Ych)); cn.append(nm)
cz=np.array(cz); cr=np.array(cr); b,e,rms=ols(np.c_[np.ones(len(cz)),np.log10(1+cz)],cr)
print(f"\n{len(cz)} clusters with both XMM (gas mass x temperature, universal shape) and Chandra (measured profile) predictions for the same column:")
print(f"   XMM/Chandra median {10**np.median(cr):.2f}, scatter {cr.std(ddof=1):.3f} dex; trend with redshift (1+z)^{b[1]:+.2f}+-{e[1]:.2f}")
out["common"]=dict(n=int(len(cz)),median=float(10**np.median(cr)),scatter=float(cr.std(ddof=1)),p=float(b[1]),p_se=float(e[1]),pts=[[float(a_),float(b_),c_] for a_,b_,c_ in zip(cz,cr,cn)])
# ---- joint with Entry 28: SPT+Chandra (73), Planck+XMM (114), and the SPT+XMM clusters that are not in the Chandra set
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import io, contextlib, importlib.util
J=json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"entry28_points.json")))
zA=np.array([p[0] for p in J["pointsA"]]); yA=np.array([p[1] for p in J["pointsA"]]); zB=np.array([p[0] for p in J["pointsB"]]); yB=np.array([p[1] for p in J["pointsB"]])
new=ok&np.array([nm not in S for nm in name]); zC=z[new]; yC=lr[new]
print(f"\njoint: SPT+Chandra {len(zA)}, Planck+XMM {len(zB)}, SPT+XMM not in the Chandra set {len(zC)}; {len(zA)+len(zB)+len(zC)} clusters, each counted once, a separate level for each pair")
zj=np.r_[zA,zB,zC]; yj=np.r_[yA,yB,yC]; D=np.zeros((len(zj),3)); D[:len(zA),0]=1; D[len(zA):len(zA)+len(zB),1]=1; D[len(zA)+len(zB):,2]=1; nj=len(zj); resJ={}
cH7=299792.458/70.; E7=lambda t: np.sqrt(0.3*(1+t)**3+0.7); DA7=np.vectorize(lambda t: cH7*integrate.quad(lambda x:1/E7(x),0,t)[0]/(1+t))
Fj={"k":lambda t: np.sqrt(cH7*np.log(1+t)/DA7(t))*(1+t)**-1.5,"klin":lambda t: np.sqrt(cH7*t/DA7(t))*(1+t)**-1.5,"ktd":lambda t: np.sqrt(cH7*np.log(1+t)/DA7(t))*(1+t)**-1.0}
for nm,lab in (("k","static, tired light"),("klin","static, Lerner's distance law"),("ktd","static, photons also stretched in time")):
    x=np.log10(Fj[nm](zj)); b,e,_=ols(np.c_[D,x],yj); bs=np.array([ols(np.c_[D[j],x[j]],yj[j])[0][3] for j in (rng.integers(0,nj,nj) for _ in range(4000))])
    resJ[nm]=dict(k=float(b[3]),k_se=float(e[3]),boot=float(bs.std()),ge1=float(np.mean(bs>=1)),sigma=float((1-b[3])/max(e[3],bs.std())))
    print(f"  {lab:40s} k = {b[3]:+.2f}+-{e[3]:.2f} (boot {bs.std():.2f}); {100*np.mean(bs>=1):.3f}% of resamples reach 1; {(1-b[3])/max(e[3],bs.std()):.1f} sigma from static, {abs(b[3])/max(e[3],bs.std()):.1f} from expanding")
b,e,_=ols(np.c_[D,np.log10(1+zj)],yj); print(f"  as a power of (1+z): {b[3]:+.2f}+-{e[3]:.2f}"); resJ["p"]=[float(b[3]),float(e[3])]
# beyond z = 0.55 only, the two SPT pairs (Chandra, and XMM for clusters not in the Chandra set)
far=zj>=0.55; Df=D[far][:,[0,2]]; keep=Df.sum(0)>0
for nm in ("k","klin","ktd"):
    x=np.log10(Fj[nm](zj[far])); b,e,_=ols(np.c_[Df,x],yj[far]); resJ[nm+"_far"]=dict(k=float(b[2]),k_se=float(e[2]),n=int(far.sum()))
print(f"  beyond z = 0.55 only ({int(far.sum())} clusters, z up to {zj[far].max():.2f}): k = {resJ['k_far']['k']:+.2f}+-{resJ['k_far']['k_se']:.2f}; Lerner {resJ['klin_far']['k']:+.2f}+-{resJ['klin_far']['k_se']:.2f}; stretched {resJ['ktd_far']['k']:+.2f}+-{resJ['ktd_far']['k_se']:.2f}")
# checks on the joint fit: only the stronger SPT+XMM detections; and the XMM clusters beyond z = 0.55 on their own
for lab,cut in (("xi8",8.0),("xi6",6.0)):
    nw=new&(xi>=cut); zq=np.r_[zA,zB,z[nw]]; yq=np.r_[yA,yB,lr[nw]]; Dq=np.zeros((len(zq),3)); Dq[:len(zA),0]=1; Dq[len(zA):len(zA)+len(zB),1]=1; Dq[len(zA)+len(zB):,2]=1
    b,e,_=ols(np.c_[Dq,np.log10(Fj["k"](zq))],yq); bl,el,_=ols(np.c_[Dq,np.log10(Fj["klin"](zq))],yq); bt,et,_=ols(np.c_[Dq,np.log10(Fj["ktd"](zq))],yq)
    resJ[lab]=dict(n=int(len(zq)),k=float(b[3]),k_se=float(e[3]),klin=float(bl[3]),klin_se=float(el[3]),ktd=float(bt[3]),ktd_se=float(et[3]))
    print(f"  joint with only SPT+XMM clusters of significance {cut:.0f} or more ({int(nw.sum())} of them; {len(zq)} in all): k = {b[3]:+.2f}+-{e[3]:.2f}; Lerner {bl[3]:+.2f}+-{el[3]:.2f}; stretched {bt[3]:+.2f}+-{et[3]:.2f}")
farC=ok&(z>=0.55); res["far"]=kfit(z,lr,"SPT+XMM alone, z >= 0.55",farC,nboot=2000)
out["joint"]=resJ
out["points"]=[[float(a_),float(b_),n_,float(c_),bool(n_ in S)] for a_,b_,n_,c_ in zip(z[m],lr[m],name[m],xi[m])]
zz=np.linspace(0.15,1.5,55); out["curve"]={"z":zz.tolist(),"static_log":np.log10(F_log(zz)).tolist(),"static_lin":np.log10(F_lin(zz)).tolist(),"static_td":np.log10(F_td(zz)).tolist()}
json.dump(out,open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"spt_xmm_distance.json"),"w"),indent=1)
