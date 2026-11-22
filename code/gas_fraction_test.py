"""Workbench Entry 34: a second test of the one static model left by Entry 33, with cluster gas fractions.
Entry 33 left a static universe one model: the distance law measured by the clustering ruler (equal to the standard comoving distance
D_M = (1+z) D_A) and photons that lose energy and arrive stretched in time. Its test there used cluster shadows against X-ray glow.
This uses X-rays alone. The hot gas of a relaxed cluster gives two masses: the gas mass, from its brightness, and the total mass, from the
pressure needed to hold the gas up. Their ratio, the gas fraction, depends on the assumed distance as d_L d_A^(1/2).
Mantz et al. (2014, 2022) measured it in a shell at 0.8-1.2 r2500 for 40 and 44 of the hottest, most relaxed clusters, using a reference
standard cosmology (h = 0.7, matter 0.3, flat). If clusters are the same at all redshifts (in static space they must be):
   expanding, reference cosmology right:   f_ref(z) = constant
   static (measured law, stretched time):  f_ref(z) = constant x (1+z)^(-1/2) x [(1+z)/E(z)]^eta
The first factor is the distance: static brightness distance equals the reference one, static size distance is (1+z) larger.
The second is the aperture: the reference shell sits at a mean density of 2500 E^2 rho_c0 in its own distances, which in static distances
is 2500 E^2/(1+z)^2, so it lands at a radius (1+z)/E times farther out in an unchanging cluster, where the shell gas fraction is higher as r^eta.
k = fitted power of the static prediction: 0 expanding, 1 static."""
import numpy as np, json, sys, os
HERE=os.path.dirname(os.path.abspath(__file__)); rng=np.random.default_rng(34)
E=lambda z: np.sqrt(0.3*(1+z)**3+0.7)
def load(fn,cols):
    R=[l.split() for l in open(os.path.join(HERE,fn)) if not l.startswith("#")]
    return np.array([r[0] for r in R]),*(np.array([float(r[c]) for r in R]) for c in cols)
def tmpl(z,eta): return -0.5*np.log(1+z)+eta*np.log((1+z)/E(z))
def fit(x,y,s,sint):
    w=1/(s**2+sint**2); X=np.c_[np.ones(len(x)),x]; cov=np.linalg.inv((X*w[:,None]).T@X); b=cov@(X*w[:,None]).T@y; r=y-X@b; return b,np.sqrt(np.diag(cov)),float(np.sum(w*r*r))
def scatter(x,y,s):
    """intrinsic scatter (in ln f) that makes chi2 equal to the degrees of freedom; 0 if the measurement errors already suffice"""
    f=lambda si: fit(x,y,s,si)[2]-(len(x)-2)
    return 0.0 if f(0.0)<=0 else float(__import__("scipy.optimize").optimize.brentq(f,0,1))
def run(name,z,f,e,eta,label,sel=None,quiet=False,nboot=4000,sint_fix=None):
    m=np.ones(len(z),bool) if sel is None else sel; zz=z[m]; y=np.log(f[m]); s=e[m]/f[m]; x=tmpl(zz,eta)
    si=scatter(x,y,s) if sint_fix is None else sint_fix; b,se,chi=fit(x,y,s,si)
    bs=[]
    for _ in range(nboot):
        j=rng.integers(0,len(zz),len(zz)); bs.append(fit(x[j],y[j],s[j],si)[0][1])
    bs=np.array(bs) if nboot else np.array([0.0]); err=max(se[1],bs.std())
    bp,sp,_=fit(np.log(1+zz),y,s,si)
    o=dict(n=int(m.sum()),k=float(b[1]),k_se=float(se[1]),k_boot=float(bs.std()),err=float(err),sint=float(si),chi=chi,p=float(bp[1]),p_se=float(sp[1]),ge1=float(np.mean(bs>=1)) if nboot else None,
           chi_exp=float(fit(np.zeros(len(zz)),y,s,si)[2]) if False else None)
    # chi2 of each fixed model (free level only)
    w=1/(s**2+si**2)
    for key,pred in (("expanding",np.zeros(len(zz))),("static",x)):
        d=y-pred; a=np.sum(w*d)/w.sum(); o["chi_"+key]=float(np.sum(w*(d-a)**2))
    if not quiet: print(f"  {label:58s} n={o['n']:2d} | k = {o['k']:+.2f} +- {err:.2f} (lsq {se[1]:.2f}, boot {bs.std():.2f}) | scatter {100*si:.1f}% | (1+z)^({bp[1]:+.2f}+-{sp[1]:.2f}) | chi2: expanding {o['chi_expanding']:.1f}, static {o['chi_static']:.1f} | static {(1-o['k'])/err:.1f} sigma away, expanding {abs(o['k'])/err:.1f}")
    return o
out={}
for tag,fn,cols,eta,title in (("2022","mantz2022_fgas.txt",(1,6,7,8,2,4),0.390,"Mantz et al. 2022: 44 clusters"),("2014","mantz2014_fgas.txt",(1,5,6,3),0.442,"Mantz et al. 2014: 40 clusters")):
    d=load(fn,cols); name,z,f,e,M=d[0],d[1],d[2],d[3],d[4]
    print(f"\n{title}, z = {z.min():.3f} to {z.max():.3f}; median error {100*np.median(e/f):.0f}%; {int((z>0.55).sum())} beyond z = 0.55, {int((z>0.9).sum())} beyond 0.9")
    if tag=="2022":
        chk=d[5]/(10*d[6])/f; print(f"   check: gas mass / total mass against the printed fraction: median {np.median(chk):.3f}, range {chk.min():.3f}-{chk.max():.3f}")
    print(f"   static prediction relative to z = 0: "+", ".join(f"z={q}: {np.exp(tmpl(np.array([q]),eta))[0]:.3f}" for q in (0.1,0.3,0.5,0.7,1.0,1.16)))
    r={}
    r["main"]=run(name,z,f,e,eta,"all clusters")
    r["noscatter"]=run(name,z,f,e,eta,"measurement errors only (no intrinsic scatter)",nboot=0,sint_fix=0.0)
    r["s074"]=run(name,z,f,e,eta,"intrinsic scatter fixed at 7.4%",nboot=0,sint_fix=0.074)
    r["eta0"]=run(name,z,f,e,0.0,"no aperture term (distance factor only)",nboot=0)
    r["eta_hi"]=run(name,z,f,e,eta+0.1,"aperture slope raised by 0.1",nboot=0)
    r["z>0.07"]=run(name,z,f,e,eta,"without the nearest cluster",z>0.05,nboot=0) if tag=="2022" else None
    r["z>0.2"]=run(name,z,f,e,eta,"z > 0.2 only",z>0.2,nboot=0)
    r["z<1.1"]=run(name,z,f,e,eta,"without the most distant cluster",z<z.max(),nboot=0)
    r["hiM"]=run(name,z,f,e,eta,"M2500 above the median",M>=np.median(M),nboot=0)
    r["loM"]=run(name,z,f,e,eta,"M2500 below the median",M<np.median(M),nboot=0)
    r["unw"]=run(name,z,f,np.full(len(z),0.01)*f,eta,"every cluster given the same weight",nboot=2000)
    # with a mass term
    y=np.log(f); s=e/f; x=tmpl(z,eta); si=r["main"]["sint"]; w=1/(s**2+si**2); X=np.c_[np.ones(len(z)),x,np.log(M/3)]; cov=np.linalg.inv((X*w[:,None]).T@X); b=cov@(X*w[:,None]).T@y
    print(f"   with a mass term: k = {b[1]:+.2f} +- {np.sqrt(cov[1,1]):.2f}; f goes as M^({b[2]:+.3f} +- {np.sqrt(cov[2,2]):.3f})"); r["withM"]=dict(k=float(b[1]),k_se=float(np.sqrt(cov[1,1])),a=float(b[2]),a_se=float(np.sqrt(cov[2,2])))
    # binned
    edges=[0,0.16,0.26,0.37,0.47,0.65,1.2]; bins=[]
    for lo,hi in zip(edges,edges[1:]):
        q=(z>=lo)&(z<hi); ww=w[q]; mu=np.sum(ww*y[q])/ww.sum(); bins.append([float(np.median(z[q])),int(q.sum()),float(np.exp(mu)),float(np.exp(mu)/np.sqrt(ww.sum())),float(np.exp(np.mean(x[q])))])
    print("   binned f_gas:",", ".join(f"z~{b_[0]:.2f} (n={b_[1]}) {b_[2]:.3f}+-{b_[3]:.3f}" for b_ in bins)); r["bins"]=bins
    hi=z>=0.55; lo=~hi; mh=np.sum(w[hi]*y[hi])/w[hi].sum(); ml=np.sum(w[lo]*y[lo])/w[lo].sum(); ed=np.hypot(1/np.sqrt(w[hi].sum()),1/np.sqrt(w[lo].sum())); pr=np.mean(x[hi])-np.mean(x[lo])
    print(f"   {int(hi.sum())} clusters at z >= 0.55 (median {np.median(z[hi]):.2f}) against {int(lo.sum())} nearer (median {np.median(z[lo]):.2f}): ratio {np.exp(mh-ml):.3f} +- {np.exp(mh-ml)*ed:.3f}; static predicts {np.exp(pr):.3f}; expanding 1")
    r["split"]=dict(n_hi=int(hi.sum()),n_lo=int(lo.sum()),z_hi=float(np.median(z[hi])),z_lo=float(np.median(z[lo])),ratio=float(np.exp(mh-ml)),err=float(np.exp(mh-ml)*ed),static=float(np.exp(pr)))
    r["points"]=[[float(a),float(b_),float(c),n_] for a,b_,c,n_ in zip(z,f,e,name)]
    out[tag]=r
# ---- gas fractions and shadows together. Both use the gas density read from the X-ray glow.
# Static space (Entry 33's last model) needs both ratios to fall: k = 1 for each. Expanding space with unchanging clusters needs k = 0 for each.
# Shadow results from Entry 33 (static_distance_law.json): power p of (1+z) in measured/predicted shadow; static predicts -1/2, so k = -2p.
from scipy.stats import chi2 as _c2, norm as _nm
sig=lambda c,dof: float(_nm.isf(_c2.sf(c,dof)/2))
S=json.load(open(os.path.join(HERE,"static_distance_law.json")))["clusters"]; kf=out["2022"]["main"]["k"]; ef=out["2022"]["main"]["err"]; pf=out["2022"]["main"]["p"]; epf=out["2022"]["main"]["p_se"]
print(f"\nGas fractions (2022 data) with the shadows of Entries 28 and 32:  k_gas = {kf:+.2f} +- {ef:.2f}")
out["joint"]={}
for key,lab in (("AB","pairs A and B, 187"),("ABCs","A, B and strong C, 203"),("ABC","A, B and C, 227")):
    p=S[key]["p"]; e=S[key]["e"]; kr=-2*p; er=2*e
    cs=((1-kf)/ef)**2+((1-kr)/er)**2; ce=(kf/ef)**2+(kr/er)**2
    g1,e1=-pf,epf; g2,e2=p,e; w1,w2=1/e1**2,1/e2**2; g=(g1*w1+g2*w2)/(w1+w2); eg=1/np.sqrt(w1+w2); cg=(g1-g)**2*w1+(g2-g)**2*w2
    print(f"   shadows ({lab}): k_shadow = {kr:+.2f} +- {er:.2f} | static (1,1): chi2 {cs:.1f}, {sig(cs,2):.1f} sigma | expanding, no change (0,0): chi2 {ce:.1f}, {sig(ce,2):.1f} sigma | the two k differ by {(kf-kr)/np.hypot(ef,er):.1f} sigma")
    print(f"        expanding with X-ray gas density changing as (1+z)^-g: g = {g:.2f} +- {eg:.2f} (gas fractions {g1:.2f}+-{e1:.2f}, shadows {g2:.2f}+-{e2:.2f}; chi2 {cg:.1f} for 1)")
    out["joint"][key]=dict(k_shadow=float(kr),e_shadow=float(er),chi_static=float(cs),sig_static=sig(cs,2),chi_exp=float(ce),sig_exp=sig(ce,2),g=float(g),g_err=float(eg),chi_g=float(cg),apart=float((kf-kr)/np.hypot(ef,er)))
zz=np.linspace(0.0,1.2,49); out["curve"]={"z":zz.tolist(),"static_2022":np.exp(tmpl(zz,0.390)).tolist(),"static_2014":np.exp(tmpl(zz,0.442)).tolist(),"distance_only":np.exp(tmpl(zz,0.0)).tolist()}
json.dump(out,open(os.path.join(HERE,"gas_fraction_test.json"),"w"),indent=1)
