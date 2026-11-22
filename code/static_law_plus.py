"""Entry 39, part 3: Entry 31's test of the static Hubble laws, rerun on Pantheon+ (diagonal errors, repeats averaged).
Every law gets a free additive constant. Same laws and cuts as sn/pantheon_check.py."""
import sys, os, json
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
import numpy as np
from scipy import integrate, optimize
d=np.loadtxt(os.path.join(HERE,"pantheon_plus_rows.txt")); zr,mr,er=d.T
zu=np.unique(zr); z=zu; m=np.array([np.sum(mr[zr==a]/er[zr==a]**2)/np.sum(1/er[zr==a]**2) for a in zu]); e=np.array([er[zr==a].min() for a in zu]); w=1/e**2; N=len(z)
def chi2(mu_model,sel=None,step=None):
    s=np.ones(N,bool) if sel is None else sel; dd=(m-mu_model)[s]; ww=w[s]
    if step is None:
        a=np.sum(ww*dd)/np.sum(ww); return float(np.sum(ww*(dd-a)**2)),float(a),0.0
    lo=z[s]<step; a1=np.sum(ww[lo]*dd[lo])/np.sum(ww[lo]); a2=np.sum(ww[~lo]*dd[~lo])/np.sum(ww[~lo]); return float(np.sum(ww[lo]*(dd[lo]-a1)**2)+np.sum(ww[~lo]*(dd[~lo]-a2)**2)),float(a2),float(a1-a2)
mu=lambda dl: 5*np.log10(dl); E=lambda zz,om: np.sqrt(om*(1+zz)**3+1-om); zg=np.linspace(0,2.4,2401)
def dl_lcdm(om,zz=None):
    zz=z if zz is None else zz; dc=integrate.cumulative_trapezoid(1/E(zg,om),zg,initial=0); return (1+zz)*np.interp(zz,zg,dc)
LAWS={"lcdm":None,"lerner":lambda t: t*np.sqrt(1+t),"cz_stretched":lambda t: t*(1+t),"tired":lambda t: np.log(1+t)*np.sqrt(1+t),"log_stretched":lambda t: np.log(1+t)*(1+t),"coasting":lambda t: t*(1+t/2),"matter_only":lambda t: 2*(1+t)*(1-1/np.sqrt(1+t))}
def fit_om(sel=None,step=None):
    f=lambda o: chi2(mu(dl_lcdm(o)),sel,step)[0]; g=optimize.minimize_scalar(f,bounds=(0.05,0.9),method="bounded"); lo=optimize.brentq(lambda o:f(o)-g.fun-1,0.05,g.x); hi=optimize.brentq(lambda o:f(o)-g.fun-1,g.x,0.9); return g.x,lo,hi,g.fun
def fit_q(sel=None,step=None):
    f=lambda q: chi2(mu(z*(1+z)**q),sel,step)[0]; g=optimize.minimize_scalar(f,bounds=(-1,3),method="bounded"); lo=optimize.brentq(lambda q:f(q)-g.fun-1,-1,g.x); hi=optimize.brentq(lambda q:f(q)-g.fun-1,g.x,3); return float(g.x),float((hi-lo)/2),float(g.fun)
out={"N":N}
om,lo,hi,c0=fit_om(); print(f"{N} supernovae. Flat LCDM: matter share {om:.3f} +- {(hi-lo)/2:.3f}, chi2 {c0:.0f}"); out["om"]=[om,(hi-lo)/2,c0]
def table(sel=None,step=None,label="all"):
    o,_,_,c=fit_om(sel,step); res={}; print(f"--- {label}: {N if sel is None else int(sel.sum())} supernovae; LCDM matter {o:.3f}, chi2 {c:.0f}")
    for k,fn in LAWS.items():
        cc,a,st=chi2(mu(dl_lcdm(o)) if fn is None else mu(fn(z)),sel,step); res[k]=[cc,cc-c,st]; print(f"   {k:14s} chi2 {cc:8.0f}  delta {cc-c:+8.0f}"+(f"  step {st:+.3f}" if step else ""))
    return res
out["all"]=table(); out["step"]=table(step=0.1,label="free step at z = 0.1"); out["far"]=table(sel=z>0.1,label="beyond z = 0.1")
out["q"]={}
for lab,sel in (("all",None),("z<0.3",z<0.3),("z>0.1",z>0.1),("z>0.3",z>0.3),("z<0.7",z<0.7)):
    q=fit_q(sel); out["q"][lab]=q; print(f"   q, {lab:6s}: {q[0]:.3f} +- {q[1]:.3f}")
a,b=out["q"]["z<0.3"],out["q"]["z>0.3"]; print(f"   halves disagree by {a[0]-b[0]:.3f} +- {np.hypot(a[1],b[1]):.3f} ({(a[0]-b[0])/np.hypot(a[1],b[1]):.1f} sigma)")
o1=fit_om(z<0.3); o2=fit_om(z>0.3); print(f"   LCDM matter share: z<0.3 {o1[0]:.3f} +- {(o1[2]-o1[1])/2:.3f}; z>0.3 {o2[0]:.3f} +- {(o2[2]-o2[1])/2:.3f}"); out["om_halves"]=[o1[0],(o1[2]-o1[1])/2,o2[0],(o2[2]-o2[1])/2]
sel=z>0.3; near=z<0.1; mu_s=mu(LAWS["lerner"](z)); a_far=chi2(mu_s,sel)[1]; dn=np.sum(w[near]*(m-mu_s-a_far)[near])/w[near].sum(); en=np.hypot(1/np.sqrt(w[near].sum()),1/np.sqrt(w[sel].sum()))
print(f"   Lerner's law levelled on z > 0.3: the nearby supernovae lie {dn:+.3f} +- {en:.3f} mag from it"); out["anchor_lerner"]=[float(dn),float(en)]
mu_t=mu(LAWS["log_stretched"](z)); a_far=chi2(mu_t,sel)[1]; dn=np.sum(w[near]*(m-mu_t-a_far)[near])/w[near].sum(); print(f"   ln(1+z) with stretched time, levelled on z > 0.3: nearby lie {dn:+.3f} +- {en:.3f}"); out["anchor_log"]=[float(dn),float(en)]
json.dump(out,open(os.path.join(HERE,"static_law_plus.json"),"w"),indent=1)
