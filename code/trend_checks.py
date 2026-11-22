"""Entry 45, checks on the trend inside the distant sample: what drives it, and does galaxy size account for it?
Run dynamical_mass_test.py first; it writes the json file this script reads. Both read vds2013_table4.txt from their own folder."""
import os, json, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); D=json.load(open(os.path.join(HERE,"dynamical_mass_test.json")))
R=[l.split() for l in open(os.path.join(HERE,"vds2013_table4.txt")) if not l.startswith("#")]
R=[r for r in R if not (r[1]=="GN5" or (r[0]=="7" and r[1]=="19627"))]
ref=np.array([int(r[0]) for r in R]); z,re,n,sig,kn,lMd,lMs=np.array([[float(x) for x in r[2:]] for r in R]).T
y0=lMd-lMs; lx=np.log10(1+z); rng=np.random.default_rng(451)
Fi={k:np.interp(z,[a for a,b in v],[b for a,b in v]) for k,v in D["F"].items()}; Fi["expanding"]=0*z
def slope(q,y,extra=None):
    X=np.c_[np.ones(q.sum()),lx[q]] if extra is None else np.c_[np.ones(q.sum()),lx[q],extra[q]]
    b=np.linalg.lstsq(X,y[q],rcond=None)[0]; bs=[]
    idx=np.where(q)[0]
    for _ in range(3000):
        j=rng.choice(idx,len(idx)); Xj=np.c_[np.ones(len(j)),lx[j]] if extra is None else np.c_[np.ones(len(j)),lx[j],extra[j]]
        bs.append(np.linalg.lstsq(Xj,y[j],rcond=None)[0])
    return b,np.std(bs,axis=0)
sel={"all 71":z>0,"z >= 0.9":z>=0.9,"z >= 1.0":z>=1.0,"without van der Wel 2008":ref!=6,"van der Wel 2008 only":ref==6,"z < 0.9 left out of van der Wel":(ref!=6)|(z>=0.9)}
for key in ("expanding","lerner","tired"):
    y=y0+Fi[key]; print(key)
    for nm,q in sel.items():
        b,e=slope(q,y); print(f"   {nm:34s} n={q.sum():2d}  power of (1+z) inside the sample: {b[1]:+.2f} +- {e[1]:.2f}")
q9=z<0.9
print(f"\nthe group below z = 0.9 (n={q9.sum()}, all from one source) against the rest (n={(~q9).sum()}): {y0[q9].mean():+.3f} and {y0[~q9].mean():+.3f} dex; difference {y0[q9].mean()-y0[~q9].mean():+.2f} +- {np.hypot(y0[q9].std(ddof=1)/q9.sum()**.5,y0[~q9].std(ddof=1)/(~q9).sum()**.5):.2f}")
print(f"median radius {np.median(re[q9]):.1f} and {np.median(re[~q9]):.1f} kpc; median dispersion {np.median(sig[q9]):.0f} and {np.median(sig[~q9]):.0f} km/s; median log stellar mass {np.median(lMs[q9]):.2f} and {np.median(lMs[~q9]):.2f}")
# does size account for the trend?  (errors in radius enter the dynamical mass, so some correlation is built in)
for key in ("expanding","lerner","tired"):
    y=y0+Fi[key]; s={"expanding":0*z,"lerner":np.log10(z*(1+z)),"tired":np.log10(np.log(1+z)*(1+z))}[key]
    # size under each law: published radius times d(1+z)/D_M; D_M from F: F = size/lum; use direct formula instead
    from scipy import integrate
    DM=np.array([integrate.quad(lambda t:1/np.sqrt(0.3*(1+t)**3+0.7),0,a)[0] for a in z])
    lr=np.log10(re)+(0 if key=="expanding" else s-np.log10(DM))
    b,e=slope(z>0,y,lr); print(f"{key:10s} with radius as a second term: power of (1+z) {b[1]:+.2f} +- {e[1]:.2f}; power of radius {b[2]:+.2f} +- {e[2]:.2f}; median radius {np.median(10**lr):.1f} kpc")
# sigma at fixed stellar mass, no sizes used
print(f"\nmedian dispersion {np.median(sig):.0f} km/s at median log stellar mass {np.median(lMs):.2f} (published) ")
for a,b_ in ((0.8,1.0),(1.0,1.5),(1.5,2.3)):
    q=(z>=a)&(z<b_); print(f"  z {a}-{b_}: sigma {np.median(sig[q]):.0f}, radius {np.median(re[q]):.1f} kpc, log M* {np.median(lMs[q]):.2f}, log Mdyn {np.median(lMd[q]):.2f}")

# ---- tied to the nearby value, with and without the group at z = 0.83; and the published BOSS trend at 0.1 < z < 0.6 ----
LOCAL=-np.log10(0.59); zl=0.06; lxa=lx-np.log10(1+zl); out={}
def F_at(key,x): v=D["F"][key]; return float(np.interp(x,[a for a,b in v],[b for a,b in v]))
print()
for nm,q in (("all 71",z>0),("without the 27 below z = 0.9",z>=0.9),("only the 27 below z = 0.9",z<0.9)):
    ex=(y0-LOCAL)[q]; s=np.sum(lxa[q]*ex)/np.sum(lxa[q]**2); bs=[]; idx=np.arange(q.sum())
    for _ in range(4000):
        j=rng.choice(idx,len(idx)); bs.append(np.sum(lxa[q][j]*ex[j])/np.sum(lxa[q][j]**2))
    print(f"published frame, tied to nearby, {nm:30s}: mean {10**ex.mean():.2f} x nearby; power of (1+z) {s:+.2f} +- {np.std(bs):.2f}"); out[nm]=[float(s),float(np.std(bs)),float(ex.mean())]
print("what each static law predicts for that power if nothing evolves (from nearby to the sample):")
for key in ("lerner","tired","lerner_t","tired_t"):
    pred=-(np.array([F_at(key,a) for a in z])-F_at(key,zl)); s=np.sum(lxa*pred)/np.sum(lxa**2); print(f"   {key:9s} {s:+.2f}"); out["pred_"+key]=float(s)
print("\nBOSS (Beifiori et al. 2014): dynamical/stellar changes as (1+z)^(-0.30 +- 0.12) over 0.1 < z < 0.6 with standard distances")
for key in ("lerner","tired","tired_t"):
    f1,f6=F_at(key,0.1),F_at(key,0.6); p=(f6-f1)/np.log10(1.6/1.1)
    print(f"   {key:8s}: the law's factor grows by {10**(f6-f1):.2f} from z = 0.1 to 0.6, a power {p:+.2f}; read with this law the BOSS ratio changes as (1+z)^({-0.30+p:+.2f} +- 0.12), and at z = 0.6 stands at {10**(-0.30*np.log10(1.6/1.06)+f6-F_at(key,zl)):.2f} x nearby")
print("\nsize factor the galaxies ask for if nothing evolves and luminosities stay as published, against what flat static space gives:")
from scipy import integrate
for a,b_ in ((0.8,1.0),(1.0,1.5),(1.5,2.3)):
    q=(z>=a)&(z<b_); m=(y0-LOCAL)[q].mean(); e=(y0-LOCAL)[q].std(ddof=1)/q.sum()**.5; zm=z[q].mean(); dm=integrate.quad(lambda t:1/np.sqrt(0.3*(1+t)**3+0.7),0,zm)[0]
    print(f"   z {zm:.2f}: asked {10**(-m):.2f} (+- {100*(10**e-1):.0f}%); Lerner {zm*(1+zm)/dm:.1f}; tired light {np.log(1+zm)*(1+zm)/dm:.1f}")
json.dump(out,open(os.path.join(HERE,"trend_checks.json"),"w"),indent=1)
