"""Workbench Entry 42: the shadow-against-glow distance test with a fifth pair of instruments, and clusters not chosen by their shadows.
Data: Czakon et al. (2015), table 3. 45 clusters at z = 0.15 to 0.89 picked from X-ray catalogues, observed with Bolocam at 140 GHz
(the Sunyaev-Zel'dovich shadow, Y inside r2500, in steradians) and Chandra (gas mass inside r2500, temperature).
Measured shadow: Y2500 as an angle on the sky. Predicted from X-rays: gas mass x temperature, turned into a shadow with the standard distance.
As in Entries 27, 28 and 32, a static universe needs measured/predicted to fall as F(z); an expanding one needs it level. k = power of F.
Both sides use the same aperture r2500, so no pressure-profile shape is assumed."""
import numpy as np, json, sys, os
from scipy import integrate
HERE=os.path.dirname(os.path.abspath(__file__))
H0,Om=70.0,0.3                                                # the cosmology of the published table
cH=299792.458/H0; Ez=lambda z: np.sqrt(Om*(1+z)**3+1-Om)
DA=np.vectorize(lambda z: cH*integrate.quad(lambda x:1/Ez(x),0,z)[0]/(1+z))
F_log=lambda z: np.sqrt(cH*np.log(1+z)/DA(z))*(1+z)**-1.5      # static, tired light, d = (c/H0) ln(1+z)
F_lin=lambda z: np.sqrt(cH*z/DA(z))*(1+z)**-1.5                # static, Lerner's d = cz/H0
F_td =lambda z: np.sqrt(cH*np.log(1+z)/DA(z))*(1+z)**-1.0      # static, photons also stretched in time
FS=(("k",F_log),("klin",F_lin),("ktd",F_td))
def wls(X,y,w):
    W=np.sqrt(w)[:,None]; b=np.linalg.lstsq(X*W,y*W[:,0],rcond=None)[0]; r=y-X@b; n,p=X.shape; s2=(w*r*r).sum()/(n-p); cov=s2*np.linalg.inv((X*W).T@(X*W)); return b,np.sqrt(np.diag(cov)),r
rng=np.random.default_rng(42)
def kfit(z,y,label,w=None,groups=None,nboot=4000,quiet=False):
    """y = log10(measured/predicted). groups: integer label per cluster, each with its own level."""
    n=len(z); w=np.ones(n) if w is None else w; g=np.zeros(n,int) if groups is None else groups; D=np.stack([(g==k).astype(float) for k in np.unique(g)],1); o=dict(n=n)
    b,e,r=wls(np.c_[D,np.log10(1+z)],y,w); o.update(p=float(b[-1]),p_se=float(e[-1]),rms=float(np.sqrt(np.mean(r*r))))
    for name,F in FS:
        x=np.log10(F(z)); b,e,_=wls(np.c_[D,x],y,w); o[name]=float(b[-1]); o[name+"_se"]=float(e[-1])
        if nboot:
            bs=[]
            for _ in range(nboot):
                j=np.concatenate([rng.choice(np.where(g==k)[0],int((g==k).sum())) for k in np.unique(g)]); bs.append(wls(np.c_[D[j],x[j]],y[j],w[j])[0][-1])
            bs=np.array(bs); o[name+"_boot"]=float(bs.std()); o[name+"_ge1"]=float(np.mean(bs>=1))
    if not quiet:
        s=lambda nm: f"{o[nm]:+.2f} +- {o.get(nm+'_boot',o[nm+'_se']):.2f} ({(1-o[nm])/o.get(nm+'_boot',o[nm+'_se']):.1f} sigma from 1)"
        print(f"  {label:52s} n={n:3d} scatter {o['rms']:.3f} dex | (1+z)^{o['p']:+.2f}+-{o['p_se']:.2f} | tired light {s('k')} | Lerner {s('klin')} | stretched {s('ktd')}")
    return o
R=[l.split() for l in open(os.path.join(HERE,"boxsz_table3.txt")) if not l.startswith("#")]; name=np.array([r[0] for r in R]); d=np.array([[float(x) for x in r[1:]] for r in R])
z,r25,Mg,Mgp,Mgm,Mt,T,eT,Y,Yp,Ym,dB,dR=d.T
Ymeas=Y*1e-10*DA(z)**2                                          # Mpc^2 with the standard distance
Ypred=1.4262e-5*Mg*T                                            # Mpc^2: (sigma_T / m_e c^2) x k T x number of electrons in the gas mass
lr=np.log10(Ymeas/Ypred); el=np.sqrt(((Yp+Ym)/2/Y)**2+((Mgp+Mgm)/2/Mg)**2+(eT/T)**2)/np.log(10)
print(f"{len(z)} clusters, z {z.min():.2f}-{z.max():.2f}; measured/predicted: median {10**np.median(lr):.2f}, scatter {lr.std(ddof=1):.3f} dex; median measurement error {np.median(el):.3f} dex")
out={"n":int(len(z)),"median":float(10**np.median(lr)),"points":[[float(a),float(b),float(c)] for a,b,c in zip(z,lr,el)]}; res={}
print("1. Bolocam with Chandra, on its own")
res["all"]=kfit(z,lr,"all 45, equal weights")
sint=np.sqrt(max(lr.std(ddof=1)**2-np.median(el)**2,0.0)); w=1/(el**2+sint**2); print(f"   intrinsic scatter beyond the errors: {sint:.3f} dex")
res["weighted"]=kfit(z,lr,"all 45, weighted by errors and scatter",w=w)
for lab,m in (("beam correction below 20% (41)",dB<0.2),("beyond z = 0.3 (33)",z>0.3),("hotter than 6.5 keV",T>6.5),("without the 3 lowest and 3 highest ratios",(lr>np.sort(lr)[2])&(lr<np.sort(lr)[-3])),("gas mass above 0.2e14",Mg>0.2)):
    res[lab]=kfit(z[m],lr[m],lab,nboot=1000)
# split by redshift
lo=z<np.median(z); print(f"   mean log ratio below and above z = {np.median(z):.2f}: {lr[lo].mean():+.3f} +- {lr[lo].std(ddof=1)/lo.sum()**0.5:.3f}, {lr[~lo].mean():+.3f} +- {lr[~lo].std(ddof=1)/(~lo).sum()**0.5:.3f}; static tired light predicts a fall of {np.mean(np.log10(F_log(z[~lo])))-np.mean(np.log10(F_log(z[lo]))):+.3f}, Lerner {np.mean(np.log10(F_lin(z[~lo])))-np.mean(np.log10(F_lin(z[lo]))):+.3f}, stretched {np.mean(np.log10(F_td(z[~lo])))-np.mean(np.log10(F_td(z[lo]))):+.3f}")
out["halves"]=[float(lr[lo].mean()),float(lr[lo].std(ddof=1)/lo.sum()**0.5),float(lr[~lo].mean()),float(lr[~lo].std(ddof=1)/(~lo).sum()**0.5),float(np.median(z))]
# 2. with the earlier pairs
P28=json.load(open(os.path.join(HERE,"entry28_points.json"))); P32=json.load(open(os.path.join(HERE,"entry32_points.json")))
print("2. With the earlier pairs, each with its own level and weighted by its own scatter")
A=[[float(q[0]),float(q[1])] for q in P28["pointsA"]]; Bp=[[float(q[0]),float(q[1])] for q in P28["pointsB"]]; namesA={q[2] for q in P28["pointsA"]}
C=[[float(q[0]),float(q[1])] for q in P32["pointsC"] if q[2] not in namesA]; D=[[a,b] for a,b in zip(z,lr)]
print(f"   pair C: {len(P32['pointsC'])} clusters, {len(C)} not already in pair A")
def stack(parts):
    zz=np.concatenate([np.array(p)[:,0] for p in parts]); yy=np.concatenate([np.array(p)[:,1] for p in parts]); gg=np.concatenate([np.full(len(p),i) for i,p in enumerate(parts)])
    ww=np.concatenate([np.full(len(p),1/np.var(np.array(p)[:,1],ddof=1)) for p in parts]); return zz,yy,gg,ww
for lab,parts in (("A + B (measured profiles, Entry 28)",[A,Bp]),("A + B + Bolocam (no profile shape assumed)",[A,Bp,D]),("A + B + C (Entry 32)",[A,Bp,C]),("A + B + C + Bolocam (five pairs of instruments)",[A,Bp,C,D])):
    zz,yy,gg,ww=stack(parts); res[lab]=kfit(zz,yy,lab,w=ww,groups=gg,nboot=2000)
# binned points for the figure
edges=[0.14,0.22,0.30,0.40,0.46,0.55,0.60,0.90]; bins=[]
for a,b in zip(edges[:-1],edges[1:]):
    q=(z>=a)&(z<b); bins.append([float(z[q].mean()),float(lr[q].mean()-lr.mean()),float(lr[q].std(ddof=1)/q.sum()**0.5),int(q.sum())])
out["bins"]=bins; zz=np.linspace(0.15,0.9,16); zm=float(np.mean(z))
out["curves"]={nm:[[float(a),float(np.log10(F(a))-np.mean(np.log10(F(z))))] for a in zz] for nm,F in (("tired",F_log),("lerner",F_lin),("stretched",F_td))}
print("   binned:",[(round(b[0],2),round(b[1],3),round(b[2],3),b[3]) for b in bins])
out["fits"]=res
json.dump(out,open(os.path.join(HERE,"boxsz_test.json"),"w"),indent=1)
