"""Workbench Entry 34, Part 3: compact radio cores as a ruler (120 quasars, Cao et al. 2017).
theta x D_M(z) should go as (1+z)^x: x = 1 in expanding space if the cores have one physical size, x = 0 in static space with the measured
distance law. D_M is the ruler distance of Entry 33 (standard comoving distance, matter share 0.297); only its shape matters."""
import numpy as np, json, os
from scipy import integrate, optimize
HERE=os.path.dirname(os.path.abspath(__file__)); rng=np.random.default_rng(341)
R=[l.split() for l in open(os.path.join(HERE,"cao2017_quasars.txt")) if not l.startswith("#")]; assert len(R)==120
z=np.array([float(r[1]) for r in R]); th=np.array([float(r[2]) for r in R]); e=np.array([float(r[3]) for r in R]); S=np.array([float(r[4]) for r in R])
zg=np.linspace(0,3,3001); dm=integrate.cumulative_trapezoid(1/np.sqrt(0.297*(1+zg)**3+0.703),zg,initial=0); DM=np.interp(z,zg,dm)
y=np.log(th*DM); x=np.log(1+z); s=np.hypot(e/th,0.10)             # the catalogue error and the 10% the authors add
def fit(x,y,s,si):
    w=1/(s**2+si**2); X=np.c_[np.ones(len(x)),x]; cov=np.linalg.inv((X*w[:,None]).T@X); b=cov@(X*w[:,None]).T@y; r=y-X@b; return b,np.sqrt(np.diag(cov)),float(np.sum(w*r*r))
si=optimize.brentq(lambda v: fit(x,y,s,v)[2]-(len(x)-2),0,3); b,se,chi=fit(x,y,s,si)
bs=np.array([fit(x[j],y[j],s[j],si)[0][1] for j in (rng.integers(0,120,120) for _ in range(4000))])
print(f"120 quasars, z {z.min():.2f}-{z.max():.2f}; weighted fit: x = {b[1]:.2f} +- {max(se[1],bs.std()):.2f} (lsq {se[1]:.2f}, boot {bs.std():.2f}); intrinsic scatter {si:.2f} in ln size")
out={"n":120,"x":float(b[1]),"x_err":float(max(se[1],bs.std())),"sint":float(si)}
# unweighted, and median in bins (robust to the few very small, very uncertain sizes)
X=np.c_[np.ones(120),x]; bu=np.linalg.lstsq(X,y,rcond=None)[0]; bsu=np.array([np.linalg.lstsq(X[j],y[j],rcond=None)[0][1] for j in (rng.integers(0,120,120) for _ in range(4000))])
print(f"   unweighted: x = {bu[1]:.2f} +- {bsu.std():.2f}"); out["x_unw"]=[float(bu[1]),float(bsu.std())]
ok=e/th<0.3; b2,se2,_=fit(x[ok],y[ok],s[ok],optimize.brentq(lambda v: fit(x[ok],y[ok],s[ok],v)[2]-(ok.sum()-2),0,3)); print(f"   sizes known to 30% or better ({int(ok.sum())}): x = {b2[1]:.2f} +- {se2[1]:.2f}"); out["x_good"]=[float(b2[1]),float(se2[1]),int(ok.sum())]
order=np.argsort(z); bins=[]
for grp in np.array_split(order,8):
    m=np.median(th[grp]); em=1.2533*np.std(np.log(th[grp]),ddof=1)/np.sqrt(len(grp))*m; bins.append([float(np.median(z[grp])),float(m),float(em),int(len(grp))])
print("   median size (mas) in 8 bins: "+", ".join(f"z~{a:.2f}: {m_:.2f}+-{e_:.2f}" for a,m_,e_,n_ in bins)); out["bins"]=bins
zb=np.array([b_[0] for b_ in bins]); mb=np.array([b_[1] for b_ in bins]); eb=np.array([b_[2] for b_ in bins]); yb=np.log(mb*np.interp(zb,zg,dm)); b3,se3,c3=fit(np.log(1+zb),yb,eb/mb,0.0)
print(f"   from the bin medians: x = {b3[1]:.2f} +- {se3[1]:.2f} (chi2 {c3:.1f} for 6)"); out["x_bins"]=[float(b3[1]),float(se3[1])]
# flux density falls with redshift in this sample; does size follow flux at fixed redshift?
Xf=np.c_[np.ones(120),x,np.log(S)]; w=1/(s**2+si**2); cov=np.linalg.inv((Xf*w[:,None]).T@Xf); bf=cov@(Xf*w[:,None]).T@y
print(f"   with a flux term: x = {bf[1]:.2f} +- {np.sqrt(cov[1,1]):.2f}; size goes as flux^({bf[2]:+.2f} +- {np.sqrt(cov[2,2]):.2f})"); out["x_flux"]=[float(bf[1]),float(np.sqrt(cov[1,1])),float(bf[2]),float(np.sqrt(cov[2,2]))]
json.dump(out,open(os.path.join(HERE,"radio_core_ruler.json"),"w"),indent=1)
