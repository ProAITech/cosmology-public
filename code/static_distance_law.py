"""Workbench Entry 33: the static universe's distance law, measured instead of assumed.
Earlier entries tested particular static distance laws (d = cz/H0, d = (c/H0) ln(1+z)). A defender can always propose another law.
This removes that freedom in three steps, each using published measurements only.
 1. In a static space a clustering feature of fixed length r subtends r/d(z) across the line of sight and spans r/(dd/dz) in redshift along it.
    The DESI DR2 baryon-acoustic measurements (Table IV of arXiv:2503.14738) are exactly these two numbers at six redshifts, so they measure
    d(z)/r and its slope with no model for d.
 2. Supernovae measure d(z) (1+z)^n, where n counts the dimming beyond the inverse-square law: n = 1/2 if photons only lose energy,
    n = 1 if they also arrive stretched in time. Dividing supernova distances by the ruler distances measures n with no model for d.
 3. With d and n measured nothing is left to adjust, and the cluster test of Entries 27, 28 and 32 becomes a prediction:
    measured/predicted shadow must go as (1+z)^(n - 3/2) in static space, against (1+z)^0 in expanding space.
Supernovae: JLA binned distances with their full covariance (static_hubble.py); Pantheon with statistical errors (pantheon_rows.py)."""
import sys, os, io, json, contextlib
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
import numpy as np
from scipy import integrate, optimize
_cwd=os.getcwd(); os.chdir(HERE)
with contextlib.redirect_stdout(io.StringIO()): import static_hubble as sh
os.chdir(_cwd)
from pantheon_rows import z as zp, m as mp
out={}
# ---------------- DESI DR2 BAO: z, D_M/r, err, D_H/r, err, correlation
BAO=[(0.510,13.588,0.167,21.863,0.425,-0.459),(0.706,17.351,0.177,19.455,0.330,-0.404),(0.934,21.576,0.152,17.641,0.193,-0.416),
     (1.321,27.601,0.318,14.176,0.221,-0.434),(1.484,30.512,0.760,12.817,0.516,-0.500),(2.330,38.988,0.531,8.632,0.101,-0.431)]
zV,DV,eDV=0.295,7.942,0.075
zb=np.array([b[0] for b in BAO]); dM=np.array([b[1] for b in BAO]); eM=np.array([b[2] for b in BAO]); dH=np.array([b[3] for b in BAO]); eH=np.array([b[4] for b in BAO]); rho=np.array([b[5] for b in BAO])
y=np.r_[DV,dM,dH]; C=np.zeros((13,13)); C[0,0]=eDV**2
for i in range(6):
    C[1+i,1+i]=eM[i]**2; C[7+i,7+i]=eH[i]**2; C[1+i,7+i]=C[7+i,1+i]=rho[i]*eM[i]*eH[i]
Ci=np.linalg.inv(C)
zg=np.linspace(0,2.6,5201)
def curve(slope):
    """slope(z) = dd/dz in units of c/H0; returns d(z) and dd/dz as functions"""
    s=slope(zg); d=integrate.cumulative_trapezoid(s,zg,initial=0); return (lambda t: np.interp(t,zg,d)),(lambda t: np.interp(t,zg,s))
def template(slope):
    d,s=curve(slope); return np.r_[(zV*d(zV)**2*s(zV))**(1/3),d(zb),s(zb)]
def chi_bao(slope):
    t=template(slope); A=(t@Ci@y)/(t@Ci@t); r=y-A*t; return float(r@Ci@r),float(A),r
h_poly=lambda a: (lambda t: 1/np.sqrt(np.clip(1+a[0]*t+a[1]*t*t+a[2]*t**3,1e-6,None)))
h_std=lambda om: (lambda t: 1/np.sqrt(om*(1+t)**3+1-om))
h_pow=lambda p: (lambda t: (1+t)**-p)
print("Step 1. The ruler distance law, from 13 DESI numbers (one free scale for each law):")
laws={}
c,A,r=chi_bao(h_pow(0.0)); laws["straight line, d = cz/H0"]=(c,0); print(f"   straight line, d proportional to z:            chi2 = {c:8.1f} (12 d.o.f.)")
c,A,r=chi_bao(h_pow(1.0)); laws["logarithmic, d = (c/H0) ln(1+z)"]=(c,0); print(f"   logarithmic, d proportional to ln(1+z):        chi2 = {c:8.1f} (12 d.o.f.)")
g=optimize.minimize_scalar(lambda p: chi_bao(h_pow(p))[0],bounds=(0,3),method="bounded"); plo=optimize.brentq(lambda p: chi_bao(h_pow(p))[0]-g.fun-1,0,g.x); phi=optimize.brentq(lambda p: chi_bao(h_pow(p))[0]-g.fun-1,g.x,3)
laws["power law"]=(g.fun,1); print(f"   slope falling as (1+z)^-p, p free:             chi2 = {g.fun:8.1f} (11 d.o.f.), p = {g.x:.3f} +{phi-g.x:.3f} -{g.x-plo:.3f}")
out["p_pow"]=[float(g.x),float(g.x-plo),float(phi-g.x),float(g.fun)]
gs=optimize.minimize_scalar(lambda o: chi_bao(h_std(o))[0],bounds=(0.05,0.9),method="bounded"); olo=optimize.brentq(lambda o: chi_bao(h_std(o))[0]-gs.fun-1,0.05,gs.x); ohi=optimize.brentq(lambda o: chi_bao(h_std(o))[0]-gs.fun-1,gs.x,0.9)
laws["standard comoving distance"]=(gs.fun,1); print(f"   the standard model's comoving distance:        chi2 = {gs.fun:8.1f} (11 d.o.f.), matter share {gs.x:.3f} +{ohi-gs.x:.3f} -{gs.x-olo:.3f}   (DESI's own fit: 0.297 +- 0.009, chi2 10.2)")
out["om_bao"]=[float(gs.x),float(gs.x-olo),float(ohi-gs.x),float(gs.fun)]
gp=optimize.minimize(lambda a: chi_bao(h_poly(a))[0],[0.9,0.9,0.3],method="Nelder-Mead",options=dict(xatol=1e-5,fatol=1e-6,maxiter=4000))
laws["free curve, three shape numbers"]=(gp.fun,3); print(f"   free curve, slope 1/sqrt(1 + a z + b z2 + c z3): chi2 = {gp.fun:8.1f} (9 d.o.f.), a,b,c = {gp.x[0]:.2f}, {gp.x[1]:.2f}, {gp.x[2]:.2f}  (standard with matter 0.3: 0.90, 0.90, 0.30)")
out["laws"]={k:[float(v[0]),v[1]] for k,v in laws.items()}; out["poly_bao"]=[float(v) for v in gp.x]
# ruler-free ratios, with errors, against each law
dstd,sstd=curve(h_std(gs.x))
def ratio(a,ea,b,eb): return a/b, a/b*np.hypot(ea/a,eb/b)
r1=ratio(dM[5],eM[5],dM[0],eM[0]); r2=ratio(dH[0],eH[0],dH[5],eH[5])
print(f"   transverse: d(2.33)/d(0.51) measured {r1[0]:.3f} +- {r1[1]:.3f}; straight line {2.33/0.51:.3f}; logarithmic {np.log(3.33)/np.log(1.51):.3f}; standard {dstd(2.33)/dstd(0.51):.3f}")
print(f"   radial: slope(0.51)/slope(2.33) measured {r2[0]:.3f} +- {r2[1]:.3f}; straight line 1.000; logarithmic {3.33/1.51:.3f}; standard {sstd(0.51)/sstd(2.33):.3f}")
out["ratios"]=dict(trans=[float(r1[0]),float(r1[1]),2.33/0.51,float(np.log(3.33)/np.log(1.51)),float(dstd(2.33)/dstd(0.51))],radial=[float(r2[0]),float(r2[1]),1.0,3.33/1.51,float(sstd(0.51)/sstd(2.33))])
# how far the measured law is from the standard comoving distance: per-point ratios
A_std=chi_bao(h_std(gs.x))[1]
print("   measured d/r over the standard curve at its best scale: "+", ".join(f"z={a:.2f}: {b/(A_std*dstd(a)):.3f}+-{c/b:.3f}" for a,b,c in zip(zb,dM,eM)))
# ---------------- Step 2. supernova distance over ruler distance
print("\nStep 2. The dimming exponent n in  d_L = d (1+z)^n  (1/2: energy loss only; 1: energy loss and stretched time, or expansion)")
# (a) JLA, point by point at the BAO redshifts. JLA's binned distances are a piecewise-linear function of log z between its 31 nodes, so the
#     value at any redshift is a fixed linear combination of two nodes, and the covariance carries over exactly.
def W_at(zq):
    W=np.zeros((len(zq),31)); lz=np.log10(sh.zb)
    for i,q in enumerate(zq):
        j=min(max(np.searchsorted(sh.zb,q)-1,0),29); f=(np.log10(q)-lz[j])/(lz[j+1]-lz[j]); W[i,j]=1-f; W[i,j+1]=f
    return W
def n_point(zq,mu,Cmu,dq,edq,label):
    yv=mu-5*np.log10(dq); Cy=Cmu+np.diag((5/np.log(10)*edq/dq)**2); X=np.c_[np.ones(len(zq)),5*np.log10(1+zq)]; Ciy=np.linalg.inv(Cy)
    cov=np.linalg.inv(X.T@Ciy@X); b=cov@X.T@Ciy@yv; r=yv-X@b; chi=float(r@Ciy@r)
    print(f"   {label:58s} n = {b[1]:.3f} +- {np.sqrt(cov[1,1]):.3f}   (chi2 {chi:.1f} for {len(zq)-2} d.o.f.)")
    return float(b[1]),float(np.sqrt(cov[1,1])),chi,yv-b[0],np.sqrt(np.diag(Cy))
sel=[0,1,2,3]; W=W_at(zb[sel]); mu_j=W@sh.mub; C_j=W@sh.C@W.T
nj=n_point(zb[sel],mu_j,C_j,dM[sel],eM[sel],"JLA (full errors) at z = 0.51, 0.71, 0.93, 1.32")
sel3=[0,1,2]; W3=W_at(zb[sel3]); nj3=n_point(zb[sel3],W3@sh.mub,W3@sh.C@W3.T,dM[sel3],eM[sel3],"JLA (full errors) at z = 0.51, 0.71, 0.93 only")
# (b) Pantheon, statistical errors: mean magnitude in a window around each BAO redshift, moved to that redshift with the local slope of a smooth curve
x=np.log10(1+zp); base=5*np.log10(zp); Avd=np.vander(x,5)
edges=np.quantile(zp,np.linspace(0,1,15)); mid=np.array([np.median(zp[(zp>=a)&(zp<=b)]) for a,b in zip(edges[:-1],edges[1:])]); sig=np.full(len(zp),0.15)
for it in range(4):
    cf=np.linalg.lstsq(Avd/sig[:,None],(mp-base)/sig,rcond=None)[0]; res=mp-base-Avd@cf
    sig=np.interp(np.log(zp),np.log(mid),np.array([res[(zp>=a)&(zp<=b)].std(ddof=1) for a,b in zip(edges[:-1],edges[1:])]))
smooth=lambda t: 5*np.log10(t)+np.vander(np.log10(1+t),5)@cf
wins=[(0.455,0.571),(0.630,0.790),(0.834,1.046),(1.18,1.42),(1.42,1.75),(2.0,2.6)]; mu_p=[];e_p=[];n_in=[]
for q,(a,b) in zip(zb,wins):
    s=(zp>=a)&(zp<b); w=1/sig[s]**2; mu_p.append(np.sum(w*(mp[s]-(smooth(zp[s])-smooth(np.array([q])))))/w.sum()); e_p.append(1/np.sqrt(w.sum())); n_in.append(int(s.sum()))
mu_p=np.array(mu_p); e_p=np.array(e_p); print("   Pantheon supernovae in each window:",n_in,"; errors (mag):"," ".join(f"{v:.3f}" for v in e_p))
npn=n_point(zb,mu_p,np.diag(e_p**2),dM,eM,"Pantheon (statistical errors) at all six redshifts")
npn4=n_point(zb[:4],mu_p[:4],np.diag(e_p[:4]**2),dM[:4],eM[:4],"Pantheon (statistical errors) at z = 0.51 to 1.32")
out["n_point"]=dict(jla=nj[:3],jla3=nj3[:3],pantheon=npn[:3],pantheon4=npn4[:3])
out["fig_n"]=dict(z=zb.tolist(),jla=[[float(a),float(b)] for a,b in zip(nj[3],nj[4])],pantheon=[[float(a),float(b)] for a,b in zip(npn[3],npn[4])],n_jla=nj[0],n_pan=npn[0])
# (c) joint fits with a free curve for d(z): all supernovae and all 13 ruler numbers together
one=np.ones(31)
def chi_jla(dfun,n):
    mod=5*np.log10(dfun(sh.zb)*(1+sh.zb)**n); d=sh.mub-mod; a=(one@sh.Ci@d)/(one@sh.Ci@one); r=d-a; return float(r@sh.Ci@r)
wp=1/sig**2
def chi_pan(dfun,n):
    mod=5*np.log10(dfun(zp)*(1+zp)**n); d=mp-mod; a=np.sum(wp*d)/wp.sum(); return float(np.sum(wp*(d-a)**2))
def joint(chi_sn,shape,x0,label,nfix=None):
    def f(par):
        n=par[-1] if nfix is None else nfix; sl=shape(par[:-1] if nfix is None else par); d,_=curve(sl); return chi_bao(sl)[0]+chi_sn(d,n)
    g=optimize.minimize(f,x0 if nfix is None else x0[:-1],method="Nelder-Mead",options=dict(xatol=1e-5,fatol=1e-6,maxiter=8000,maxfev=8000)); return g
def profile_n(chi_sn,shape,x0,label):
    g=joint(chi_sn,shape,x0,label); nb=g.x[-1]
    prof=lambda n: joint(chi_sn,shape,list(g.x),label,nfix=n).fun-g.fun-1
    lo=optimize.brentq(prof,nb-0.6,nb); hi=optimize.brentq(prof,nb,nb+0.6)
    c05=joint(chi_sn,shape,list(g.x),label,nfix=0.5).fun; c10=joint(chi_sn,shape,list(g.x),label,nfix=1.0).fun
    print(f"   {label:58s} n = {nb:.3f} +{hi-nb:.3f} -{nb-lo:.3f}   chi2 {g.fun:.1f}; with n fixed at 1: {c10:.1f}; at 1/2: {c05:.1f} (worse by {c05-g.fun:.1f})")
    return dict(n=float(nb),lo=float(nb-lo),hi=float(hi-nb),chi=float(g.fun),chi_n1=float(c10),chi_n05=float(c05),par=[float(v) for v in g.x])
J={}
J["jla_free"]=profile_n(chi_jla,h_poly,[0.9,0.9,0.3,1.0],"JLA (full errors) + ruler, free curve")
J["jla_std"]=profile_n(chi_jla,lambda a: h_std(a[0]),[0.3,1.0],"JLA (full errors) + ruler, standard-shaped curve")
J["pan_free"]=profile_n(chi_pan,h_poly,[0.9,0.9,0.3,1.0],"Pantheon (statistical) + ruler, free curve")
J["pan_std"]=profile_n(chi_pan,lambda a: h_std(a[0]),[0.3,1.0],"Pantheon (statistical) + ruler, standard-shaped curve")
out["joint"]=J
# a static model with its own named law and its own n, against supernovae and ruler together (JLA, full errors)
print("\n   Named static models against ruler and supernovae together (JLA with full errors; best free-curve fit for comparison):")
best=J["jla_free"]["chi"]; out["named"]={}
for lab,sl,n in (("Lerner: d = cz/H0, energy loss only",h_pow(0.0),0.5),("d = cz/H0 with stretched time",h_pow(0.0),1.0),("tired light: d = ln(1+z), energy loss only",h_pow(1.0),0.5),("d = ln(1+z) with stretched time",h_pow(1.0),1.0),("standard comoving distance law with energy loss only",h_std(gs.x),0.5),("standard comoving distance law with stretched time",h_std(gs.x),1.0)):
    d,_=curve(sl); cb=chi_bao(sl)[0]; cs=chi_jla(d,n); out["named"][lab]=[cb,cs,cb+cs]; print(f"     {lab:52s} ruler {cb:8.1f} + supernovae {cs:7.1f} = {cb+cs:8.1f}   (free curve: {best:.1f})")
# ---------------- Step 3. the cluster test as a prediction
print("\nStep 3. Clusters. With d equal to the measured ruler distance (the standard comoving distance D_M = (1+z) D_A) and dimming n, static space")
print("        predicts measured/predicted shadow going as (1+z)^(n - 3/2): -1/2 for n = 1, -1 for n = 1/2. Expanding space predicts 0.")
P28=json.load(open(os.path.join(HERE,"entry28_points.json"))); P32=json.load(open(os.path.join(HERE,"entry32_points.json")))["pointsC"]
zA=np.array([p[0] for p in P28["pointsA"]]); yA=np.array([p[1] for p in P28["pointsA"]]); zB=np.array([p[0] for p in P28["pointsB"]]); yB=np.array([p[1] for p in P28["pointsB"]])
newC=[p for p in P32 if not p[4]]; rng=np.random.default_rng(33)
def power(sets,label):
    zz=np.concatenate([s[0] for s in sets]); yy=np.concatenate([s[1] for s in sets]); D=np.zeros((len(zz),len(sets))); k=0
    for i,s in enumerate(sets): D[k:k+len(s[0]),i]=1; k+=len(s[0])
    X=np.c_[D,np.log10(1+zz)]; b=np.linalg.lstsq(X,yy,rcond=None)[0]; r=yy-X@b; s2=(r*r).sum()/(len(yy)-X.shape[1]); se=np.sqrt(s2*np.linalg.inv(X.T@X)[-1,-1])
    bs=np.array([np.linalg.lstsq(X[j],yy[j],rcond=None)[0][-1] for j in (rng.integers(0,len(yy),len(yy)) for _ in range(3000))]); e=max(se,bs.std())
    print(f"   {label:46s} {len(zz)} clusters: (1+z)^({b[-1]:+.2f} +- {e:.2f});  static with n = 1 is {(b[-1]+0.5)/e:.1f} sigma away, with n = 1/2 {(b[-1]+1.0)/e:.1f}")
    return dict(n=int(len(zz)),p=float(b[-1]),e=float(e))
cl={}
cl["AB"]=power([(zA,yA),(zB,yB)],"Entry 28 (pairs A and B)")
cl["ABCs"]=power([(zA,yA),(zB,yB),(np.array([p[0] for p in newC if p[3]>=8]),np.array([p[1] for p in newC if p[3]>=8]))],"with the strong detections of Entry 32")
cl["ABC"]=power([(zA,yA),(zB,yB),(np.array([p[0] for p in newC]),np.array([p[1] for p in newC]))],"with all of Entry 32")
out["clusters"]=cl
# the two measurements of the dimming, on one scale: n from supernovae over ruler, and n needed by the clusters (n_cl = p + 3/2)
n_sn=J["jla_free"]["n"]; e_sn=max(J["jla_free"]["lo"],J["jla_free"]["hi"])
print("\n   Dimming exponent that a static universe needs, from two sides:")
for k,lab in (("AB","pairs A and B"),("ABCs","with strong Entry 32"),("ABC","with all of Entry 32")):
    ncl=cl[k]["p"]+1.5; print(f"     supernovae over ruler: {n_sn:.2f} +- {e_sn:.2f};  clusters ({lab}): {ncl:.2f} +- {cl[k]['e']:.2f};  apart by {(ncl-n_sn)/np.hypot(e_sn,cl[k]['e']):.1f} sigma")
    cl[k]["n_needed"]=float(ncl); cl[k]["apart"]=float((ncl-n_sn)/np.hypot(e_sn,cl[k]['e']))
# ---------------- figure data
zz=np.linspace(0.05,2.5,50); fig={"z":zz.tolist(),"dM":[[float(a),float(b),float(c)] for a,b,c in zip(zb,dM,eM)],"dH":[[float(a),float(b),float(c)] for a,b,c in zip(zb,dH,eH)]}
for key,sl in (("lin",h_pow(0.0)),("log",h_pow(1.0)),("std",h_std(gs.x))):
    A=chi_bao(sl)[1]; d,s=curve(sl); fig[key]={"d":(A*d(zz)).tolist(),"s":(A*s(zz)).tolist()}
out["fig"]=fig
# ratios for the power law at its best p, and figure arrays relative to the standard curve
pb=out["p_pow"][0]; dpw,spw=curve(h_pow(pb)); out["ratios"]["pow"]=[float(dpw(2.33)/dpw(0.51)),float(spw(0.51)/spw(2.33))]
print(f"\n   power law at p = {pb:.3f}: transverse ratio {dpw(2.33)/dpw(0.51):.3f}, radial ratio {spw(0.51)/spw(2.33):.3f}")
A_s=chi_bao(h_std(gs.x))[1]; zz2=np.linspace(0.2,2.5,47)
rel={"z":[round(float(v),3) for v in zz2],"dM":[[float(a),round(float(b/(A_s*dstd(a))),4),round(float(c/b),4)] for a,b,c in zip(zb,dM,eM)],"dH":[[float(a),round(float(b/(A_s*sstd(a))),4),round(float(c/b),4)] for a,b,c in zip(zb,dH,eH)],
     "dV":[zV,round(float(DV/(A_s*(zV*dstd(zV)**2*sstd(zV))**(1/3))),4),round(eDV/DV,4)]}
for key,sl in (("lin",h_pow(0.0)),("log",h_pow(1.0)),("pow",h_pow(pb))):
    A=chi_bao(sl)[1]; d,sv=curve(sl); rel[key]={"d":[round(float(v),4) for v in A*d(zz2)/(A_s*dstd(zz2))],"s":[round(float(v),4) for v in A*sv(zz2)/(A_s*sstd(zz2))]}
out["rel"]=rel
# supernova over ruler, relative to n = 1, each set at its own level; and the n = 1/2 line at its best level through the Pantheon points
fn={}
for key,(nres,zq) in (("jla",(nj,zb[[0,1,2,3]])),("pantheon",(npn,zb))):
    v=np.array(nres[3])-5*np.log10(1+zq)*1.0; e=np.array(nres[4]); v=v+5*np.log10(1+zq)*0  # nres[3] is y minus fitted intercept
    w=1/e**2; v=v-np.sum(w*v)/w.sum(); fn[key]=[[float(a),round(float(b),3),round(float(c),3)] for a,b,c in zip(zq,v,e)]
zq=zb; e=np.array(npn[4]); w=1/e**2; v=np.array([p[1] for p in fn["pantheon"]]); line=-2.5*np.log10(1+zz2); c0=np.sum(w*(v-(-2.5*np.log10(1+zq))))/w.sum()
fn["half"]=[[round(float(a),3),round(float(b+c0),3)] for a,b in zip(zz2,line)]
out["fig_n2"]=fn
print("   supernova over ruler relative to n = 1:", fn["jla"], fn["pantheon"]); print("   n = 1/2 line at z = 0.5, 1.3, 2.3:",[p for p in fn["half"] if p[0] in (0.5,1.3,2.3)])
json.dump(out,open(os.path.join(HERE,"static_distance_law.json"),"w"),indent=1)
