"""Workbench Entry 46: which frame is at rest?  The Sun's motion measured from nearby supernovae alone, against the CMB dipole.

Input: 630 rows of the Pantheon+ table (Brout et al. 2022; VizieR J/ApJ/938/110/table7) with 0.01 < zHD < 0.1, typed in with
heliocentric and CMB-frame redshifts and positions.  For a trial velocity v of the Sun, each heliocentric redshift is moved to the
frame in which the Sun moves at v, and the supernova magnitudes are compared with a smooth Hubble curve at those redshifts.
The velocity that leaves the least scatter is the Sun's motion relative to the supernova hosts.  No CMB input is used.
The expected spread from large-scale flows in the standard model is computed from a linear power spectrum (Eisenstein & Hu 1998, no wiggles)."""
import os, json, sys
import numpy as np
from scipy import integrate, optimize, interpolate, special, stats
HERE=os.path.dirname(os.path.abspath(__file__)); C=299792.458
rows=[l.split() for l in open(os.path.join(HERE,"pantheon_plus_nearby_coords.txt"))]
name=np.array([r[0] for r in rows]); zhd,zc,zh,mb,em,ra,de=np.array([[float(x) for x in r[1:]] for r in rows]).T
def unit(lon,lat): lon,lat=np.radians(lon),np.radians(lat); return np.array([np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)])
E2G=np.array([[-0.0548755604,-0.8734370902,-0.4838350155],[0.4941094279,-0.4448296300,0.7469822445],[-0.8676661490,-0.1980763734,0.4559837762]])
def lb(v): g=E2G@v; a=np.linalg.norm(g); return a, np.degrees(np.arctan2(g[1],g[0]))%360, np.degrees(np.arcsin(g[2]/a))
VCMB=E2G.T@(369.82*unit(264.021,48.253))                       # Planck 2018 solar dipole, equatorial components, km/s
# ---- one entry per host: rows that share a redshift and lie within 0.1 degree are one host; magnitudes averaged, smallest error kept ----
N=unit(ra,de).T; grp=-np.ones(len(zh),int); g=0
for i in range(len(zh)):
    if grp[i]>=0: continue
    same=(abs(zhd-zhd[i])<1e-6)&(N@N[i]>np.cos(np.radians(0.1))); grp[same]=g; g+=1
H=[]
for k in range(g):
    q=grp==k; w=1/em[q]**2; H.append([zhd[q][0],zc[q][0],zh[q][0],np.sum(w*mb[q])/w.sum(),em[q].min(),*N[q][0],q.sum()])
H=np.array(H); zhd,zc,zh,mb,em=H[:,:5].T; N=H[:,5:8]; print(f"{len(rows)} rows, {g} hosts, {int((H[:,8]>1).sum())} of them with more than one light curve")
OM=0.334
zz=np.linspace(0,0.13,1400); DMt=integrate.cumulative_trapezoid(1/np.sqrt(OM*(1+zz)**3+1-OM),zz,initial=0); DM=interpolate.interp1d(zz,DMt,kind="cubic")
def model(p,zhel,n):                                           # p = M, vx, vy, vz
    zbar=(1+zhel)*(1+(n@p[1:])/C)-1; return p[0]+5*np.log10((1+zhel)*DM(zbar))
def fit(q,zhel=None,label="",quiet=False):
    zhel=zh if zhel is None else zhel
    f=lambda p:(mb[q]-model(p,zhel[q],N[q]))/em[q]
    r=optimize.least_squares(f,[np.median(mb[q]-5*np.log10(DM(zhel[q]))),0,0,0],x_scale=[0.01,100,100,100])
    chi=float(np.sum(r.fun**2)); dof=q.sum()-4; cov=np.linalg.inv(r.jac.T@r.jac)*max(1,chi/dof)
    return r.x,cov,chi,dof
def describe(v,cov,tag):
    a,l,b=lb(v); ea=np.sqrt(v@cov@v)/a; rng=np.random.default_rng(1); S=rng.multivariate_normal(v,cov,20000)
    ang=np.degrees(np.arccos(np.clip((S@VCMB)/np.linalg.norm(S,axis=1)/369.82,-1,1))); a0=np.degrees(np.arccos(np.clip(v@VCMB/a/369.82,-1,1)))
    d=v-VCMB; chi=float(d@np.linalg.solve(cov,d))
    print(f"{tag}: {a:4.0f} +- {ea:3.0f} km/s toward l = {l:5.1f}, b = {b:+5.1f}; {a0:4.1f} deg from the CMB apex (68% of draws within {np.percentile(ang,68):.0f} deg); difference from the CMB velocity {np.linalg.norm(d):4.0f} km/s, chi2 {chi:.1f} for 3 (noise only)")
    return dict(amp=float(a),err=float(ea),l=float(l),b=float(b),angle=float(a0),diff=float(np.linalg.norm(d)),chi_noise=chi)
out={"hosts":int(g)}
allq=np.ones(g,bool)
print("\n--- the Sun's motion from supernovae alone (heliocentric redshifts; no CMB input) ---")
p,cov,chi,dof=fit(allq); print(f"all hosts 0.01-0.1: chi2 {chi:.0f} for {dof}; errors scaled by {max(1,chi/dof)**0.5:.2f}")
out["all"]=describe(p[1:],cov[1:,1:],"all hosts, 0.01 < z < 0.1"); out["all"]["chi"]=[chi,int(dof)]; out["all"]["v_eq"]=[float(x) for x in p[1:]]; out["all"]["cov"]=cov[1:,1:].tolist()
# how strongly do the supernovae choose a frame?
def chi_at(v,q=allq,zhel=None):
    zhel=zh if zhel is None else zhel
    d=mb[q]-model(np.r_[0.0,v],zhel[q],N[q]); w=1/em[q]**2; M=np.sum(w*d)/w.sum(); return float(np.sum(w*(d-M)**2))
c_hel=chi_at(np.zeros(3)); c_cmb=chi_at(VCMB); c_best=chi
print(f"chi2 with the Sun at rest {c_hel:.1f}; with the CMB velocity {c_cmb:.1f}; best {c_best:.1f}.  The CMB velocity gains {c_hel-c_cmb:.1f}; the best fit gains {c_cmb-c_best:.1f} more for 3 numbers")
# errors from the data: resample hosts
rng0=np.random.default_rng(46); B=[]
for _ in range(600):
    j=rng0.integers(0,g,g); f_=lambda p:(mb[j]-model(p,zh[j],N[j]))/em[j]; B.append(optimize.least_squares(f_,p,x_scale=[0.01,100,100,100]).x[1:])
B=np.array(B); covB=np.cov(B.T); ampB=np.linalg.norm(B,axis=1); angB=np.degrees(np.arccos(np.clip((B@VCMB)/ampB/369.82,-1,1))); dB=np.linalg.norm(B-VCMB,axis=1)
print(f"resampling the hosts 600 times: speed {ampB.mean():.0f} +- {ampB.std():.0f} km/s; angle from the CMB apex 68% within {np.percentile(angB,68):.0f} deg, 95% within {np.percentile(angB,95):.0f}; difference from the CMB velocity {dB.mean():.0f} +- {dB.std():.0f}; chi2 of the difference with these errors {(p[1:]-VCMB)@np.linalg.solve(covB,p[1:]-VCMB):.1f} for 3; speed is zero at {np.linalg.norm(p[1:])/ampB.std():.1f} sigma, chi2 {p[1:]@np.linalg.solve(covB,p[1:]):.1f} for 3")
out["boot"]=dict(amp_err=float(ampB.std()),ang68=float(np.percentile(angB,68)),ang95=float(np.percentile(angB,95)),cov=covB.tolist(),chi_diff=float((p[1:]-VCMB)@np.linalg.solve(covB,p[1:]-VCMB)),chi_zero=float(p[1:]@np.linalg.solve(covB,p[1:])))
out["chi"]=dict(sun_at_rest=float(c_hel),cmb=float(c_cmb),best=float(c_best))
# the same with the published flow corrections (zHD: CMB frame plus a model of local flows): what is left over should be zero
zhd_as_hel=zhd.copy(); pH,covH,chiH,dofH=fit(allq,zhel=zhd_as_hel)
a,l,b=lb(pH[1:]); print(f"check with the published redshifts already in the CMB frame and corrected for modelled flows: left-over velocity {a:.0f} km/s (chi2 {pH[1:]@np.linalg.solve(covH[1:,1:],pH[1:]):.1f} for 3); chi2 {chiH:.0f}")
out["zhd_left"]=float(a)
print("\n--- by depth ---"); out["shells"]=[]
for lo,hi in ((0.01,0.02),(0.02,0.03),(0.03,0.05),(0.05,0.1),(0.01,0.03),(0.03,0.1),(0.023,0.1),(0.015,0.1)):
    q=(zhd>=lo)&(zhd<hi); p2,cov2,chi2,dof2=fit(q); d=describe(p2[1:],cov2[1:,1:],f"z {lo}-{hi} ({q.sum():3d} hosts, median z {np.median(zhd[q]):.3f})"); d.update(lo=lo,hi=hi,n=int(q.sum()),zmed=float(np.median(zhd[q])),v_eq=[float(x) for x in p2[1:]],cov=cov2[1:,1:].tolist()); out["shells"].append(d)

# ---- what the standard model expects: covariance of the fitted velocity from large-scale flows ----
h,Om,Ob,ns,s8=0.674,0.315,0.049,0.965,0.811
def T_nw(k):                                                   # k in h/Mpc; Eisenstein & Hu 1998 no-wiggle transfer function
    om,ob=Om*h*h,Ob*h*h; fb=Ob/Om; th=2.7255/2.7; s=44.5*np.log(9.83/om)/np.sqrt(1+10*ob**0.75)
    al=1-0.328*np.log(431*om)*fb+0.38*np.log(22.3*om)*fb**2; kk=k*h; G=Om*h*(al+(1-al)/(1+(0.43*kk*s)**4)); q=k*th**2/G
    L=np.log(2*np.e+1.8*q); Cq=14.2+731/(1+62.5*q); return L/(L+Cq*q*q)
k=np.logspace(-4,1.3,6000); P=k**ns*T_nw(k)**2
W=lambda x:3*(np.sin(x)-x*np.cos(x))/x**3; P*=s8**2/np.trapezoid(P*k**3*W(8*k)**2/(2*np.pi**2),np.log(k))
f=Om**0.55; pref=(100*f)**2/(2*np.pi**2)
rg=np.r_[0,np.logspace(-1,3.2,500)]; j0=lambda x:np.sinc(x/np.pi); j1x=lambda x:np.where(x<1e-3,1/3-x*x/30,(np.sin(x)-x*np.cos(x))/np.where(x<1e-3,1,x)**3)
psi_perp=np.array([pref*np.trapezoid(P*k*j1x(k*r),np.log(k)) for r in rg]); psi_par=np.array([pref*np.trapezoid(P*k*(j0(k*r)-2*j1x(k*r)),np.log(k)) for r in rg])
print(f"\nlinear theory: one-dimensional velocity spread {psi_perp[0]**0.5:.0f} km/s; correlation falls to half by {rg[np.argmax(psi_perp<psi_perp[0]/2)]:.0f} Mpc/h (across) ")
Pp=interpolate.interp1d(rg,psi_perp,bounds_error=False,fill_value=0.0); Pl=interpolate.interp1d(rg,psi_par,bounds_error=False,fill_value=0.0)
def cosmic_cov(q,extra_sigma=0.0):
    n=N[q]; X=(2997.92458*DM(zc[q]))[:,None]*n; d=X[:,None,:]-X[None,:,:]; r=np.linalg.norm(d,axis=2); rh=d/np.where(r>0,r,1)[:,:,None]
    a=np.einsum("ik,ijk->ij",n,rh); b_=np.einsum("jk,ijk->ij",n,rh); Cs=Pp(r)*(n@n.T)+(Pl(r)-Pp(r))*a*b_; np.fill_diagonal(Cs,0.0)   # cross terms only: each host's own random velocity is already in its error
    zq=zc[q]; gmag=(5/np.log(10))*(1/DM(zq))*np.gradient(DMt,zz)[np.searchsorted(zz,zq)]*(1+zq)/C       # magnitudes per km/s of line-of-sight velocity
    A=np.c_[np.ones(q.sum()),gmag[:,None]*n]; Ninv=1/em[q]**2; F=np.linalg.inv(A.T@(Ninv[:,None]*A)); Pm=F@(A.T*Ninv)
    Cm=(gmag[:,None]*gmag[None,:])*Cs; return (Pm@Cm@Pm.T)[1:,1:], F[1:,1:]
print("\n--- is the supernova frame the CMB frame, allowing for flows of the size the standard model expects? ---"); out["flow"]=[]
for lo,hi in ((0.01,0.1),(0.015,0.1),(0.023,0.1),(0.01,0.03),(0.03,0.1),(0.05,0.1)):
    q=(zhd>=lo)&(zhd<hi); p2,cov2,chi2,dof2=fit(q); u=VCMB-p2[1:]; Cc,Fn=cosmic_cov(q); tot=cov2[1:,1:]+Cc
    a,l,b=lb(u); chi_u=float(u@np.linalg.solve(tot,u)); pv=float(stats.chi2.sf(chi_u,3)); exp=float(np.sqrt(np.trace(Cc))); noi=float(np.sqrt(np.trace(cov2[1:,1:])))
    print(f"z {lo}-{hi} ({q.sum():3d}): flow of the hosts in the CMB frame {a:4.0f} km/s toward l = {l:5.1f}, b = {b:+5.1f}; measurement error {noi:3.0f}, expected flow {exp:3.0f} (three-dimensional rms); chi2 {chi_u:4.1f} for 3, p = {pv:.2f}")
    out["flow"].append(dict(lo=lo,hi=hi,n=int(q.sum()),amp=float(a),l=float(l),b=float(b),noise=noi,expected=exp,chi=chi_u,p=pv))
    if (lo,hi)==(0.01,0.1):
        c2=float(u@np.linalg.solve(covB+Cc,u)); print(f"      with the resampled measurement errors instead of the table's: chi2 {c2:.1f} for 3, p = {stats.chi2.sf(c2,3):.2f}; measurement error {np.sqrt(np.trace(covB)):.0f} km/s"); out["flow_boot"]=[c2,float(stats.chi2.sf(c2,3))]
# ---- could the Sun be moving faster than the CMB dipole says, as the source-count dipoles would need? ----
print("\n--- the Sun's speed along the CMB direction, in units of the CMB speed ---"); out["k"]=[]
nc=VCMB/369.82
for lo,hi in ((0.01,0.1),(0.023,0.1),(0.03,0.1)):
    q=(zhd>=lo)&(zhd<hi); p2,cov2,chi2,dof2=fit(q); Cc,Fn=cosmic_cov(q); tot=cov2[1:,1:]+Cc
    kk=float(p2[1:]@nc/369.82); en=float(np.sqrt(nc@cov2[1:,1:]@nc)/369.82); et=float(np.sqrt(nc@tot@nc)/369.82); line=f"z {lo}-{hi}: {kk:.2f} +- {en:.2f} (measurement) +- {et:.2f} (with flows)"
    res=dict(lo=lo,hi=hi,k=kk,e_noise=en,e_total=et,tests=[])
    for kt,lab in ((1.0,"CMB"),(2.0,"twice"),(2.2,"quasar counts"),(3.67,"radio counts")):
        d=p2[1:]-kt*VCMB; c=float(d@np.linalg.solve(tot,d)); sg=float(stats.norm.isf(stats.chi2.sf(c,3)/2)) if c<1400 else float("inf"); line+=f" | {lab} ({kt}): chi2 {c:.1f}, {sg:.1f} sigma"; res["tests"].append([kt,c,sg])
    print(line); out["k"].append(res)
# where the hosts are heading, and the Shapley concentration (Abell 3558: RA 201.99, Dec -31.50)
sh=unit(201.99,-31.50); a_,l_,b_=lb(sh); u=VCMB-np.array(out["all"]["v_eq"]); print(f"Shapley concentration at l = {l_:.1f}, b = {b_:+.1f}; the hosts' flow points {np.degrees(np.arccos(u@sh/np.linalg.norm(u))):.0f} deg from it")
out["shapley"]=[float(l_),float(b_),float(np.degrees(np.arccos(u@sh/np.linalg.norm(u))))]
# sky coverage of the hosts
g_=(E2G@N.T).T; bl=np.degrees(np.arcsin(g_[:,2])); print(f"hosts north of the Galactic plane: {np.mean(bl>0)*100:.0f}%; median distance {np.median(2997.92458*DM(zc)):.0f} Mpc/h; within 20 deg of the plane: {np.mean(abs(bl)<20)*100:.0f}%; hosts in the hemisphere of the CMB apex: {np.mean(N@nc>0)*100:.0f}%")
json.dump(out,open(os.path.join(HERE,"sn_frame.json"),"w"),indent=1)
