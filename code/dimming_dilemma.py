"""Workbench Entry 31, part 3: is there any dimming law that saves a static distance?
In static flat space with a distance law d(z), let photons lose energy as (1+z)^-1 and let their arrival rate fall as (1+z)^-j
(j = 0: nothing but energy loss; j = 1: arrival stretched in time as supernova light curves are).
 - Supernovae: flux ~ 1/(d^2 (1+z)^(1+j)), so brightness distance = d (1+z)^((1+j)/2).
 - Clusters: X-ray photons per area of sky fall as (1+z)^-j, so the shadow predicted from the X-rays with standard assumptions is off by
   (d/D_A)^(1/2) (1+z)^((j-3)/2), and the ratio measured/predicted follows that (Entries 27, 28).
Each data set gives its own j. A consistent static model needs them to agree."""
import numpy as np, json, sys, os
from scipy import integrate, optimize
here=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,here)
import static_hubble as sh                      # runs the supernova fits and exposes zb, mub, chi2, mu
cH=299792.458/70.; Ez=lambda z: np.sqrt(0.3*(1+z)**3+0.7)
DA=np.vectorize(lambda z: cH*integrate.quad(lambda x:1/Ez(x),0,z)[0]/(1+z))
laws={"d = cz/H0":lambda z: cH*z,"d = (c/H0) ln(1+z)":lambda z: cH*np.log(1+z)}
# cluster ratios of Entry 28 (two pairs of telescopes), recomputed from the published data files
X2=here if os.path.exists(os.path.join(here,"distance_test_repeat.json")) else os.path.join(here,"..","xsz2")     # written by distance_test_repeat.py (Entry 28)
J=json.load(open(os.path.join(X2,"distance_test_repeat.json")))
zA=np.array([p[0] for p in J["pointsA"]]); rA=np.array([p[1] for p in J["pointsA"]]); zB=np.array([p[0] for p in J["pointsB"]]); rB=np.array([p[1] for p in J["pointsB"]])
z=np.r_[zA,zB]; r=np.r_[rA,rB]; s1=np.r_[np.ones(len(zA)),np.zeros(len(zB))]; s2=1-s1
rng=np.random.default_rng(3); out={}
print(f"clusters: {len(z)} (South Pole Telescope with Chandra, Planck with XMM); supernovae: 31 binned distances (740 supernovae)\n")
for name,d in laws.items():
    # supernovae: d_L = d (1+z)^q, q = (1+j)/2
    f=lambda q: sh.chi2(sh.mu(d(sh.zb)*(1+sh.zb)**q))[0]; g=optimize.minimize_scalar(f,bounds=(0,3),method="bounded"); q=g.x
    qlo=optimize.brentq(lambda x:f(x)-g.fun-1,0,q); qhi=optimize.brentq(lambda x:f(x)-g.fun-1,q,3); j_sn=2*q-1; e_sn=(qhi-qlo)
    # clusters: log R = level + 0.5 log(d/D_A) + (j-3)/2 log(1+z)
    y=r-0.5*np.log10(d(z)/DA(z)); x=np.log10(1+z); A=np.c_[s1,s2,x]; b=np.linalg.lstsq(A,y,rcond=None)[0]; res=y-A@b; cov=res@res/(len(y)-3)*np.linalg.inv(A.T@A)
    bs=[]
    for _ in range(4000):
        k=rng.integers(0,len(z),len(z)); bs.append(np.linalg.lstsq(A[k],y[k],rcond=None)[0][2])
    e_slope=max(np.sqrt(cov[2,2]),np.std(bs)); j_cl=3+2*b[2]; e_cl=2*e_slope
    sep=(j_cl-j_sn)/np.hypot(e_cl,e_sn)
    print(f"{name:22s} supernovae: j = {j_sn:+.2f} +- {e_sn:.2f}  (chi2 {g.fun:.1f});  clusters: j = {j_cl:+.2f} +- {e_cl:.2f};  apart by {sep:.1f} sigma")
    for jj,lab in ((0.0,"energy loss only"),(1.0,"stretched time")):
        c_sn=f((1+jj)/2); pull=(j_cl-jj)/e_cl; print(f"      {lab:18s} (j = {jj:.0f}): supernova chi2 {c_sn:6.1f} (+{c_sn-33.6:.1f} over standard), clusters {abs(pull):.1f} sigma away")
    # one j for both: weighted mean and its chi2
    w1,w2=1/e_sn**2,1/e_cl**2; jbar=(w1*j_sn+w2*j_cl)/(w1+w2); chi=w1*(j_sn-jbar)**2+w2*(j_cl-jbar)**2
    print(f"      best single j = {jbar:.2f}; misfit chi2 = {chi:.1f} for one degree of freedom")
    out[name]=dict(j_sn=float(j_sn),e_sn=float(e_sn),j_cl=float(j_cl),e_cl=float(e_cl),sep=float(sep),jbar=float(jbar),chi=float(chi))
# the same bookkeeping for the standard model: photon counts per area fall as (1+z)^-3, distances standard: the cluster ratio should be level
x=np.log10(1+z); A=np.c_[s1,s2,x]; b=np.linalg.lstsq(A,r,rcond=None)[0]; res=r-A@b; e=np.sqrt(res@res/(len(r)-3)*np.linalg.inv(A.T@A)[2,2]); print(f"\nstandard model, same clusters: photon surface brightness falls as (1+z)^-n with n = {3+2*b[2]:.2f} +- {2*e:.2f} (expected 3); in energy, (1+z)^-{4+2*b[2]:.2f}")
out["tolman"]=dict(n_photon=float(3+2*b[2]),n_energy=float(4+2*b[2]),err=float(2*e))
json.dump(out,open(os.path.join(here,"dimming_dilemma.json"),"w"),indent=1)
