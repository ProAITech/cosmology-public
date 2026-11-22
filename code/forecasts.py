"""Workbench Entry 38: what to measure next. Forecasts from the notebook's own fits.
Each block takes a result of Entries 31-37, adds hypothetical data, and asks how much is needed to decide the question it left open."""
import sys, os, json
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
import numpy as np
from scipy import optimize
import cosmo as c
from dark_energy_fit import fit, sig, make
out={}
# ---------- 1. Gas fractions (Entry 34): how many more distant relaxed clusters?
R=[l.split() for l in open(os.path.join(HERE,"mantz2022_fgas.txt")) if not l.startswith("#")]; z=np.array([float(r[1]) for r in R]); f=np.array([float(r[6]) for r in R]); e=np.array([float(r[7]) for r in R])
E=lambda t: np.sqrt(0.3*(1+t)**3+0.7); tm=lambda t: -0.5*np.log(1+t)+0.390*np.log((1+t)/E(t)); sint=0.0613
def sig_k(zs,ss):
    X=np.c_[np.ones(len(zs)),tm(zs)]; w=1/(ss**2+sint**2); return float(np.sqrt(np.linalg.inv((X*w[:,None]).T@X)[1,1]))
s0=sig_k(z,e/f); print(f"1. Gas fractions: 44 clusters give an error of {s0:.3f} on k (Entry 34: 0.21)")
rng=np.random.default_rng(38); rows=[]
for N in (0,6,12,20,30,50,80):
    zn=np.linspace(0.75,1.2,N) if N else np.array([]); s=sig_k(np.r_[z,zn],np.r_[e/f,np.full(N,0.12)]); rows.append([N,s,0.5/s,1.0/s])
    print(f"   + {N:2d} clusters at z = 0.75-1.2 measured to 12%: error {s:.3f}; 0 and 0.5 apart by {0.5/s:.1f} sigma; 0 and 1 by {1/s:.1f}")
out["fgas"]=rows
need=lambda target: next((N for N in range(0,400) if sig_k(np.r_[z,np.linspace(0.75,1.2,N)] if N else z,np.r_[e/f,np.full(N,0.12)])<=target),None)
out["fgas_need"]=dict(three=need(0.5/3),five=need(0.5/5)); print(f"   to tell 0 from 0.5 at 3 sigma: {out['fgas_need']['three']} more; at 5 sigma: {out['fgas_need']['five']} more")
# ---------- 2. Shadows with measured profiles (Entries 28, 32, 33)
P28=json.load(open(os.path.join(HERE,"entry28_points.json"))); zA=np.array([p[0] for p in P28["pointsA"]]); yA=np.array([p[1] for p in P28["pointsA"]]); zB=np.array([p[0] for p in P28["pointsB"]]); yB=np.array([p[1] for p in P28["pointsB"]])
def sig_p(extra_z,scat=0.10):
    zz=np.r_[zA,zB,extra_z]; D=np.zeros((len(zz),3)); D[:len(zA),0]=1; D[len(zA):len(zA)+len(zB),1]=1; D[len(zA)+len(zB):,2]=1
    if len(extra_z)==0: D=D[:,:2]
    X=np.c_[D,np.log10(1+zz)]; s=np.r_[np.full(len(zA),yA.std(ddof=1)),np.full(len(zB),yB.std(ddof=1)),np.full(len(extra_z),scat)]; w=1/s**2
    return float(np.sqrt(np.linalg.inv((X*w[:,None]).T@X)[-1,-1]))
print(f"\n2. Shadows: pairs A and B give an error of {sig_p(np.array([])):.3f} on the power of (1+z) (Entry 33: 0.18 with resampling)")
rows=[]
for N in (0,57,100,150,200,300):
    zn=np.linspace(0.3,1.4,N) if N else np.array([]); s=sig_p(zn); rows.append([N,s,0.5/s]); print(f"   + {N:3d} clusters with measured profiles, z = 0.3-1.4, scatter 0.10 dex: error {s:.3f}; the last static model ({-0.5}) is {0.5/s:.1f} sigma from level")
out["shadows"]=rows
# ---------- 3. The Hubble constant (Entry 37)
print("\n3. Hubble constant: precision a direct measurement needs to tell the shorter sound horizon (69.28 +- 0.53) from changing dark energy")
rows=[]
for name,h,e_ in (("+ Pantheon+",67.51,0.59),("+ DES 5-year",66.74,0.56),("+ Union3",65.91,0.84)):
    d=69.28-h; r=[]
    for ns in (2,3):
        v=(d/ns)**2-0.53**2-e_**2; r.append(float(np.sqrt(v)) if v>0 else None)
    rows.append([name,d,r[0],r[1]]); print(f"   changing dark energy {name:14s}: apart by {d:.2f}; needs +-{r[0] if r[0] else float('nan'):.2f} for 2 sigma, "+(f"+-{r[1]:.2f}" if r[1] else "not reachable")+" for 3 sigma")
out["h0"]=rows
# ---------- 4. Clumping (Entry 36, Part 4)
C=json.load(open(os.path.join(HERE,"..","de","clumping_fit.json"))) if os.path.exists(os.path.join(HERE,"..","de","clumping_fit.json")) else None
if C:
    b,lo,hi,_=C["b"]; print(f"\n4. Clumping: ruler + CMB ask for b = {b:.2f} +{hi:.2f} -{lo:.2f}. A direct measurement with error s would detect it at {b:.2f}/s sigma: s = 0.08 gives {b/0.08:.1f} sigma; a null with s = 0.08 leaves b < 0.16 at 2 sigma, {(b-0.0)/np.hypot(lo,0.08):.1f} sigma below what the ruler asks")
    out["clump"]=dict(b=b,lo=lo,hi=hi)
# ---------- 5. Supernova calibration between nearby and distant (Entries 31, 36)
print("\n5. Supernovae: the dark-energy preference (ruler + CMB + 1,048 supernovae) as the nearby ones (z < 0.1) are shifted in brightness")
from cosmo import ZP, MP, WP
def chi_pan_shift(s):
    mag=MP+np.where(ZP<0.1,s,0.0)
    def f(m):
        d=mag-5*np.log10((1+ZP)*m.DM(ZP)); a=np.sum(WP*d)/WP.sum(); return float(np.sum(WP*(d-a)**2))
    return f
bao=lambda m: c.chi_bao(m); cmb=lambda m: c.chi_cmb(m); rows=[]
for s in (-0.02,0.0,0.02,0.04,0.06,0.08):
    chi=make([bao,cmb,chi_pan_shift(s)]); L=fit(chi,'lcdm'); W=fit(chi,'w0wa'); g=L['chi']-W['chi']; rows.append([s,g,sig(g),W['w0'],W['wa']])
    print(f"   nearby supernovae made {s:+.2f} mag fainter: gain {g:5.2f} = {sig(g):.1f} sigma, w0 {W['w0']:+.2f}, wa {W['wa']:+.2f}")
out["sn_shift"]=rows
json.dump(out,open(os.path.join(HERE,"forecasts.json"),"w"),indent=1)
