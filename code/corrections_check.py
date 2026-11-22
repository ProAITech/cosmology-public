"""Entry 41: could the corrections applied to the supernova magnitudes create the same-camera step that sinks the static law?
Subsamples of the Pantheon+ table with the light-curve parameters: Foundation at z = 0.03-0.04 and 0.05-0.10 (77), Pan-STARRS Medium Deep
at 0.20-0.25, 0.35-0.45 and 0.45-0.65 (121). Columns: zHD, corrected magnitude, raw peak magnitude mB, stretch x1, colour c.
The published corrected magnitude is mB + alpha x1 - beta c minus corrections for selection, scatter model and host. A 'plain' magnitude
keeps only the first three terms, with one alpha and beta for everyone. The difference is everything a model put in."""
import os, json, sys
HERE=os.path.dirname(os.path.abspath(__file__))
import numpy as np
from scipy import integrate
ld=lambda f: np.loadtxt(os.path.join(HERE,"corr",f))
F=np.vstack([ld("found_c.txt"),ld("found_a.txt"),ld("found_d.txt"),ld("found_b.txt")]); P=np.vstack([ld("ps1_a.txt"),ld("ps1_b.txt"),ld("ps1_e.txt"),ld("ps1_c.txt"),ld("ps1_d.txt")]); A=np.loadtxt(os.path.join(HERE,"pantheon_plus_rows.txt"))
key={(round(r[0],5),round(r[1],3)):r[2] for r in A}
miss=[tuple(r[:2]) for r in np.vstack([F,P]) if (round(r[0],5),round(r[1],3)) not in key]; print("rows not found in the main table:",miss)
eF=np.array([key[(round(r[0],5),round(r[1],3))] for r in F]); eP=np.array([key[(round(r[0],5),round(r[1],3))] for r in P]); out={"n":[len(F),len(P)]}
zg=np.linspace(0,2.4,2401)
def mu_l(z,om=0.30): dc=integrate.cumulative_trapezoid(1/np.sqrt(om*(1+zg)**3+1-om),zg,initial=0); return 5*np.log10((1+z)*np.interp(z,zg,dc))
CUR={"standard":mu_l,"Lerner":lambda z:5*np.log10(z*np.sqrt(1+z)),"ln(1+z) stretched":lambda z:5*np.log10(np.log(1+z)*(1+z)),"coasting":lambda z:5*np.log10(z*(1+z/2))}
def seam(mF,mP,wF,wP,curve):
    rF=mF-curve(F[:,0]); rP=mP-curve(P[:,0]); return float(np.sum(wF*rF)/wF.sum()-np.sum(wP*rP)/wP.sum())
# 1. the corrections themselves
print("1. Corrected minus plain magnitude (alpha 0.148)")
out["corr"]={}
for be in (3.09,2.8,2.5):
    dF=F[:,1]-(F[:,2]+0.148*F[:,3]-be*F[:,4]); dP=P[:,1]-(P[:,2]+0.148*P[:,3]-be*P[:,4])
    print(f"   beta {be}: Foundation {dF.mean():+.4f} +- {dF.std(ddof=1)/len(dF)**0.5:.4f}; Pan-STARRS {dP.mean():+.4f} +- {dP.std(ddof=1)/len(dP)**0.5:.4f}; Pan-STARRS minus Foundation {dP.mean()-dF.mean():+.4f} +- {np.hypot(dF.std(ddof=1)/len(dF)**0.5,dP.std(ddof=1)/len(dP)**0.5):.4f}")
    out["corr"][str(be)]=[float(dF.mean()),float(dP.mean()),float(dP.mean()-dF.mean()),float(np.hypot(dF.std(ddof=1)/len(dF)**0.5,dP.std(ddof=1)/len(dP)**0.5))]
    for lab,sl in (("0.10-0.25",P[:,0]<0.25),("0.25-0.35",(P[:,0]>0.25)&(P[:,0]<0.35)),("0.35-0.65",P[:,0]>0.35)): print(f"        Pan-STARRS {lab}: {dP[sl].mean():+.4f} +- {dP[sl].std(ddof=1)/sl.sum()**0.5:.4f} ({int(sl.sum())}); mean colour {P[sl,4].mean():+.3f}, stretch {P[sl,3].mean():+.3f}")
print(f"   mean colour: Foundation {F[:,4].mean():+.3f}, Pan-STARRS {P[:,4].mean():+.3f}; mean stretch: {F[:,3].mean():+.3f}, {P[:,3].mean():+.3f}")
out["colour"]=[float(F[:,4].mean()),float(P[:,4].mean())]; out["stretch"]=[float(F[:,3].mean()),float(P[:,3].mean())]
# 2. the seam of these subsamples, corrected and plain
print("2. Foundation minus Pan-STARRS, these subsamples (negative = nearby too bright for the curve)")
wF=1/eF**2; wP=1/eP**2; err=float((1/wF.sum()+1/wP.sum())**0.5); out["seam"]={}; out["err"]=err
for nm,cv in CUR.items():
    s0=seam(F[:,1],P[:,1],wF,wP,cv); row=[s0]
    txt=f"   {nm:18s} corrected {s0:+.3f}"
    for be in (3.09,2.8,2.5):
        s=seam(F[:,2]+0.148*F[:,3]-be*F[:,4],P[:,2]+0.148*P[:,3]-be*P[:,4],wF,wP,cv); row.append(s); txt+=f" | plain, beta {be}: {s:+.3f}"
    # no colour or stretch correction at all
    s=seam(F[:,2],P[:,2],wF,wP,cv); row.append(s); txt+=f" | raw peak magnitude: {s:+.3f}"
    print(txt+f"   (+- {err:.3f})"); out["seam"][nm]=row
# 3. plain magnitudes with alpha, beta fitted to these 198 supernovae themselves (one pair for all), each curve
print("3. alpha and beta fitted on these supernovae, with a free step between the two surveys")
out["fit"]={}
for nm,cv in CUR.items():
    z=np.r_[F[:,0],P[:,0]]; y=np.r_[F[:,2],P[:,2]]-cv(z); x1=np.r_[F[:,3],P[:,3]]; c=np.r_[F[:,4],P[:,4]]; isF=np.r_[np.ones(len(F)),np.zeros(len(P))]
    X=np.c_[np.ones(len(z)),isF,-x1,c]; co,res,rk,sv=np.linalg.lstsq(X,y,rcond=None); r=y-X@co; s2=r.var(ddof=4); cov=s2*np.linalg.inv(X.T@X)
    print(f"   {nm:18s} step {co[1]:+.3f} +- {cov[1,1]**0.5:.3f}; alpha {co[2]:.3f}, beta {co[3]:.2f}; scatter {s2**0.5:.3f}"); out["fit"][nm]=[float(co[1]),float(cov[1,1]**0.5),float(co[2]),float(co[3]),float(s2**0.5)]
json.dump(out,open(os.path.join(HERE,"corrections_check.json"),"w"),indent=1)
