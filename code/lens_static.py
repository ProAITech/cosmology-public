"""Workbench Entry 112: lens time delays against static distance laws (open question from Entry 34).

What a lens measures. Imaging gives the Fermat-potential difference between images, in angle squared; the delay gives c*dt. Their ratio is the
time-delay distance D_dt = c dt / d(phi), in Mpc, whatever the cosmology. Stellar kinematics give a second, dimensionless number,
R = D_s/D_ds, because the velocity dispersion goes as theta_E * R in any geometry where light bends as in general relativity.

Predictions.
  Flat expanding space:              D_dt = (1+z_l) D_l D_s / D_ls = D_M,l D_M,s / (D_M,s - D_M,l);  R = D_M,s / (D_M,s - D_M,l)
  Static Euclidean space, distance d: the delay is made at the lens and travels to us. By the crest rule of Entry 111 it arrives stretched by
    S = 1 + z_l if the model stretches supernova light curves (a growing travel time), and S = 1 if photons only lose energy.
    D_dt = S d_l d_s / (d_s - d_l);  R = d_s / (d_s - d_l)
  So a static model with no stretch and d proportional to the comoving distance D_M predicts exactly what flat expanding space predicts, for
  both numbers. A static model that stretches time predicts D_dt larger by 1 + z_l.
Every law is cz/H0 at low redshift, so the local Hubble constant fixes each law's scale; the lenses then give each law's own H0, which must agree.

Data: per-lens D_dt from H0LiCOW XIII (Wong et al. 2020, Table 2, power-law or composite mass models) and DES J0408-5354 (Shajib et al. 2020).
TDCOSMO 2026 (eleven lenses, mass-sheet transform fixed by resolved kinematics, lambda_int = 0.936 +- 0.050) gets 71.1 (+2.4, -2.3) for flat
LCDM; each law is rescaled to that sample by its own ratio to flat space, averaged over the eleven lenses' redshifts. Three of the new lenses
have per-lens H0 values in the paper and join the spread test.
Local H0: SH0ES 73.0 +- 1.0; CCHP 70.4 +- 1.9."""
import json, os
import numpy as np
from scipy import integrate, optimize
HERE=os.path.dirname(os.path.abspath(__file__)); C=299792.458
# name, z_l, z_s, D_dt median, +err, -err (Mpc)
L=[("B1608+656",0.6304,1.394,5156,296,236),("RXJ1131−1231",0.295,0.654,2096,98,83),("HE0435−1223",0.4546,1.693,2707,183,168),
   ("SDSS1206+4332",0.745,1.789,5769,589,471),("WFI2033−4723",0.6575,1.662,4784,399,248),("PG1115+080",0.311,1.722,1470,137,127),
   ("DESJ0408−5354",0.597,2.375,3382,146,115)]
names=[a[0] for a in L]; zl=np.array([a[1] for a in L]); zs=np.array([a[2] for a in L]); D=np.array([a[3] for a in L],float)
sl=np.array([(a[4]+a[5])/2/a[3] for a in L])       # fractional (log) error
zg=np.linspace(0,2.5,25001)
def dm_shape(om=0.297):                              # comoving distance shape that fits DESI DR2's ruler (Entry 33), units c/H0
    E=np.sqrt(om*(1+zg)**3+1-om); return integrate.cumulative_trapezoid(1/E,zg,initial=0)
DMS=dm_shape()
LAWS={"ruler":lambda z: np.interp(z,zg,DMS),"line":lambda z: z,"log":lambda z: np.log1p(z)}
LAWNAME={"ruler":"d follows the clustering ruler (the measured static law, Entry 33)","line":"d = cz/H0 (Lerner)","log":"d = (c/H0) ln(1+z) (tired light)"}
def Ddt_unit(law,stretch):                           # D_dt in units of c/H0
    d=LAWS[law]; dl,ds=d(zl),d(zs); S=(1+zl) if stretch else 1.0
    return S*dl*ds/(ds-dl)
def R_of(law):
    d=LAWS[law]; dl,ds=d(zl),d(zs); return ds/(ds-dl)
out={"lenses":[dict(name=n,zl=float(a),zs=float(b),D=float(c),s=float(e)) for n,a,b,c,e in zip(names,zl,zs,D,sl)]}
print("Check: flat expanding space with the ruler's matter share 0.297, per-lens H0 and combined")
def fitH(u):
    h=C*u/D; w=1/sl**2; Hm=np.exp(np.sum(w*np.log(h))/w.sum()); e=Hm/np.sqrt(w.sum()); chi=float(np.sum(w*np.log(h/Hm)**2)); return h,Hm,e,chi
h,Hm,e,chi=fitH(Ddt_unit("ruler",False)); out["flat"]=dict(h=h.tolist(),H=Hm,e=e,chi=chi)
print("   "+", ".join(f"{n} {v:.1f}" for n,v in zip(names,h))+f"\n   combined {Hm:.1f} +- {e:.1f}, chi2 {chi:.1f} for 6 d.o.f.  (H0LiCOW+DES0408 published: 73.3 and 74.2)")
# trend of per-lens H0 with lens redshift, as published concerns noted
w=1/sl**2; X=np.c_[np.ones(7),np.log(1+zl)]; cov=np.linalg.inv(X.T@(X*w[:,None])); b=cov@X.T@(w*np.log(h))
out["flat"]["trend"]=[float(b[1]),float(np.sqrt(cov[1,1]))]; print(f"   per-lens H0 against lens redshift: H0 ∝ (1+z_l)^({b[1]:+.2f} ± {np.sqrt(cov[1,1]):.2f})")
print("\nStatic laws. Each law's H0 from the lenses must equal the local value (73.0 ± 1.0 SH0ES; 70.4 ± 1.9 CCHP).")
out["models"]={}
for law in ("ruler","line","log"):
    for st in (False,True):
        h,Hm,e,chi=fitH(Ddt_unit(law,st)); Rr=R_of(law)
        key=f"{law}_{'stretch' if st else 'loss'}"
        # rescale to the TDCOSMO 2025 mass-sheet-free scale: same ratio to flat space as above
        Htd=71.6*Hm/out["flat"]["H"]; etd=Htd*np.hypot(3.6/71.6,e/Hm)
        n73=(Hm-73.0)/np.hypot(e,1.0); n70=(Hm-70.4)/np.hypot(e,1.9); m73=(Htd-73.0)/np.hypot(etd,1.0)
        out["models"][key]=dict(h=h.tolist(),H=Hm,e=e,chi=chi,Htd=Htd,etd=etd,sig73=n73,sig70=n70,sig73_td=m73,R=Rr.tolist())
        print(f"   {LAWNAME[law]:62s} {'stretched' if st else 'energy loss only':16s} H0 = {Hm:6.1f} ± {e:3.1f}; spread chi2 {chi:5.1f}/6;"
              f" vs SH0ES {n73:+5.1f}σ, vs CCHP {n70:+5.1f}σ; on the TDCOSMO 2025 scale {Htd:5.1f} ± {etd:3.1f} ({m73:+.1f}σ vs SH0ES)")
# kinematic ratio R = d_s/d_ls: H0-free; compare the laws to flat space
print("\nKinematic ratio R = d_s/(d_s - d_l), which the velocity dispersion measures (no H0 in it):")
Rf=R_of("ruler"); out["R"]={"flat":Rf.tolist()}
for law in ("line","log"):
    r=R_of(law)/Rf; out["R"][law]=R_of(law).tolist(); print(f"   {LAWNAME[law]:62s} R relative to flat space: "+", ".join(f"{v:.2f}" for v in r))
# the kinematic ratio measured: R = D_dt / ((1+z_l) D_d), from the published D_d (five lenses)
DD={"B1608+656":(1228,177,151),"RXJ1131−1231":(804,141,112),"SDSS1206+4332":(1805,555,398),"PG1115+080":(697,186,144),"DESJ0408−5354":(1711,376,280)}
idx=[names.index(k) for k in DD]; dd=np.array([DD[k][0] for k in DD],float); sd=np.array([(DD[k][1]+DD[k][2])/2/DD[k][0] for k in DD])
Robs=D[idx]/((1+zl[idx])*dd); sR=np.hypot(sd,sl[idx])        # covariance between D_dt and D_d ignored: errors overstated
out["Rtest"]=dict(names=list(DD),obs=Robs.tolist(),err=sR.tolist())
print("\nMeasured kinematic ratio against each law (log errors; D_d dominates):")
for law in ("ruler","line","log"):
    r=R_of(law)[idx]; chi=float(np.sum((np.log(Robs/r)/sR)**2)); mean=float(np.exp(np.sum(np.log(Robs/r)/sR**2)/np.sum(1/sR**2))); em=float(mean/np.sqrt(np.sum(1/sR**2)))
    out["Rtest"][law]=dict(pred=r.tolist(),chi=chi,ratio=mean,err=em); print(f"   {LAWNAME[law]:62s} measured/predicted {mean:.2f} ± {em:.2f}; chi2 {chi:.1f} for 5")
# what mass-sheet factor would rescue each stretched law at the local H0
print("\nFactor by which every lens's D_dt would have to be wrong (a common mass sheet) to bring each law to 73.0:")
out["lam"]={}
for key,v in out["models"].items():
    lam=v["H"]/73.0; out["lam"][key]=lam; print(f"   {key:16s} lambda = {lam:.2f}")
# ---- TDCOSMO 2026: eleven lenses, scale with the mass sheet free
T26=[("DESJ0408",0.597,2.375),("HE0435",0.455,1.693),("HE1104",0.729,2.319),("PG1115",0.310,1.727),("RXJ1131",0.295,0.657),("SDSS1206",0.745,1.789),
     ("SDSS1433",0.407,2.737),("J1537",0.590,1.721),("B1608",0.630,1.394),("WFI2033",0.657,1.662),("WGD2038",0.228,0.777)]
zl2=np.array([a[1] for a in T26]); zs2=np.array([a[2] for a in T26])
def unit2(law,stretch):
    d=LAWS[law]; dl,ds=d(zl2),d(zs2); return ((1+zl2) if stretch else 1.0)*dl*ds/(ds-dl)
uf=unit2("ruler",False); out["t26"]={}
print("\nOn the TDCOSMO 2026 sample (71.1 +2.4 -2.3 for flat space, mass sheet free):")
for law in ("ruler","line","log"):
    for st in (False,True):
        r=float(np.exp(np.mean(np.log(unit2(law,st)/uf)))); H=71.1*r; e=2.35*r; key=f"{law}_{'stretch' if st else 'loss'}"
        s73=(H-73.0)/np.hypot(e,1.0); s70=(H-70.4)/np.hypot(e,1.9); out["t26"][key]=dict(ratio=r,H=H,e=e,sig73=s73,sig70=s70)
        print(f"   {key:16s} ratio {r:.3f}: H0 = {H:6.1f} ± {e:3.1f}; vs SH0ES {s73:+5.1f}σ; vs CCHP {s70:+5.1f}σ")
# kinematics set the mass sheet through sigma^2 ~ lambda * R * theta_E. A law with a different R needs lambda changed by R_flat/R_law,
# and the inferred D_dt moves by R_law/R_flat, so the law's H0 moves by R_flat/R_law. Estimate on the eleven-lens redshifts.
def R2(law):
    d=LAWS[law]; dl,ds=d(zl2),d(zs2); return ds/(ds-dl)
out["t26_kin"]={}
for law in ("line","log"):
    f=float(np.exp(np.mean(np.log(R2("ruler")/R2(law))))); out["t26_kin"][law]=f
    for st in (False,True):
        key=f"{law}_{'stretch' if st else 'loss'}"; H=out["t26"][key]["H"]*f; e=out["t26"][key]["e"]*f
        out["t26_kin"][key]=dict(H=H,e=e,sig73=(H-73.0)/np.hypot(e,1.0),sig70=(H-70.4)/np.hypot(e,1.9))
        print(f"   with its own kinematic ratio (mean shift x{f:.3f}): {key:16s} H0 = {H:6.1f} ± {e:3.1f}; vs SH0ES {(H-73)/np.hypot(e,1):+5.1f}σ; vs CCHP {(H-70.4)/np.hypot(e,1.9):+5.1f}σ")
# ---- per-lens spread with ten lenses (three new per-lens values from TDCOSMO 2026, flat-space H0 with errors)
NEW=[("HE1104−1805",0.729,2.319,64.2,5.4),("SDSS1433+6007",0.407,2.737,73.2,4.75),("J1537−3010",0.590,1.721,75.5,7.55)]
zl3=np.r_[zl,[a[1] for a in NEW]]; zs3=np.r_[zs,[a[2] for a in NEW]]
hf=np.r_[np.array(out["flat"]["h"]),[a[3] for a in NEW]]; s3=np.r_[sl,[a[4]/a[3] for a in NEW]]; w3=1/s3**2
def unit3(law,stretch):
    d=LAWS[law]; dl,ds=d(zl3),d(zs3); return ((1+zl3) if stretch else 1.0)*dl*ds/(ds-dl)
uf3=unit3("ruler",False); out["spread10"]={}
X=np.c_[np.ones(10),np.log(1+zl3)]; cv=np.linalg.inv(X.T@(X*w3[:,None])); bb=cv@X.T@(w3*np.log(hf)); out["spread10"]["trend"]=[float(bb[1]),float(np.sqrt(cv[1,1]))]
print(f"\nTen lenses, flat-space per-lens H0 against lens redshift: (1+z_l)^({bb[1]:+.2f} ± {np.sqrt(cv[1,1]):.2f}); static with stretch and the ruler's law predicts −1 here")
for law in ("ruler","line","log"):
    for st in (False,True):
        h=hf*unit3(law,st)/uf3; m=np.exp(np.sum(w3*np.log(h))/w3.sum()); chi=float(np.sum(w3*np.log(h/m)**2)); key=f"{law}_{'stretch' if st else 'loss'}"
        out["spread10"][key]=dict(h=h.tolist(),chi=chi); print(f"   {key:16s} spread chi2 {chi:5.1f} for 9")
out["spread10"]["zl"]=zl3.tolist(); out["spread10"]["err"]=s3.tolist(); out["spread10"]["names"]=names+[a[0] for a in NEW]
json.dump(out,open(os.path.join(HERE,"lens_static.json"),"w"),indent=1,ensure_ascii=False)
