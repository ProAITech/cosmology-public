"""Workbench Entry 111: the four DemystiCon 2025 redshift posters (Coste, Booth, Kelly, Graneau), sorted into the table of Entry 21 and,
where a formula exists, tested.

Only Booth gives formulas (Preprints.org 202412.0837 v2, not peer reviewed):
  Eq. 32  d_A = (c/H0) * integral_0^z [1 + log(1+z)] / E(z') dz'   with E(z) = 1 + z for his non-expanding universe
  Eq. 35  T(z) = (1/H0) * integral_z^inf [1 + log(1+z')] / [(1+z') E(z')] dz'
  Eq. 2   d_L = (1+z) d_A ;  H0 = 70.74 km/s/Mpc from his contraction rate.
Two things are not fixed by the text we could read: whether the bracket is evaluated at z (outside the integral) or z' (inside), and the base
of the logarithm. All four readings are tested. His own timing statements (excess age over the standard model of more than 6 Gyr at its
largest, more than 2 Gyr extra at z = 25) are used to pick the reading he computed with.

Key point used throughout: d_A in Eq. 32 has the form of a comoving distance, integral of c dz / H_eff with H_eff = H0 (1+z) / [1 + log(1+z)].
So for supernovae (d_L = (1+z) d_A) and for the clustering ruler (transverse distance d_A, radial distance its slope) the model makes the same
predictions as an expanding universe with that H_eff, whichever way one interprets it. That is the Entry 21 rule at work.

Data: DESI DR2 D_M/D_H at six redshifts (as Entry 21) and all 13 ruler numbers (as Entry 33); Pantheon+ (diagonal errors, repeats averaged,
as Entry 39); JLA binned distances with full covariance (as Entry 31)."""
import sys, os, io, json, contextlib
HERE=os.path.dirname(os.path.abspath(__file__))
CODE=HERE if os.path.exists(os.path.join(HERE,"pantheon_plus_rows.txt")) else os.path.abspath(os.path.join(HERE,"..","artifact-files","4379c2ac-1195-4a41-978e-eb44b945a265","code"))
sys.path.insert(0,CODE)
import numpy as np
from scipy import integrate, optimize
_cwd=os.getcwd(); os.chdir(CODE)
with contextlib.redirect_stdout(io.StringIO()): import static_hubble as sh
os.chdir(_cwd)
out={}
LN10=np.log(10.0)
# ---------------- Booth's distance law, four readings. d in units of c/H0; s = dd/dz.
def booth(reading):
    k=1.0 if reading.endswith("ln") else 1/LN10          # log base: natural (k = 1) or base 10 (k = 1/ln 10)
    inside=reading.startswith("a")                        # bracket inside the integral (a) or evaluated at z (b)
    if inside:
        d=lambda z: np.log1p(z)+k*np.log1p(z)**2/2; s=lambda z: (1+k*np.log1p(z))/(1+z)
    else:
        d=lambda z: np.log1p(z)*(1+k*np.log1p(z)); s=lambda z: (1+2*k*np.log1p(z))/(1+z)
    return d,s
READ={"a_ln":"bracket inside the integral, natural log","b_ln":"bracket at z, natural log","a_10":"bracket inside the integral, base-10 log","b_10":"bracket at z, base-10 log"}
# ---------------- 1. which reading did Booth compute with? His ages against the standard model
H0B=70.74; GYR=977.792/H0B                               # 1/H0 in Gyr for H0 = 70.74
def T_booth(z,reading):
    k=1.0 if reading.endswith("ln") else 1/LN10
    if reading.startswith("a"):                           # integral_z^inf (1 + k ln u)/u^2 du with u = 1+z' : [1 + k (1 + ln u)]/u
        return GYR*(1+k*(1+np.log1p(z)))/(1+z)
    return GYR*(1+k*np.log1p(z))/(1+z)                    # bracket at z: (1 + k ln(1+z)) * integral_z^inf du/u^2
def T_lcdm(z,H0=67.36,om=0.3153):
    orad=4.18e-5/(H0/100)**2; ol=1-om-orad
    f=lambda u: 1/(u*np.sqrt(om*u**3+orad*u**4+ol))      # u = 1+z
    return 977.792/H0*integrate.quad(f,1+z,np.inf,limit=200)[0]
zz=np.r_[np.linspace(0,5,501)[1:],np.linspace(5.5,30,50)]; TL=np.array([T_lcdm(z) for z in zz]); out["ages"]={}
print(f"Booth's own timing statements (H0 = {H0B}; standard model: Planck 2018, age today {T_lcdm(0):.2f} Gyr)")
for r in READ:
    ex=T_booth(zz,r)-TL; i=int(np.argmax(ex)); e25=float(T_booth(25,r)-T_lcdm(25))
    out["ages"][r]=dict(age0=float(T_booth(0,r)),max_excess=float(ex[i]),z_max=float(zz[i]),excess_z25=e25,T_z2=float(T_booth(2,r)))
    print(f"   {READ[r]:42s} age today {T_booth(0,r):5.1f} Gyr; largest excess {ex[i]:5.1f} Gyr at z = {zz[i]:.2f}; extra at z = 25: {e25:4.2f} Gyr")
print("   paper: largest excess > 6 Gyr; extra at z = 25 > 2 Gyr; standard model at z = 25: %.0f Myr" % (1e3*T_lcdm(25)))
# ---------------- 2. the shape test (scale-free): D_M/D_H = d / s
SHAPE=[(0.510,0.622,0.017),(0.706,0.892,0.021),(0.934,1.223,0.019),(1.321,1.948,0.045),(1.484,2.386,0.136),(2.330,4.518,0.097)]
zs=np.array([a for a,_,_ in SHAPE]); ys=np.array([b for _,b,_ in SHAPE]); es=np.array([c for _,_,c in SHAPE])
def F_lcdm(z,om=0.315):
    E=lambda x: np.sqrt(om*(1+x)**3+1-om); return integrate.quad(lambda x:1/E(x),0,z)[0]*E(z)
out["shape"]={"z":zs.tolist(),"obs":ys.tolist(),"err":es.tolist(),"lcdm":[F_lcdm(z) for z in zs]}
c_l=float(np.sum(((ys-np.array(out["shape"]["lcdm"]))/es)**2)); out["shape"]["chi_lcdm"]=c_l
print(f"\nShape test, D_M/D_H at six DESI redshifts (no ruler length, no galaxy evolution). Standard model (matter 0.315): chi2 = {c_l:.1f}")
for r in READ:
    d,s=booth(r); F=d(zs)/s(zs); sig=(ys-F)/es; c=float(np.sum(sig**2)); out["shape"][r]=dict(pred=F.tolist(),chi=c,sig=sig.tolist())
    print(f"   {READ[r]:42s} " + " ".join(f"{v:5.3f}" for v in F) + f"   chi2 = {c:7.0f}; worst point {np.max(np.abs(sig)):.0f} sigma")
print("   DESI                                       " + " ".join(f"{v:5.3f}" for v in ys))
zf=np.linspace(0.3,2.6,47); out["shape"]["curve"]={"z":zf.tolist(),"lcdm":[F_lcdm(z) for z in zf]}
for r in READ:
    d,s=booth(r); out["shape"]["curve"][r]=(d(zf)/s(zf)).tolist()
out["shape"]["curve"]["line"]=zf.tolist()                           # d = cz/H0 in static space: D_M/D_H = z
out["shape"]["curve"]["log"]=((1+zf)*np.log1p(zf)).tolist()          # d = (c/H0) ln(1+z)
# ---------------- 3. all 13 ruler numbers, one free scale (Entry 33's method)
BAO=[(0.510,13.588,0.167,21.863,0.425,-0.459),(0.706,17.351,0.177,19.455,0.330,-0.404),(0.934,21.576,0.152,17.641,0.193,-0.416),
     (1.321,27.601,0.318,14.176,0.221,-0.434),(1.484,30.512,0.760,12.817,0.516,-0.500),(2.330,38.988,0.531,8.632,0.101,-0.431)]
zV,DV,eDV=0.295,7.942,0.075; zb=np.array([b[0] for b in BAO])
y=np.r_[DV,[b[1] for b in BAO],[b[3] for b in BAO]]; C=np.zeros((13,13)); C[0,0]=eDV**2
for i,b in enumerate(BAO): C[1+i,1+i]=b[2]**2; C[7+i,7+i]=b[4]**2; C[1+i,7+i]=C[7+i,1+i]=b[5]*b[2]*b[4]
Ci=np.linalg.inv(C)
def chi_ruler(d,s):
    t=np.r_[(zV*d(zV)**2*s(zV))**(1/3),d(zb),s(zb)]; A=(t@Ci@y)/(t@Ci@t); r=y-A*t; return float(r@Ci@r),float(A)
def lcdm_ds(om):
    zg=np.linspace(0,2.6,2601); E=np.sqrt(om*(1+zg)**3+1-om); dc=integrate.cumulative_trapezoid(1/E,zg,initial=0)
    return (lambda z: np.interp(z,zg,dc)),(lambda z: 1/np.sqrt(om*(1+z)**3+1-om))
g=optimize.minimize_scalar(lambda o: chi_ruler(*lcdm_ds(o))[0],bounds=(0.1,0.6),method="bounded"); out["ruler"]={"lcdm":[float(g.fun),float(g.x)]}
print(f"\nAll 13 DESI ruler numbers, one free scale (standard model, matter share fitted: {g.x:.3f}, chi2 = {g.fun:.1f} for 11 d.o.f.)")
for r in READ:
    c,A=chi_ruler(*booth(r)); out["ruler"][r]=[c,A]; rl=A*0+A                    # A = (c/H0)/r in Booth's units
    print(f"   {READ[r]:42s} chi2 = {c:8.0f} (12 d.o.f.); implied ruler length r = (c/H0)/{A:.1f} = {299792.458/H0B/A:.0f} Mpc for H0 = {H0B}")
# ---------------- 4. supernovae: Pantheon+ (diagonal, repeats averaged) and JLA (full covariance). d_L = (1+z) d. Free offset.
dd=np.loadtxt(os.path.join(CODE,"pantheon_plus_rows.txt")); zr,mr,er=dd.T
zu=np.unique(zr); zp=zu; mp=np.array([np.sum(mr[zr==a]/er[zr==a]**2)/np.sum(1/er[zr==a]**2) for a in zu]); ep=np.array([er[zr==a].min() for a in zu]); wp=1/ep**2
def chi_pp(mu):
    dl=mp-mu; a=np.sum(wp*dl)/wp.sum(); return float(np.sum(wp*(dl-a)**2))
def chi_jla(mu):
    one=np.ones(31); dl=sh.mub-mu; a=(one@sh.Ci@dl)/(one@sh.Ci@one); r=dl-a; return float(r@sh.Ci@r)
def mu_lcdm(z,om):
    zg=np.linspace(0,2.4,2401); dc=integrate.cumulative_trapezoid(1/np.sqrt(om*(1+zg)**3+1-om),zg,initial=0); return 5*np.log10((1+z)*np.interp(z,zg,dc))
gp=optimize.minimize_scalar(lambda o: chi_pp(mu_lcdm(zp,o)),bounds=(0.05,0.9),method="bounded"); gj=optimize.minimize_scalar(lambda o: chi_jla(mu_lcdm(sh.zb,o)),bounds=(0.05,0.9),method="bounded")
out["sn"]={"N_pp":int(len(zp)),"lcdm_pp":[float(gp.fun),float(gp.x)],"lcdm_jla":[float(gj.fun),float(gj.x)]}
print(f"\nSupernovae, d_L = (1+z) d, free offset. Pantheon+ {len(zp)} redshifts: standard model chi2 {gp.fun:.0f} (matter {gp.x:.3f}); JLA 31 bins: {gj.fun:.1f} (matter {gj.x:.3f})")
for r in READ:
    d,_=booth(r); cp=chi_pp(5*np.log10((1+zp)*d(zp))); cj=chi_jla(5*np.log10((1+sh.zb)*d(sh.zb))); out["sn"][r]=[cp,cp-gp.fun,cj,cj-gj.fun]
    print(f"   {READ[r]:42s} Pantheon+ chi2 {cp:7.0f} (worse by {cp-gp.fun:6.0f});  JLA {cj:7.1f} (worse by {cj-gj.fun:6.1f})")
# how far off in magnitudes at z = 1 once both are levelled on the nearby supernovae (z < 0.1)
near=zp<0.1
for r in ("a_ln","b_ln"):
    d,_=booth(r); mB=5*np.log10((1+zp)*d(zp)); mL=mu_lcdm(zp,gp.x); off=np.sum(wp[near]*(mB-mL)[near])/wp[near].sum()
    z1=np.array([1.0]); dm=float(5*np.log10((1+z1)*d(z1))[0]-mu_lcdm(z1,gp.x)[0]-off); out["sn"][r+"_dm_z1"]=dm
    print(f"   {READ[r]:42s} levelled on z < 0.1: supernovae at z = 1 predicted {dm:+.2f} mag fainter than the standard fit")
# binned Pantheon+ residuals for the figure: relative to the best standard fit, in 12 bins
edges=np.quantile(zp[zp>0.02],np.linspace(0,1,13)); res_l=mp-mu_lcdm(zp,gp.x); res_l-=np.sum(wp*res_l)/wp.sum()
out["sn"]["bins"]=[[float(np.median(zp[(zp>=a)&(zp<=b)])),float(np.sum(wp[(zp>=a)&(zp<=b)]*res_l[(zp>=a)&(zp<=b)])/wp[(zp>=a)&(zp<=b)].sum()),float(1/np.sqrt(wp[(zp>=a)&(zp<=b)].sum()))] for a,b in zip(edges[:-1],edges[1:])]
# ---------------- 5. low-redshift expansion: d_L = (c/H0)[z + (1 - q0) z^2 / 2 + ...]
print("\nDeceleration parameter implied at low redshift (standard model with matter 0.315: q0 = %.2f)" % (1.5*0.315-1))
out["q0"]={}
for r in READ:
    d,_=booth(r); h=1e-4; f=lambda z: (1+z)*d(z); q0=1-(f(2*h)-2*f(h)+f(0))/h**2; out["q0"][r]=float(q0); print(f"   {READ[r]:42s} q0 = {q0:+.2f}")
# ---------------- 6. Coste: gravity between photons. Size of the effect over a Hubble distance.
OMG=2.47e-5/0.674**2                                    # photon share of the critical density today (CMB; starlight adds ~10%)
pot=3*OMG/(8*np.pi)*(4*np.pi/3)/2                        # G rho_gamma R^2 / c^2 with R = c/H0, times a sphere's (4 pi/3)/2 potential factor
need=np.expm1(1.0)                                       # redshift at a light-travel distance of c/H0 in the standard picture is ~1.7; any law needs O(1)
out["coste"]=dict(omega_gamma=float(OMG),max_shift=float(pot),shortfall=float(np.log1p(need)/pot))
print(f"\nCoste: photon share of the cosmic density {OMG:.1e}; the whole photon population's potential across c/H0 shifts a photon's energy by at most {pot:.1e}")
print(f"   against a needed ln(1+z) of about 1 over that distance: short by a factor of {1/pot:.0e}")
# ---------------- 7. Kelly's second mechanism: photons crossing deepening potential wells. Redshift = rate of growth of the light-travel time.
# z = d(tau)/dt_emit, where tau is the travel time. Through wells the Shapiro delay is (2/c^3) int |Phi| dl. To give z = H0 d / c at distance d the
# path-averaged 2|Phi|/c^2 must grow by H0 per unit time, i.e. by 1/(age-like 1/H0) ~ 7% per Gyr, everywhere along every line of sight.
H0s=70.0/977.792                                          # per Gyr
phi_cl=(1000/299792.458)**2; phi_lss=1.0e-5             # deepest typical wells: rich clusters, 1000 km/s; large-scale structure ~1e-5
need_rate=H0s                                             # d<2|Phi|/c^2>/dt needed, per Gyr
best_rate=2*phi_lss/1.0                                   # every well along the path deepening by its full depth each Gyr (generous)
out["kelly"]=dict(need_per_gyr=float(need_rate),max_per_gyr=float(best_rate),shortfall=float(need_rate/best_rate),phi_cluster=float(phi_cl))
print(f"\nKelly (evolving wells): needed growth of path-averaged 2|Phi|/c^2: {need_rate:.3f} per Gyr; wells are {phi_cl:.1e} deep (clusters)")
print(f"   even if every well on every path deepened by its full depth each Gyr: {best_rate:.1e} per Gyr; short by {need_rate/best_rate:.0f}x")
# Kelly's first mechanism: a statistical Doppler bias. Peculiar speeds ~300 km/s do not grow with distance.
vpec=300.0; out["kelly"]["doppler"]=[vpec/299792.458,(vpec/299792.458)**2/2]
print(f"   Doppler from random motions: |z| ~ {vpec/299792.458:.4f} at any distance; mean bias (transverse Doppler) {(vpec/299792.458)**2/2:.1e}; observed z reaches 14")
json.dump(out,open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"poster_models.json"),"w"),indent=1)
