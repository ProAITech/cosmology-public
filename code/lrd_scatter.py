"""Workbench Entry 119: little red dots - plasma cocoons or overmassive black holes?

JWST's 'little red dots' (LRDs) at z ~ 4-9 show broad Balmer lines. Read as gas orbiting a black hole, the widths give masses of 1e7-1e8 Msun,
often more than the host's stars, which has been cited as a problem for early-universe models. Rusakov et al. (2026) argue the broad wings come
from Thomson scattering by free electrons in a dense ionized cocoon, so the masses are 100 times smaller. Others find counter-examples.
This is a plasma-physics question, so we compute it:
 1. Monte Carlo of line photons scattering off thermal electrons in a uniform sphere: the wing shape and width for each optical depth.
 2. The electron column that the measured wings need, against the X-ray weakness of LRDs (Compton-thick above 1.5e24 cm^-2).
 3. A test that does not depend on the model: Thomson scattering is the same at every wavelength, so all hydrogen lines formed together must
    show the same wings (Brazzini et al. 2025 measured three lines in one LRD).
 4. Whatever the masses, does the LRD population strain the expanding model? Their black-hole mass density against the baryons, and their
    number density against halo abundance (Sheth-Tormen, Planck 2018, the power spectrum of Entry 101)."""
import os, sys, json, math
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
CODE=os.path.abspath(os.path.join(HERE,"..","artifact-files","4379c2ac-1195-4a41-978e-eb44b945a265","code")) if not os.path.exists(os.path.join(HERE,"progenitors.py")) else HERE
sys.path.insert(0,CODE)
rng=np.random.default_rng(119)
C=299792.458; KB=1.380649e-23; ME=9.1093837e-31; SIGT=6.6524587e-29
def iso(n):
    v=rng.normal(size=(n,3)); return v/np.linalg.norm(v,axis=1)[:,None]
def run(tauR,Te,n=60000,maxs=4000):
    """photons from the centre of a uniform sphere of radial Thomson depth tauR; returns velocity shift (km/s) and number of scatterings"""
    sig=math.sqrt(KB*Te/ME)/1e3                                   # km/s per component of electron velocity
    pos=np.zeros((n,3)); d=iso(n); shift=np.zeros(n); ns=np.zeros(n,int); alive=np.ones(n,bool)
    for it in range(maxs):
        idx=np.where(alive)[0]
        if not len(idx): break
        L=rng.exponential(size=len(idx)); p=pos[idx]+d[idx]*L[:,None]
        out=np.linalg.norm(p,axis=1)>tauR; alive[idx[out]]=False
        sc=idx[~out]
        if not len(sc): break
        pos[sc]=p[~out]; ve=rng.normal(size=(len(sc),3))*sig; dn=iso(len(sc))
        shift[sc]+=np.sum(ve*(dn-d[sc]),axis=1); d[sc]=dn; ns[sc]+=1
    return shift,ns
def wing(shift,ns):
    """exponential scale W of the scattered photons' |shift| between its 50th and 95th percentiles, and the FWHM of all photons"""
    s=np.abs(shift[ns>0]); lo,hi=np.percentile(s,[50,95]); m=(s>=lo)&(s<=hi)
    h,e=np.histogram(s[m],bins=30); c=0.5*(e[1:]+e[:-1]); k=h>0
    W=-1/np.polyfit(c[k],np.log(h[k]),1)[0]
    return float(W)
out={"runs":[]}
print("Thomson-scattered line from the centre of a uniform sphere (isotropic scattering):")
for Te in (1e4,2e4,5e4):
    for tau in (0.05,0.1,0.25,0.5,1,2,3,5,8):
        sh,ns=run(tau,Te); W=wing(sh,ns); f=float(np.mean(ns>0))
        out["runs"].append(dict(Te=Te,tau=tau,W=W,f_scat=f,mean_N=float(ns.mean())))
        print(f"   T_e = {Te:.0e} K, tau = {tau:4.2f}: scattered share {f:.2f}, mean scatterings {ns.mean():5.1f}, exponential wing W = {W:5.0f} km/s")
# what the measured wings need
R=[r for r in out["runs"] if r["Te"]==1e4]; Ws=np.array([r["W"] for r in R]); ts=np.array([r["tau"] for r in R])
for Wobs in (1000,1500,2000):
    t=float(np.interp(Wobs,Ws,ts)); Ncol=t/(SIGT*1e4)
    out[f"tau_for_W{Wobs}"]=t; out[f"col_for_W{Wobs}"]=Ncol
    print(f"   a wing of W = {Wobs} km/s at 1e4 K needs tau ≈ {t:.1f}: an electron column of {Ncol:.1e} cm^-2 (Compton-thick above 1.5e24)")
# black-hole masses: virial M ∝ FWHM^2; if the broad width is scattering, the orbital width is the narrow core
for core,broad in ((300,3000),(500,3000),(1000,3000)):
    print(f"   if the orbiting gas has FWHM {core} km/s and scattering makes it look {broad}: masses overestimated {broad**2/core**2:.0f}x")
# the three-line test (Brazzini et al. 2025, GN-28074, z = 2.26): exponential widths and scattered shares
lines={"Halpha":(1110,0.74),"Hbeta":(4688,0.25),"Pabeta":(1424,0.07)}; out["three_lines"]=lines
print("   GN-28074: W = 1110 (Hα), 4688 (Hβ), 1424 (Paβ) km/s; scattered share 0.74, 0.25, 0.07. Thomson scattering would give one W and one share for all three")
fs=np.array([r["f_scat"] for r in R]); Wt=Ws                      # 1e4 K grid; W scales as sqrt(T_e) at fixed tau (checked: 2e4/1e4 ratio below)
rat=[out["runs"][i+9]["W"]/out["runs"][i]["W"] for i in range(9)]; out["W_ratio_2e4_1e4"]=rat
print("   check: W(2e4 K)/W(1e4 K) at each tau =", " ".join(f"{x:.2f}" for x in rat), "(sqrt 2 = 1.41)")
out["line_fit"]={}
for k,(Wobs,f) in lines.items():
    t=float(np.interp(f,fs,ts)); W1=float(np.interp(t,ts,Wt)); Te=1e4*(Wobs/W1)**2
    out["line_fit"][k]=dict(tau=t,W_1e4=W1,Te_needed=Te)
    print(f"   {k:7s}: share {f:.2f} -> tau {t:.2f}; at 1e4 K that gives W = {W1:4.0f} km/s; the observed {Wobs} km/s needs T_e = {Te:.1e} K")
# profile shapes for the figure (tau = 1 and 3 at 1e4 K) against a Gaussian of the same FWHM core
prof={}
for tau in (1,3):
    sh,ns=run(tau,1e4,n=200000); h,e=np.histogram(sh,bins=np.linspace(-6000,6000,121)); prof[str(tau)]=dict(v=(0.5*(e[1:]+e[:-1])).tolist(),f=(h/h.max()).tolist())
out["profiles"]=prof
# 4. cosmology: number density, black-hole mass density, halo abundance
import contextlib, io
with contextlib.redirect_stdout(io.StringIO()): import progenitors as P
def st_n_above(M,z):                                        # Sheth-Tormen number density above M (per Mpc^3), Planck 2018
    Ms=np.geomspace(M,1e16,300); s=np.array([P.sigM(m) for m in Ms])*P.D(z)/P.D(0)*P.s8/P.sigma(8/P.h)
    nu=1.686/s; a,p,A=0.707,0.3,0.3222; anu=math.sqrt(a)*nu
    f=A*math.sqrt(2/math.pi)*anu*(1+anu**(-2*p))*np.exp(-anu**2/2)
    dlns=np.gradient(np.log(s),np.log(Ms)); dn=P.rho_m/Ms*f*np.abs(dlns)          # dn/dlnM
    return float(np.trapezoid(dn,np.log(Ms)))
nLRD=6e-5; out["n_LRD"]=nLRD; out["halo"]={}
for z in (5,7):
    row={}
    for M in (1e10,1e11,1e12):
        row[f"{M:.0e}"]=st_n_above(M,z)
    out["halo"][str(z)]=row
    print(f"   halos at z = {z}: n(>1e10) {row['1e+10']:.1e}, n(>1e11) {row['1e+11']:.1e}, n(>1e12) {row['1e+12']:.1e} per Mpc^3 (LRDs: about {nLRD:.0e})")
rho_b=P.rho_m*P.fb; out["rho_b"]=rho_b
for Mbh in (1e5,1e7,1e8):
    r=nLRD*Mbh; print(f"   black holes of {Mbh:.0e} Msun in every LRD: {r:.1e} Msun/Mpc^3, {r/rho_b:.1e} of all baryons")
    out[f"rhoBH_{Mbh:.0e}"]=r
# matching numbers: the halo mass M with n(>M) = n_LRD
from scipy.optimize import brentq
for z in (5,7):
    Mm=10**brentq(lambda lm: math.log(st_n_above(10**lm,z)/nLRD),10,13); out[f"Mmatch_z{z}"]=Mm
    print(f"   at z = {z}, putting every LRD in its own halo means halos above {Mm:.1e} Msun, holding {Mm*P.fb:.1e} Msun of baryons each")
# growth time: Eddington-limited growth with radiative efficiency 0.1 doubles the mass every ln2 * 50 Myr
tE=0.45*0.1/0.9*1e3; out["t_efold_Myr"]=tE
for zs,Ms in ((25,100),(20,1e5)):
    for zt,Mt in ((7.04,5e7),(5,1e8)):
        need=tE*math.log(Mt/Ms); have=(P.age(zt)-P.age(zs))*1e3
        out[f"grow_{zs}_{Ms:.0e}_{zt}"]=(need,have)
        print(f"   seed {Ms:.0e} Msun at z = {zs} -> {Mt:.0e} Msun by z = {zt}: needs {need:4.0f} Myr at the Eddington limit; the expanding model gives {have:4.0f} Myr (duty cycle {need/have:.2f})")
rho_bh0=4.6e5; out["rho_bh0"]=rho_bh0
print(f"   today's black holes hold about {rho_bh0:.1e} Msun/Mpc^3 (Marconi et al. 2004): LRDs at 1e8 each would be {6e3/rho_bh0*100:.1f}% of it")
json.dump(out,open(os.path.join(HERE,"lrd_scatter.json"),"w"),indent=1)
