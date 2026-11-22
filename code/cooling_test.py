"""Workbench Entry 47: the claimed cooling of the CMB near nearby galaxies and filaments, and energy loss in plasma.

Published claim (Luparello et al. 2023; Hansen et al. 2023, 2025; Dominguez Feldman et al. 2026): the CMB is about 27 microkelvin colder
inside 0.2 Mpc of large spiral galaxies in dense filaments at z < 0.02, by the same amount from 41 to 217 GHz.
This script asks what kind of process could do that, using only the published numbers:
  A. the spectrum each kind of process would give across the bands;
  B. the size of the known processes that do shift temperature evenly;
  C. energy loss per electron crossed: the coefficient needed, Brynjolfsson's 'plasma redshift' coefficient, and the electron columns
     that fast radio bursts measure directly;
  D. what the small-scale CMB would look like if every such galaxy out to redshift 1-3 did the same."""
import os, json
import numpy as np
from scipy import integrate, special
HERE=os.path.dirname(os.path.abspath(__file__)); out={}
T0=2.7255; DT=-27.0; eps=-DT*1e-6/T0; print(f"a decrement of {DT} uK is a fractional shift of {eps:.2e}: as a redshift, cz = {eps*299792.458:.1f} km/s")
out["eps"]=eps
# ---------- A. spectra ----------
bands=[41,61,70,94,100,143,217]; x=np.array(bands)/56.79       # h nu / k T for T = 2.7255 K
shapes={"temperature shift (energy loss per photon, Doppler, gravity)":np.ones_like(x),
        "thermal SZ (hot electrons)":x/np.tanh(x/2)-4,
        "photons removed, same share at every frequency":(np.exp(x)-1)/(x*np.exp(x)),
        "photons removed, share rising as frequency (resonant conversion)":1-np.exp(-x),
        "photons removed, share rising as frequency squared":x*(1-np.exp(-x))}
print("\nA. thermodynamic temperature change in each band, scaled to -27 uK at 100 GHz"); out["bands"]=bands; out["spectra"]={}
i100=bands.index(100)
for nm,g in shapes.items():
    v=DT*g/g[i100]; out["spectra"][nm]=[round(float(a),2) for a in v]
    print(f"  {nm:66s} "+" ".join(f"{a:7.1f}" for a in v)+f"   spread {v.max()-v.min():5.1f} uK")
print("  published: the depth varies by less than 1 uK from 41 to 217 GHz")
# conversion of photons to axion-like particles in a magnetic field, away from resonance (natural units, eV)
eV_per_Mpc=1.5637e29; gauss=1.953e-2; hnu=4.1357e-15*100e9          # 1 Mpc^-1 in eV^-1 inverse; 1 gauss in eV^2; photon energy at 100 GHz
out["alp"]=[]
for g_,B_,ne_ in ((1e-11,1.0,1e-3),(6e-11,10.0,1e-4)):
    d_ag=g_*1e-9*B_*1e-6*gauss/2; wp2=(3.7e-11)**2*ne_; d_pl=wp2/(2*hnu); P=(2*d_ag/d_pl)**2
    print(f"  photon-axion mixing at 100 GHz, coupling {g_:.0e}/GeV, field {B_} uG, {ne_:.0e} electrons/cm3: mixing term {d_ag*eV_per_Mpc:.2e} per Mpc, plasma term {d_pl*eV_per_Mpc:.1e} per Mpc; largest conversion away from resonance {P:.1e}; needed about 1e-5")
    out["alp"].append([g_,B_,ne_,float(d_ag*eV_per_Mpc),float(d_pl*eV_per_Mpc),float(P)])
    Lg=0.1*eV_per_Mpc; Pres=np.pi*(2*d_ag)**2*hnu/(wp2/Lg)                         # Landau-Zener, density falling over 100 kpc
    print(f"      at a resonance (axion mass equal to the plasma frequency, {wp2**0.5:.1e} eV): conversion {Pres:.1e}, rising in proportion to frequency"); out["alp"][-1].append(float(Pres))
# ---------- B. known even shifts ----------
c=2.99792458e10; G=6.674e-8; Msun=1.989e33; Mpc=3.0857e24; sT=6.6524e-25; mp=1.6726e-24; fb=0.157; mue=1.14
print("\nB. known processes that shift temperature evenly")
out["halo"]=[]
for Mh in (1e12,1e13):
    Ne=fb*Mh*Msun/(mue*mp)/(np.pi*(0.2*Mpc)**2); tau=sT*Ne; v=eps/tau*c/1e5
    alpha=4*G*Mh*Msun/c**2/(0.2*Mpc); ml=alpha*1e-3*T0*1e6
    print(f"  halo of {Mh:.0e} Msun with every baryon inside 0.2 Mpc: column {Ne:.2e} electrons/cm2 (DM {Ne/3.0857e18:.0f} pc/cm3), optical depth {tau:.1e}; moving-electron (kinetic SZ) shift needs {v:,.0f} km/s along the line of sight; at 300 km/s it gives {tau*1e-3*T0*1e6:.2f} uK; moving-lens shift at 300 km/s: {ml:.3f} uK")
    out["halo"].append(dict(M=Mh,Ne=float(Ne),DM=float(Ne/3.0857e18),tau=float(tau),v_needed=float(v),ksz300=float(tau*1e-3*T0*1e6),lens=float(ml)))
for nm,phi,R in (("galaxy halo",1e-6,0.2),("cluster",1.1e-5,1.0),("supercluster",1e-5,30.0)):
    print(f"  changing gravity (Rees-Sciama order of magnitude) for a {nm}: potential {phi:.0e} c^2, size {R} Mpc: about {2*phi*(R/4448)*T0*1e6:.3f} uK")
# ---------- C. energy loss per electron ----------
print("\nC. energy loss per electron crossed")
kB=sT/2; print(f"  Brynjolfsson's coefficient: half the Thomson cross-section = {kB:.3e} cm2 of full energy loss per electron")
for h_ in out["halo"]: h_["kappa_needed"]=eps/h_["Ne"]; print(f"  coefficient needed for 27 uK through a {h_['M']:.0e} halo: {eps/h_['Ne']:.1e} cm2 = {eps/h_['Ne']/sT:.3f} Thomson (a lower bound if the halo holds fewer baryons)")
# fast radio bursts: name, total dispersion measure (pc/cm3), host redshift (Macquart et al. 2020, gold sample)
FRB=[("180924",361.42,0.3214),("181112",589.27,0.4755),("190102",363.6,0.291),("190608",338.7,0.1178),("190711",593.1,0.522)]
h=0.6774; Ob=0.0486; Om=0.3089; fd=0.85; chi=0.875     # as in the source paper
rho_c=1.8785e-29*h*h; ne0=Ob*rho_c/mp*chi*fd; cH=c/(h*100*1e5/Mpc)/3.0857e18          # c/H0 in pc
E=lambda z:np.sqrt(Om*(1+z)**3+1-Om); DMc=lambda z:ne0*cH*integrate.quad(lambda t:(1+t)/E(t),0,z)[0]
print(f"  mean electron column the standard baryon density gives: {DMc(0.5):.0f} pc/cm3 to z = 0.5, {DMc(1.0):.0f} to z = 1 (diffuse share {fd})")
out["frb"]=[]; out["dmc"]=[[round(z,3),round(float(DMc(z)),1)] for z in np.r_[0.02,np.arange(0.05,1.01,0.05)]]
for nm,dm,z in FRB:
    N=dm*3.0857e18; share=kB*N/np.log(1+z); need=np.log(1+z)/N; cosmic=dm-30-50-50/(1+z)
    print(f"  FRB {nm}: z = {z:.4f}, total DM {dm:6.1f} -> {N:.2e} electrons/cm2; plasma redshift from that column {np.expm1(kB*N):.1e} = {100*share:.3f}% of the redshift; coefficient needed {need/sT:.0f} Thomson; column needed {np.log(1+z)/kB/3.0857e18:,.0f} pc/cm3; beyond the Galaxy and host {cosmic:.0f}, standard prediction {DMc(z):.0f}")
    out["frb"].append(dict(name=nm,z=z,dm=dm,share=float(share),need_thomson=float(need/sT),dm_needed=float(np.log(1+z)/kB/3.0857e18),cosmic=float(cosmic),pred=float(DMc(z))))
sh=np.array([f["share"] for f in out["frb"]]); rat=np.array([f["cosmic"]/f["pred"] for f in out["frb"]]); sf=np.array([f["dm_needed"]/f["dm"] for f in out["frb"]])
print(f"  all five: plasma redshift supplies {100*sh.mean():.2f}% of the redshift (range {100*sh.min():.2f}-{100*sh.max():.2f}); the electrons fall short by a factor {sf.mean():.0f} ({sf.min():.0f}-{sf.max():.0f}); column beyond Galaxy and host over the standard prediction: {rat.mean():.2f} +- {rat.std(ddof=1)/5**.5:.2f}")
out["frb_summary"]=dict(share=float(sh.mean()),short=float(sf.mean()),short_range=[float(sf.min()),float(sf.max())],ratio=[float(rat.mean()),float(rat.std(ddof=1)/5**.5)])
# ---------- D. if every such galaxy did it ----------
print("\nD. small-scale CMB if every large spiral out to redshift zmax cooled the CMB by 27 uK inside 0.2 Mpc")
hh=0.674; cH0=2997.92458/hh; V=4*np.pi/3*((0.02*cH0)**3-(0.004*cH0)**3)*(1-np.sin(np.radians(10))); 
for Ngal,lab in ((2500,"all large spirals in the claim's sample"),(1250,"the half in dense filaments")):
    n=Ngal/V; print(f"  {lab}: {Ngal} in {V:.2e} Mpc3 -> {n:.2e} per Mpc3")
    E2=lambda z:np.sqrt(0.315*(1+z)**3+0.685); DMz=lambda z:cH0*integrate.quad(lambda t:1/E2(t),0,z)[0]
    def Dl(l,zmax,zmin=0.05,rd=0.2,dT=27.0):
        def f(z):
            dm=DMz(z); th=rd*(1+z)/dm; a=l*th; W=2*special.j1(a)/a if a>1e-6 else 1.0
            return n*(cH0/E2(z))*dm*dm*(dT*np.pi*th*th*W)**2
        return l*(l+1)/(2*np.pi)*integrate.quad(f,zmin,zmax,limit=200)[0]
    ncross=lambda zmax:integrate.quad(lambda z:n*np.pi*0.2**2*(1+z)**2*cH0/E2(z),0,zmax)[0]
    res=[]
    for zmax in (0.5,1.0,2.0,3.0):
        d=Dl(3000,zmax); res.append([zmax,float(d),float(ncross(zmax))]); print(f"    out to z = {zmax}: {ncross(zmax):.2f} such discs on an average sightline; added blackbody power at l = 3000: {d:6.1f} uK2; measured, all of it: 3.0 +- 1.0; largest decrement that fits under 5 uK2: {27*np.sqrt(5/d):.1f} uK")
    out.setdefault("power",{})[lab]=res
json.dump(out,open(os.path.join(HERE,"cooling_test.json"),"w"),indent=1)
