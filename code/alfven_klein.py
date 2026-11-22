"""Workbench Entry 35: Alfven and Klein's antimatter universe, in numbers.
The model (Klein 1971; Alfven 1966): equal amounts of matter and antimatter; a huge dilute cloud of both contracts under gravity; annihilation
radiation turns the contraction into the expansion we see; redshift is a Doppler shift of real motion.
Five checks, each an order-of-magnitude account or a fit to published numbers:
 1. the cloud as described, against its own gravity;
 2. the energy the observed speeds need, and what share of the matter must annihilate to supply it;
 3. where annihilation energy goes (gamma rays, electrons, neutrinos), and what the gamma-ray sky allows today;
 4. the distance law of a free explosion (the coasting, or Milne, law) against the clustering ruler and the supernovae;
 5. how bright a collision between a matter cluster and an antimatter cluster would be."""
import sys, os, io, json, contextlib
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
import numpy as np
from scipy import integrate, optimize
out={}
G=6.674e-8; c=2.998e10; mp=1.6726e-24; ly=9.461e17; Msun=1.989e33; eV=1.602e-12
# ---- 1. the starting cloud, as quoted from Alfven by Stepanian & Kohandel (arXiv:1308.1849): radius 1e12 light years, one particle per 100 m^3
R=1e12*ly; n=1e-8; rho=n*mp; M=4/3*np.pi*R**3*rho; Rs=2*G*M/c**2; RH=c/np.sqrt(8*np.pi*G*rho/3)
print(f"1. Starting cloud: mass {M/Msun:.1e} Msun; Schwarzschild radius {Rs/ly:.1e} light years = {Rs/R:.1f} x its own radius; free-fall time {np.sqrt(3*np.pi/(32*G*rho))/3.156e7:.1e} yr")
print(f"   the same in other words: the cloud is {R/RH:.1f} times wider than c / sqrt(8 pi G rho / 3)")
out["cloud"]=dict(M_sun=float(M/Msun),Rs_over_R=float(Rs/R),R_over_RH=float(R/RH))
# density at which annihilation matters ('about one particle per cubic unit' in the same source: read as 1 per cm3 and 1 per m3)
for nn,lab in ((1e-6,"1 per m3"),(1.0,"1 per cm3")):
    Rn=R*(n/nn)**(1/3); print(f"   compressed to {lab}: radius {Rn/ly:.1e} light years, {Rs/Rn:.0f} x inside its Schwarzschild radius")
    out["cloud"][lab]=dict(R_ly=float(Rn/ly),Rs_over_R=float(Rs/Rn))
# ---- 2. energy for the observed speeds, redshift read as a Doppler shift
beta=lambda z: ((1+z)**2-1)/((1+z)**2+1); gam=lambda z: ((1+z)**2+1)/(2*(1+z))
# share of energy per annihilation: 5 pions on average, 3.05 charged and 1.93 neutral, with equal mean energies
f_g=1.93/4.98; f_ch=3.05/4.98; f_nu=f_ch*(0.213+0.787*0.65); f_e=f_ch*0.787*0.35
print(f"\n2. Annihilation energy goes: gamma rays {100*f_g:.0f}%, electrons and positrons {100*f_e:.0f}%, neutrinos {100*f_nu:.0f}%")
out["shares"]=dict(gamma=float(f_g),e=float(f_e),nu=float(f_nu))
eps=f_g+f_e; rows=[]
print("   redshift | speed | kinetic energy / rest energy | share of the original matter that must annihilate (all energy used; neutrinos lost)")
for z in (0.5,1,2,3,6,10,14):
    k=gam(z)-1; x1=k/(k+1.0); x2=k/(k+eps); rows.append([z,float(beta(z)),float(k),float(x1),float(x2)]); print(f"   {z:5.1f}    | {beta(z):.3f}c | {k:6.3f} | {100*x1:5.1f}% ; {100*x2:5.1f}%")
out["energy"]=rows
# ---- 3. what annihilation is allowed today
u_cmb=0.2606; nb=0.02237*1.878e-29/mp; ub=nb*938.272e6
I100=7.2e-6; ng=4*np.pi*I100/c; per_ann=2*1.93*0.7       # photons above 100 MeV per annihilation (0.7 of the photons from neutral pions)
x_gamma=2*ng/(per_ann*nb); x_cmb=u_cmb/(eps*ub)
print(f"\n3. Baryons {nb:.2e} per cm3, rest energy {ub:.0f} eV per cm3; microwave background {u_cmb} eV per cm3")
print(f"   gamma rays above 100 MeV: {ng:.2e} per cm3 = {ng/nb:.1e} per baryon; if every one came from annihilation, the share of matter annihilated within the gamma-ray horizon is below {x_gamma:.1e}")
print(f"   if the whole microwave background were thermalised annihilation energy, the share is below {x_cmb:.1e}")
H0=2.2e-18; print(f"   as a rate today: fewer than {x_gamma:.1e} of the baryons per {1/H0/3.156e7/1e9:.1f} billion years")
out["today"]=dict(nb=float(nb),ub=float(ub),ngamma=float(ng),x_gamma=float(x_gamma),x_cmb=float(x_cmb))
# ---- 4. the free-explosion (coasting, Milne) distance law against ruler and supernovae
_cwd=os.getcwd(); os.chdir(HERE)
with contextlib.redirect_stdout(io.StringIO()): import static_hubble as sh
os.chdir(_cwd)
BAO=[(0.510,13.588,0.167,21.863,0.425,-0.459),(0.706,17.351,0.177,19.455,0.330,-0.404),(0.934,21.576,0.152,17.641,0.193,-0.416),(1.321,27.601,0.318,14.176,0.221,-0.434),(1.484,30.512,0.760,12.817,0.516,-0.500),(2.330,38.988,0.531,8.632,0.101,-0.431)]
zV,DV,eDV=0.295,7.942,0.075; zb=np.array([b[0] for b in BAO]); y=np.r_[DV,[b[1] for b in BAO],[b[3] for b in BAO]]; C=np.zeros((13,13)); C[0,0]=eDV**2
for i,b in enumerate(BAO): C[1+i,1+i]=b[2]**2; C[7+i,7+i]=b[4]**2; C[1+i,7+i]=C[7+i,1+i]=b[5]*b[2]*b[4]
Ci=np.linalg.inv(C)
def chi_bao(dM,dH):
    t=np.r_[(zV*dM(zV)**2*dH(zV))**(1/3),dM(zb),dH(zb)]; A=(t@Ci@y)/(t@Ci@t); r=y-A*t; return float(r@Ci@r),r
one=np.ones(31)
def chi_sn(dL):
    d=sh.mub-5*np.log10(dL(sh.zb)); a=(one@sh.Ci@d)/(one@sh.Ci@one); r=d-a; return float(r@sh.Ci@r)
milne_M=lambda z: np.sinh(np.log(1+z)); milne_H=lambda z: 1/(1+z); milne_L=lambda z: (1+z)*np.sinh(np.log(1+z))
zg=np.linspace(0,2.6,5201)
def std(om):
    s=1/np.sqrt(om*(1+zg)**3+1-om); d=integrate.cumulative_trapezoid(s,zg,initial=0); return (lambda t: np.interp(t,zg,d)),(lambda t: np.interp(t,zg,s))
g=optimize.minimize_scalar(lambda o: chi_bao(*std(o))[0]+chi_sn(lambda t: (1+t)*std(o)[0](t)),bounds=(0.1,0.6),method="bounded"); dM,dH=std(g.x)
cb_m,rm=chi_bao(milne_M,milne_H); cs_m=chi_sn(milne_L); cb_s=chi_bao(dM,dH)[0]; cs_s=chi_sn(lambda t:(1+t)*dM(t))
print(f"\n4. Free explosion (coasting) law: ruler chi2 {cb_m:.1f} (13 numbers), supernovae {cs_m:.1f} (31), total {cb_m+cs_m:.1f}")
print(f"   standard law, matter share {g.x:.3f}: ruler {cb_s:.1f}, supernovae {cs_s:.1f}, total {cb_s+cs_s:.1f}")
print(f"   ruler-free ratios: d(2.33)/d(0.51) coasting {milne_M(2.33)/milne_M(0.51):.3f} (measured 2.869 +- 0.053); slope(0.51)/slope(2.33) coasting {milne_H(0.51)/milne_H(2.33):.3f} (measured 2.533 +- 0.057)")
sig=np.sqrt(np.diag(C)); pulls=rm/sig; print("   coasting, pulls on the 13 numbers:"," ".join(f"{p:+.1f}" for p in pulls))
out["milne"]=dict(chi_bao=cb_m,chi_sn=cs_m,chi_bao_std=cb_s,chi_sn_std=cs_s,om=float(g.x),trans=float(milne_M(2.33)/milne_M(0.51)),radial=float(milne_H(0.51)/milne_H(2.33)))
# a free explosion slowed by its own gravity: matter only, open, density parameter Om (0 = coasting)
def open_matter(om):
    Ez=lambda t: np.sqrt(om*(1+t)**3+(1-om)*(1+t)**2); s=1/Ez(zg); d=integrate.cumulative_trapezoid(s,zg,initial=0); k=np.sqrt(1-om)
    return (lambda t: np.sinh(k*np.interp(t,zg,d))/k),(lambda t: np.interp(t,zg,s))
go=optimize.minimize_scalar(lambda o: chi_bao(*open_matter(o))[0]+chi_sn(lambda t:(1+t)*open_matter(o)[0](t)),bounds=(0.001,0.99),method="bounded"); a,b=open_matter(go.x)
print(f"   explosion slowed by gravity (matter only, no dark energy): best density parameter {go.x:.3f}; ruler {chi_bao(a,b)[0]:.1f}, supernovae {chi_sn(lambda t:(1+t)*a(t)):.1f}, total {go.fun:.1f}")
out["open"]=dict(om=float(go.x),chi_bao=chi_bao(a,b)[0],chi_sn=chi_sn(lambda t:(1+t)*a(t)),total=float(go.fun))
# figure arrays: ruler numbers and the two explosion laws, each at its own best scale, divided by the standard curve
def scale(dMf,dHf):
    t=np.r_[(zV*dMf(zV)**2*dHf(zV))**(1/3),dMf(zb),dHf(zb)]; return (t@Ci@y)/(t@Ci@t)
dS,sS=std(0.297); A_s=scale(dS,sS); zz=np.linspace(0.2,2.5,24); fig={"z":[round(float(v),2) for v in zz]}
fig["dM"]=[[float(b[0]),round(float(b[1]/(A_s*dS(b[0]))),4),round(float(b[2]/b[1]),4)] for b in BAO]; fig["dH"]=[[float(b[0]),round(float(b[3]/(A_s*sS(b[0]))),4),round(float(b[4]/b[3]),4)] for b in BAO]
for key,(fm,fh) in (("coast",(milne_M,milne_H)),("open",open_matter(go.x))):
    A=scale(fm,fh); fig[key]={"d":[round(float(v),4) for v in A*fm(zz)/(A_s*dS(zz))],"s":[round(float(v),4) for v in A*fh(zz)/(A_s*sS(zz))]}
out["fig"]=fig; print("   coasting over standard, distance:",fig["coast"]["d"][::6],"slope:",fig["coast"]["s"][::6]); print("   slowed explosion over standard, distance:",fig["open"]["d"][::6],"slope:",fig["open"]["s"][::6])
# ---- 5. a matter cluster meeting an antimatter cluster (Steigman 2008, eqs 1-3)
FX,TkeV=4.7e-12,14.0; T8=TkeV/8.6; Fg=lambda f: 5.4e-14*f/np.sqrt(T8)*FX/(1.4e-23*np.sqrt(T8))
print(f"\n5. Bullet Cluster with a share f of antimatter in the hot gas: gamma rays above 100 MeV {Fg(1):.1e} f photons per cm2 per s")
print(f"   f = 0.1: {Fg(0.1):.1e}; the brightest steady source in the sky (Vela pulsar) is about 1e-5; EGRET limit used by Steigman 3.5e-8 -> f < {3.5e-8/Fg(1):.1e}")
out["bullet"]=dict(coef=float(Fg(1)),f_limit=float(3.5e-8/Fg(1)))
json.dump(out,open(os.path.join(HERE,"alfven_klein.json"),"w"),indent=1)
