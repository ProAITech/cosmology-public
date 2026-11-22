"""Workbench Entry 53, checks.
(1) How thermal is the plasma during nucleosynthesis?  Coulomb collision rates against nuclear reaction rates and the expansion rate.
(2) A magnetic field holding 8 to 13% of the energy (the published magnetic route): what the extra energy alone does to the yields,
    and what happens when that energy is returned to the photons afterwards.
(3) The baryon density track, for the deuterium-lithium plane."""
import os, json, sys
import numpy as np
from scipy.integrate import quad
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import bbn_net as bn
HERE=os.path.dirname(os.path.abspath(__file__)); out={}
HBARC=1.973270e-11; E2=1.439964e-13                   # MeV cm
base=bn.run(dense=True); sol=base['sol']; Tt=sol.y[0]
def at(T9):
    i=np.argmin(np.abs(Tt-T9)); return sol.y[1:,i]
print("(1) collisions against reactions")
print("   T9    kT keV  baryons/cm3  pairs/cm3   lnL  Coulomb rate of a nucleus (1/s)   fastest reaction per nucleus (1/s)   expansion (1/s)   reactions per collision")
out['therm']=[]
for T9 in (1.0,0.8,0.6,0.4,0.3,0.2,0.1):
    kT=T9/bn.MEV9; x=bn.ME/kT; re,S,dS=bn.plasma(T9); nb=6.10e-10*bn.NGAM*T9**3*S; Y=at(T9); R=nb/bn.NAV; f=bn.rates(T9)
    I=quad(lambda y: y*y*np.exp(-np.sqrt(y*y+x*x))/(1+np.exp(-np.sqrt(y*y+x*x))),0,x+60)[0]; npair=2/np.pi**2*(kT/HBARC)**3*I   # electrons plus positrons from pairs
    ne=npair+nb*(Y[bn.P]+2*Y[bn.A]); TeV=kT*1e6
    lamD=np.sqrt(kT/(4*np.pi*ne*E2)); bmin=max(E2/kT,HBARC/np.sqrt(938.272*kT)); lnL=np.log(lamD/bmin)
    nu_ie=3.25e-9*ne*lnL/TeV**1.5                       # energy exchange of a proton with the electron-positron gas (non-relativistic formula)
    nu_ii=4.80e-8*nb*Y[bn.P]*lnL/TeV**1.5               # proton-proton
    nu=nu_ie+nu_ii
    rx={'d':R*((f['ddn']+f['ddp'])*Y[bn.D]+f['dpg']*Y[bn.P]+f['tdn']*Y[bn.T]+f['hdp']*Y[bn.H]),'3He':R*(f['hdp']*Y[bn.D]+f['hent']*Y[bn.N]+f['hag']*Y[bn.A]),'t':R*(f['tdn']*Y[bn.D]+f['tag']*Y[bn.A])}
    k=max(rx,key=rx.get); Tn=T9*(4/11*S)**(1/3); rho=bn.RHOG*T9**4*(1+re)+3*7/8*bn.RHOG*Tn**4; H=np.sqrt(8*np.pi*bn.GN*rho/3)
    print(f"  {T9:4.1f}  {kT*1e3:6.1f}   {nb:9.2e}   {npair:9.2e}  {lnL:4.1f}   {nu:9.2e} (pairs {nu_ie:8.1e}, ions {nu_ii:8.1e})   {rx[k]:9.2e} ({k})   {H:9.2e}   {rx[k]/nu:8.1e}")
    out['therm'].append(dict(T9=T9,kT_keV=kT*1e3,nb=nb,npair=npair,lnL=lnL,nu=nu,nu_ie=nu_ie,nu_ii=nu_ii,rx=rx[k],rx_sp=k,H=H,ratio=rx[k]/nu))
print("   the published tail (q = 1.07) removes 70% of the nuclei at 5 kT; reactions can remove at most the last column per collision time")
print("\n(2) a magnetic field holding a share fB of the early energy density (scaling as radiation)")
b=dict(Yp=base['Yp'],D=base['D'],Li7=base['Li7']); out['fB']=[]
for fB in (0.0,0.04,0.08,0.11,0.13):
    r=b if fB==0 else bn.run(fB=fB); dN=6.14*fB/(1-fB); B0=np.sqrt(8*np.pi*fB/(1-fB)*5.375*(4/11)**(4/3)*4.175e-13)*1e6
    print(f"  fB = {fB:4.2f} (extra neutrino species {dN:4.2f}; {B0:4.2f} microgauss today if it survived): Yp {r['Yp']:.4f} ({r['Yp']-b['Yp']:+.4f})  D/H {r['D']*1e5:.3f}e-5 (x{r['D']/b['D']:.3f})  7Li/H {r['Li7']*1e10:.2f}e-10 (x{r['Li7']/b['Li7']:.3f})")
    out['fB'].append([fB,dN,B0,float(r['Yp']),float(r['D']),float(r['Li7'])])
r=bn.run(fB=0.11,eta10=6.10*1.135); print(f"  fB = 0.11 with the field's energy later given to the photons (13.5% more photons, so 13.5% more baryons per photon during BBN):")
print(f"      Yp {r['Yp']:.4f} ({r['Yp']-b['Yp']:+.4f})  D/H {r['D']*1e5:.3f}e-5 (x{r['D']/b['D']:.3f})  7Li/H {r['Li7']*1e10:.2f}e-10 (x{r['Li7']/b['Li7']:.3f})")
out['fB_dissipated']=[0.11,float(r['Yp']),float(r['D']),float(r['Li7'])]
print("  the field that eases the Hubble tension in Entry 48 is 5 to 10 picogauss: "+f"{1.18e-6/1e-11:.1e} to {1.51e-6/5e-12:.1e} times weaker in strength than 1.18 to 1.51 microgauss")
print("\n(3) baryon density track"); out['eta']=[]
for e10 in (3.5,4.0,4.5,5.0,5.5,6.1,7.0,8.0):
    r=bn.run(eta10=e10); print(f"  eta10 = {e10:3.1f}: Yp {r['Yp']:.4f}  D/H {r['D']*1e5:.3f}e-5  7Li/H {r['Li7']*1e10:.2f}e-10"); out['eta'].append([e10,float(r['Yp']),float(r['D']),float(r['Li7'])])
json.dump(out,open(os.path.join(HERE,'plasma_checks.json'),'w'),indent=1)
