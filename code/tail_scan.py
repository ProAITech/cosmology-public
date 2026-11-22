"""Workbench Entry 53: can a non-thermal tail in the primordial plasma fix lithium?
Charged nuclei must tunnel, so their reaction rates are set by the few nuclei several kT up the tail.  Change the tail and the rates
change, most for the reactions with the highest barrier.  3He + 4He -> 7Be (charges 2 x 2) makes the lithium; d + d and d + p (1 x 1)
burn the deuterium.  The tail is given the Tsallis form used in the published proposal (Hou et al. 2017):
     f(E) ~ [1 - (q-1) E/kT]^(1/(q-1)),   cut off at E = kT/(q-1) for q > 1;   q < 1 gives the power-law ('kappa') tails of space plasmas.
Rates for charged pairs are rescaled by the ratio of the tunnelling integral with this tail to the thermal one; neutron reactions
(cross-section ~ 1/v) do not depend on the shape; photodisintegration keeps the thermal photon spectrum."""
import os, json, sys
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicSpline
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import bbn_net as bn
HERE=os.path.dirname(os.path.abspath(__file__))
# charge product, reduced mass (u), and a shape for the astrophysical S-factor (E in MeV); only the shape matters for the ratio
PAIR={'dpg':(1,0.6718,lambda E:0.214+5.56*E+9.3*E*E),'ddn':(1,1.00705,lambda E:1+4.0*E),'ddp':(1,1.00705,lambda E:1+2.4*E),
      'tdn':(1,1.2076,lambda E:1/((E-0.048)**2+0.040**2)),'hdp':(2,1.2076,lambda E:1/((E-0.21)**2+0.13**2)),
      'hag':(4,1.7200,lambda E:np.exp(-0.6*E)),'tag':(2,1.7200,lambda E:np.exp(-0.8*E)),'lipa':(3,0.8812,lambda E:1+2.0*E)}
def weight(qq):
    if abs(qq-1)<1e-9: return (lambda x: np.exp(-x)), np.inf
    if qq>1:
        xm=1/(qq-1); w=lambda x: (1-(qq-1)*x)**(1/(qq-1)) if x<xm else 0.0
    else:
        xm=np.inf; w=lambda x: (1+(1-qq)*x)**(-1/(1-qq))
    norm=quad(lambda x: np.sqrt(x)*w(x),0,min(xm,400),limit=200)[0]/(np.sqrt(np.pi)/2)       # same number of particles as the thermal case
    return (lambda x: w(x)/norm), xm
TG=np.exp(np.linspace(np.log(0.02),np.log(12.0),34))
def ratio_table(qq):
    w,xm=weight(qq); tab={}
    for k,(zz,mu,S) in PAIR.items():
        b=0.98951*zz*np.sqrt(mu); r=[]
        for T9 in TG:
            kT=T9/bn.MEV9; E0=(b*kT/2)**(2/3); hi=min(xm*kT,E0+12*np.sqrt(E0*kT)+30*kT)
            pts=[p for p in (E0,0.048,0.21) if 0<p<hi]
            num=quad(lambda E: S(E)*np.exp(-b/np.sqrt(E))*w(E/kT),1e-9,hi,points=pts,limit=300)[0]
            den=quad(lambda E: S(E)*np.exp(-b/np.sqrt(E)-E/kT),1e-9,E0+12*np.sqrt(E0*kT)+30*kT,points=pts,limit=300)[0]
            r.append(num/den)
        tab[k]=CubicSpline(np.log(TG),np.log(np.maximum(r,1e-300)))
    return lambda k,T9: float(np.exp(tab[k](np.log(min(max(T9,0.02),12.0))))), tab
base=bn.run(); out={'base':{k:float(base[k]) for k in ('Yp','D','He3','Li7')},'scan':[]}
print(f"thermal plasma: Yp {base['Yp']:.4f}  D/H {base['D']*1e5:.3f}e-5  3He/H {base['He3']*1e5:.3f}e-5  7Li/H {base['Li7']*1e10:.2f}e-10")
f,tab=ratio_table(1.0755); print("\nrate with the tail of q = 1.0755 over the thermal rate, at T9 = 1.0, 0.7, 0.5, 0.3:")
for k in PAIR: print(f"  {k:5s} (Z1Z2 = {PAIR[k][0]}): "+"  ".join(f"{f(k,T):.2f}" for T in (1.0,0.7,0.5,0.3)))
out['ratios_q10755']={k:[f(k,T) for T in (1.0,0.7,0.5,0.3)] for k in PAIR}
print("\n   q      Yp     D/H e-5  x base   3He/H e-5   7Li/H e-10  x base")
for qq in (0.92,0.95,0.97,0.98,0.99,1.0,1.005,1.01,1.02,1.03,1.04,1.05,1.06,1.069,1.0755,1.082,1.09,1.10):
    r=base if qq==1.0 else bn.run(tail=ratio_table(qq)[0])
    print(f"  {qq:6.4f}  {r['Yp']:.4f}  {r['D']*1e5:6.3f}  {r['D']/base['D']:5.3f}   {r['He3']*1e5:6.3f}     {r['Li7']*1e10:6.2f}    {r['Li7']/base['Li7']:5.3f}")
    out['scan'].append([qq,float(r['Yp']),float(r['D']),float(r['He3']),float(r['Li7'])])
print("\npublished for q = 1.069-1.082 (Hou et al. 2017): D/H 3.14-3.25 (x1.22-1.26 of their thermal 2.57), 3He/H 1.46-1.50 (x1.40-1.44), 7Li/H 1.62-1.90 (x0.31-0.36)")
json.dump(out,open(os.path.join(HERE,'tail_scan.json'),'w'),indent=1)
