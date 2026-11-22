"""Entry 38: two checks on the supernova scan.
1. The matter share a supernova-only fit returns for each shift, to compare with the published supernova-only values.
2. The scan carried on to +0.12 mag, for Union3."""
import sys, os, json
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
import numpy as np
from scipy.optimize import minimize_scalar
import cosmo as c
from dark_energy_fit import fit, sig, make
from cosmo import ZP, MP, WP, Model
print("supernovae nearer than z = 0.1:", int((ZP<0.1).sum()), "of", len(ZP))
def chi_sn(mag):
    def f(m):
        d=mag-5*np.log10((1+ZP)*m.DM(ZP)); a=np.sum(WP*d)/WP.sum(); return float(np.sum(WP*(d-a)**2))
    return f
def om_only(s):
    mag=MP+np.where(ZP<0.1,s,0.0); f=chi_sn(mag); h=0.68; ob=0.02237
    g=lambda om: f(Model(ob,om*h*h-0.00064,h))          # 0.00064 = neutrino density, so om is the full matter share
    r=minimize_scalar(g,bounds=(0.15,0.55),method='bounded',options={'xatol':1e-5}); c0=r.fun
    lo=minimize_scalar(lambda x:(g(x)-c0-1)**2,bounds=(0.15,r.x),method='bounded').x; hi=minimize_scalar(lambda x:(g(x)-c0-1)**2,bounds=(r.x,0.55),method='bounded').x
    return float(r.x), float((hi-lo)/2)
out={"om":[],"scan":[]}
for s in np.round(np.arange(-0.02,0.1201,0.01),2):
    om,e=om_only(s); out["om"].append([float(s),om,e]); print(f"{s:+.2f}: supernova-only matter share {om:.3f} +- {e:.3f}")
bao=lambda m: c.chi_bao(m); cmb=lambda m: c.chi_cmb(m)
for s in (0.10,0.11,0.12):
    mag=MP+np.where(ZP<0.1,s,0.0); chi=make([bao,cmb,chi_sn(mag)]); L=fit(chi,'lcdm'); W=fit(chi,'w0wa'); g=L['chi']-W['chi']
    out["scan"].append([s,float(g),sig(g),float(W['w0']),float(W['wa'])]); print(f"{s:+.2f}: gain {g:6.2f} = {sig(g):.2f} sigma; w0 {W['w0']:+.3f}, wa {W['wa']:+.3f}")
json.dump(out,open(os.path.join(HERE,"sn_shift_check.json"),"w"))
