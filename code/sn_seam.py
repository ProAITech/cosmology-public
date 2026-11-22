"""Entry 38: does the scan depend on where "nearby" ends? The same shift applied below z = 0.05 and below z = 0.2,
and the offset the nearby supernovae show against the curve the distant ones prefer."""
import sys, os, json
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
import numpy as np
import cosmo as c
from dark_energy_fit import fit, sig, make
from cosmo import ZP, MP, WP
bao=lambda m: c.chi_bao(m); cmb=lambda m: c.chi_cmb(m); out={"breaks":[]}
def chi_sn(mag):
    def f(m):
        d=mag-5*np.log10((1+ZP)*m.DM(ZP)); a=np.sum(WP*d)/WP.sum(); return float(np.sum(WP*(d-a)**2))
    return f
# the offset of the nearby supernovae from the curve the distant ones prefer (with ruler and CMB)
far=ZP>0.1
for kind in ('lcdm','w0wa'):
    W=fit(make([bao,cmb,lambda m: c.chi_pan(m,sel=far)]),kind); m=c.Model(W['ob'],W['obc'],W['h'],W['w0'],W['wa'])
    d=MP-5*np.log10((1+ZP)*m.DM(ZP)); off=np.sum(WP[~far]*d[~far])/WP[~far].sum()-np.sum(WP[far]*d[far])/WP[far].sum(); e=(1/WP[~far].sum()+1/WP[far].sum())**0.5
    out["offset_"+kind]=[float(off),float(e),W['w0'],W['wa']]; print(f"{kind}: distant-only fit w0 {W['w0']:+.3f} wa {W['wa']:+.3f}; nearby minus distant offset {off:+.4f} +- {e:.4f} mag")
for zb in (0.05,0.2):
    n=int((ZP<zb).sum())
    for s in (0.02,0.05,0.08):
        chi=make([bao,cmb,chi_sn(MP+np.where(ZP<zb,s,0.0))]); L=fit(chi,'lcdm'); W=fit(chi,'w0wa'); g=L['chi']-W['chi']
        out["breaks"].append([zb,n,s,float(g),sig(g),W['w0'],W['wa']]); print(f"break {zb} ({n} supernovae) shift {s:+.2f}: gain {g:6.2f} = {sig(g):.2f} sigma; w0 {W['w0']:+.3f}, wa {W['wa']:+.3f}")
json.dump(out,open(os.path.join(HERE,"sn_seam.json"),"w"))
