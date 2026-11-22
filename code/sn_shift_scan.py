"""Entry 38, block 5 in finer steps: dark-energy preference against a brightness shift of the supernovae nearer than z = 0.1."""
import sys, os, json
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
import numpy as np
import cosmo as c
from dark_energy_fit import fit, sig, make
from cosmo import ZP, MP, WP
def chi_pan_shift(s):
    mag=MP+np.where(ZP<0.1,s,0.0)
    def f(m):
        d=mag-5*np.log10((1+ZP)*m.DM(ZP)); a=np.sum(WP*d)/WP.sum(); return float(np.sum(WP*(d-a)**2))
    return f
bao=lambda m: c.chi_bao(m); cmb=lambda m: c.chi_cmb(m); rows=[]
for s in np.round(np.arange(-0.03,0.0901,0.01),2):
    chi=make([bao,cmb,chi_pan_shift(s)]); L=fit(chi,'lcdm'); W=fit(chi,'w0wa'); g=L['chi']-W['chi']; rows.append([float(s),float(g),sig(g),float(W['w0']),float(W['wa'])]); print(f"{s:+.2f}: gain {g:6.2f} = {sig(g):.2f} sigma; w0 {W['w0']:+.3f}, wa {W['wa']:+.3f}")
json.dump(rows,open(os.path.join(HERE,"sn_shift_scan.json"),"w"))
