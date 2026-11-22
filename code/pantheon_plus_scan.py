"""Entry 39: the nearby-supernova scan of Entry 38 run on Pantheon+ in finer steps, and both samples binned against one curve."""
import sys, os, json
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
import numpy as np
import cosmo as c
from dark_energy_fit import fit, sig, make
from cosmo import Model, ZP, MP, WP
d=np.loadtxt(os.path.join(HERE,"pantheon_plus_rows.txt")); zr,mr,er=d.T
zu=np.unique(zr); Z=np.array(zu); M=np.array([np.sum(mr[zr==z0]/er[zr==z0]**2)/np.sum(1/er[zr==z0]**2) for z0 in zu]); E=np.array([er[zr==z0].min() for z0 in zu]); W=1/E**2
out={}
# both samples binned against the constant curve fitted to ruler + CMB (matter share 0.3008, Entry 36), each relative to its own mean beyond z = 0.1
ref=Model(0.02240485553264522,0.13984741416087412,0.6834652312221218)
edges=[0.01,0.02,0.03,0.05,0.1,0.15,0.2,0.3,0.4,0.6,1.0,2.3]
def binned(z,m,w):
    r=m-5*np.log10((1+z)*ref.DM(z)); far=z>0.1; r=r-np.sum(w[far]*r[far])/w[far].sum(); rows=[]
    for a,b in zip(edges[:-1],edges[1:]):
        q=(z>=a)&(z<b)
        if q.sum()>2: rows.append([float(np.exp(np.mean(np.log(z[q])))),float(np.sum(w[q]*r[q])/w[q].sum()),float(w[q].sum()**-0.5),int(q.sum())])
    return rows
out["bins_plus"]=binned(Z,M,W); out["bins_2018"]=binned(ZP,MP,WP)
for k in ("bins_2018","bins_plus"):
    print(k); [print(f"   z {r[0]:.3f}: {r[1]:+.3f} +- {r[2]:.3f}  ({r[3]})") for r in out[k]]
json.dump(out,open(os.path.join(HERE,"pantheon_plus_scan.json"),"w"),indent=1)
bao=lambda m: c.chi_bao(m); cmb=lambda m: c.chi_cmb(m)
def chi_sn(z,m,w,shift):
    mm=m+np.where(z<0.1,shift,0.0)
    def f(mod):
        r=mm-5*np.log10((1+z)*mod.DM(z)); a=np.sum(w*r)/w.sum(); return float(np.sum(w*(r-a)**2))
    return f
rows=[]
for s in (-0.10,-0.07,-0.05,-0.03,-0.01,0.0,0.01,0.03):
    chi=make([bao,cmb,chi_sn(Z,M,W,s)]); L=fit(chi,'lcdm'); Wf=fit(chi,'w0wa'); g=L['chi']-Wf['chi']; rows.append([s,float(g),sig(g),Wf['w0'],Wf['wa']]); print(f"   {s:+.2f}: gain {g:6.2f} = {sig(g):.2f} sigma; w0 {Wf['w0']:+.3f}, wa {Wf['wa']:+.3f}",flush=True)
    out["scan"]=rows; json.dump(out,open(os.path.join(HERE,"pantheon_plus_scan.json"),"w"),indent=1)
