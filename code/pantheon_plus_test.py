"""Workbench Entry 39: the supernova test of Entry 38 on the Pantheon+ table itself.
Data: 1,590 light curves with zHD > 0.01 from Brout et al. (2022), table 7 (VizieR J/ApJ/938/110): zHD, corrected magnitude, and the
diagonal error. The full error matrix is not used. Light curves of one supernova (same zHD) are averaged and keep the smallest error.
Fits are those of Entries 36 and 38: DESI's second-release ruler, the three-number CMB summary, supernovae with a free brightness level."""
import sys, os, json
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
import numpy as np
from scipy.optimize import minimize_scalar
import cosmo as c
from dark_energy_fit import fit, sig, make
from cosmo import Model, ZP, MP, WP
out={}
d=np.loadtxt(os.path.join(HERE,"pantheon_plus_rows.txt")); zr,mr,er=d.T
zu=np.unique(zr); Z=[];M=[];E=[]
for z0 in zu:
    q=zr==z0; w=1/er[q]**2; Z.append(z0); M.append(np.sum(w*mr[q])/w.sum()); E.append(er[q].min())
Z=np.array(Z);M=np.array(M);E=np.array(E);W=1/E**2
print(f"{len(zr)} light curves, {len(Z)} after averaging repeats; {int((Z<0.1).sum())} nearer than z = 0.1"); out["n"]=[len(zr),len(Z),int((Z<0.1).sum())]
def chi_sn(z,m,w,shift=0.0,zb=0.1):
    mm=m+np.where(z<zb,shift,0.0)
    def f(mod):
        r=mm-5*np.log10((1+z)*mod.DM(z)); a=np.sum(w*r)/w.sum(); return float(np.sum(w*(r-a)**2))
    return f
def om_only(f):
    h=0.68; g=lambda om: f(Model(0.02237,om*h*h-0.00064,h)); r=minimize_scalar(g,bounds=(0.15,0.6),method='bounded',options={'xatol':1e-5}); c0=r.fun
    lo=minimize_scalar(lambda x:(g(x)-c0-1)**2,bounds=(0.15,r.x),method='bounded').x; hi=minimize_scalar(lambda x:(g(x)-c0-1)**2,bounds=(r.x,0.6),method='bounded').x
    return float(r.x),float((hi-lo)/2),float(c0)
# 1. supernovae alone
om,e,c0=om_only(chi_sn(Z,M,W)); print(f"1. Supernovae alone, constant dark energy: matter share {om:.3f} +- {e:.3f}, chi2 {c0:.0f} for {len(Z)} (published with the full error matrix: 0.334 +- 0.018)"); out["om_all"]=[om,e,c0]
far=Z>0.1; om2,e2,_=om_only(chi_sn(Z[far],M[far],W[far])); print(f"   beyond z = 0.1 only: {om2:.3f} +- {e2:.3f}"); out["om_far"]=[om2,e2]
omr,er_,_=om_only(chi_sn(zr,mr,1/er**2)); print(f"   every light curve counted separately: {omr:.3f} +- {er_:.3f}"); out["om_raw"]=[omr,er_]
# 2. with ruler and CMB
bao=lambda m: c.chi_bao(m); cmb=lambda m: c.chi_cmb(m)
def cmp(label,snchi,key):
    chi=make([bao,cmb,snchi]); L=fit(chi,'lcdm'); Wf=fit(chi,'w0wa'); g=L['chi']-Wf['chi']
    print(f"   {label:58s} gain {g:6.2f} = {sig(g):.2f} sigma; w0 {Wf['w0']:+.3f}, wa {Wf['wa']:+.3f}; LCDM matter share {L['Om']:.4f}"); out[key]=dict(gain=float(g),sigma=sig(g),w0=Wf['w0'],wa=Wf['wa'],L=L,W=Wf); return L,Wf
print("2. Ruler + CMB + Pantheon+ (DESI, full CMB and error matrix: gain 10.7 = 2.8 sigma, w0 -0.838, wa -0.62)")
La,Wa=cmp("all supernovae beyond z = 0.01",chi_sn(Z,M,W),"all")
Lf,Wf_=cmp("supernovae beyond z = 0.1 only",chi_sn(Z[far],M[far],W[far]),"far")
def chi_step(mod):
    r=M-5*np.log10((1+Z)*mod.DM(Z)); return float(sum(np.sum(W[q]*(r[q]-np.sum(W[q]*r[q])/W[q].sum())**2) for q in (far,~far)))
cmp("all supernovae, free step at z = 0.1",chi_step,"step")
z25=Z>0.025; cmp("supernovae beyond z = 0.025 (the DES sample's lower limit)",chi_sn(Z[z25],M[z25],W[z25]),"z025")
# 3. the seam: nearby minus distant against the curves the distant supernovae prefer
def seam(z,m,w,mod,zb=0.1):
    r=m-5*np.log10((1+z)*mod.DM(z)); n=z<zb; return float(np.sum(w[n]*r[n])/w[n].sum()-np.sum(w[~n]*r[~n])/w[~n].sum()), float((1/w[n].sum()+1/w[~n].sum())**0.5)
print("3. Nearby minus distant supernovae (negative = nearby too bright)")
for nm,P in (("constant, from ruler + CMB + distant Pantheon+",Lf),("changing, from ruler + CMB + distant Pantheon+",Wf_)):
    mod=Model(P['ob'],P['obc'],P['h'],P['w0'],P['wa']); a,ea=seam(Z,M,W,mod); b,eb=seam(ZP,MP,WP,mod)
    print(f"   against {nm}: Pantheon+ {a:+.4f} +- {ea:.4f}; 2018 sample {b:+.4f} +- {eb:.4f}; difference {a-b:+.4f}"); out["seam_"+nm.split(',')[0]]=[a,ea,b,eb]
ref=Model(0.02237,0.1424,0.6736); a,ea=seam(Z,M,W,ref); b,eb=seam(ZP,MP,WP,ref); print(f"   against the Planck constant model: Pantheon+ {a:+.4f} +- {ea:.4f}; 2018 {b:+.4f} +- {eb:.4f}; difference {a-b:+.4f} (Entry 38 predicted +0.064)"); out["seam_planck"]=[a,ea,b,eb]
for zb in (0.05,0.2):
    a,ea=seam(Z,M,W,ref,zb); b,eb=seam(ZP,MP,WP,ref,zb); print(f"      boundary at z = {zb}: Pantheon+ {a:+.4f} +- {ea:.4f}; 2018 {b:+.4f} +- {eb:.4f}; difference {a-b:+.4f}"); out[f"seam_planck_{zb}"]=[a,ea,b,eb]
json.dump(out,open(os.path.join(HERE,"pantheon_plus_test.json"),"w"),indent=1)
# 4. the scan of Entry 38 on Pantheon+ itself
print("4. Pantheon+ with its nearby supernovae shifted")
rows=[]
for s in (-0.08,-0.06,-0.04,-0.02,0.02,0.04):
    chi=make([bao,cmb,chi_sn(Z,M,W,shift=s)]); L=fit(chi,'lcdm'); Wf2=fit(chi,'w0wa'); g=L['chi']-Wf2['chi']; rows.append([s,float(g),sig(g),Wf2['w0'],Wf2['wa']]); print(f"   {s:+.2f}: gain {g:6.2f} = {sig(g):.2f} sigma; w0 {Wf2['w0']:+.3f}, wa {Wf2['wa']:+.3f}")
out["scan"]=rows
json.dump(out,open(os.path.join(HERE,"pantheon_plus_test.json"),"w"),indent=1)
