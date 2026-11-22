"""Workbench Entry 40: Pantheon+ with no survey's brightness level trusted against another's.
Each survey (or each camera) gets its own free level, so only the shape of the Hubble diagram inside a survey is used.
That removes every calibration difference between surveys, at the cost of most of the lever arm.
Rows: pantheon_plus_rows.txt (zHD, corrected magnitude, diagonal error); survey of each row: survey_labels.txt."""
import sys, os, json
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
import numpy as np
from scipy import integrate, optimize
A=np.loadtxt(os.path.join(HERE,"pantheon_plus_rows.txt")); z,m,e=A.T; w=1/e**2
lab=np.array(list("".join("".join(l.split()[1:]) for l in open(os.path.join(HERE,"survey_labels.txt")) if l.startswith("c")))); assert len(lab)==len(z)
NAMES={"F":"Foundation","P":"Pan-STARRS Medium Deep","S":"SDSS","N":"SNLS","D":"DES","H":"Hubble","C":"CSP","1":"CfA1","2":"CfA2","K":"CfA3K","3":"CfA3S","4":"CfA4p1","5":"CfA4p2","L":"LOSS1","M":"LOSS2","U":"SOUSA","Z":"LOWZ","I":"CNIa0.02"}
def groups(mode):
    """group index per row. 'one': one level for all; 'survey': one per survey; 'camera': one per survey, Foundation and Pan-STARRS sharing; 'seam': nearby/distant at z = 0.1"""
    if mode=="one": g=np.zeros(len(z),int)
    elif mode=="seam": g=(z>0.1).astype(int)
    else:
        keys=sorted(set(lab)); idx={k:i for i,k in enumerate(keys)}
        if mode=="camera": idx["F"]=idx["P"]
        g=np.array([idx[k] for k in lab])
    return g
def chi_resid(r,g,sel=None):
    s=np.ones(len(z),bool) if sel is None else sel; c=0.0
    for k in np.unique(g[s]):
        q=s&(g==k); a=np.sum(w[q]*r[q])/w[q].sum(); c+=np.sum(w[q]*(r[q]-a)**2)
    return float(c)
mu=lambda dl:5*np.log10(dl); zg=np.linspace(0,2.4,2401)
def dl_lcdm(om): dc=integrate.cumulative_trapezoid(1/np.sqrt(om*(1+zg)**3+1-om),zg,initial=0); return (1+z)*np.interp(z,zg,dc)
def fit1(f,lo,hi):
    g=optimize.minimize_scalar(f,bounds=(lo,hi),method="bounded",options={"xatol":1e-6})
    try: a=optimize.brentq(lambda x:f(x)-g.fun-1,lo,g.x)
    except ValueError: a=lo
    try: b=optimize.brentq(lambda x:f(x)-g.fun-1,g.x,hi)
    except ValueError: b=hi
    return float(g.x),float((b-a)/2),float(g.fun)
out={"counts":{NAMES[k]:[int((lab==k).sum()),float(z[lab==k].min()),float(z[lab==k].max())] for k in NAMES}}
print("1. The static law, with brightness levels freed step by step")
LAWS={"lerner":lambda t:t*np.sqrt(1+t),"log_stretched":lambda t:np.log(1+t)*(1+t),"coasting":lambda t:t*(1+t/2)}
out["static"]={}
for mode,desc in (("one","one level for everything"),("seam","two levels: nearer and farther than z = 0.1"),("camera","one level per camera"),("survey","one level per survey")):
    g=groups(mode); q=fit1(lambda qq:chi_resid(m-mu(z*(1+z)**qq),g),-1,3); om=fit1(lambda o:chi_resid(m-mu(dl_lcdm(o)),g),0.02,0.95); d={k:chi_resid(m-mu(fn(z)),g)-om[2] for k,fn in LAWS.items()}
    print(f"   {desc:44s} q = {q[0]:.3f} +- {q[1]:.3f} | LCDM matter {om[0]:.3f} +- {om[1]:.3f}, chi2 {om[2]:.0f} | above LCDM: Lerner {d['lerner']:+.1f}, ln+stretch {d['log_stretched']:+.1f}, coasting {d['coasting']:+.1f}")
    out["static"][mode]=dict(q=q,om=om,d=d,ngroups=int(len(np.unique(g))))
print("2. Inside single surveys")
out["single"]={}
for k in ("S","P","D","N","F","H"):
    s=lab==k; g=np.zeros(len(z),int); q=fit1(lambda qq:chi_resid(m-mu(z*(1+z)**qq),g,s),-3,5); om=fit1(lambda o:chi_resid(m-mu(dl_lcdm(o)),g,s),0.02,0.95); dL=chi_resid(m-mu(z*np.sqrt(1+z)),g,s)-om[2]
    print(f"   {NAMES[k]:24s} {int(s.sum()):4d} light curves, z {z[s].min():.2f}-{z[s].max():.2f}: q = {q[0]:.3f} +- {q[1]:.3f}; LCDM matter {om[0]:.3f} +- {om[1]:.3f}; Lerner above LCDM {dL:+.1f}")
    out["single"][NAMES[k]]=dict(n=int(s.sum()),zmin=float(z[s].min()),zmax=float(z[s].max()),q=q,om=om,dlerner=dL)
big=np.isin(lab,list("SPDN")); g=groups("survey"); q=fit1(lambda qq:chi_resid(m-mu(z*(1+z)**qq),g,big),-1,3); print(f"   the four large distant surveys together, each with its own level: q = {q[0]:.3f} +- {q[1]:.3f}"); out["big4"]=q
json.dump(out,open(os.path.join(HERE,"pantheon_plus_surveys.json"),"w"),indent=1)
print("3. Dark energy with ruler and CMB, supernova levels freed step by step")
import cosmo as c
from dark_energy_fit import fit, sig, make
bao=lambda mod: c.chi_bao(mod); cmb=lambda mod: c.chi_cmb(mod); out["de"]={}
for mode,desc in (("one","one level for everything"),("seam","two levels: nearer and farther than z = 0.1"),("camera","one level per camera"),("survey","one level per survey")):
    g=groups(mode); sn=lambda mod: chi_resid(m-5*np.log10((1+z)*mod.DM(z)),g); chi=make([bao,cmb,sn]); L=fit(chi,'lcdm'); Wf=fit(chi,'w0wa'); gn=L['chi']-Wf['chi']
    print(f"   {desc:44s} gain {gn:6.2f} = {sig(gn):.2f} sigma; w0 {Wf['w0']:+.3f}, wa {Wf['wa']:+.3f}",flush=True); out["de"][mode]=dict(gain=float(gn),sigma=sig(gn),w0=Wf['w0'],wa=Wf['wa'])
    json.dump(out,open(os.path.join(HERE,"pantheon_plus_surveys.json"),"w"),indent=1)
