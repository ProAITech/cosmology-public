"""Entry 41, second part: the same-camera pair with and without the model corrections, in the two fits of Entry 40.
(a) the power q of the pair alone; (b) the dark-energy preference with one level per camera, the pair's magnitudes swapped for plain ones."""
import sys, os, json
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
import numpy as np
from scipy import optimize
ld=lambda f: np.loadtxt(os.path.join(HERE,"corr",f))
F=np.vstack([ld("found_c.txt"),ld("found_a.txt"),ld("found_d.txt"),ld("found_b.txt")]); P=np.vstack([ld("ps1_a.txt"),ld("ps1_b.txt"),ld("ps1_e.txt"),ld("ps1_c.txt"),ld("ps1_d.txt")])
A=np.loadtxt(os.path.join(HERE,"pantheon_plus_rows.txt")); z,m,e=A.T; w=1/e**2
lab=np.array(list("".join("".join(l.split()[1:]) for l in open(os.path.join(HERE,"survey_labels.txt")) if l.startswith("c"))))
key={(round(r[0],5),round(r[1],3)):i for i,r in enumerate(A)}
def swapped(al,be,raw=False):
    mm=m.copy()
    for r in np.vstack([F,P]):
        i=key[(round(r[0],5),round(r[1],3))]; mm[i]=r[2] if raw else r[2]+al*r[3]-be*r[4]
    return mm
VAR={"corrected, as published":m,"plain, published coefficients":swapped(0.148,3.09),"plain, coefficients fitted on the pair":swapped(0.129,2.51),"raw peak magnitude":swapped(0,0,True)}
pair=(lab=="F")|((lab=="P")&(z>0.1)); out={"q":{},"de":{}}
print("(a) power q of the same-camera pair alone, one level")
for nm,mm in VAR.items():
    def f(q):
        r=(mm-5*np.log10(z*(1+z)**q))[pair]; a=np.sum(w[pair]*r)/w[pair].sum(); return float(np.sum(w[pair]*(r-a)**2))
    g=optimize.minimize_scalar(f,bounds=(-1,3),method="bounded",options={"xatol":1e-6}); a=optimize.brentq(lambda x:f(x)-g.fun-1,-1,g.x); b=optimize.brentq(lambda x:f(x)-g.fun-1,g.x,3)
    print(f"   {nm:40s} q = {g.x:.3f} +- {(b-a)/2:.3f}"); out["q"][nm]=[float(g.x),float((b-a)/2)]
json.dump(out,open(os.path.join(HERE,"corrections_de.json"),"w"),indent=1)
import cosmo as c
from dark_energy_fit import fit, sig, make
keys=sorted(set(lab)); idx={k:i for i,k in enumerate(keys)}; idx["F"]=idx["P"]; g=np.array([idx[k] for k in lab])
def chi_resid(r):
    cc=0.0
    for k in np.unique(g):
        q=g==k; a=np.sum(w[q]*r[q])/w[q].sum(); cc+=np.sum(w[q]*(r[q]-a)**2)
    return float(cc)
bao=lambda mod: c.chi_bao(mod); cmb=lambda mod: c.chi_cmb(mod)
print("(b) ruler + CMB + Pantheon+, one level per camera, the pair's magnitudes as stated")
for nm,mm in VAR.items():
    sn=lambda mod,mm=mm: chi_resid(mm-5*np.log10((1+z)*mod.DM(z))); chi=make([bao,cmb,sn]); L=fit(chi,'lcdm'); Wf=fit(chi,'w0wa'); gn=L['chi']-Wf['chi']
    print(f"   {nm:40s} gain {gn:6.2f} = {sig(gn):.2f} sigma; w0 {Wf['w0']:+.3f}, wa {Wf['wa']:+.3f}",flush=True); out["de"][nm]=dict(gain=float(gn),sigma=sig(gn),w0=Wf['w0'],wa=Wf['wa'])
    json.dump(out,open(os.path.join(HERE,"corrections_de.json"),"w"),indent=1)
