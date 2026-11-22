"""Entry 39, part 3: the seam inside one telescope. Foundation (nearby) and the Pan-STARRS Medium Deep supernovae (distant) were both
observed with the Pan-STARRS1 camera. Neither was in the 2018 sample's nearby anchor. Survey membership from the Survey column of the
Pantheon+ table (FOUND, PS1MD), matched to the main rows by redshift and magnitude."""
import sys, os, json
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
import numpy as np
import cosmo as c
from dark_energy_fit import fit, sig, make
from cosmo import Model, ZP, MP, WP
A=np.loadtxt(os.path.join(HERE,"pantheon_plus_rows.txt")); key={(round(r[0],5),round(r[1],3)):i for i,r in enumerate(A)}
tag=np.zeros(len(A),int)
for fn,t in (("pantheon_plus_foundation.txt",1),("pantheon_plus_ps1md.txt",2)):
    for r in np.loadtxt(os.path.join(HERE,fn)): tag[key[(round(r[0],5),round(r[1],3))]]=t
z,m,e=A.T; w=1/e**2; out={}
def avg(q):
    zz=z[q]; mm=m[q]; ee=e[q]; zu=np.unique(zz); return zu,np.array([np.sum(mm[zz==a]/ee[zz==a]**2)/np.sum(1/ee[zz==a]**2) for a in zu]),np.array([ee[zz==a].min() for a in zu])
def seam(near,far,mod):
    zn,mn,en=avg(near); zf,mf,ef=avg(far); rn=mn-5*np.log10((1+zn)*mod.DM(zn)); rf=mf-5*np.log10((1+zf)*mod.DM(zf)); wn=1/en**2; wf=1/ef**2
    return float(np.sum(wn*rn)/wn.sum()-np.sum(wf*rf)/wf.sum()), float((1/wn.sum()+1/wf.sum())**0.5), len(zn), len(zf)
ref=Model(0.02237,0.1424,0.6736)
near=z<0.1; far=~near; F=tag==1; P=tag==2; old=near&(tag==0)
print("Nearby minus distant against the Planck constant model (negative = nearby too bright)")
cases=[("all nearby, all distant",near,far),("Foundation, Pan-STARRS distant (one camera)",F&near,P&far),("Foundation, all distant",F&near,far),("Foundation, distant other than Pan-STARRS",F&near,far&~P),
       ("older nearby surveys, all distant",old,far),("older nearby surveys, Pan-STARRS distant",old,P&far),("all nearby, Pan-STARRS distant",near,P&far)]
for nm,a,b in cases:
    s,es,n1,n2=seam(a,b,ref); print(f"   {nm:46s} {s:+.4f} +- {es:.4f}   ({n1} nearby, {n2} distant)"); out[nm]=[s,es,n1,n2]
zn,mn,en=ZP[ZP<0.1],MP[ZP<0.1],WP[ZP<0.1]; r=MP-5*np.log10((1+ZP)*ref.DM(ZP)); n=ZP<0.1
s18=float(np.sum(WP[n]*r[n])/WP[n].sum()-np.sum(WP[~n]*r[~n])/WP[~n].sum()); e18=float((1/WP[n].sum()+1/WP[~n].sum())**0.5); print(f"   2018 sample                                    {s18:+.4f} +- {e18:.4f}"); out["2018"]=[s18,e18]
# Foundation against the older nearby surveys directly, same redshift range
zf,mf,ef=avg(F&near); zo,mo,eo=avg(old); rf=mf-5*np.log10((1+zf)*ref.DM(zf)); ro=mo-5*np.log10((1+zo)*ref.DM(zo)); d=np.sum(rf/ef**2)/np.sum(1/ef**2)-np.sum(ro/eo**2)/np.sum(1/eo**2); ed=(1/np.sum(1/ef**2)+1/np.sum(1/eo**2))**0.5
print(f"Foundation minus older nearby surveys: {d:+.4f} +- {ed:.4f}"); out["found_minus_old"]=[float(d),float(ed)]
json.dump(out,open(os.path.join(HERE,"pantheon_plus_ps1.json"),"w"),indent=1)
# fits with ruler and CMB
bao=lambda mod: c.chi_bao(mod); cmb=lambda mod: c.chi_cmb(mod)
def chi_sub(q):
    zz,mm,ee=avg(q); ww=1/ee**2
    def f(mod):
        r=mm-5*np.log10((1+zz)*mod.DM(zz)); a=np.sum(ww*r)/ww.sum(); return float(np.sum(ww*(r-a)**2))
    return f,len(zz)
print("Ruler + CMB + part of Pantheon+")
for nm,q in (("Pan-STARRS camera only (Foundation + Medium Deep)",F|P),("without Foundation",~F),("older nearby surveys + all distant",old|far),("Foundation + all distant",(F&near)|far)):
    f,nn=chi_sub(q); chi=make([bao,cmb,f]); L=fit(chi,'lcdm'); Wf=fit(chi,'w0wa'); g=L['chi']-Wf['chi']; print(f"   {nm:50s} ({nn}) gain {g:6.2f} = {sig(g):.2f} sigma; w0 {Wf['w0']:+.3f}, wa {Wf['wa']:+.3f}",flush=True); out["fit_"+nm]=[nn,float(g),sig(g),Wf['w0'],Wf['wa']]
    json.dump(out,open(os.path.join(HERE,"pantheon_plus_ps1.json"),"w"),indent=1)
