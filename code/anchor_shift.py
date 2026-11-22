"""Workbench Entry 43: what moved the nearby anchor between the 2018 supernova compilation and Pantheon+?
Supernovae in both tables are matched by name. For each: the change in corrected magnitude, the change in the redshift used for the
Hubble diagram, and what that redshift change does to the predicted magnitude. Nearby: CfA and CSP light curves below z = 0.1.
Distant reference: every supernova in two redshift slices (0.20-0.215 and 0.40-0.43) present in both tables.
2018: Scolnic et al., VizieR J/ApJ/859/101/fullz (Name, zcmb, corrected magnitude). 2022: Brout et al., J/ApJ/938/110/table7 (Name, Survey, zHD, corrected magnitude)."""
import os, json, sys
HERE=os.path.dirname(os.path.abspath(__file__))
import numpy as np
from scipy import integrate
rd=lambda f: [l.split() for l in open(os.path.join(HERE,f)) if l.strip()]
norm=lambda s: s.lstrip("0") if s.isdigit() else s
n18={norm(r[0]):(float(r[1]),float(r[2])) for r in rd("near2018.txt")}; f18={norm(r[0]):(float(r[1]),float(r[2])) for r in rd("far2018.txt")}
nP={}; 
for r in rd("nearplus.txt"): nP.setdefault(norm(r[0]),[]).append((r[1],float(r[2]),float(r[3])))
fP={norm(r[0]):(r[1],float(r[2]),float(r[3])) for r in rd("farplus.txt")}
zg=np.linspace(0,1.0,100001); dc=integrate.cumulative_trapezoid(1/np.sqrt(0.3*(1+zg)**3+0.7),zg,initial=0)
mu=lambda z: 5*np.log10((1+z)*np.interp(z,zg,dc))
out={}
# nearby
rows=[]
for k,(z0,m0) in n18.items():
    if k in nP:
        zs=np.mean([a[1] for a in nP[k]]); ms=np.mean([a[2] for a in nP[k]]); rows.append((k,z0,zs,m0,ms,"+".join(sorted({a[0] for a in nP[k]}))))
N=np.array([[r[1],r[2],r[3],r[4]] for r in rows]); print(f"nearby: {len(n18)} in the 2018 table, {len(nP)} supernovae with CfA or CSP light curves in Pantheon+, {len(rows)} in both")
rowsF=[(k,f18[k][0],fP[k][1],f18[k][1],fP[k][2],fP[k][0]) for k in f18 if k in fP]; Fa=np.array([[r[1],r[2],r[3],r[4]] for r in rowsF]); print(f"distant: {len(f18)} and {len(fP)} in the slices, {len(rowsF)} in both: "+", ".join(f"{s} {sum(1 for r in rowsF if r[5]==s)}" for s in sorted({r[5] for r in rowsF})))
def parts(A,label):
    z0,z1,m0,m1=A.T; dm=m1-m0; dmu=mu(z1)-mu(z0); dr=dm-dmu; se=lambda x: x.std(ddof=1)/len(x)**0.5
    print(f"   {label:34s} n={len(dm):3d}: magnitude {dm.mean():+.4f} +- {se(dm):.4f} (scatter {dm.std(ddof=1):.3f}) | redshift {np.mean((z1-z0)/z0)*100:+.2f}% -> predicted magnitude {dmu.mean():+.4f} +- {se(dmu):.4f} | residual {dr.mean():+.4f} +- {se(dr):.4f}")
    return dict(n=int(len(dm)),dm=[float(dm.mean()),float(se(dm))],dmu=[float(dmu.mean()),float(se(dmu))],dr=[float(dr.mean()),float(se(dr))],dz=float(np.mean((z1-z0)/z0)),dm_sd=float(dm.std(ddof=1)),dmu_sd=float(dmu.std(ddof=1)))
print("Change from 2018 to Pantheon+, same supernovae")
out["near"]=parts(N,"nearby, CfA + CSP"); out["far"]=parts(Fa,"distant slices")
for lab,sel in (("   nearby z < 0.02",N[:,0]<0.02),("   nearby 0.02 - 0.03",(N[:,0]>=0.02)&(N[:,0]<0.03)),("   nearby 0.03 - 0.05",(N[:,0]>=0.03)&(N[:,0]<0.05)),("   nearby beyond 0.05",N[:,0]>=0.05)): out[lab.strip()]=parts(N[sel],lab)
for lab,sel in (("   distant z ~ 0.2",Fa[:,0]<0.3),("   distant z ~ 0.4",Fa[:,0]>0.3)): out[lab.strip()]=parts(Fa[sel],lab)
for s in ("SDSS","PS1MD","SNLS"):
    q=np.array([r[5]==s for r in rowsF]); out["far_"+s]=parts(Fa[q],"   distant, "+s)
for s in ("CSP","CfA"):
    q=np.array([s in r[5] for r in rows]); out["near_"+s]=parts(N[q],"   nearby with a "+s+" light curve")
a,b=out["near"],out["far"]; e=lambda k: float(np.hypot(a[k][1],b[k][1]))
print(f"\nSeam shift on the same supernovae (nearby minus distant): {a['dr'][0]-b['dr'][0]:+.4f} +- {e('dr'):.4f}")
print(f"   from the magnitudes: {a['dm'][0]-b['dm'][0]:+.4f} +- {e('dm'):.4f};  from the redshifts: {-(a['dmu'][0]-b['dmu'][0]):+.4f} +- {e('dmu'):.4f}")
out["shift"]=dict(total=[a['dr'][0]-b['dr'][0],e('dr')],mag=[a['dm'][0]-b['dm'][0],e('dm')],z=[-(a['dmu'][0]-b['dmu'][0]),e('dmu')])
out["rows_near"]=[[r[0],r[1],r[2],r[3],round(r[4],4),r[5]] for r in rows]; out["rows_far"]=[[r[0],r[1],r[2],r[3],r[4],r[5]] for r in rowsF]
# the largest redshift changes
dz=(N[:,1]-N[:,0])/N[:,0]; j=np.argsort(-abs(dz))[:8]; print("   largest redshift changes nearby:",[(rows[i][0],rows[i][1],round(rows[i][2],5)) for i in j])
print(f"   nearby redshift change: mean {dz.mean()*100:+.2f}%, median {np.median(dz)*100:+.2f}%, scatter {dz.std()*100:.2f}%; share raised {np.mean(dz>0)*100:.0f}%")
out["dz_near"]=dict(mean=float(dz.mean()),median=float(np.median(dz)),sd=float(dz.std()),up=float(np.mean(dz>0)))
json.dump(out,open(os.path.join(HERE,"anchor_shift.json"),"w"),indent=1)
# ---- is it the supernovae or their treatment? The same 126, placed in each compilation's own Hubble diagram
A=np.loadtxt(os.path.join(HERE,"pantheon_plus_rows.txt")); zA,mA,eA=A.T; wA=1/eA**2; keyP={(round(a,5),round(b,3)):1/c**2 for a,b,c in A}
sys.path.insert(0,HERE)
from cosmo import ZP,MP,WP
k18={(round(float(a),5),round(float(b),3)):float(c) for a,b,c in zip(ZP,MP,WP)}
farP=zA>0.1; refP=np.sum(wA[farP]*(mA-mu(np.minimum(zA,0.99)))[farP])/wA[farP].sum() if False else None
zg2=np.linspace(0,2.4,240001); dc2=integrate.cumulative_trapezoid(1/np.sqrt(0.3*(1+zg2)**3+0.7),zg2,initial=0); mu2=lambda z: 5*np.log10((1+z)*np.interp(z,zg2,dc2))
refP=np.sum(wA[farP]*(mA-mu2(zA))[farP])/wA[farP].sum(); far18=ZP>0.1; ref18=np.sum(WP[far18]*(MP-mu2(ZP))[far18])/WP[far18].sum()
common={r[0] for r in rows}
def seamP(names):
    num=den=0.0
    for k in names:
        for s,zz,mm in nP[k]:
            w_=keyP[(round(zz,5),round(mm,3))]; num+=w_*(mm-mu2(zz)-refP); den+=w_
    return num/den, den**-0.5
def seam18(names):
    num=den=0.0
    for k in names:
        zz,mm=n18[k]; w_=k18.get((round(zz,5),round(mm,3)))
        if w_ is None: w_=min(k18.items(),key=lambda t: abs(t[0][0]-zz)+abs(t[0][1]-mm))[1]
        num+=w_*(mm-mu2(zz)-ref18); den+=w_
    return num/den, den**-0.5
a=seamP(common); b=seamP(set(nP)-common); c18=seam18(common); d18=seam18(set(n18)-common)
print(f"\nSeams against one curve (matter share 0.30), each relative to its own compilation's distant supernovae:")
print(f"   Pantheon+: the {len(common)} also in the 2018 table {a[0]:+.3f} +- {a[1]:.3f}; the {len(set(nP)-common)} CfA and CSP supernovae new to Pantheon+ {b[0]:+.3f} +- {b[1]:.3f}")
print(f"   2018 table: the same {len(common)} {c18[0]:+.3f} +- {c18[1]:.3f}; the other {len(set(n18)-common)} nearby {d18[0]:+.3f} +- {d18[1]:.3f}")
print(f"   the same {len(common)} supernovae moved by {a[0]-c18[0]:+.3f}")
out["seams"]=dict(plus_common=a,plus_new=b,old_common=c18,old_other=d18)
# figure data
out["fig"]=dict(near=[[round(r[1],4),round((r[4]-r[3])-(mu(r[2])-mu(r[1])),3)] for r in rows],far=[[round(r[1],4),round((r[4]-r[3])-(mu(r[2])-mu(r[1])),3)] for r in rowsF])
json.dump(out,open(os.path.join(HERE,"anchor_shift.json"),"w"),indent=1,default=float)
