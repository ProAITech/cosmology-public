"""Workbench Entry 37: two follow-ups on the dark-energy signal of Entry 36.
 1. The galaxies at redshift 0.71 carry half of it. How does that measurement compare with the same redshift in DESI's first release and in
    SDSS, and what happens to the signal when one is swapped for another?
 2. The two readings of the signal (changing dark energy, or a shorter sound horizon) predict different Hubble constants.
    Which do the direct measurements favour?"""
import sys, os, json
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
import numpy as np
from scipy import optimize
import cosmo as c
from dark_energy_fit import fit, sig, make   # reruns Entry 36's fits on import (prints suppressed below)
out={}
# ---------- ruler data sets: ('V', z, value, err) or ('MH', z, dM, eM, dH, eH, rho)
DR2=[('V',0.295,7.942,0.075),('MH',0.510,13.588,0.167,21.863,0.425,-0.459),('MH',0.706,17.351,0.177,19.455,0.330,-0.404),('MH',0.934,21.576,0.152,17.641,0.193,-0.416),('MH',1.321,27.601,0.318,14.176,0.221,-0.434),('MH',1.484,30.512,0.760,12.817,0.516,-0.500),('MH',2.330,38.988,0.531,8.632,0.101,-0.431)]
DR1=[('V',0.295,7.93,0.15),('MH',0.510,13.62,0.25,20.98,0.61,-0.445),('MH',0.706,16.85,0.32,20.08,0.60,-0.420),('MH',0.930,21.71,0.28,17.88,0.35,-0.389),('MH',1.317,27.79,0.69,13.82,0.42,-0.444),('V',1.491,26.07,0.67),('MH',2.330,39.71,0.94,8.52,0.17,-0.477)]
def sdss(rho): return [('V',0.15,4.47,0.17),('MH',0.38,10.23,0.17,25.00,0.76,rho),('MH',0.51,13.36,0.21,22.33,0.58,rho),('MH',0.698,17.86,0.33,19.33,0.53,rho),('V',0.85,18.33,0.60),('MH',1.48,30.69,0.80,13.26,0.55,rho),('MH',2.33,37.6,1.9,8.93,0.28,rho),('MH',2.33,37.3,1.7,9.08,0.34,rho)]
def chi_set(data):
    def f(m):
        rd=m.rd(); tot=0.0
        for d in data:
            if d[0]=='V': tot+=(((d[1]*m.DM(d[1])**2*m.DH(d[1]))**(1/3)/rd-d[2])/d[3])**2
            else:
                _,z,dM,eM,dH,eH,rho=d; a=(m.DM(z)/rd-dM)/eM; b=(m.DH(z)/rd-dH)/eH; tot+=(a*a+b*b-2*rho*a*b)/(1-rho*rho)
        return float(tot)
    return f
cmb=lambda m: c.chi_cmb(m)
def compare(label,parts):
    chi=make(parts); L=fit(chi,'lcdm'); W=fit(chi,'w0wa'); d=L['chi']-W['chi']
    print(f"  {label:66s} LCDM chi2 {L['chi']:6.2f} (Om {L['Om']:.4f}, H0 {100*L['h']:.2f}) | changing: chi2 {W['chi']:6.2f}, w0 {W['w0']:+.2f}, wa {W['wa']:+.2f} | gain {d:5.2f} = {sig(d):.1f} sigma")
    return dict(lcdm=L,w0wa=W,gain=float(d),sigma=sig(d))
print("1. The measurement near redshift 0.7 in three data sets, against the best constant-dark-energy fit to DESI's second release with the CMB")
ref=fit(make([chi_set(DR2),cmb]),'lcdm'); m=c.Model(ref['ob'],ref['obc'],ref['h']); rd=m.rd(); rows=[]
for name,(z,dM,eM,dH,eH) in (("SDSS (2020), z = 0.698",(0.698,17.86,0.33,19.33,0.53)),("DESI first release (2024), z = 0.706",(0.706,16.85,0.32,20.08,0.60)),("DESI second release (2025), z = 0.706",(0.706,17.351,0.177,19.455,0.330))):
    pm=m.DM(z)/rd; ph=m.DH(z)/rd; rows.append([name,float(dM/pm),float(eM/pm),float(dH/ph),float(eH/ph)])
    print(f"   {name:38s} distance {dM/pm:.3f} +- {eM/pm:.3f} ({(dM-pm)/eM:+.1f} sigma) | slope {dH/ph:.3f} +- {eH/ph:.3f} ({(dH-ph)/eH:+.1f} sigma) of the fitted curve")
out["z07"]=rows
# per-tracer chi2 under the two best fits (DR2 + CMB)
W=fit(make([chi_set(DR2),cmb]),'w0wa'); mw=c.Model(W['ob'],W['obc'],W['h'],W['w0'],W['wa']); per=[]
for d,nm in zip(DR2,c.BAO_NAMES):
    a=chi_set([d])(m); b=chi_set([d])(mw); per.append([nm,float(a),float(b)])
print("   chi2 by tracer, constant -> changing: "+"; ".join(f"{nm} {a:.1f} -> {b:.1f}" for nm,a,b in per)+f"; CMB {cmb(m):.1f} -> {cmb(mw):.1f}")
out["per_tracer"]=per; out["cmb_split"]=[float(cmb(m)),float(cmb(mw))]
print("\n2. The signal with other rulers (each with the CMB summary)")
R={}
R["dr2"]=compare("DESI second release",[chi_set(DR2),cmb])
swap=[d for d in DR2 if not (d[0]=='MH' and d[1]==0.706)]+[('MH',0.698,17.86,0.33,19.33,0.53,-0.4)]
R["dr2_sdss07"]=compare("DESI second release, its z = 0.71 point replaced by SDSS's",[chi_set(swap),cmb])
swap1=[d for d in DR2 if not (d[0]=='MH' and d[1]==0.706)]+[DR1[2]]
R["dr2_dr107"]=compare("DESI second release, its z = 0.71 point replaced by the first release's",[chi_set(swap1),cmb])
R["dr1"]=compare("DESI first release",[chi_set(DR1),cmb])
R["sdss"]=compare("SDSS, correlations taken as -0.4",[chi_set(sdss(-0.4)),cmb])
R["sdss0"]=compare("SDSS, correlations taken as 0",[chi_set(sdss(0.0)),cmb])
out["rulers"]=R
print("\n3. The Hubble constant as referee")
# our three readings from DESI second release + CMB summary; error on the shorter-sound-horizon value by profile
A=fit(make([chi_set(DR2),cmb]),'lcdm',free_alpha=True)
def prof_h(hh):
    chi=make([chi_set(DR2),cmb])
    def f(x):
        ob,obc,al=x[0]*0.02223,x[1]*0.142,x[2]
        if not (0.01<ob<0.04 and 0.05<obc<0.3 and 0.8<al<1.2): return 1e9
        mm=c.Model(ob,obc,hh,-1,0,al); return chi(mm) if mm.ode>0 else 1e9
    return optimize.minimize(f,[A['ob']/0.02223,A['obc']/0.142,A['alpha']],method='Nelder-Mead',options=dict(xatol=1e-6,fatol=1e-6,maxiter=4000)).fun
lo=optimize.brentq(lambda v: prof_h(v)-A['chi']-1,A['h']-0.04,A['h']); hi=optimize.brentq(lambda v: prof_h(v)-A['chi']-1,A['h'],A['h']+0.04)
print(f"   constant dark energy, shorter sound horizon: H0 = {100*A['h']:.2f} +{100*(hi-A['h']):.2f} -{100*(A['h']-lo):.2f}")
readings=[["constant, standard sound horizon",68.17,0.28,"DESI + CMB"],["constant, shorter sound horizon",100*A['h'],100*max(hi-A['h'],A['h']-lo),"our fit"],
          ["changing, ruler + CMB",63.6,2.1,"DESI + CMB"],["changing, + Pantheon+",67.51,0.59,"DESI"],["changing, + Union3",65.91,0.84,"DESI"],["changing, + DES 5-year",66.74,0.56,"DESI"]]
direct=[["SH0ES Cepheids and supernovae",73.04,1.04],["Chicago-Carnegie red giants",70.39,float(np.sqrt(1.22**2+1.33**2+0.70**2))],["Lens time delays (TDCOSMO)",71.6,3.6],["Megamasers",73.9,3.0]]
def mean(rows):
    w=np.array([1/r[2]**2 for r in rows]); v=np.array([r[1] for r in rows]); return float(np.sum(w*v)/w.sum()),float(1/np.sqrt(w.sum()))
m4=mean(direct); m3=mean(direct[1:]); print(f"   direct measurements: all four {m4[0]:.2f} +- {m4[1]:.2f}; without SH0ES {m3[0]:.2f} +- {m3[1]:.2f}")
tab=[]
for name,h,e,src in readings:
    c4=sum(((d[1]-h)**2/(d[2]**2+e**2)) for d in direct); c3=sum(((d[1]-h)**2/(d[2]**2+e**2)) for d in direct[1:])
    s4=(m4[0]-h)/np.hypot(m4[1],e); s3=(m3[0]-h)/np.hypot(m3[1],e); tab.append([name,float(h),float(e),float(c4),float(c3),float(s4),float(s3)])
    print(f"   {name:36s} H0 {h:5.2f} +- {e:.2f} | against all four: chi2 {c4:5.1f}, mean is {s4:.1f} sigma above | without SH0ES: chi2 {c3:4.1f}, {s3:.1f} sigma above")
out["readings"]=tab; out["direct"]=direct; out["means"]=dict(all=m4,noshoes=m3); out["alpha_h"]=[100*A['h'],100*(A['h']-lo),100*(hi-A['h'])]
json.dump(out,open(os.path.join(HERE,"followup.json"),"w"),indent=1)
