"""Workbench Entry 55: the count dipole after one erratum and one new measurement, and what redshift slices would show.
(A) The Quaia quasar dipole after its authors' erratum.
(B) The DESI 'redshift dipole' of Wu & Xia (2026): four tracers, four velocity vectors.  Do they agree with the CMB, with each other,
    and with a Sun moving 2.2 or 3.67 times faster?
(C) The kinematic count dipole in redshift slices.  It is NOT the same in every slice: sources are shifted across the slice edges,
    so the dipole follows the slope of the redshift distribution (von Hausegger & Dalang 2024).  A forecast for telling a faster Sun
    from an intrinsic dipole."""
import os, json
import numpy as np
from scipy.optimize import minimize, brentq
from scipy.stats import chi2, norm
from scipy.integrate import quad
HERE=os.path.dirname(os.path.abspath(__file__)); out={}
C=299792.458; VC=369.82; BETA=VC/C
def unit(l,b): l,b=np.radians(l),np.radians(b); return np.array([np.cos(b)*np.cos(l),np.cos(b)*np.sin(l),np.sin(b)])
CMB=unit(264.02,48.25)
# ---------------- (A)
print("(A) Quaia, G < 20.0, 40-degree mask (Mittal, Oayda & Lewis 2024, with erratum)")
meas,lo95,hi95=0.011,0.005,0.006; s1=(lo95+hi95)/2/1.96
for nm,exp_ in (("as first published",0.0080),("after the erratum",0.0048)):
    print(f"   expected {exp_:.4f} ({nm}): measured/expected = {meas/exp_:.2f} +- {s1/exp_:.2f}; expectation is {(meas-exp_)/(lo95/1.96):.1f} sigma below the measurement; excess as an added dipole {1e2*(meas-exp_):.2f} +- {1e2*s1:.2f}%")
out['quaia']=dict(meas=meas,sig=s1,ratio_old=meas/0.0080,ratio_new=meas/0.0048,ratio_err=s1/0.0048,z=(meas-0.0048)/(lo95/1.96),added=meas-0.0048)
add={'CatWISE quasars':(0.0085,0.24*0.0072),'radio galaxies':(0.0115,0.49*0.0115/2.67),'Quaia quasars':(meas-0.0048,s1)}
w=np.array([1/v[1]**2 for v in add.values()]); m=np.array([v[0] for v in add.values()]); mean=(w*m).sum()/w.sum(); ch=(w*(m-mean)**2).sum()
print("   excess as an added dipole: "+"; ".join(f"{k} {1e2*v[0]:.2f} +- {1e2*v[1]:.2f}%" for k,v in add.items())+f"; mean {1e2*mean:.2f} +- {1e2/np.sqrt(w.sum()):.2f}%, chi2 {ch:.1f} for 2")
mult={'CatWISE quasars':(2.22,0.24),'radio galaxies':(3.67,0.49),'Quaia quasars':(meas/0.0048,s1/0.0048)}
w2=np.array([1/v[1]**2 for v in mult.values()]); m2=np.array([v[0] for v in mult.values()]); mean2=(w2*m2).sum()/w2.sum(); ch2=(w2*(m2-mean2)**2).sum()
print(f"   excess as a multiple of the expectation: 2.22 +- 0.24, 3.67 +- 0.49, {meas/0.0048:.2f} +- {s1/0.0048:.2f}; mean {mean2:.2f}, chi2 {ch2:.1f} for 2 (chance {chi2.sf(ch2,2):.3f}); as an added dipole chance {chi2.sf(ch,2):.2f}")
out['added']=dict(mean=float(mean),err=float(1/np.sqrt(w.sum())),chi2=float(ch),chi2_mult=float(ch2),p_mult=float(chi2.sf(ch2,2)),p_add=float(chi2.sf(ch,2)))
# ---------------- (B)
T={'BGS, z 0.1-0.4':(335.37,20.57,60.88,118.41,104.10,47.99),'LRG, z 0.4-1.1':(158.21,52.79,18.54,262.61,34.12,31.09),'ELG, z 0.8-1.6':(337.03,-73.09,20.61,220.05,33.13,28.39),'QSO, z 0.8-2.1':(203.19,17.03,32.13,357.95,55.05,48.47)}
print("\n(B) DESI redshift dipole (Wu & Xia 2026): velocity of the Sun from each tracer")
print("   tracer            speed km/s        angle from CMB  along CMB (x CMB speed)  vs CMB vector: amp, dir -> chance   speed vs 2.2x   vs 3.67x")
out['desi']=[]
def chi_vec(V,t):
    l,b,d68,v,ep,em=t; n=unit(l,b); a=np.linalg.norm(V); ang=np.degrees(np.arccos(np.clip(V@n/a,-1,1))); s=ep if a>v else em
    return ((a-v)/s)**2+(ang/(d68/1.515))**2
for k,t in T.items():
    l,b,d68,v,ep,em=t; n=unit(l,b); ang=np.degrees(np.arccos(n@CMB)); za=(VC-v)/ep; zd=ang/(d68/1.515); c2=za**2+zd**2; p=chi2.sf(c2,3)
    z22=(2.2*VC-v)/ep; z367=(3.67*VC-v)/ep; along=v*np.cos(np.radians(ang))/VC
    print(f"   {k:16s} {v:6.1f} +{ep:5.1f} -{em:5.1f}   {ang:5.1f} deg      {along:+.2f}                  {za:+.1f}, {zd:.1f} sigma -> {p:.1e}        {z22:4.1f} sigma   {z367:4.1f} sigma")
    out['desi'].append(dict(name=k,l=l,b=b,d68=d68,v=v,ep=ep,em=em,angle=float(ang),along=float(along),z_amp=float(za),z_dir=float(zd),p_cmb=float(p),z22=float(z22),z367=float(z367)))
f=lambda V: sum(chi_vec(V,t) for t in T.values()); best=min((minimize(f,x0,method='Nelder-Mead',options=dict(xatol=1e-3,fatol=1e-6,maxiter=4000)) for x0 in [VC*CMB]+[t[3]*unit(t[0],t[1]) for t in T.values()]),key=lambda r:r.fun)
V=best.x; a=np.linalg.norm(V); lb=(np.degrees(np.arctan2(V[1],V[0]))%360,np.degrees(np.arcsin(V[2]/a)))
print(f"   one velocity for all four: best {a:.0f} km/s toward l = {lb[0]:.0f}, b = {lb[1]:+.0f} ({np.degrees(np.arccos(V@CMB/a)):.0f} deg from the CMB); chi2 = {best.fun:.1f} for 9, chance {chi2.sf(best.fun,9):.1e}")
print(f"   all four against the CMB vector: chi2 = {f(VC*CMB):.1f} for 12, chance {chi2.sf(f(VC*CMB),12):.1e};  against 2.2 x CMB: {f(2.2*VC*CMB):.0f};  3.67 x: {f(3.67*VC*CMB):.0f}")
q=T['QSO, z 0.8-2.1']; print(f"   quasars alone: speed / CMB speed = {q[3]/VC:.2f} +{q[4]/VC:.2f} -{q[5]/VC:.2f}")
out['desi_common']=dict(v=float(a),l=float(lb[0]),b=float(lb[1]),chi2=float(best.fun),p=float(chi2.sf(best.fun,9)),chi2_cmb=float(f(VC*CMB)),p_cmb=float(chi2.sf(f(VC*CMB),12)),chi2_22=float(f(2.2*VC*CMB)),chi2_367=float(f(3.67*VC*CMB)))
# ---------------- (C)
print("\n(C) the kinematic count dipole in redshift slices")
def make_n(a,b,zmed=1.47):
    g=lambda z0: quad(lambda z: z**a*np.exp(-(z/z0)**b),0,zmed)[0]/quad(lambda z: z**a*np.exp(-(z/z0)**b),0,8)[0]-0.5
    z0=brentq(g,0.2,5); nrm=quad(lambda z: z**a*np.exp(-(z/z0)**b),0,8)[0]; return (lambda z: z**a*np.exp(-(z/z0)**b)/nrm), z0
def slices(n,nb,base):
    zz=np.linspace(0,8,8001); cdf=np.cumsum(n(zz))*(zz[1]-zz[0]); cdf/=cdf[-1]; edges=np.interp(np.linspace(0.005,0.995,nb+1),cdf,zz); A=[]
    for z1,z2 in zip(edges[:-1],edges[1:]):
        Nb=quad(n,z1,z2)[0]; A.append(base+((1+z2)*n(z2)-(1+z1)*n(z1))/Nb)
    return edges,np.array(A)
# error on a dipole component along the CMB direction, with a Galactic cut: Fisher matrix for (mean, three dipole components)
def penalty(bcut):
    rng=np.random.default_rng(1); v=rng.standard_normal((400000,3)); v/=np.linalg.norm(v,axis=1)[:,None]; v=v[np.abs(v[:,2])>np.sin(np.radians(bcut))]
    f=np.hstack([np.ones((len(v),1)),v]); F=f.T@f/len(f); cov=np.linalg.inv(F)[1:,1:]; return float(CMB@cov@CMB), len(v)/400000
out['pen']={}
for bc in (0,30,40): var,fs=penalty(bc); out['pen'][bc]=[var,fs]; print(f"   Galactic cut {bc:2d} deg: sky kept {fs:.2f}; variance of the dipole along the CMB direction = {var:.2f} / N (3 / N for the whole sky)")
base_q=0.0043/BETA; print(f"   Quaia G < 20.5: expected whole-sample dipole 0.0043 = {base_q:.2f} beta")
out['fore']=[]
for (a_,b_) in ((2.0,1.6),(1.5,2.0),(3.0,1.2)):
    n,z0=make_n(a_,b_); e,A=slices(n,8,base_q)
    print(f"   n(z) ~ z^{a_} exp(-(z/{z0:.2f})^{b_}), median 1.47; eight equal slices, edges "+" ".join(f"{x:.2f}" for x in e))
    print("      kinematic dipole per slice, in units of beta: "+" ".join(f"{x:+.1f}" for x in A)+f"   (whole sample {base_q:.1f})")
    row=dict(a=a_,b=b_,edges=[float(x) for x in e],A=[float(x) for x in A],cases=[])
    for name,N,bc in (("Quaia G<20.5, 30-degree cut",1295502,30),("ten million sources, 30-degree cut",1.0e7,30),("a hundred million, 30-degree cut",1.0e8,30)):
        var,fs=out['pen'][bc]; Nuse=N*fs if N<2e6 else N; sb=np.sqrt(var/(Nuse/8)); lss=3e-4; sb=np.sqrt(sb**2+lss**2)
        # two readings of a whole-sample dipole 2.3 times the expectation: a faster Sun, or the CMB speed plus an intrinsic dipole
        Dv=2.3*BETA*A; Di=BETA*A+1.3*BETA*base_q; dchi=np.sum(((Dv-Di)/sb)**2)
        Fm=np.array([[np.sum((BETA*A)**2),np.sum(BETA*A)],[np.sum(BETA*A),8.0]])/sb**2; cov=np.linalg.inv(Fm)
        print(f"      {name:36s}: error per slice {1e2*sb:.2f}%; speed multiple to +-{np.sqrt(cov[0,0]):.2f}, intrinsic dipole to +-{1e2*np.sqrt(cov[1,1]):.2f}%; faster Sun against intrinsic dipole: {np.sqrt(dchi):.1f} sigma")
        row['cases'].append(dict(name=name,N=float(N),sb=float(sb),sk=float(np.sqrt(cov[0,0])),sd=float(np.sqrt(cov[1,1])),sep=float(np.sqrt(dchi))))
    out['fore'].append(row)
json.dump(out,open(os.path.join(HERE,'slices.json'),'w'),indent=1)
