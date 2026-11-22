"""Workbench Entry 50: two new forecasts and the scoreboard.
1. How many galaxies a stack needs to find, or rule out, a 27 uK cooling of the CMB around spirals at redshift 0.05 to 0.3.
   Per-galaxy error = primary CMB left inside a disc of 0.2 Mpc after removing the largest scales, plus map noise.
   The CMB spectrum is an approximate table of the measured one (good to about 20%), not a calculation.
2. How many sources per redshift slice a count dipole needs to tell our motion alone from the published excess."""
import os, json, re
import numpy as np
from scipy import special, integrate
HERE=os.path.dirname(os.path.abspath(__file__)); out={}
# approximate D_l = l(l+1)C_l/2pi of the CMB temperature, uK^2 (read off published spectra)
LT=np.array([2,10,30,100,220,320,412,540,670,810,1000,1130,1300,1450,1700,2000,2500,3000,4000,6000]); DT=np.array([1000,900,1000,2500,5700,2400,1700,2600,1800,2500,1050,1240,800,800,400,250,100,45,15,3.0])
Dl=lambda l: np.exp(np.interp(np.log(l),np.log(LT),np.log(DT)))
Om=0.315; h=0.674; cH0=2997.92458/h; DA=lambda z: cH0*integrate.quad(lambda t:1/np.sqrt(Om*(1+t)**3+1-Om),0,z)[0]/(1+z)
def per_galaxy(z,noise,fwhm,lmin,rd=0.2):
    th=rd/DA(z); sb=np.radians(fwhm/60)/2.3548; l=np.logspace(np.log10(lmin),np.log10(2e4),4000); a=l*th
    W=2*special.j1(a)/a; cmb=np.trapezoid(Dl(l)*np.exp(-l*l*sb*sb)*W*W,np.log(l)); area=np.pi*(np.degrees(th)*60)**2+2*np.pi*(fwhm/2.3548)**2
    return float(np.sqrt(cmb+noise**2/area)), float(np.degrees(th)*60), float(np.sqrt(cmb)), float(noise/np.sqrt(area))
print("1. cooling stack: error on the mean temperature inside 0.2 Mpc of one galaxy")
out["stack"]=[]
for nm,noise,fwhm in (("ACT-like, 15 uK-arcmin, 1.4 arcmin beam",15,1.4),("Simons-Observatory-like, 6 uK-arcmin, 1.4 arcmin beam",6,1.4),("Planck-like, 45 uK-arcmin, 5 arcmin beam",45,5.0)):
    for z in (0.05,0.1,0.2,0.3):
        th=0.2/DA(z); lmin=max(300,int(np.pi/(2.0/DA(z))/1))          # remove scales larger than about 2 Mpc at the galaxy
        s,tha,sc,sn=per_galaxy(z,noise,fwhm,lmin); keep=0.6           # share of a 27 uK profile fading over 1-2 Mpc that survives the filter (assumed)
        dil=np.pi*tha**2/(np.pi*tha**2+2*np.pi*(fwhm/2.3548)**2)       # a disc smaller than the beam is washed out
        n5=(5*s/(27*keep*dil))**2; n3=(3*s/(5*dil))**2
        print(f"   {nm:52s} z = {z}: disc radius {tha:4.1f} arcmin, scales below multipole {lmin:5d} removed: {s:4.1f} uK per galaxy (CMB {sc:4.1f}, noise {sn:4.1f}); to find 27 uK at 5 sigma: {n5:5.0f} galaxies; to rule out 5 uK at 3 sigma: {n3:6.0f}")
        out["stack"].append(dict(inst=nm,z=z,theta=tha,lmin=lmin,per=s,cmb=sc,noise=sn,dil=float(dil),n_find=float(n5),n_rule=float(n3)))
print("\n2. count dipole in redshift slices")
out["dipole"]=[]
for nm,dk,dx in (("quasars: expected 0.007, published 0.0155",0.007,0.0155),("radio: expected 0.0043, published 3.67 times",0.0043,0.0158)):
    for fsky in (1.0,0.5,0.34):
        s3=(dx-dk)/3; N=3/s3**2/fsky
        print(f"   {nm:44s} sky share {fsky:.2f}: to tell the two apart at 3 sigma in one slice needs an error of {s3:.4f}, or {N:,.0f} sources in the slice")
        out["dipole"].append(dict(sample=nm,fsky=fsky,err=s3,N=float(N)))
print(f"   check against a published result: 1.18 million quasars over about a fifth of the sky give {np.sqrt(3/1.18e6/0.2):.4f}; the paper's 444 +- 204 km/s is {204/370*0.007:.4f} on an expected 0.007")
# 3. scoreboard from the Gap Map
g=open(os.path.join(HERE,"..","gapmap.html")).read(); sc={}
for m in re.finditer(r"\{id:'([a-z0-9_]+)',theme:'([a-z0-9]+)'(?:,hub:true)?,label:'((?:[^'\\]|\\.)*)',status:'([a-z]+)',aff:'([A-Z])'",g):
    sc.setdefault(m.group(5),{}).setdefault(m.group(4),[]).append(m.group(3))
names={"P":"plasma-cosmology claims","S":"standard-model results, with the plasma physics inside it","L":"tensions in the standard model","A":"changed-gravity alternatives"}
print("\n3. scoreboard, from the Gap Map")
for a in "PSLA":
    print("  ",names[a]+":",{k:len(v) for k,v in sc[a].items()})
    for k,v in sc[a].items(): print("       ",k,"-",'; '.join(v))
out["score"]={a:{k:v for k,v in sc[a].items()} for a in sc}
json.dump(out,open(os.path.join(HERE,"forecasts2.json"),"w"),indent=1,ensure_ascii=False)
