"""Workbench Entry 45: stellar against dynamical mass in 71 quiescent galaxies at z = 0.8 to 2.2, under expanding and static distances.

A galaxy's velocity dispersion comes straight from its spectrum and needs no distance. Its radius (and so its dynamical mass,
k sigma^2 r / G) scales with the angular-size distance. Its luminosity (and so its stellar mass) scales with the luminosity distance squared.
Changing the distance law therefore multiplies dynamical/stellar mass by a known factor F(z) and leaves everything else alone.
Entry 12 did this for 13 galaxies. This is the homogeneous compilation of van de Sande et al. (2013, table 4), anchored on the nearby
value from the same paper (SDSS quiescent galaxies above 10^11 solar masses: stellar/dynamical = 0.59, scatter 0.12 dex)."""
import os, json
import numpy as np
from scipy import integrate
HERE=os.path.dirname(os.path.abspath(__file__))
SRC={0:"van de Sande 2013",1:"Bezanson 2013",2:"van Dokkum 2009",3:"Onodera 2012",4:"Cappellari 2009",5:"Newman 2010",6:"van der Wel 2008",7:"Toft 2012"}   # identified from names and redshifts
R=[l.split() for l in open(os.path.join(HERE,"vds2013_table4.txt")) if not l.startswith("#")]
ref=np.array([int(r[0]) for r in R]); ids=np.array([r[1] for r in R]); z,re,n,sig,kn,lMd,lMs=np.array([[float(x) for x in r[2:]] for r in R]).T
G=4.3009e-6                                                   # kpc (km/s)^2 / Msun
chk=np.log10(kn*sig**2*re/G)-lMd
print(f"{len(z)} rows. Check of k sigma^2 r/G against the listed dynamical mass: {np.sum(abs(chk)<=0.03)} rows within 0.03 dex; the other {np.sum(abs(chk)>0.03)} are all source {set(ref[abs(chk)>0.03])}, listed {abs(chk[abs(chk)>0.03]).min():.2f} to {abs(chk).max():.2f} dex lower. Listed values are used.")
bad=(lMs-lMd)>0.6; print("dropped, stellar mass 13 times the dynamical mass (a likely misprint in the source table):",[(ids[i],lMd[i],lMs[i]) for i in np.where(bad)[0]])
dup=(ref==7)&(ids=="19627"); a=np.where((ref==0)&(ids=="19627"))[0][0]; b=np.where(dup)[0][0]
print(f"one galaxy measured twice (UDS-19627): radius {re[a]} and {re[b]} kpc, dispersion {sig[a]:.0f} and {sig[b]:.0f} km/s, log dynamical mass {lMd[a]} and {lMd[b]}: a factor {10**(lMd[b]-lMd[a]):.1f} apart. The compilers' own measurement is kept.")
ok=~bad&~dup; z,re,n,sig,lMd,lMs,ref,ids=[v[ok] for v in (z,re,n,sig,lMd,lMs,ref,ids)]; N=len(z); print(f"{N} galaxies kept, redshift {z.min():.2f} to {z.max():.2f}, median {np.median(z):.2f}")
Om=0.3; E=lambda x: np.sqrt(Om*(1+x)**3+1-Om); DM=np.vectorize(lambda x: integrate.quad(lambda t:1/E(t),0,x)[0])   # units of c/H0
# name: (distance d, power q in d_L = d (1+z)^q).  Size distance is d itself in flat static space.
STATIC={"lerner":("Lerner: d = cz/H0",lambda x:x,0.5),"tired":("tired light: d = (c/H0) ln(1+z)",lambda x:np.log(1+x),0.5),
        "lerner_t":("d = cz/H0, time also stretched",lambda x:x,1.0),"tired_t":("ln(1+z), time also stretched",lambda x:np.log(1+x),1.0)}
def parts(x,d,q):                                              # size and luminosity relative to the published (expanding) values
    x=np.asarray(x,float); return d(x)*(1+x)/DM(x), (d(x)*(1+x)**q/((1+x)*DM(x)))**2
def F(key,x):
    if key=="expanding": return np.ones_like(np.asarray(x,float))
    s,l=parts(x,*STATIC[key][1:]); return s/l
LOCAL=-np.log10(0.59); zl=0.06
y0=lMd-lMs; rng=np.random.default_rng(45); out={"n":int(N),"local":float(LOCAL),"laws":{}}
bins=((0.8,1.0),(1.0,1.5),(1.5,2.3)); lx=np.log10(1+z)-np.log10(1+zl)
print(f"nearby anchor: dynamical/stellar = {10**LOCAL:.2f} (stars {100*0.59:.0f}% of the mass inside the half-light radius)")
print(f"scatter of the {N} about their mean: {y0.std(ddof=1):.2f} dex; standard error of the mean {y0.std(ddof=1)/N**0.5:.3f} dex")
for key in ("expanding",)+tuple(STATIC):
    nm="expanding (as published)" if key=="expanding" else "static, "+STATIC[key][0]
    y=y0+np.log10(F(key,z)); loc=LOCAL+np.log10(F(key,[zl]))[0]; res={"name":nm,"local":float(loc),"bins":[]}
    print(f"\n{nm}\n  nearby {10**loc:.2f}")
    for a_,b_ in bins:
        q=(z>=a_)&(z<b_); m=y[q].mean(); e=y[q].std(ddof=1)/q.sum()**0.5; res["bins"].append([float(z[q].mean()),float(m),float(e),int(q.sum())])
        print(f"  z {a_}-{b_} (n={q.sum():2d}, mean z {z[q].mean():.2f}): {10**m:5.2f}  = {10**(m-loc):.2f} x nearby (log {m-loc:+.2f} +- {e:.2f}); stars {100*10**(-m):.0f}% of the mass")
    ex=y-loc; m=ex.mean(); e=ex.std(ddof=1)/N**0.5
    s_a=float(np.sum(lx*ex)/np.sum(lx*lx)); X=np.c_[np.ones(N),lx]; s_f=float(np.linalg.lstsq(X,ex,rcond=None)[0][1]); ba=[];bf=[]
    for _ in range(4000):
        j=rng.integers(0,N,N); ba.append(np.sum(lx[j]*ex[j])/np.sum(lx[j]**2)); bf.append(np.linalg.lstsq(X[j],ex[j],rcond=None)[0][1])
    cut=lMs>=11.0; mc=ex[cut].mean(); ec=ex[cut].std(ddof=1)/cut.sum()**0.5
    lo=[]; 
    for r_ in sorted(set(ref)):
        k=ref!=r_; lo.append(ex[k].mean())
    srcs=[(SRC[r_],int((ref==r_).sum()),float(ex[ref==r_].mean())) for r_ in sorted(set(ref)) if (ref==r_).sum()>=5]
    over=int(np.sum(y<0))
    print(f"  all {N}: {10**m:.2f} x nearby (log {m:+.3f} +- {e:.3f}, {m/e:+.1f} standard errors); median ratio {10**np.median(y):.1f}, stars {100*10**(-np.median(y)):.0f}%")
    print(f"  power of (1+z): {s_a:+.2f} +- {np.std(ba):.2f} tied to the nearby value; {s_f:+.2f} +- {np.std(bf):.2f} inside the sample alone")
    print(f"  stellar mass above 10^11 only (n={cut.sum()}): {10**mc:.2f} x nearby (log {mc:+.3f} +- {ec:.3f}); leaving out one source at a time: {10**min(lo):.2f} to {10**max(lo):.2f}")
    print("  by source: "+"; ".join(f"{a_} ({b_}) {10**c_:.2f} x" for a_,b_,c_ in srcs)+f"; galaxies with more stars than dynamical mass: {over}")
    res.update(excess=[float(m),float(e)],power_tied=[s_a,float(np.std(ba))],power_free=[s_f,float(np.std(bf))],cut11=[float(mc),float(ec),int(cut.sum())],leave_one=[float(min(lo)),float(max(lo))],sources=srcs,over=over,median=float(np.median(y)))
    out["laws"][key]=res
print("\nwhat the static laws change, relative to the published values (size, luminosity):")
for key,(nm,d,q) in STATIC.items():
    s1,l1=parts([1.0],d,q); s2,l2=parts([2.0],d,q); sm,lm=parts(z,d,q)
    print(f"  {nm:34s} z=1: size x{s1[0]:.2f}, luminosity x{l1[0]:.2f} | z=2: size x{s2[0]:.2f}, luminosity x{l2[0]:.2f} | sample medians: radius {np.median(re):.1f} -> {np.median(re*sm):.1f} kpc, log stellar mass {np.median(lMs):.2f} -> {np.median(lMs+np.log10(lm)):.2f}")
    out["laws"][key]["z1"]=[float(s1[0]),float(l1[0])]; out["laws"][key]["z2"]=[float(s2[0]),float(l2[0])]
print(f"median dispersion {np.median(sig):.0f} km/s (all); {np.median(sig[z<1.0]):.0f} at 0.8-1.0; {np.median(sig[z>1.5]):.0f} beyond 1.5")
print(f"light-travel time on Lerner's scale (z/H0, H0 = 70): {0.83*13.97:.1f} Gyr at z = 0.83, {1.3*13.97:.0f} at 1.3, {2*13.97:.0f} at 2")
need=10**(-out['laws']['lerner']['excess'][0]); print(f"to bring Lerner's law to the nearby value at fixed luminosity, sizes must be {1/need:.1f} times smaller: a size distance of {need*np.median(z):.2f} c/H0 at the median redshift {np.median(z):.2f}, against {np.median(z):.2f} c/H0; the expanding value there is {DM(np.median(z))/(1+np.median(z)):.2f}")
out["points"]=[[float(a_),float(b_),int(c_),str(d_)] for a_,b_,c_,d_ in zip(z,y0,ref,ids)]; out["sources"]=SRC
zz=np.r_[0.03,0.06,np.arange(0.1,2.31,0.05)]; out["F"]={k:[[round(float(a_),3),round(float(np.log10(F(k,[a_])[0])),4)] for a_ in zz] for k in STATIC}
json.dump(out,open(os.path.join(HERE,"dynamical_mass_test.json"),"w"),indent=1)
