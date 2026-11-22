"""Workbench Entry 54: the Milky Way's vertical magnetic flux, and what could supply it.
Entry 18 showed that the halo's field rings are the Galaxy's own vertical field wound up by rotation, and left the source of that
vertical field open.  A field frozen into the gas keeps its flux.  So: how much flux threads the disc, and how much would a field
left over from the early universe have threaded through the gas that became the Galaxy?"""
import os, json
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); out={}
KPC=3.0857e21; MPC=3.0857e24; MSUN=1.98847e33; G=6.6743e-8; PC=3.0857e18
# ---- 1. the flux through the disc, four ways
est={}
est['local caps: 0.30 uG over 2-15 kpc (Entry 18)']=0.30e-6*np.pi*(15**2-2**2)*KPC**2
est['Unger & Farrar base model: 0.98 uG out to 7.29 kpc']=0.98e-6*np.pi*(7.29*KPC)**2
est['Unger & Farrar exponential: 5.8 uG, scale 2.5 kpc']=2*np.pi*5.8e-6*(2.5*KPC)**2
est['Jansson & Farrar 2012: 4.6 uG, scale 2.9 kpc, tilted 49 deg']=2*np.pi*4.6e-6*(2.9*KPC)**2*np.sin(np.radians(49))
print("1. vertical flux through the disc (maxwell = gauss cm2)")
for k,v in est.items(): print(f"   {k:62s} {v:.2e}")
lo,hi=min(est.values()),max(est.values()); print(f"   range {lo:.1e} to {hi:.1e}"); out['flux']={k:float(v) for k,v in est.items()}
# ---- 2. does the vertical field follow the baryons?
Rd=2.6; R0=8.2
print("\n2. profile: a frozen-in field follows the surface density of baryons, an exponential of scale 2.6 +- 0.5 kpc")
for Phi in (lo,hi):
    Bc=Phi/(2*np.pi*(Rd*KPC)**2); print(f"   flux {Phi:.1e}: centre {Bc*1e6:.1f} uG, at the Sun {Bc*np.exp(-R0/Rd)*1e6:.2f} uG   (scale 2.1: {Phi/(2*np.pi*(2.1*KPC)**2)*np.exp(-R0/2.1)*1e6:.2f}; scale 3.1: {Phi/(2*np.pi*(3.1*KPC)**2)*np.exp(-R0/3.1)*1e6:.2f})")
print(f"   fitted exponentials: scale 2.5 +- 0.1 kpc (Unger & Farrar), 2.9 +- 0.1 (Jansson & Farrar); at the Sun {5.8*np.exp(-R0/2.5):.2f} and {4.6*np.sin(np.radians(49))*np.exp(-R0/2.9):.2f} uG; measured south cap 0.30 +- 0.03")
Mb=(54.3e9+1.1e10+1.2e9)*MSUN; lam_crit=1/(2*np.pi*np.sqrt(G))
print(f"   baryons {Mb/MSUN:.2e} Msun; mass-to-flux ratio {Mb/hi:.1e} to {Mb/lo:.1e} g/Mx = {Mb/hi/lam_crit:.0f} to {Mb/lo/lam_crit:.0f} times the value at which a field can resist gravity")
out['profile']=dict(Rd=Rd,sun=[float(lo/(2*np.pi*(Rd*KPC)**2)*np.exp(-R0/Rd)*1e6),float(hi/(2*np.pi*(Rd*KPC)**2)*np.exp(-R0/Rd)*1e6)],mtf=[float(Mb/hi/lam_crit),float(Mb/lo/lam_crit)])
# ---- 3. the patch of the early universe that became the Galaxy
h=0.674; rho_m=0.315*2.775e11*h*h; fb=0.157
print("\n3. the patch that became the Galaxy (mean matter density today %.2e Msun/Mpc3)"%rho_m); out['patch']=[]
for M in (1.0e12,1.3e12,1.6e12):
    RL=(3*M/(4*np.pi*rho_m))**(1/3); A=np.pi*(RL*MPC)**2; cap=(Mb/MSUN)/(fb*M); f_sph=cap**(2/3); f_cyl=1-(1-cap)**(2/3)
    print(f"   mass {M:.1e}: radius {RL:.2f} Mpc (comoving), area {A:.2e} cm2; the disc holds {100*cap:.0f}% of its baryons, which filled {100*f_sph:.0f}% of the area if they came from the centre, {100*f_cyl:.0f}% if from a column along the axis")
    row=[M,RL,A,cap,f_sph,f_cyl]
    for nm,f in (('all of the patch',1.0),('central sphere',f_sph),('axial column',f_cyl)):
        print(f"        field needed along the spin axis, {nm:16s}: {lo/(f*A)*1e12:5.1f} to {hi/(f*A)*1e12:5.1f} pG"); row+=[lo/(f*A)*1e12,hi/(f*A)*1e12]
    out['patch'].append([float(x) for x in row])
# ---- 4. what each candidate supplies (mass 1.3e12, central-sphere capture; full capture in brackets)
M=1.3e12; RL=(3*M/(4*np.pi*rho_m))**(1/3); A=np.pi*(RL*MPC)**2; cap=(Mb/MSUN)/(fb*M); fs=cap**(2/3); Reff=RL*np.sqrt(fs)
L=json.load(open(os.path.join(HERE,'loop_flux.json'))); c=L['tangled']['c']
print("\n4. supply, against %.1e to %.1e needed"%(lo,hi)); out['supply']=[]
def show(name,phi,phi_full=None):
    s=f"   {name:66s} {phi:.1e} Mx = {phi/hi:.1e} to {phi/lo:.1e} of the need"
    if phi_full: s+=f"   (all of the patch: {phi_full:.1e} = {phi_full/hi:.2f} to {phi_full/lo:.2f})"
    print(s); out['supply'].append([name,float(phi),float(phi/hi),float(phi/lo),float(phi_full) if phi_full else None])
for B in (5e-12,10e-12): show(f"tangled field, {B*1e12:.0f} pG, coherence 1 kpc (causal, eases the Hubble tension)",c*B*KPC**1.5*np.sqrt(Reff*MPC),c*B*KPC**1.5*np.sqrt(RL*MPC))
for B in (5e-12,10e-12,20e-12,45e-12): show(f"equal power on all scales, {B*1e12:.0f} pG, typical galaxy (rms)",0.4*B*fs*A,0.4*B*A)
for B in (20e-12,45e-12): show(f"the same, field along the spin axis (best case), {B*1e12:.0f} pG",B/np.sqrt(3)*np.sqrt(0.5)*np.sqrt(3)*fs*A,B*np.sqrt(0.5)*A)
# battery at the central black hole: Poynting-Robertson drag on electrons drives a ring voltage; d(flux)/dt = L sigma_T v / (2 r c e)
Mbh=4.3e6*MSUN; rg=G*Mbh/9e20; Ledd=1.26e38*4.3e6; sT=6.652e-25; e=4.803e-10; c_=2.998e10
dphi=Ledd*sT*(0.4*c_)/(2*6*rg*c_*e); t=1e8*3.156e7; show("battery at the central black hole, full power for 100 million years",dphi*t)
print(f"      (black hole 4.3e6 Msun, Eddington power {Ledd:.1e} erg/s, ring at 6 gravitational radii: {dphi:.1e} Mx/s)")
# ---- 5. a neighbour: would Andromeda's vertical field run with or against its spin?
def radec(a,d): a,d=np.radians(a),np.radians(d); return np.array([np.cos(d)*np.cos(a),np.cos(d)*np.sin(a),np.sin(d)])
def basis(a,d):
    a,d=np.radians(a),np.radians(d); return np.array([-np.sin(a),np.cos(a),0]), np.array([-np.sin(d)*np.cos(a),-np.sin(d)*np.sin(a),np.cos(d)])
ngp=radec(192.85948,27.12825)                                 # Galactic north pole; the Milky Way's spin points the other way
los=radec(10.6847,41.2687); E,Nn=basis(10.6847,41.2687); pa=np.radians(38.0); inc=np.radians(77.0)
m=np.cos(pa)*Nn+np.sin(pa)*E                                   # north-east end of the major axis, which is receding
p=np.cos(pa-np.pi/2)*Nn+np.sin(pa-np.pi/2)*E                   # north-west minor axis, the near side
u=p*np.cos(inc)-los*np.sin(inc)                                # in the disc plane, toward the near side
s=np.cross(m,-u); s/=np.linalg.norm(s)                         # spin: position x velocity at the receding end
ang=np.degrees(np.arccos(np.dot(s,ngp)))
print(f"\n5. Andromeda: spin axis {ang:.0f} degrees from the Milky Way's north pole (cosine {np.dot(s,ngp):+.2f}); the Milky Way's own spin is at 180 degrees")
print(f"   the Milky Way's vertical field points to its north pole, against its spin.  A shared primordial direction within {abs(90-ang):.0f} degrees of that pole would thread Andromeda {'along' if np.dot(s,ngp)>0 else 'against'} its spin")
for B0 in (45,60,80): 
    for need in (28,39): 
        if need<=B0: print(f"   limit {B0} pG, need {need} pG along the Milky Way's axis: the direction lies within {np.degrees(np.arccos(need/B0)):.0f} degrees of the pole")
out['m31']=dict(angle=float(ang),cos=float(np.dot(s,ngp)),spin_radec=[float(np.degrees(np.arctan2(s[1],s[0]))%360),float(np.degrees(np.arcsin(s[2])))])
json.dump(out,open(os.path.join(HERE,'flux_budget.json'),'w'),indent=1)
# share of allowed directions that thread Andromeda against its spin, for cones of different half-angle about the Milky Way's north pole
rng=np.random.default_rng(5); v=rng.standard_normal((400000,3)); v/=np.linalg.norm(v,axis=1)[:,None]; cz=v@ngp; cs=v@s; res=[]
for half in (30,40,52,62,70,90):
    q=cz>=np.cos(np.radians(half)); res.append([half,float(np.mean(cs[q]<0))]); print(f"   primordial direction within {half} degrees of the pole: threads Andromeda against its spin in {100*np.mean(cs[q]<0):.0f}% of cases")
out['m31']['against']=res; json.dump(out,open(os.path.join(HERE,'flux_budget.json'),'w'),indent=1)
