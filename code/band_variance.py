"""Workbench Entry 66: how much should the density of a volume the size of the Giant Arc's band differ from the mean,
in the standard model, and how does that compare with the counting noise of the absorber catalogue?
Band: 1000 x 200 x 338 Mpc (and 1000 x 300 x 338) at redshift 0.8.  Linear theory, Eisenstein-Hu transfer function without wiggles,
sigma8 = 0.81, matter 0.31, h = 0.68, n_s = 0.965."""
import json, os
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
trap=getattr(np,'trapezoid',None) or getattr(np,'trapz')
h,Om,Ob,ns,s8=0.68,0.31,0.049,0.965,0.81
def T_eh(k):                      # k in 1/Mpc; Eisenstein & Hu 1998, zero-baryon-wiggle form
    th=2.7255/2.7; omh2=Om*h*h; obh2=Ob*h*h; fb=Ob/Om
    s=44.5*np.log(9.83/omh2)/np.sqrt(1+10*obh2**0.75)
    a=1-0.328*np.log(431*omh2)*fb+0.38*np.log(22.3*omh2)*fb**2
    G=Om*h*(a+(1-a)/(1+(0.43*k*s)**4)); q=(k/h)*th*th/G       # k/h is the wavenumber in h/Mpc
    L0=np.log(2*np.e+1.8*q); C0=14.2+731/(1+62.5*q); return L0/(L0+C0*q*q)
k=np.logspace(-5,1.5,6000); P0=k**ns*T_eh(k)**2
def sig_sphere(R,P):
    x=k*R; W=3*(np.sin(x)-x*np.cos(x))/x**3; return np.sqrt(trap(k**3*P*W**2/(2*np.pi**2),np.log(k)))
A=(s8/sig_sphere(8/h,P0))**2; P=A*P0
def growth(z):
    E=lambda a:np.sqrt(Om/a**3+1-Om)
    def g(zz):
        a=np.linspace(1e-4,1/(1+zz),4000); return E(a[-1])*trap(1/(a*E(a))**3,a)
    return g(z)/g(0)
def sig_box(Lx,Ly,Lz):
    mu,wm=np.polynomial.legendre.leggauss(96); ph=np.linspace(0,np.pi/2,97); ph=(ph[1:]+ph[:-1])/2          # one octant is enough by symmetry
    kk=np.logspace(-5,0.3,1500); out=np.zeros_like(kk)
    for i,kv in enumerate(kk):
        m=(mu[:,None]+1)/2; st=np.sqrt(1-m*m); kx=kv*st*np.cos(ph)[None,:]; ky=kv*st*np.sin(ph)[None,:]; kz=kv*m*np.ones_like(kx)
        W=np.sinc(kx*Lx/2/np.pi)*np.sinc(ky*Ly/2/np.pi)*np.sinc(kz*Lz/2/np.pi); out[i]=(W**2*wm[:,None]/2).sum()/len(ph)
    Pk=np.interp(kk,k,P); return np.sqrt(trap(kk**3*Pk*out/(2*np.pi**2),np.log(kk)))
D=growth(0.8); out={"growth_z0.8":float(D),"sigma8_check":float(sig_sphere(8/h,P))}
print(f"growth to z = 0.8: {D:.3f}; sigma in 8/h spheres today: {sig_sphere(8/h,P):.3f}")
for name,(Lx,Ly,Lz),n in (("band 1000 x 200 x 338",(1000,200,338.2),26.8),("band 1000 x 300 x 338",(1000,300,338.2),40.1)):
    s0=sig_box(Lx,Ly,Lz); s=s0*D; shot=1/np.sqrt(n)
    print(f"{name} Mpc: matter differs from the mean by {100*s0:.1f}% today, {100*s:.1f}% at z = 0.8; with a bias of 1.5, {150*s:.1f}%. Counting noise for {n} expected absorbers: {100*shot:.0f}%, {shot/(1.5*s):.0f} times larger")
    out[name]=dict(sigma_today=float(s0),sigma_z=float(s),biased=float(1.5*s),shot=float(shot),ratio=float(shot/(1.5*s)))
R=(3*1000*200*338.2/(4*np.pi))**(1/3); print(f"sphere of equal volume: radius {R:.0f} Mpc, sigma at z = 0.8: {100*sig_sphere(R,P)*D:.1f}%")
for mult in (7.6,):
    n=26.8*mult; print(f"with {mult} times the absorbers (DESI first release against this catalogue, 270,529 / 35,629): counting noise {100/np.sqrt(n):.1f}%; an excess of 34% would be {0.34*np.sqrt(n):.1f} errors")
    out["denser"]=dict(mult=mult,shot=float(1/np.sqrt(n)),n_sigma_34=float(0.34*np.sqrt(n)))
json.dump(out,open(os.path.join(HERE,"band_variance.json"),"w"),indent=1)
