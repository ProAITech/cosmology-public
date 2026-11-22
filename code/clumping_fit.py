"""Workbench Entry 36, Part 4: the free factor on the sound horizon replaced by the clumping model of Entry 13.
Entry 13 ran our recombination code with the early plasma clumped by primordial magnetic fields (clumping strength b) and recorded, at
Planck's densities, the sound horizon at last scattering r_s, the redshift of last scattering z*, and the sound horizon at the ruler r_d
(clump_track.json). Here those three are applied as ratios to the standard values, so that b replaces the single factor alpha of Part 3."""
import sys, os, json
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
import numpy as np
from scipy import optimize, integrate
import cosmo as c
from dark_energy_fit import sig
T=json.load(open(os.path.join(HERE,"clump_track.json"))); B=np.array([t["b"] for t in T]); RS=np.array([t["rs"] for t in T]); RD=np.array([t["rd"] for t in T]); ZS=np.array([t["zstar"] for t in T])
f_s=lambda b: np.interp(b,B,RS/RS[0]); f_d=lambda b: np.interp(b,B,RD/RD[0]); dz=lambda b: np.interp(b,B,ZS-ZS[0])
print("clumping b -> sound horizon at the ruler: "+", ".join(f"b={b_:.1f}: {100*(f_d(b_)-1):+.2f}%" for b_ in B[1:]))
class ClumpModel(c.Model):
    def __init__(s,ob,obc,h,w0=-1.0,wa=0.0,b=0.0):
        super().__init__(ob,obc,h,w0,wa); s.b=b
    def rd(s): return super().rd()*f_d(s.b)
    def theta_star(s):
        th=super().theta_star(); zs=s.zstar(); a=1/(1+zs)                      # correct the distance for the later/earlier last scattering
        dDM=dz(s.b)*c.C_KMS/100/np.sqrt(s.obc/a**3+(c.OG_H2*(1+0.2271*c.NEFF))/a**4)
        rs_over_dm=th; dm=1.0                                                 # theta = rs/DM ; DM -> DM + dDM
        DM=np.interp(min(zs,c.ZG[-1]),c.ZG,s.dc); return th*f_s(s.b)/(1+dDM/DM)
def chi(m): return c.chi_bao(m)+c.chi_cmb(m)
def fit(kind,free_b=True,extra=None):
    def unpack(x):
        ob,obc,h=x[0]*0.02223,x[1]*0.142,x[2]*0.68; i=3; w0,wa=-1.0,0.0
        if kind=='w0wa': w0,wa=x[3],x[4]; i=5
        b=x[i] if free_b else 0.0; return ob,obc,h,w0,wa,b
    def f(x):
        ob,obc,h,w0,wa,b=unpack(x)
        if not (0.01<ob<0.04 and 0.05<obc<0.3 and 0.4<h<1.0 and -3<w0<1 and -3<wa<2 and w0+wa<0 and 0<=b<=0.8): return 1e9
        m=ClumpModel(ob,obc,h,w0,wa,b)
        if m.ode<=0: return 1e9
        return chi(m)+(extra(m) if extra else 0.0)
    best=None
    for s in ((1.0,1.0),(0.98,1.02),(1.02,0.99)):
        for ws in ([(-1.0,0.0),(-0.5,-1.5),(-0.8,-0.6)] if kind=='w0wa' else [None]):
            for b0 in ([0.05,0.3] if free_b else [None]):
                x0=[1.0,s[0],s[1]]+(list(ws) if ws else [])+([b0] if free_b else [])
                r=optimize.minimize(f,x0,method='Nelder-Mead',options=dict(xatol=1e-5,fatol=1e-5,maxiter=6000,maxfev=6000)); r=optimize.minimize(f,r.x,method='Nelder-Mead',options=dict(xatol=1e-6,fatol=1e-6,maxiter=6000,maxfev=6000))
                if best is None or r.fun<best.fun: best=r
    ob,obc,h,w0,wa,b=unpack(best.x); return dict(chi=float(best.fun),ob=ob,obc=obc,h=h,w0=w0,wa=wa,b=b,Om=(obc+c.OMNU0)/h**2)
out={}
# check of the mapping against Entry 13's own track: Planck densities, b = 0.3, H0 = 69.65 should return Planck's acoustic angle
m=ClumpModel(0.02237,0.14307-c.OMNU0,0.6965,b=0.3); print(f"check: Entry 13's b = 0.3 solution gives 100 theta* = {100*c._KAPPA*m.theta_star():.4f} (Planck 1.0411); H0 x r_d = {m.h*m.rd():.2f} (track 101.68)")
L0=fit('lcdm',free_b=False); Lb=fit('lcdm'); W0=fit('w0wa',free_b=False); Wb=fit('w0wa')
print(f"constant dark energy, no clumping:   chi2 {L0['chi']:.2f}, H0 {100*L0['h']:.2f}, matter share {L0['Om']:.4f}")
print(f"constant dark energy, clumping free: chi2 {Lb['chi']:.2f}, b = {Lb['b']:.3f}, H0 {100*Lb['h']:.2f}, matter share {Lb['Om']:.4f}; gain {L0['chi']-Lb['chi']:.2f} for one number = {sig(L0['chi']-Lb['chi'],1):.1f} sigma")
print(f"changing dark energy, no clumping:   chi2 {W0['chi']:.2f}; gain over constant {L0['chi']-W0['chi']:.2f} = {sig(L0['chi']-W0['chi']):.1f} sigma")
print(f"changing dark energy, clumping free: chi2 {Wb['chi']:.2f}, b = {Wb['b']:.3f}, w0 {Wb['w0']:+.2f}, wa {Wb['wa']:+.2f}; gain over constant with clumping {Lb['chi']-Wb['chi']:.2f} = {sig(Lb['chi']-Wb['chi']):.1f} sigma")
# profile of b
def prof(b):
    def f(x):
        ob,obc,h=x[0]*0.02223,x[1]*0.142,x[2]*0.68
        if not (0.01<ob<0.04 and 0.05<obc<0.3 and 0.4<h<1.0): return 1e9
        m=ClumpModel(ob,obc,h,-1,0,b); return chi(m) if m.ode>0 else 1e9
    r=optimize.minimize(f,[Lb['ob']/0.02223,Lb['obc']/0.142,Lb['h']/0.68],method='Nelder-Mead',options=dict(xatol=1e-6,fatol=1e-6,maxiter=4000)); return r.fun,r.x[2]*0.68
bb=np.linspace(0,0.8,17); pr=[prof(b) for b in bb]
print("profile: "+", ".join(f"b={b:.2f}: +{p[0]-Lb['chi']:.2f} (H0 {100*p[1]:.1f})" for b,p in zip(bb,pr)))
g=lambda b: prof(b)[0]-Lb['chi']-1
lo=optimize.brentq(g,0.0,Lb['b']) if g(0.0)>0 else 0.0; hi=optimize.brentq(g,Lb['b'],0.8)
g4=lambda b: prof(b)[0]-Lb['chi']-4; hi2=optimize.brentq(g4,Lb['b'],0.8) if g4(0.8)>0 else 0.8
print(f"b = {Lb['b']:.2f} +{hi-Lb['b']:.2f} -{Lb['b']-lo:.2f}; below {hi2:.2f} at 2 sigma; b = 0 is worse by {prof(0.0)[0]-Lb['chi']:.2f}")
out=dict(L0=L0,Lb=Lb,W0=W0,Wb=Wb,b=[Lb['b'],Lb['b']-lo,hi-Lb['b'],hi2],profile=[[float(b),float(p[0]-Lb['chi']),float(100*p[1])] for b,p in zip(bb,pr)],track=[[float(b_),float(100*(f_d(b_)-1))] for b_ in B])
json.dump(out,open(os.path.join(HERE,"clumping_fit.json"),"w"),indent=1)
