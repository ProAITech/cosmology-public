"""A small Big Bang nucleosynthesis network (Workbench Entry 53).
Eight species (n, p, d, t, 3He, 4He, 7Li, 7Be), the neutron-proton weak rates, and the eleven reactions that carry the light elements.
Reaction rates are the fits in Kawano's public code (Smith, Kawano & Malaney 1993), read from a public copy of that code.
Weak rates are the Born integrals normalised to the neutron lifetime.  The photon-electron plasma is followed through pair annihilation
with entropy conserved; neutrinos are decoupled.  The purpose is to measure how the yields RESPOND to changes in the plasma, so
agreement with full codes to a few percent is enough."""
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.interpolate import CubicSpline
ME=0.51099895; QNP=1.29333236; MEV9=11.6045; TAUN=878.4; GN=6.6743e-8; RHOG=8.4180; NGAM=2.02868e28; NAV=6.02214076e23
q=QNP/ME
# ---- electron-positron plasma: energy and pressure relative to the photons, as functions of x = m_e/T
def _re(x): return 30/np.pi**4*quad(lambda y: y*y*np.sqrt(y*y+x*x)*np.exp(-np.sqrt(y*y+x*x))/(1+np.exp(-np.sqrt(y*y+x*x))),0,x+60,limit=200)[0]
def _pe(x): return 30/np.pi**4*quad(lambda y: y**4/(3*np.sqrt(y*y+x*x))*np.exp(-np.sqrt(y*y+x*x))/(1+np.exp(-np.sqrt(y*y+x*x))),0,x+60,limit=200)[0]
_lt=np.linspace(np.log(2e-3),np.log(3e3),260); _T=np.exp(_lt); _x=ME*MEV9/_T
_rev=np.array([_re(x) if x<200 else 0.0 for x in _x]); _pev=np.array([_pe(x) if x<200 else 0.0 for x in _x])
_S=1+0.75*(_rev+_pev)                                   # entropy of photons + pairs over entropy of photons
sp_re=CubicSpline(_lt,_rev); sp_lnS=CubicSpline(_lt,np.log(_S))
def plasma(T9):
    lt=np.log(T9); S=np.exp(sp_lnS(lt)); return float(sp_re(lt)), float(S), float(sp_lnS(lt,1))   # rho_e/rho_gamma, S, dlnS/dlnT
# ---- weak rates
lam0=quad(lambda e: e*np.sqrt(e*e-1)*(q-e)**2,1,q)[0]; KW=1/(TAUN*lam0)
def _fd(u): return 1/(1+np.exp(np.clip(u,-700,700)))
def weak(T9,Tnu9):
    z=ME*MEV9/T9; zn=ME*MEV9/Tnu9; ph=lambda e: e*np.sqrt(max(e*e-1,0.0)); w=1/min(z,zn)
    def tailint(fn,lo):                                                             # integral from lo upward, in pieces a few T wide
        edges=lo+w*np.array([0,2,5,9,14,20,28,40,60]); return sum(quad(fn,edges[i],edges[i+1],limit=100)[0] for i in range(8))
    a1=tailint(lambda e: ph(e)*(e-q)**2*_fd((e-q)*zn)*_fd(-e*z),q)                   # n + nu -> p + e-
    a2=tailint(lambda e: ph(e)*(e+q)**2*_fd(e*z)*_fd(-(e+q)*zn),1)                   # n + e+ -> p + nubar
    a3=quad(lambda e: ph(e)*(q-e)**2*_fd(-e*z)*_fd(-(q-e)*zn),1,q,limit=200)[0]      # n -> p + e- + nubar
    b1=tailint(lambda e: ph(e)*(e-q)**2*_fd(e*z)*_fd(-(e-q)*zn),q)                   # p + e- -> n + nu
    b2=tailint(lambda e: ph(e)*(e+q)**2*_fd((e+q)*zn)*_fd(-e*z),1)                   # p + nubar -> n + e+
    b3=quad(lambda e: ph(e)*(q-e)**2*_fd(e*z)*_fd((q-e)*zn),1,q,limit=200)[0]        # p + e- + nubar -> n
    return KW*(a1+a2+a3), KW*(b1+b2+b3)
_Tw=np.exp(np.linspace(np.log(5e-3),np.log(2e2),160)); _wr=[]
for T in _Tw:
    S=np.exp(sp_lnS(np.log(T))); _wr.append(weak(T,T*(4/11*S)**(1/3)))
_wr=np.array(_wr); sp_np=CubicSpline(np.log(_Tw),np.log(_wr[:,0])); sp_pn=CubicSpline(np.log(_Tw),np.log(np.maximum(_wr[:,1],1e-300)))
# ---- thermonuclear rates, N_A <sigma v> in cm3 s-1 mol-1 (Kawano's code, subroutine RATE2)
def rates(t9):
    ex=np.exp; t913=t9**(1/3); t923=t913*t913; t943=t923*t923; t953=t9*t923; t912=np.sqrt(t9); t932=t9*t912; m23=1/t923; m32=1/t932
    t9a=t9/(1+13.076*t9); t9d=t9/(1+0.759*t9); t9e=t9/(1+0.1378*t9); t9f=t9/(1+0.1071*t9); f={}
    f['npdg']=4.742e4*(1-.8504*t912+.4895*t9-.09623*t932+8.471e-3*t9*t9-2.80e-4*t9*t932)
    f['hent']=7.21e8*(1-.508*t912+.228*t9)
    f['benp']=2.675e9*(1-.560*t912+.179*t9-.0283*t932+2.214e-3*t9*t9-6.851e-5*t9*t932)+9.391e8*t9a**1.5*m32+4.467e7*m32*ex(-0.07486/t9)
    f['dpg']=2.65e3*m23*ex(-3.720/t913)*(1+.112*t913+1.99*t923+1.56*t9+.162*t943+.324*t953)
    f['lipa']=(1.096e9*m23*ex(-8.472/t913)-4.830e8*t9d**(5/6)*m32*ex(-8.472/t9d**(1/3))+1.06e10*m32*ex(-30.442/t9)
               +1.56e5*m23*ex((-8.472/t913)-(t9/1.696)**2)*(1+.049*t913-2.498*t923+.860*t9+3.518*t943+3.08*t953)+1.55e6*m32*ex(-4.478/t9))
    f['tag']=3.032e5*m23*ex(-8.090/t913)*(1+.0516*t913+.0229*t923+8.28e-3*t9-3.28e-4*t943-3.01e-4*t953)+5.109e5*t9e**(5/6)*m32*ex(-8.068/t9e**(1/3))
    f['hag']=4.817e6*m23*ex(-14.964/t913)*(1+.0325*t913-1.04e-3*t923-2.37e-4*t9-8.11e-5*t943-4.69e-5*t953)+5.938e6*t9f**(5/6)*m32*ex(-12.859/t9f**(1/3))
    f['ddn']=3.95e8*m23*ex(-4.259/t913)*(1+.098*t913+.765*t923+.525*t9+9.61e-3*t943+.0167*t953)
    f['ddp']=4.17e8*m23*ex(-4.258/t913)*(1+.098*t913+.518*t923+.355*t9-.010*t943-.018*t953)
    f['tdn']=1.063e11*m23*ex(-4.559/t913-(t9/.0754)**2)*(1+.092*t913-.375*t923-.242*t9+33.82*t943+55.42*t953)+8.047e8*m23*ex(-0.4857/t9)
    f['hdp']=5.021e10*m23*ex(-7.144/t913-(t9/.270)**2)*(1+.058*t913+.603*t923+.245*t9+6.97*t943+7.19*t953)+5.212e8/t912*ex(-1.762/t9)
    return f
CHARGED=('dpg','ddn','ddp','tdn','hdp','hag','tag','lipa')
N,P,D,T,H,A,L,B=range(8); AMASS=np.array([1,1,2,3,3,4,7,7.])
def run(eta10=6.10, nnu=3.0, tail=None, fB=0.0, scale=None, T9i=100.0, T9f=0.01, tau=TAUN, dense=False):
    """tail: function (name, T9) -> factor multiplying the forward rate of a charged-particle reaction (non-thermal tail).
       fB: extra energy density scaling as radiation (a magnetic field), as a fraction of the standard total at early times.
       scale: dict of constant factors on forward rates."""
    eta=eta10*1e-10; scale=scale or {}; rhoB0=fB/(1-fB)*RHOG*(1+7/4+nnu*7/8) if fB>0 else 0.0   # per T_nu^4, set where T_nu = T
    def rhs(t,y):
        T9=y[0]; Y=np.maximum(y[1:],0.0); re,S,dS=plasma(T9); Tn=T9*(4/11*S)**(1/3)
        nb=eta*NGAM*T9**3*S; R=nb/NAV
        rho=RHOG*T9**4*(1+re)+nnu*7/8*RHOG*Tn**4+rhoB0*Tn**4+R
        Hub=np.sqrt(8*np.pi*GN*rho/3); dT=-Hub*T9/(1+dS/3)
        lt=np.log(min(max(T9,5.1e-3),199.)); lnp=np.exp(sp_np(lt))*TAUN/tau; lpn=np.exp(sp_pn(lt))*TAUN/tau
        f0=rates(T9); f=dict(f0)                      # f: forward rates with any tail change; f0: thermal, used for the reverse of captures
        for k in CHARGED:
            if tail is not None: f[k]=f0[k]*tail(k,T9)
        for k,v in scale.items(): f[k]=f[k]*v; f0[k]=f0[k]*v
        e=np.exp; t32=T9**1.5
        r1=R*f['npdg']*Y[N]*Y[P]-f0['npdg']*0.471e10*t32*e(-25.82/T9)*Y[D]
        r2=R*f['dpg']*Y[D]*Y[P]-f0['dpg']*1.63e10*t32*e(-63.750/T9)*Y[H]
        r3=R*f['ddn']*(Y[D]**2/2-1.73*e(-37.935/T9)*Y[N]*Y[H])
        r4=R*f['ddp']*(Y[D]**2/2-1.73*e(-46.798/T9)*Y[P]*Y[T])
        r5=R*f['hent']*(Y[H]*Y[N]-1.002*e(-8.863/T9)*Y[P]*Y[T])
        r6=R*f['tdn']*(Y[T]*Y[D]-5.54*e(-204.117/T9)*Y[N]*Y[A])
        r7=R*f['hdp']*(Y[H]*Y[D]-5.55*e(-212.980/T9)*Y[P]*Y[A])
        r8=R*f['hag']*Y[H]*Y[A]-f0['hag']*1.11e10*t32*e(-18.423/T9)*Y[B]
        r9=R*f['tag']*Y[T]*Y[A]-f0['tag']*1.11e10*t32*e(-28.640/T9)*Y[L]
        r10=R*f['benp']*(Y[B]*Y[N]-0.998*e(-19.081/T9)*Y[P]*Y[L])
        r11=R*f['lipa']*(Y[L]*Y[P]-4.69*e(-201.291/T9)*Y[A]**2/2)
        w=lnp*Y[N]-lpn*Y[P]; d=np.zeros(9); d[0]=dT
        d[1+N]=-w-r1+r3-r5+r6-r10
        d[1+P]= w-r1-r2+r4+r5+r7+r10-r11
        d[1+D]= r1-r2-2*r3-2*r4-r6-r7
        d[1+T]= r4+r5-r6-r9
        d[1+H]= r2+r3-r5-r7-r8
        d[1+A]= r6+r7-r8-r9+2*r11
        d[1+L]= r9+r10-r11
        d[1+B]= r8-r10
        return d
    xn=1/(1+np.exp(QNP*MEV9/T9i)); y0=np.array([T9i,xn,1-xn,1e-30,1e-30,1e-30,1e-30,1e-30,1e-30])
    stop=lambda t,y: y[0]-T9f; stop.terminal=True
    sol=solve_ivp(rhs,(0,1e9),y0,method='BDF',rtol=1e-7,atol=1e-22,events=stop,dense_output=dense,first_step=1e-9)
    Y=sol.y[1:,-1]; out=dict(Yp=4*Y[A]/np.sum(AMASS*Y),D=Y[D]/Y[P],He3=(Y[H]+Y[T])/Y[P],Li7=(Y[L]+Y[B])/Y[P],Be_frac=Y[B]/(Y[L]+Y[B]),n=Y[N],sumA=float(np.sum(AMASS*Y)),t_end=float(sol.t[-1]))
    if dense: out['sol']=sol
    return out
if __name__=="__main__":
    import time; t0=time.time(); r=run()
    print(f"standard run, eta10 = 6.10, three neutrinos, neutron lifetime {TAUN} s  ({time.time()-t0:.1f} s)")
    print(f"  Yp = {r['Yp']:.4f}   D/H = {r['D']*1e5:.3f}e-5   3He/H = {r['He3']*1e5:.3f}e-5   7Li/H = {r['Li7']*1e10:.2f}e-10   (7Be share {100*r['Be_frac']:.0f}%)   baryons conserved to {abs(r['sumA']-1):.1e}")
    print("  full codes with these rates: Yp 0.246-0.247, D/H 2.5-2.6e-5, 3He/H 1.0e-5, 7Li/H 4.4-5.2e-10")
    for e10 in (3.0,10.0):
        r=run(eta10=e10); print(f"  eta10 = {e10}: Yp {r['Yp']:.4f}  D/H {r['D']*1e5:.2f}e-5  3He/H {r['He3']*1e5:.2f}e-5  7Li/H {r['Li7']*1e10:.2f}e-10")
