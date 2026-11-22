"""Distances and likelihoods for Workbench Entry 36: flat models with dark energy w(a) = w0 + wa (1 - a), photons, and neutrinos
(N_eff = 3.044, one of 0.06 eV). Data: DESI DR2 clustering ruler (13 numbers), DESI's compressed CMB prior on (theta*, omega_b, omega_bc),
JLA binned supernovae with full covariance, Pantheon supernovae with statistical errors."""
import sys, os, io, contextlib
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
import numpy as np
from scipy import integrate
C_KMS=299792.458; OG_H2=2.4728e-5; TNU_EV=1.6765e-4; NEFF=3.044; MNU=0.06
def nu_h2(a):
    """neutrino density x h^2 / a^-4 scaling handled inside: returns omega_nu(a) such that rho ~ omega_nu(a)/a^4"""
    f=lambda y: (1+(0.3173*y)**1.83)**(1/1.83)
    return 0.2271*OG_H2*(NEFF/3)*(2+f(MNU*a/TNU_EV))
OMNU0=nu_h2(1.0)                       # today, as a matter-like density x h^2
ZG=np.concatenate([np.linspace(0,3,3001)[:-1],np.geomspace(3,1200,1500)])
class Model:
    def __init__(s,ob,obc,h,w0=-1.0,wa=0.0,alpha=1.0):
        s.ob,s.obc,s.h,s.w0,s.wa,s.alpha=ob,obc,h,w0,wa,alpha
        s.ode=h*h-obc-OG_H2-OMNU0     # flat
        a=1/(1+ZG); de=s.ode*a**(-3*(1+w0+wa))*np.exp(-3*wa*(1-a))
        s.E=np.sqrt(np.clip(obc/a**3+(OG_H2+nu_h2(a))/a**4+de,1e-30,None))        # H(z)/(100 km/s/Mpc)
        s.dc=integrate.cumulative_trapezoid(C_KMS/100/s.E,ZG,initial=0)           # comoving distance, Mpc
    def DM(s,z): return np.interp(z,ZG,s.dc)
    def DH(s,z): return C_KMS/100/np.interp(z,ZG,s.E)
    def rd(s): return s.alpha*147.05*(s.ob/0.02236)**-0.13*(s.obc/0.1432)**-0.23*(NEFF/3.04)**-0.1
    def zstar(s):
        ob,om=s.ob,s.obc; g1=0.0783*ob**-0.238/(1+39.5*ob**0.763); g2=0.560/(1+21.1*ob**1.81); return 1048*(1+0.00124*ob**-0.738)*(1+g1*om**g2)
    def theta_star(s):
        zs=s.zstar(); zz=np.geomspace(zs,1e7,4000); a=1/(1+zz)
        E=np.sqrt(s.obc/a**3+(OG_H2*(1+0.2271*NEFF))/a**4)                         # early times: matter + radiation
        R=31500*s.ob*(2.7255/2.7)**-4*a; cs=C_KMS/np.sqrt(3*(1+R)); rs=integrate.trapezoid(cs/(100*E),zz)
        # distance to z*: grid to 1200, then matter + radiation beyond
        if zs<=ZG[-1]: dm=np.interp(zs,ZG,s.dc)
        else:
            z2=np.linspace(ZG[-1],zs,200); a2=1/(1+z2); dm=s.dc[-1]+integrate.trapezoid(C_KMS/100/np.sqrt(s.obc/a2**3+(OG_H2+nu_h2(a2))/a2**4),z2)
        return s.alpha*rs/dm
# calibrate theta* to Planck 2018 (TT,TE,EE+lowE+lensing): omega_b 0.02237, omega_c 0.1200, H0 67.36 -> 100 theta* = 1.04110
_KAPPA=0.0104110/Model(0.02237,0.02237+0.1200,0.6736).theta_star()
# ---------------- data
BAO=[(0.510,13.588,0.167,21.863,0.425,-0.459),(0.706,17.351,0.177,19.455,0.330,-0.404),(0.934,21.576,0.152,17.641,0.193,-0.416),(1.321,27.601,0.318,14.176,0.221,-0.434),(1.484,30.512,0.760,12.817,0.516,-0.500),(2.330,38.988,0.531,8.632,0.101,-0.431)]
ZV,DV,EDV=0.295,7.942,0.075; ZB=np.array([b[0] for b in BAO]); YB=np.r_[DV,[b[1] for b in BAO],[b[3] for b in BAO]]; CB=np.zeros((13,13)); CB[0,0]=EDV**2
for i,b in enumerate(BAO): CB[1+i,1+i]=b[2]**2; CB[7+i,7+i]=b[4]**2; CB[1+i,7+i]=CB[7+i,1+i]=b[5]*b[2]*b[4]
CBI=np.linalg.inv(CB)
BAO_NAMES=["BGS 0.30","LRG1 0.51","LRG2 0.71","LRG3+ELG1 0.93","ELG2 1.32","QSO 1.48","Lya 2.33"]
def bao_mask(drop=None):
    """indices of the 13 numbers kept when tracer number `drop` (0-6) is removed"""
    keep=np.ones(13,bool)
    if drop is not None:
        if drop==0: keep[0]=False
        else: keep[drop]=False; keep[6+drop]=False
    return keep
def chi_bao(m,drop=None):
    rd=m.rd(); t=np.r_[(ZV*m.DM(ZV)**2*m.DH(ZV))**(1/3),m.DM(ZB),m.DH(ZB)]/rd; r=YB-t
    if drop is None: return float(r@CBI@r)
    k=bao_mask(drop); return float(r[k]@np.linalg.inv(CB[np.ix_(k,k)])@r[k])
MU_CMB=np.array([0.01041,0.02223,0.14208]); C_CMB=1e-9*np.array([[0.006621,0.12444,-1.1929],[0.12444,21.344,-94.001],[-1.1929,-94.001,1488.4]]); CI_CMB=np.linalg.inv(C_CMB)
def chi_cmb(m,which=(0,1,2)):
    v=np.array([_KAPPA*m.theta_star(),m.ob,m.obc])-MU_CMB; w=list(which); return float(v[w]@np.linalg.inv(C_CMB[np.ix_(w,w)])@v[w])
_cwd=os.getcwd(); os.chdir(HERE)
with contextlib.redirect_stdout(io.StringIO()): import static_hubble as sh
os.chdir(_cwd)
def chi_jla(m,sel=None):
    s=np.ones(31,bool) if sel is None else sel; z=sh.zb[s]; Ci=np.linalg.inv(sh.C[np.ix_(s,s)]); one=np.ones(s.sum())
    d=sh.mub[s]-5*np.log10((1+z)*m.DM(z)); a=(one@Ci@d)/(one@Ci@one); r=d-a; return float(r@Ci@r)
from pantheon_rows import z as ZP, m as MP
def _scatter():
    x=np.log10(1+ZP); base=5*np.log10(ZP); A=np.vander(x,5); edges=np.quantile(ZP,np.linspace(0,1,15)); mid=np.array([np.median(ZP[(ZP>=a)&(ZP<=b)]) for a,b in zip(edges[:-1],edges[1:])]); sig=np.full(len(ZP),0.15)
    for it in range(4):
        cf=np.linalg.lstsq(A/sig[:,None],(MP-base)/sig,rcond=None)[0]; res=MP-base-A@cf
        sig=np.interp(np.log(ZP),np.log(mid),np.array([res[(ZP>=a)&(ZP<=b)].std(ddof=1) for a,b in zip(edges[:-1],edges[1:])]))
    return sig
SIGP=_scatter(); WP=1/SIGP**2
def chi_pan(m,sel=None,step=None):
    s=np.ones(len(ZP),bool) if sel is None else sel; z=ZP[s]; w=WP[s]; d=MP[s]-5*np.log10((1+z)*m.DM(z))
    if step is None:
        a=np.sum(w*d)/w.sum(); return float(np.sum(w*(d-a)**2))
    lo=z<step; return float(sum(np.sum(w[q]*(d[q]-np.sum(w[q]*d[q])/w[q].sum())**2) for q in (lo,~lo)))
