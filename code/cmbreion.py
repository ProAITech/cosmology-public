"""Semi-analytic CMB temperature AND E-mode polarization spectra from primordial-plasma physics
(Entry 92: extends the Entry 12 code with the polarization source of Zaldarriaga & Seljak 1997).
Entry 94: adds reionization (tanh in (1+z)^1.5, width 0.5 in z, as in CAMB; helium II at z = 3.5),
with the redshift chosen to give a target optical depth TAU_RE, and samples the late visibility.
Tight-coupled photon-baryon fluid + CDM + neutrino fluid in conformal Newtonian gauge,
Peebles/RECFAST-style hydrogen recombination, Silk damping, line-of-sight projection
with visibility function and early ISW. Accuracy target ~10% on peak heights.
Usage: TAU_RE=0.054 python3 cmbreion.py 0.0222 0.2975 0.6851 0.965 tau054   (about 8 minutes; also tau0 with TAU_RE=0
and tau09 with TAU_RE=0.09; writes cmbreion_<label>.json, which reion_compare.py reads)
Entry 95: TAU_RE=0.09 ZRE=6.95 FLASH=15,30 python3 cmbreion.py 0.0222 0.2975 0.6851 0.965 flash09  adds an early burst
between z = 15 and 30 to a late step at z = 6.95, with its strength set to reach the total TAU_RE
"""
import sys, json, math, os
HERE=os.path.dirname(os.path.abspath(__file__))
import numpy as np
from scipy.integrate import solve_ivp, cumulative_trapezoid
from scipy.special import spherical_jn

wb, Om, h, ns = [float(x) for x in sys.argv[1:5]]
label = sys.argv[5] if len(sys.argv) > 5 else 'run'

# ---------- background ----------
T0 = 2.7255; Yp = 0.247; Neff = 3.046
wg = 2.469e-5; wnu = 0.2271 * Neff * wg
wm = Om * h * h; wc = wm - wb
Og, On, Ob, Oc = wg/h**2, wnu/h**2, wb/h**2, wc/h**2
Or = Og + On; OL = 1 - Om - Or
H0 = h / 2997.92458  # Mpc^-1
def Ha(a): return H0 * np.sqrt(Or/a**4 + Om/a**3 + OL)
lna = np.linspace(np.log(1e-8), 0, 40000); a_g = np.exp(lna)
deta = 1/(a_g*Ha(a_g))                      # d eta / d ln a
eta_g = cumulative_trapezoid(deta, lna, initial=0) + 1e-8/(H0*np.sqrt(Or))  # analytic start
eta0 = eta_g[-1]
def eta_of_a(a): return np.interp(np.log(a), lna, eta_g)

# ---------- recombination (H only, Peebles 3-level, RECFAST fudge 1.14) ----------
kB=1.380649e-23; hP=6.62607015e-34; me=9.1093837e-31; c=2.99792458e8; sigT=6.6524587e-29
Eion=2.179872e-18; E2=Eion/4; Ealpha=0.75*Eion; lam_a=121.5668e-9; L2s=8.2245809
mH=1.6735575e-27; Mpc=3.0857e22
nH0 = (1-Yp) * wb * 1.87847e-26 / mH
def alphaB(T):
    t=T/1e4; return 1.14*1e-19*4.309*t**-0.6166/(1+0.6703*t**0.5300)
def rhs_x(z, x):
    T=T0*(1+z); nH=nH0*(1+z)**3; xe=min(max(x[0],1e-12),1.0)
    al=alphaB(T); beta=al*(2*np.pi*me*kB*T/hP**2)**1.5*np.exp(-E2/(kB*T))
    Hs=Ha(1/(1+z))*c/Mpc
    La=8*np.pi*Hs/(lam_a**3*nH*max(1-xe,1e-12))
    Cr=(L2s+La)/(L2s+La+beta)
    return [Cr*(nH*xe*xe*al - beta*(1-xe)*np.exp(-Ealpha/(kB*T)))/(Hs*(1+z))]
zs_=np.linspace(1700, 0, 3401)
import os
B_CLUMP=float(os.environ.get('CLUMP_B','0')); MU=float(os.environ.get('ME_RATIO','1')); CDM_CS2=float(os.environ.get('CDM_CS2','0'))
def zones_M1(b):
    # Jedamzik & Pogosian (2020) model M1: f2=1/3, D1=0.1, D2=1; f1,f3,D3 from mean density and clumping b
    from scipy.optimize import brentq
    D1,D2,f2=0.1,1.0,1/3
    def resid(f3):
        f1=1-f2-f3; D3=(1-f1*D1-f2*D2)/f3
        return f1*D1**2+f2*D2**2+f3*D3**2-(1+b)
    f3=brentq(resid,1e-4,2/3-1e-6); f1=1-f2-f3; D3=(1-f1*D1-f2*D2)/f3
    return [(f1,D1),(f2,D2),(f3,D3)]
ZONES=zones_M1(B_CLUMP) if B_CLUMP>0 else [(1.0,1.0)]
Eion*=MU; E2*=MU; Ealpha*=MU; lam_a/=MU; L2s*=MU; sigT/=MU**2; me*=MU
_alphaB0=alphaB
def alphaB(T): return MU**-2*_alphaB0(T/MU)
def solve_zone(D):
    def rhs_z(z,x):
        T=T0*(1+z); nH=D*nH0*(1+z)**3; xe=min(max(x[0],1e-12),1.0)
        al=alphaB(T); beta=al*(2*np.pi*me*kB*T/hP**2)**1.5*np.exp(-E2/(kB*T))
        Hs=Ha(1/(1+z))*c/Mpc
        La=8*np.pi*Hs/(lam_a**3*nH*max(1-xe,1e-12))
        Cr=(L2s+La)/(L2s+La+beta)
        return [Cr*(nH*xe*xe*al - beta*(1-xe)*np.exp(-Ealpha/(kB*T)))/(Hs*(1+z))]
    T=T0*1701; nH=D*nH0*1701**3; S=(2*np.pi*me*kB*T/hP**2)**1.5*np.exp(-Eion/(kB*T))/nH
    x0=(-S+np.sqrt(S*S+4*S))/2
    return solve_ivp(rhs_z,(1700,0),[x0],t_eval=zs_,method='LSODA',rtol=1e-6,atol=1e-10).y[0]
xe_z=sum(f*D*solve_zone(D) for f,D in ZONES)   # <x_e> = sum f_V Delta x_e (Eq. A2)
def xe_rec(z):
    z=np.asarray(z); out=np.interp(-z, -zs_, xe_z)
    return np.where(z>1700, 1.0, out)
TAU_RE=float(os.environ.get('TAU_RE','0.054')); fHe=Yp/(4*(1-Yp))
# Entry 95: an optional early burst, FLASH="zlo,zhi": a fraction F_FL of hydrogen (and helium's first electron) ionized
# between zlo and zhi, edges 1 in z wide. With ZRE set, the late step is fixed there and F_FL is chosen to reach TAU_RE.
FLASH=os.environ.get('FLASH'); ZRE_FIX=os.environ.get('ZRE'); F_FL=0.0
FL_LO,FL_HI=(float(v) for v in FLASH.split(',')) if FLASH else (0.0,0.0)
def xe_re(z, zre):
    z=np.asarray(z,float); y=(1+z)**1.5; yre=(1+zre)**1.5; dy=1.5*np.sqrt(1+zre)*0.5
    fl=F_FL*(1+fHe)*0.5*(np.tanh(z-FL_LO)-np.tanh(z-FL_HI)) if FLASH else 0.0
    return (1+fHe)/2*(1+np.tanh((yre-y)/dy)) + fHe/2*(1+np.tanh((3.5-z)/0.5)) + fl
def tau_reion(zre):
    zz=np.linspace(0,40,8001); a=1/(1+zz); x=np.maximum(xe_rec(zz),xe_re(zz,zre))-xe_rec(zz)
    dtdz=1/((1+zz)*Ha(a)*c/Mpc)            # dt/dz in s (Ha in Mpc^-1 -> s^-1 via c/Mpc)
    return float(np.trapezoid(x*nH0*(1+zz)**3*sigT*c*dtdz, zz))
if TAU_RE>0 and FLASH and ZRE_FIX:
    Z_RE=float(ZRE_FIX); lo,hi=0.0,1.0
    for _ in range(60):
        F_FL=0.5*(lo+hi)
        if tau_reion(Z_RE)<TAU_RE: lo=F_FL
        else: hi=F_FL
    F_FL=0.5*(lo+hi); print('flash fraction',F_FL,'between z',FL_LO,FL_HI,flush=True)
elif TAU_RE>0:
    lo,hi=3.0,25.0
    for _ in range(60):
        mid=0.5*(lo+hi)
        if tau_reion(mid)<TAU_RE: lo=mid
        else: hi=mid
    Z_RE=0.5*(lo+hi)
else: Z_RE=None
print('z_re',Z_RE,'tau_re',tau_reion(Z_RE) if Z_RE else 0,flush=True)
def xe_of_z(z):
    r=xe_rec(z)
    return np.maximum(r, xe_re(z,Z_RE)) if Z_RE else r
# optical depth along eta
a_r = a_g; z_r = 1/a_r - 1
ne = xe_of_z(np.minimum(z_r, 1e6)) * nH0 * (1+z_r)**3
dtau = ne*sigT*a_r*Mpc                         # d tau / d eta  [Mpc^-1]
tau_from_end = -cumulative_trapezoid((dtau*deta)[::-1], lna[::-1], initial=0)[::-1]
tau = tau_from_end
g = dtau*np.exp(-tau)                          # visibility per unit eta
R_g = 3*wb/(4*wg)*a_r
kD2inv = cumulative_trapezoid(deta/(6*(1+R_g)*dtau)*(R_g**2/(1+R_g)+16/15), lna, initial=0)
imax=np.argmax(g); z_star=z_r[imax]; eta_star=eta_g[imax]

# eta sample points: visibility window + early ISW tail
gnorm=g/np.max(g)
win=np.where((gnorm>1e-3)&(z_r>200))[0]
a_lo, a_hi = a_r[win[0]], a_r[win[-1]]
a_vis=np.exp(np.linspace(np.log(a_lo),np.log(a_hi),28))
a_isw=np.exp(np.linspace(np.log(a_hi),np.log(0.15),31))[1:]
a_re=np.exp(np.linspace(np.log(1/41),np.log(1/3.0),90 if FLASH else 60))
a_pts=np.unique(np.concatenate([a_vis,a_isw,a_re]))
eta_pts=eta_of_a(a_pts)
g_pts=np.interp(np.log(a_pts),lna,g); expmt=np.interp(np.log(a_pts),lna,np.exp(-tau))
kD2_pts=np.interp(np.log(a_pts),lna,kD2inv)
w_eta=np.gradient(eta_pts)                     # quadrature weights in eta


# ---------- perturbations: CDM + baryons + photons + free-streaming neutrinos (Ma & Bertschinger 1995, Newtonian gauge) ----------
Rnu=On/(Og+On); LN=10; LG=8
lna_tau=lna; dtau_g=dtau
def taudot(a): return np.interp(np.log(a),lna_tau,dtau_g)
a_dec=np.exp(np.interp(-0.0,-np.log(1.0)+0*lna,lna)) if False else None
# switch off tight coupling just after the visibility peak (z ~ 1000), when photons start to free-stream
i_dec=np.searchsorted(-z_r, -1000.0)  # z_r decreasing with index
Z_TC=float(sys.argv[7]) if len(sys.argv)>7 else 1600.0
x_dec=float(np.log(1/(1+Z_TC)))
kD2_dec=float(np.interp(x_dec,lna,kD2inv))
LP=6
NV=8+(LN-1)+(LG-1)+(LP+1)+1   # dc thc dg thg db thb dn thn | Fn2..FnLN | Fg2..FgLG | G0..GLP | phi
iFn=8; iFg=8+(LN-1); iG=iFg+(LG-1); iphi=NV-1
def solve_k(k, check=False):
    ai=1e-7; eta_i=eta_of_a(ai)
    psi0=20/(15+4*Rnu); phi0=(1+0.4*Rnu)*psi0; th=0.5*k*k*eta_i*psi0
    y0=np.zeros(NV); y0[:8]=[-1.5*psi0,th,-2*psi0,th,-1.5*psi0,th,-2*psi0,th]
    y0[iFn]=(2/15)*(k*eta_i)**2*psi0; y0[iphi]=phi0
    def parts(a,y,tca,evolve):
        dc,thc,dg,thg,db,thb,dn,thn=y[:8]
        Fn=y[iFn:iFn+LN-1]; Fg=y[iFg:iFg+LG-1]
        fm=1.5*H0**2/a; fr=1.5*H0**2/a**2; Hc=a*Ha(a)
        sig=0.5*(On*Fn[0]+(0 if tca else Og*Fg[0]))
        Sd=fm*(Oc*dc+Ob*db)+fr*(Og*dg+On*dn)
        St=fm*(Oc*thc+Ob*thb)+fr*(4/3)*(Og*thg+On*thn)
        if evolve:
            phi=y[iphi]; psi=phi-4*fr*sig/k**2; dphi=-Hc*psi-(k*k*phi+Sd)/(3*Hc)
        else:
            phi=-(Sd+3*Hc*St/k**2)/k**2; psi=phi-4*fr*sig/k**2; dphi=St/k**2-Hc*psi
        return Hc,fr,phi,psi,dphi
    def deriv(x,y,tca,evolve):
        a=np.exp(x); Hc,fr,phi,psi,dphi=parts(a,y,tca,evolve)
        dc,thc,dg,thg,db,thb,dn,thn=y[:8]
        Fn=y[iFn:iFn+LN-1]; Fg=y[iFg:iFg+LG-1]
        R=3*wb/(4*wg)*a; td=taudot(a); eta=eta_of_a(a)
        d=np.zeros(NV)
        d[0]=-thc+3*dphi; d[1]=-Hc*thc+CDM_CS2*k*k*dc+k*k*psi
        d[2]=-(4/3)*thg+4*dphi
        if tca:
            d[3]=(-Hc*R*thg+0.25*k*k*dg)/(1+R)+k*k*psi
            d[4]=0.75*d[2]; d[5]=d[3]
        else:
            sg=0.5*Fg[0]
            d[3]=k*k*(0.25*dg-sg)+k*k*psi+td*(thb-thg)
            d[4]=-thb+3*dphi; d[5]=-Hc*thb+(td/R)*(thg-thb)+k*k*psi
            # photon hierarchy l=2..LG (unpolarized Thomson)
            F=np.concatenate([[dg, 4*thg/(3*k)],Fg])  # F0=delta, F1=4theta/(3k)
            dF=np.zeros(LG-1)
            G=y[iG:iG+LP+1]; Pi=F[2]+G[0]+G[2]
            dF[0]=(8/15)*thg-(3/5)*k*F[3]-0.9*td*F[2]+0.1*td*(G[0]+G[2])
            for l in range(3,LG):
                dF[l-2]=k/(2*l+1)*(l*F[l-1]-(l+1)*F[l+1])-td*F[l]
            dF[LG-2]=k*F[LG-1]-(LG+1)/eta*F[LG]-td*F[LG]
            d[iFg:iFg+LG-1]=dF
            dG=np.zeros(LP+1)
            dG[0]=-k*G[1]+td*(-G[0]+0.5*Pi)
            for l in range(1,LP):
                dG[l]=k/(2*l+1)*(l*G[l-1]-(l+1)*G[l+1])-td*G[l]
            dG[2]+=td*0.1*Pi
            dG[LP]=k*G[LP-1]-(LP+1)/eta*G[LP]-td*G[LP]
            d[iG:iG+LP+1]=dG
        sn=0.5*Fn[0]
        d[6]=-(4/3)*thn+4*dphi; d[7]=k*k*(0.25*dn-sn)+k*k*psi
        F=np.concatenate([[dn,4*thn/(3*k)],Fn]); dF=np.zeros(LN-1)
        dF[0]=(8/15)*thn-(3/5)*k*F[3]
        for l in range(3,LN):
            dF[l-2]=k/(2*l+1)*(l*F[l-1]-(l+1)*F[l+1])
        dF[LN-2]=k*F[LN-1]-(LN+1)/eta*F[LN]
        d[iFn:iFn+LN-1]=dF
        if evolve: d[iphi]=dphi
        return d/Hc, (phi,psi,dphi)
    xe=np.log(a_pts); xend=xe[-1]
    xsw=float(np.interp(2.0/k, eta_g, lna))
    bounds=sorted(set([np.log(ai),min(xsw,xend),min(x_dec,xend),xend]))
    out_y=np.zeros((NV,len(xe))); out_aux=np.zeros((4,len(xe)))
    y=y0.copy()
    for s0,s1 in zip(bounds[:-1],bounds[1:]):
        if s1<=s0: continue
        mid=0.5*(s0+s1); tca=mid<x_dec; evolve=mid<xsw
        if not tca and mid-s0<1e-12: pass
        if (not tca) and abs(s0-x_dec)<1e-12:
            a=np.exp(s0); F2=(32/45)*y[3]/taudot(a); y[iFg]=F2; y[iG]=1.25*F2; y[iG+2]=0.25*F2   # tight-coupling quadrupole and polarization at the switch
        mask=(xe>=s0)&(xe<s1) if s1<xend else (xe>=s0)&(xe<=s1)
        tev=np.append(xe[mask],s1) if (len(xe[mask])==0 or xe[mask][-1]<s1) else xe[mask]
        sol=solve_ivp(lambda x,yy: deriv(x,yy,tca,evolve)[0],(s0,s1),y,t_eval=tev,method='DOP853',rtol=1e-6,atol=1e-9)
        idx=np.where(mask)[0]
        for j,ii in enumerate(idx):
            yy=sol.y[:,j]; dd,(phi,psi,dphi)=deriv(xe[ii],yy,tca,evolve)
            a=np.exp(xe[ii]); Hc=a*Ha(a); fr=1.5*H0**2/a**2
            sig=0.5*(On*yy[iFn]+(0 if tca else Og*yy[iFg])); dsig=0.5*(On*dd[iFn]+(0 if tca else Og*dd[iFg]))*Hc
            dpsi=dphi-(4/k**2)*(-2*Hc*fr*sig+fr*dsig)
            out_y[:,ii]=yy; out_aux[:,ii]=[phi,psi,dphi,dpsi]
            out_aux[3,ii]=dpsi
            out_aux_tca[ii]=tca
        y=sol.y[:,-1]
    psi=out_aux[1]; isw=out_aux[2]+out_aux[3]
    damp=np.where(out_aux_tca,np.exp(-k*k*kD2_pts),np.exp(-k*k*min(kD2_dec,1e9)))
    SM=(0.25*out_y[2]+psi)*damp
    SD=(out_y[5]/k)*damp          # baryon velocity (Doppler)
    if check: print('k',k,'phi*',round(out_aux[0][np.argmax(g_pts)],4),'psi*',round(psi[np.argmax(g_pts)],4))
    PI=(out_y[iFg]+out_y[iG]+out_y[iG+2])*np.where(out_aux_tca,0.0,np.exp(-k*k*min(kD2_dec,1e9)))   # Pi = F2+G0+G2 (Ma & Bertschinger)
    return SM,SD,isw,PI
out_aux_tca=np.zeros(len(a_pts),dtype=bool)

if len(sys.argv)>6 and sys.argv[6]=='check':
    import time
    for k in (2e-4,0.002,0.0159,0.05,0.15):
        t=time.time(); SMk,SDk,DPk,PIk=solve_k(k,True); i=np.argmax(g_pts)
        print('   SM*',round(SMk[i],4),'SD*',round(SDk[i],4),'isw*',f'{DPk[i]:.2e}','t',round(time.time()-t,2))
    sys.exit()

kgrid=np.arange(2e-4,0.20,1.5e-4)
import time; t0=time.time()
SMs=[];SDs=[];DPs=[];PIs=[]
for k in kgrid:
    a1,a2,a3,a4=solve_k(k); SMs.append(a1);SDs.append(a2);DPs.append(a3);PIs.append(a4)
SM=np.array(SMs); SD=np.array(SDs); DP=np.array(DPs); PI=np.array(PIs)
print('ode done',round(time.time()-t0),flush=True)
np.savez(os.path.join(HERE,f'srcR_{label}.npz'),kgrid=kgrid,SM=SM,SD=SD,DP=DP,PI=PI,eta_pts=eta_pts,g_pts=g_pts,expmt=expmt,w_eta=w_eta,eta0=eta0)
dist=eta0-eta_pts
ells=np.unique(np.concatenate([np.arange(2,40,1),np.arange(40,200,8),np.arange(200,1300,20)]))
NS=sorted(set([ns,0.965]))
Dls={n:[] for n in NS}; Ees={n:[] for n in NS}; TEs={n:[] for n in NS}
for l in ells:
    x=np.outer(kgrid,dist)
    jl=spherical_jn(l,x); jlp=spherical_jn(l,x,derivative=True)
    Th=(SM*g_pts*w_eta*jl).sum(1)+(SD*g_pts*w_eta*jlp).sum(1)+(DP*expmt*w_eta*jl).sum(1)
    El=math.sqrt((l+2)*(l+1)*l*(l-1))*(0.75*PI*g_pts*w_eta*jl/x**2).sum(1)
    for n in NS:
        Dls[n].append(l*(l+1)*np.trapezoid(kgrid**(n-1)*Th**2/kgrid,kgrid))
        Ees[n].append(l*(l+1)*np.trapezoid(kgrid**(n-1)*El**2/kgrid,kgrid))
        TEs[n].append(l*(l+1)*np.trapezoid(kgrid**(n-1)*Th*El/kgrid,kgrid))
print('proj done',round(time.time()-t0),flush=True)
out={}
def extrema(y,lmin=100):
    pk=[];dp=[]
    for i in range(1,len(ells)-1):
        if ells[i]<lmin: continue
        if y[i]>y[i-1] and y[i]>=y[i+1]:
            # parabolic refinement
            l0,l1,l2=ells[i-1:i+2]; y0,y1,y2=y[i-1:i+2]; A=np.polyfit([l0,l1,l2],[y0,y1,y2],2); lp=-A[1]/(2*A[0]) if A[0]!=0 else l1
            pk.append((round(float(lp),1),float(np.polyval(A,lp))))
        if y[i]<y[i-1] and y[i]<=y[i+1]:
            l0,l1,l2=ells[i-1:i+2]; y0,y1,y2=y[i-1:i+2]; A=np.polyfit([l0,l1,l2],[y0,y1,y2],2); lp=-A[1]/(2*A[0]) if A[0]!=0 else l1
            dp.append((round(float(lp),1),float(np.polyval(A,lp))))
    return pk,dp
for n in NS:
    Dl=np.array(Dls[n]); El=np.array(Ees[n]); Xl=np.array(TEs[n])
    n1=max(Dl[(ells>150)&(ells<300)])
    tp,td=extrema(Dl,150); ep,ed=extrema(El,100); xp,xd=extrema(Xl,100)
    zc=[round(float(ells[i]-Xl[i]*(ells[i+1]-ells[i])/(Xl[i+1]-Xl[i])),1) for i in range(len(ells)-1) if ells[i]>=100 and Xl[i]*Xl[i+1]<0]
    out[str(n)]=dict(TT_peaks=[(l,round(v/n1,4)) for l,v in tp[:7]],TT_dips=[(l,round(v/n1,4)) for l,v in td[:6]],
                     EE_peaks=[(l,round(v/n1,6)) for l,v in ep[:6]],EE_dips=[(l,round(v/n1,6)) for l,v in ed[:5]],
                     TE_max=[(l,round(v/n1,5)) for l,v in xp[:7]],TE_min=[(l,round(v/n1,5)) for l,v in xd[:7]],TE_zero=zc[:12],
                     TT=(Dl/n1).round(6).tolist(),EE=(El/n1).round(8).tolist(),TE=(Xl/n1).round(7).tolist())
res=dict(z_re=Z_RE,tau_re=TAU_RE,flash=[FL_LO,FL_HI,F_FL] if FLASH else None,label=label,wb=wb,Om=Om,h=h,z_star=round(float(z_star),1),eta0=round(float(eta0),1),ells=ells.tolist(),runs=out)
json.dump(res,open(os.path.join(HERE,f'cmbreion_{label}.json'),'w'))
print(json.dumps({k:{kk:vv for kk,vv in v.items() if kk not in ('TT','EE','TE')} for k,v in out.items()}))
