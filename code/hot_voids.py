"""Entry 146: hot nearby voids and cold galaxies (Hansen et al. 2026). What can make +30-40 uK through a void at z < 0.03?
(1) Linear ISW of a 2MRS-like projected void in LCDM, and what sign flip or growth rate it would need.
(2) The most any change of the void's potential can do: 2|dPhi|/c^2.
(3) An achromatic energy loss per electron crossed (plasma redshift): coefficient needed, against the Sun and cluster SZ spectra."""
import numpy as np, json
c=299792.458; H0=67.7; Om=0.31; T0=2.7255e6  # uK
ch=c/H0                      # Mpc
z_max=0.03
def Ez(z): return np.sqrt(Om*(1+z)**3+1-Om)
zz=np.linspace(0,z_max,2001); rmax=ch*np.trapezoid(1/Ez(zz),zz)
Omz=lambda z: Om*(1+z)**3/Ez(z)**2
f0=Omz(0.015)**0.55
pref=(3*Om/(8*np.pi))/ch**2   # Phi/c^2 = -pref * int delta/|x-x'| d^3x'  (lengths in Mpc)
def cone_cells(th, r_in, r_out, Nr=70, Nt=24, Np=40, th_in=0.0):
    r=np.linspace(r_in,r_out,Nr+1); rc=0.5*(r[1:]+r[:-1]); dr=np.diff(r)
    mu=np.linspace(np.cos(th),np.cos(th_in),Nt+1); muc=0.5*(mu[1:]+mu[:-1]); dmu=np.diff(mu)
    ph=np.linspace(0,2*np.pi,Np+1)[:-1]+np.pi/Np; dph=2*np.pi/Np
    R,MU,PH=np.meshgrid(rc,muc,ph,indexing='ij'); S=np.sqrt(1-MU**2)
    xyz=np.stack([R*S*np.cos(PH),R*S*np.sin(PH),R*MU],-1).reshape(-1,3)
    dV=(rc[:,None,None]**2*dr[:,None,None]*dmu[None,:,None]*dph*np.ones_like(R)).reshape(-1)
    return xyz,dV
def phi_line(cells, dirn, s):
    xyz,w=cells; pts=s[:,None]*dirn[None,:]
    out=np.empty(len(s))
    for i,p in enumerate(pts):
        d=np.sqrt(((xyz-p)**2).sum(1)); d=np.maximum(d,0.5)
        out[i]=-pref*(w/d).sum()
    return out
def void_isw(th_deg, delta, r_in=10.0, comp=False):
    th=np.radians(th_deg)
    cells=cone_cells(th,r_in,rmax); xyz,dV=cells; w=delta*dV
    if comp:   # compensating wall in the cone shell th..1.4 th, same depth range
        x2,d2=cone_cells(1.4*th,r_in,rmax,th_in=th); w2=-delta*dV.sum()/d2.sum()*d2
        xyz=np.vstack([xyz,x2]); w=np.concatenate([w,w2])
    s=np.linspace(0,4*rmax,500)
    axis=np.array([0,0,1.0]); ref=np.array([np.sin(2.5*th),0,np.cos(2.5*th)])
    pa=phi_line((xyz,w),axis,s); pr=phi_line((xyz,w),ref,s)
    dphi=pa-pr
    # linear ISW: dT/T = 2/c^3 int dPhi/dt dl ; dPhi/dt = H0 (f-1) Phi  ->  dT/T = 2 (f-1) (H0/c) int Phi/c^2 dl
    isw_unit=2*(1/ch)*np.trapezoid(dphi,s)*T0          # uK per unit (f-1)
    return dict(th=th_deg,delta=delta,comp=comp,R_mid_Mpc=float(np.tan(th)*rmax/2),
                phi_axis_max=float(dphi.min()) if delta<0 else float(dphi.max()),
                isw_lcdm=float(isw_unit*(f0-1)), isw_per_fm1=float(isw_unit),
                max_any=float(2*np.abs(dphi).max()*T0))
res=[]
for th in (5,10,15,20):
    for dl in (-0.3,-0.5,-0.8):
        for comp in (False,True):
            r=void_isw(th,dl,comp=comp); res.append(r)
            print(f"th {th:2d} delta {dl:+.1f} comp {comp!s:5} Rmid {r['R_mid_Mpc']:5.1f} Mpc  ISW_LCDM {r['isw_lcdm']:+.3f} uK  per(f-1) {r['isw_per_fm1']:+.3f}  max-any {r['max_any']:.2f} uK  Phi/c2 {r['phi_axis_max']:.2e}")
obs=35.0
base=[r for r in res if r['th']==10 and r['delta']==-0.5 and not r['comp']][0]
need_fm1=obs/base['isw_per_fm1']
print('f0',round(f0,3),'rmax',round(rmax,1),'needed (f-1) for +35 uK, 10 deg, -0.5, uncompensated:',round(need_fm1,1))
# ---- plasma: achromatic loss kappa per electron per cm^2 ----
mp=1.67262e-24; Mpc=3.0857e24; sigT=6.6524587e-25; kB=sigT/2*1.0  # Brynjolfsson: 3.326e-25
kB=3.326e-25
rho_c=1.87847e-29*(H0/100)**2; Obh2=0.0224
ne=Obh2*1.87847e-29/mp*(0.755+0.245/2)   # cm^-3, all baryons ionized
pl={}
for fig in (1.0,0.85):
    for dl in (-0.3,-0.5,-0.8):
        N=ne*fig*abs(dl)*(rmax-10)*Mpc
        k=obs/T0/N
        pl[f'{fig}_{dl}']=dict(N=N,kappa=k,kappa_over_kB=k/kB)
        print(f'f_IGM {fig} delta {dl}: column deficit {N:.2e} cm^-2, kappa needed {k:.2e} cm^2 = {k/kB:.2f} x Brynjolfsson')
# bounds
x=353/56.78; g353=x/np.tanh(x/2)-4
theta=4.4/511.0
k_sz=g353*theta*sigT
print('g353',round(g353,3),'SZ bound kappa <',f'{k_sz:.2e}','=',round(k_sz/kB,3),'kB')
k_sun=17/531*kB
print('Sun bound kappa <',f'{k_sun:.2e}','=',round(17/531,3),'kB')
# what kappa_B would do through a 4.4 keV cluster relative to its RJ SZ decrement (2y): ratio = kB/(2 theta sigT)
print('kB through a 4.4 keV cluster = ',round(kB/(2*theta*sigT),1),'x the Rayleigh-Jeans SZ decrement, at every frequency')
# Hubble-law coefficient for comparison
kH=(H0/ (c*1e5) * 1e5/ Mpc * 1e5)  # placeholder, recomputed below
H0s=H0*1e5/Mpc  # s^-1
kH=H0s/(2.99792458e10*ne*0.85)
print('kappa for a plasma-redshift Hubble law:',f'{kH:.2e}','=',round(kH/kB),'x kB')
json.dump(dict(voids=res,f0=f0,rmax=rmax,need_fm1=need_fm1,plasma=pl,ne=ne,k_sz=k_sz,k_sun=k_sun,kB=kB,kH=kH,g353=g353,obs=obs),open('hot_voids.json','w'),indent=1)
