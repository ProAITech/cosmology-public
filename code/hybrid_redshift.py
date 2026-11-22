"""Hybrid redshift: 1+z = (1+z_c)(1+z_t). z_c from expansion (FLRW, matter + Lambda), z_t from photon energy loss.
Family A: constant log-fraction, ln(1+z_t) = f ln(1+z).
Family B: constant energy-loss rate, d ln E/dt = -f H0 along the path, expansion rate (1-f) H0 E(z_c).
Only expansion stretches time intervals and heats the CMB: stretch = 1+z_c, T = T0 (1+z_c).
DESI DR2 D_M/D_H (Table IV) tests the geometry."""
import numpy as np, json
from scipy.integrate import quad, solve_ivp
from scipy.optimize import minimize_scalar, brentq
data=[(0.510,0.622,0.017),(0.706,0.892,0.021),(0.934,1.223,0.019),(1.321,1.948,0.045),(1.484,2.386,0.136),(2.330,4.518,0.097)]
E=lambda zc,Om: np.sqrt(Om*(1+zc)**3+1-Om)
def F_A(z,f,Om):
    zc=(1+z)**(1-f)-1
    chi=quad(lambda x:1/E(x,Om),0,zc)[0]            # in c/H_c0
    return chi*E(zc,Om)*(1+z)/((1-f)*(1+zc))        # D_M/D_H
def F_B(z,f,Om):
    # d ln(1+z)/d ln(1+z_c) = 1 + f/((1-f) E(z_c));  chi in units of c/H0: d chi = d z_c /((1-f) E)
    g=lambda zc: np.log(1+zc)+ (f/(1-f))*quad(lambda x:1/((1+x)*E(x,Om)),0,zc)[0]-np.log(1+z)
    zc=brentq(g,0,z+1e-9) if f>0 else z
    chi=quad(lambda x:1/((1-f)*E(x,Om)),0,zc)[0]
    dz_dzc=(1+z)*(1/(1+zc)+ (f/(1-f))/((1+zc)*E(zc,Om)))
    return chi*(1-f)*E(zc,Om)*dz_dzc, zc
def chi2(fn,f,Om): return sum(((o-(fn(z,f,Om) if fn is F_A else fn(z,f,Om)[0]))/s)**2 for z,o,s in data)
out={}
for name,fn in (('A',F_A),('B',F_B)):
    fs=np.round(np.concatenate([np.arange(0,0.30,0.01),np.arange(0.30,0.96,0.05)]),3); prof=[]
    for f in fs:
        r=minimize_scalar(lambda Om:chi2(fn,f,Om),bounds=(0.02,0.98),method='bounded'); prof.append((float(f),float(r.x),float(r.fun)))
    c0=min(p[2] for p in prof); best=min(prof,key=lambda p:p[2])
    lim95=None
    for i in range(1,len(prof)):
        if prof[i][2]-c0>3.84 and prof[i][0]>best[0]:
            a,b=prof[i-1],prof[i]; lim95=a[0]+(3.84-(a[2]-c0))*(b[0]-a[0])/((b[2]-c0)-(a[2]-c0)); break
    out[name]=dict(best=best,chi2_f0=prof[0][2],Om_f0=prof[0][1],lim95=lim95,profile=prof)
    print(name,'best f,Om,chi2',[round(x,3) for x in best],'| f=0: Om',round(prof[0][1],3),'chi2',round(prof[0][2],2),'| 95% upper f',None if lim95 is None else round(lim95,3))
    for p in prof[::5]: print('   f=%.2f Om=%.3f chi2=%.1f'%p)
# time dilation (DES): stretch exponent b = 1.003 +- 0.005 stat +- 0.010 sys ; family A: b = 1-f
sb=np.hypot(0.005,0.010); fa=1-1.003
print('time dilation: f = %.3f +- %.3f ; 95%% upper (one-sided, f>=0): %.3f'%(fa,sb,fa+1.645*sb))
# CMB T(z) = T0 (1+z)^(1-alpha), alpha = -0.007 +- 0.013 ; family A: alpha = f
print('CMB temperature: f = %.3f +- %.3f ; 95%% upper: %.3f'%(-0.007,0.013,-0.007+1.645*0.013))
w1,w2=1/sb**2,1/0.013**2; fc=(fa*w1-0.007*w2)/(w1+w2); sc=(w1+w2)**-0.5
print('combined clocks+temperature: f = %.4f +- %.4f ; 95%% upper: %.3f ; 99.7%% upper %.3f'%(fc,sc,fc+1.645*sc,fc+3*sc))
out['clock']=dict(f=fa,s=float(sb)); out['temp']=dict(f=-0.007,s=0.013); out['comb']=dict(f=float(fc),s=float(sc),u95=float(fc+1.645*sc))
# family B: effective stretch exponent at z=0.5,1 for a given f (Om 0.3)
for f in (0.05,0.1,0.2):
    for z in (0.5,1.0):
        zc=F_B(z,f,0.3)[1]; print('  family B f=%.2f z=%.1f: effective time-dilation exponent %.3f'%(f,z,np.log(1+zc)/np.log(1+z)))
out['B_exp']={str(f):float(np.log(1+F_B(0.7,f,0.3)[1])/np.log(1.7)) for f in (0.01,0.02,0.05,0.1,0.2)}
print(out['B_exp'])
json.dump(out,open('hybrid.json','w'))   # full profiles
