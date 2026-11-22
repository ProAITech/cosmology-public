"""Workbench Entry 32: does selection bend the SPT-with-XMM test? The simulation of Entry 28, rerun for this sample.
The 57 clusters reach down to SPT significance 4.7, lower than the Chandra sample's 6.8, so noise-boosted detections matter more.
Same ingredients as ratio_selection_sim.py of Entry 28: Sheth-Tormen mass function (evolving for the expanding model, frozen for the static one),
universal pressure profile for the shadow inside 0.75 arcmin, SPT noise from the catalogue, scatter between shadow and X-ray pressure, X-ray error.
Each simulated sample has the redshifts of the real one and is passed through the same fits."""
import numpy as np, json, sys, os
from scipy import integrate
dd=sys.argv[1] if len(sys.argv)>1 else "."
cH=299792.458/70.; Ez=lambda z: np.sqrt(0.3*(1+z)**3+0.7)
DA=np.vectorize(lambda z: cH*integrate.quad(lambda x:1/Ez(x),0,z)[0]/(1+z)); d_s=lambda z: cH*np.log(1+z)
F_log=lambda z: np.sqrt(d_s(z)/DA(z))*(1+z)**-1.5
am=180*60/np.pi
c5,ga,al,be=1.177,0.3081,1.0510,5.4905
p=lambda x: 1/((c5*x)**ga*(1+(c5*x)**al)**((be-ga)/al))
los=lambda b: 2*integrate.quad(lambda l:p(np.hypot(b,l)),0,np.sqrt(64-b*b),limit=200)[0]
I0=los(1e-6); ag=np.geomspace(0.03,4,60); Kg=np.array([integrate.quad(lambda b:2*np.pi*b*los(b)/I0,0,a,limit=200)[0] for a in ag])
Kc=lambda a: np.exp(np.interp(np.log(a),np.log(ag),np.log(Kg)))
def growth(z):
    g=lambda zz: (lambda om,ol: 2.5*om/(om**(4/7)-ol+(1+om/2)*(1+ol/70))/(1+zz))(0.3*(1+zz)**3/Ez(zz)**2,0.7/Ez(zz)**2)
    return g(z)/g(0.0)
def dndlnM(M,z,frozen):
    sig=0.954*(1.5*M/1.43e14)**-0.24*(1.0 if frozen else growth(z)); a,pp,dc=0.707,0.3,1.686; nu2=a*dc*dc/sig**2
    return (1/M)*(1+nu2**-pp)*np.sqrt(nu2)*np.exp(-nu2/2)
lM=np.linspace(np.log(0.8e14),np.log(6e15),700); Mg=np.exp(lM)
def shadow075(M,z,model):
    E=1.0 if model=="static" else Ez(z); dist=d_s(z) if model=="static" else DA(z)
    y0=4.95e-5*E*E*(M/3e14)**1.08; R=(3*M/(4*np.pi*500*1.3598e11*E*E))**(1/3); th=R/dist*am
    return y0*th*th*Kc(0.75/th), th
D=[l.split() for l in open(dd+"/bulbul_spt.txt") if not l.startswith("#")]
keep=[r for r in D if r[0] not in ("0330-5228","2017-6258")]
zs=np.array([float(r[6]) for r in keep]); xi_d=np.array([float(r[2]) for r in keep]); eY=np.array([float(r[5]) for r in keep])*1e-6; Yd=np.array([float(r[4]) for r in keep])*1e-6
sigY=np.median(eY); cxi=np.median(xi_d/(Yd/eY))
print(f"data: {len(zs)} clusters; SPT noise on Y(0.75'): median {sigY*1e6:.0f}e-6 arcmin2; xi/(Y/noise) median {cxi:.2f}; lowest xi {xi_d.min():.2f}; {int((xi_d<8).sum())} below 8")
A=7.2*sigY/cxi/shadow075(4.3e14,0.6,"expanding")[0]
ols=lambda X,y: np.linalg.lstsq(X,y,rcond=None)[0]
rng=np.random.default_rng(32)
def pool(z,model,xi_min,s_gas,s_int,s_x,npool=250):
    w=dndlnM(Mg,z,model=="static"); w=w/w.sum(); out=[]
    while len(out)<npool:
        M=Mg[rng.choice(len(Mg),40000,p=w)]; S,th=shadow075(M,z,model); gas=10**rng.normal(0,s_gas,len(M))
        Ytrue=A*S*gas*10**rng.normal(0,s_int,len(M)); Yobs=Ytrue+rng.normal(0,sigY,len(M)); xi=cxi*Yobs/sigY; sel=xi>=xi_min
        Xpred=A*S*gas*10**rng.normal(0,s_x,len(M))/(F_log(z) if model=="static" else 1.0)
        for j in np.flatnonzero(sel)[:npool-len(out)]: out.append((np.log10(Yobs[j]/Xpred[j]),xi[j]))
    return np.array(out)
def run(model,xi_min=4.7,s_gas=0.07,s_int=0.05,s_x=0.10,nrep=400):
    pools=[pool(z,model,xi_min,s_gas,s_int,s_x) for z in zs]; x=np.log10(F_log(zs)); res=[]
    for _ in range(nrep):
        pick=np.array([P[rng.integers(len(P))] for P in pools]); lr,xi=pick.T; row=[]
        for cut in (0,6,8):
            m=xi>=cut; row.append(ols(np.c_[np.ones(m.sum()),x[m]],lr[m])[1] if m.sum()>8 else np.nan); row.append(m.sum())
        row.append(10**np.median(lr)); row.append(np.median(xi)); res.append(row)
    return np.array(res)
out={}
for model in ("expanding","static"):
    for lab,kw in (("baseline",{}),("more scatter",dict(s_int=0.08,s_gas=0.10)),("threshold 5.5",dict(xi_min=5.5))):
        r=run(model,**kw); mu=np.nanmean(r,0); sd=np.nanstd(r,0)
        out[f"{model}|{lab}"]=dict(k=[float(mu[0]),float(sd[0])],k6=[float(mu[2]),float(sd[2]),float(mu[3])],k8=[float(mu[4]),float(sd[4]),float(mu[5])],ratio=float(mu[6]),xi_med=float(mu[7]))
        print(f"  simulated {model:9s} {lab:14s}: k = {mu[0]:+.2f}+-{sd[0]:.2f} | xi>=6 ({mu[3]:.0f}): {mu[2]:+.2f}+-{sd[2]:.2f} | xi>=8 ({mu[5]:.0f}): {mu[4]:+.2f}+-{sd[4]:.2f} | median ratio {mu[6]:.2f}, median xi {mu[7]:.1f}")
json.dump(out,open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"selection_sim_xmm.json"),"w"),indent=1)
