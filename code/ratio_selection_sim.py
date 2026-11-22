"""Workbench Entry 28: does the way the clusters were chosen bend the shadow-against-glow test?
The SPT clusters were chosen because their shadow was strong (significance xi above about 6.5), and the shadow measured from the same map is what
goes into the ratio. Noise that helps a cluster over the threshold also raises its ratio, and it does so more where clusters are rarer.
This simulates the SPT sample under each hypothesis and passes it through the same fit as the data.
Ingredients: Sheth-Tormen mass function (evolving for the expanding model, frozen at today's for the static one); universal pressure profile for the
shadow inside 0.75 arcmin; SPT noise from the catalogue; scatter between shadow and X-ray pressure; X-ray measurement error."""
import numpy as np, json, sys
from scipy import integrate
dd=sys.argv[1] if len(sys.argv)>1 else "."
cH=299792.458/70.; Ez=lambda z: np.sqrt(0.3*(1+z)**3+0.7)
DA=np.vectorize(lambda z: cH*integrate.quad(lambda x:1/Ez(x),0,z)[0]/(1+z)); d_s=lambda z: cH*np.log(1+z)
F_log=lambda z: np.sqrt(d_s(z)/DA(z))*(1+z)**-1.5
am=180*60/np.pi
# pressure-profile column integral inside a/theta500, in units of y0 theta500^2
c5,ga,al,be=1.177,0.3081,1.0510,5.4905
p=lambda x: 1/((c5*x)**ga*(1+(c5*x)**al)**((be-ga)/al))
los=lambda b: 2*integrate.quad(lambda l:p(np.hypot(b,l)),0,np.sqrt(64-b*b),limit=200)[0]
I0=los(1e-6); ag=np.geomspace(0.03,4,60); Kg=np.array([integrate.quad(lambda b:2*np.pi*b*los(b)/I0,0,a,limit=200)[0] for a in ag])
Kc=lambda a: np.exp(np.interp(np.log(a),np.log(ag),np.log(Kg)))
# mass function shape
def growth(z):
    g=lambda zz: (lambda om,ol: 2.5*om/(om**(4/7)-ol+(1+om/2)*(1+ol/70))/(1+zz))(0.3*(1+zz)**3/Ez(zz)**2,0.7/Ez(zz)**2)
    return g(z)/g(0.0)
def dndlnM(M,z,frozen):
    sig=0.954*(1.5*M/1.43e14)**-0.24*(1.0 if frozen else growth(z)); a,pp,dc=0.707,0.3,1.686; nu2=a*dc*dc/sig**2
    return (1/M)*(1+nu2**-pp)*np.sqrt(nu2)*np.exp(-nu2/2)
lM=np.linspace(np.log(1.2e14),np.log(6e15),700); Mg=np.exp(lM)
def shadow075(M,z,model):
    E=1.0 if model=="static" else Ez(z); dist=d_s(z) if model=="static" else DA(z)
    y0=4.95e-5*E*E*(M/3e14)**1.08; R=(3*M/(4*np.pi*500*1.3598e11*E*E))**(1/3); th=R/dist*am
    return y0*th*th*Kc(0.75/th), th
D=[l.split() for l in open(dd+"/sanders_spt.txt") if not l.startswith("#")]
zs=np.array([float(r[1]) for r in D]); zb=np.array([float(r[6]) for r in D]); ok=np.abs(zs-zb)<=0.03*(1+zs); zs=zs[ok]; xi_d=np.array([float(r[2]) for r in D])[ok]
eY=np.array([float(r[5]) for r in D])[ok]*1e-6; Yd=np.array([float(r[4]) for r in D])[ok]*1e-6
print(f"data: {len(zs)} clusters; SPT noise on Y(0.75'): median {np.median(eY)*1e6:.0f}e-6 arcmin2; xi/(Y/noise) median {np.median(xi_d/(Yd/eY)):.2f}; lowest xi {xi_d.min():.1f}, 10th percentile {np.percentile(xi_d,10):.1f}")
sigY=np.median(eY); cxi=np.median(xi_d/(Yd/eY))
# one normalisation for the shadow, fixed so that an expanding-model cluster of 4.3e14 Msun at z = 0.6 sits at xi = 7.2 (SPT's own mass scale)
A=7.2*sigY/cxi/shadow075(4.3e14,0.6,"expanding")[0]
def ols(X,y):
    b=np.linalg.lstsq(X,y,rcond=None)[0]; return b
rng=np.random.default_rng(4)
def pool(z,model,xi_min,s_gas,s_int,s_x,s_T=0.0,npool=300):
    w=dndlnM(Mg,z,model=="static"); w=w/w.sum(); out=[]
    while len(out)<npool:
        M=Mg[rng.choice(len(Mg),20000,p=w)]; S,th=shadow075(M,z,model); gas=10**rng.normal(0,s_gas,len(M))
        Ytrue=A*S*gas*10**rng.normal(0,s_int,len(M)); Yobs=Ytrue+rng.normal(0,sigY,len(M)); xi=cxi*Yobs/sigY; sel=xi>=xi_min
        E=1.0 if model=="static" else Ez(z); eT=rng.normal(0,s_T,len(M)) if s_T>0 else 0.0     # eT: X-ray temperature above or below the pressure-weighted one; it enters the measured T and the X-ray prediction alike
        T=5*(M*E/4.31e14)**(1/1.53)*10**(rng.normal(0,np.sqrt(max(0.06**2-s_T**2,1e-6)),len(M))+eT)
        Xpred=A*S*gas*10**(rng.normal(0,np.sqrt(max(s_x**2-s_T**2,1e-6)),len(M))+eT)/(F_log(z) if model=="static" else 1.0)     # what a standard-distance X-ray analysis would predict
        for j in np.flatnonzero(sel)[:npool-len(out)]: out.append((np.log10(Yobs[j]/Xpred[j]),xi[j],T[j],M[j]))
    return np.array(out)
def run(model,xi_min=6.8,s_gas=0.07,s_int=0.05,s_x=0.07,s_T=0.0,nrep=400):
    pools=[pool(z,model,xi_min,s_gas,s_int,s_x,s_T) for z in zs]; x=np.log10(F_log(zs)); res=[]
    for _ in range(nrep):
        pick=np.array([P[rng.integers(len(P))] for P in pools]); lr,xi,T,M=pick.T; row=[]
        for cut in (0,8,10):
            m=xi>=cut; row.append(ols(np.c_[np.ones(m.sum()),x[m]],lr[m])[1] if m.sum()>8 else np.nan); row.append(m.sum())
        bT=ols(np.c_[np.ones(len(x)),x,np.log10(T/6)],lr); row.append(bT[1]); row.append(bT[2]); row.append(np.median(T)); row.append(np.median(M[zs>0.75])/1e14); row.append(10**np.median(lr)); res.append(row)
    return np.array(res)
out={"norm":float(A),"sigY":float(sigY),"cxi":float(cxi)}
print("simulated samples, same redshifts as the data, same fit (mean and spread over 400 samples):")
print("  case                                   k (all)        k (xi>=8)      k (xi>=10)     k with T term   median kT   median mass z>0.75")
for model in ("expanding","static"):
    for label,kw in (("baseline",{}),("more scatter (0.08 shadow, 0.10 gas)",dict(s_int=0.08,s_gas=0.10)),("less scatter (0.03 shadow)",dict(s_int=0.03)),("threshold xi 6.0",dict(xi_min=6.0)),("threshold xi 7.5",dict(xi_min=7.5)),("X-ray error 0.12 dex",dict(s_x=0.12)),("temperature error shared, 0.05 dex",dict(s_T=0.05)),("temperature error shared, 0.06 dex",dict(s_T=0.06,s_x=0.08))):
        r=run(model,**kw); mm=np.nanmean(r,0); ss=np.nanstd(r,0); key=model+"|"+label
        out[key]=dict(k=[float(mm[0]),float(ss[0])],k8=[float(mm[2]),float(ss[2]),float(mm[3])],k10=[float(mm[4]),float(ss[4]),float(mm[5])],kT=[float(mm[6]),float(ss[6])],Tslope=[float(mm[7]),float(ss[7])],Tmed=float(mm[8]),Mhi=float(mm[9]),ratio=float(mm[10]))
        print(f"  {model:9s} {label:28s} {mm[0]:+.2f}+-{ss[0]:.2f}   {mm[2]:+.2f}+-{ss[2]:.2f} (n~{mm[3]:.0f})   {mm[4]:+.2f}+-{ss[4]:.2f} (n~{mm[5]:.0f})   {mm[6]:+.2f}+-{ss[6]:.2f} (T power {mm[7]:+.2f}+-{ss[7]:.2f})    {mm[8]:4.1f} keV   {mm[9]:4.1f}e14")
json.dump(out,open("ratio_selection_sim.json","w"),indent=1)
