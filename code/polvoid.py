"""Entry 147: how much of a void-scale CMB temperature spot can E-mode polarization predict?
Uses the notebook's own semi-analytic TT/EE/TE spectra with reionization at tau = 0.054 (Entry 94, cmbreion_tau054.json).
That run stores E in the Ma-Bertschinger normalisation; as in Entry 94 (reion_compare.py) EE is divided by 16 and TE by 4,
matching the Entry 92 run (EE4/TT1 = 0.00778 vs Planck 0.00743). r_l = TE/sqrt(TT EE) does not depend on that factor.
Absolute scale: first TT peak = 5,733 uK^2 (Planck)."""
import json, numpy as np, sys
from scipy.special import eval_legendre
d=json.load(open(sys.argv[1])); run=list(d['runs'].values())[0]
L0=np.array(d['ells'],float); TT1=5733.0
TT=np.array(run['TT'])*TT1; EE=np.array(run['EE'])*TT1/16; TE=np.array(run['TE'])*TT1/4
ell=np.arange(2,601)
iD=lambda y: np.interp(ell,L0,y)
DTT,DEE,DTE=iD(TT),iD(EE),iD(TE)
pre=2*np.pi/(ell*(ell+1)); CTT,CEE,CTE=DTT*pre,DEE*pre,DTE*pre
arcmin=np.pi/180/60
def NEE(w_uKarcmin):  # white noise, beam irrelevant at l < 200
    return (w_uKarcmin*arcmin)**2*np.ones_like(ell,float)
def W_smhw(R_deg):
    x=ell*np.radians(R_deg); return x**2*np.exp(-x**2/2)
def W_disk(th_deg):
    c=np.cos(np.radians(th_deg)); return (eval_legendre(ell-1,c)-eval_legendre(ell+1,c))/((2*ell+1)*(1-c))
def W_cdisk(th_deg):  # disk minus equal-area surrounding ring (compensated top hat)
    c1=np.cos(np.radians(th_deg)); c2=2*c1-1  # equal area ring: 1-c2 = 2(1-c1)
    A1=1-c1; Wd=(eval_legendre(ell-1,c1)-eval_legendre(ell+1,c1))/((2*ell+1)*A1)
    Wb=(eval_legendre(ell-1,c2)-eval_legendre(ell+1,c2))/((2*ell+1)*(1-c2))
    Wring=(Wb*(1-c2)-Wd*A1)/((1-c2)-A1); return Wd-Wring
def stats(W, isw=0.0):
    g=(2*ell+1)/(4*np.pi)*W**2
    ctt=CTT*(1+isw*np.exp(-(ell-2)/10.0))   # optional late-ISW power, uncorrelated with E, fading by l ~ 30
    var=np.sum(g*ctt)
    out={'sigma_uK':float(np.sqrt(var))}
    for name,w in (('ideal',0.0),('LiteBIRD_2',2.0),('Planck_40',40.0),('Planck_60',60.0)):
        out['F_'+name]=float(np.sum(g*CTE**2/(CEE+NEE(w)))/var)
    # l range carrying half the filtered variance
    cum=np.cumsum(g*ctt)/var; out['l_med']=int(ell[np.searchsorted(cum,0.5)]); out['l_16_84']=[int(ell[np.searchsorted(cum,0.16)]),int(ell[np.searchsorted(cum,0.84)])]
    return out
res={}
r=lambda l: float(DTE[l-2]/np.sqrt(DTT[l-2]*DEE[l-2]))
res['r_l']={l:round(r(l),3) for l in (2,5,10,15,20,25,30,40,50,60,80,100,150)}
res['D_check']={'TT10':float(DTT[8]),'TT20':float(DTT[18]),'EE5':float(DEE[3]),'EE20':float(DEE[18]),'EE100':float(DEE[98]),'TE5':float(DTE[3]),'TE20':float(DTE[18])}
for lab,W in (('SMHW_4deg',W_smhw(4)),('SMHW_2deg',W_smhw(2)),('disk_5deg',W_disk(5)),('cdisk_5deg',W_cdisk(5)),('cdisk_2.5deg',W_cdisk(2.5))):
    res[lab]=stats(W); res[lab+'_isw15']=stats(W,isw=0.15)
S=np.array([3.5,4.0])
for lab in ('SMHW_4deg','cdisk_5deg'):
    for k in ('ideal','LiteBIRD_2','Planck_40','Planck_60'):
        F=res[lab]['F_'+k]; res[lab]['Z_'+k]=[float(np.sqrt(F)*s) for s in S]
print(json.dumps(res,indent=1))
json.dump(res,open('polvoid.json','w'),indent=1)

# ---- Part 2: what the claim's significance implies, and the next-slab replication ----
s4=res['SMHW_4deg']['sigma_uK']; s2=res['SMHW_2deg']['sigma_uK']; dT=35.0
Neff=[(S_*s4/dT)**2 for S_ in (3.5,4.0)]
Nnext=[4*n for n in Neff]
sig_next=[s2/np.sqrt(n) for n in Nnext]
rep=dict(sigma_4deg=s4,sigma_2deg=s2,N_eff=Neff,N_next=Nnext,sigma_stack_next=sig_next,Z_column=[dT/s for s in sig_next],
         upper3_next=[3*s for s in sig_next])
# Planck noise vs E signal at l = 20
for w in (40,60):
    rep[f'NEE_D20_{w}']=float(20*21/(2*np.pi)*(w*arcmin)**2)
rep['EE_D20']=float(DEE[18])
res['replication']=rep
# figure data: r_l and the 4-degree wavelet's variance weight
lf=np.arange(2,151); g=(2*lf+1)*W_smhw(4)[:149]**2*CTT[:149]; g=g/g.max()
res['fig']=dict(l=lf.tolist(),r=[round(float(DTE[l-2]/np.sqrt(DTT[l-2]*DEE[l-2])),3) for l in lf],w=[round(float(x),3) for x in g],
                bars=[dict(lab='Polarization, Planck noise',lo=round(min(res['SMHW_4deg']['Z_Planck_60']),2),hi=round(max(res['SMHW_4deg']['Z_Planck_40']),2),kind='pol'),
                      dict(lab='Polarization, noise-free E',lo=round(min(res['SMHW_4deg']['Z_ideal']),2),hi=round(max(res['SMHW_4deg']['Z_ideal']),2),kind='pol'),
                      dict(lab='Next slab, z = 0.03–0.06',lo=round(min(rep['Z_column']),1),hi=round(max(rep['Z_column']),1),kind='rep')])
print(json.dumps(rep,indent=1)); print(json.dumps(res['fig']['bars']))
json.dump(res,open('polvoid.json','w'),indent=1)
