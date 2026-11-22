"""Entry 148: a frozen void catalogue for the deeper-slab test of the hot-void claim, and its forecast.
Voids: found in the browser from the 2M++ catalogue (Lavaux & Hudson 2011; VizieR J/MNRAS/416/2840), see void_finder.js.
Slab 1: 0.002 < z < 0.03, K < 11.5, M_K < -24.09 (volume-limited), |b| >= 10 deg.
Slab 2: 0.03 <= z < 0.06, K < 12.5 regions (SDSS, 6dFGRS), M_K < -24.57, |b| >= 10 deg.
Catalogue columns: l, b, radius (deg), Delta(<radius), N, N expected, log10 Poisson P(<=N) [, inner disk overlaps a slab-1 inner disk].
CMB scatter per void: compensated top hat of radius theta/2, from the notebook's spectra (Entry 94 run, as in Entry 147)."""
import json, numpy as np, sys
from scipy.special import eval_legendre
def load(f): return np.array([list(map(float,r.split(','))) for r in open(f).read().strip().split(';')])
S1=load('voids_slab1.txt'); S2=load('voids_slab2.txt')
d=json.load(open(sys.argv[1])); run=list(d['runs'].values())[0]
L0=np.array(d['ells'],float); TT=np.array(run['TT'])*5733.0
ell=np.arange(2,1201); DTT=np.interp(ell,L0,TT); CTT=DTT*2*np.pi/(ell*(ell+1))
def W_cdisk(th):
    c1=np.cos(np.radians(th)); c2=2*c1-1; A1=1-c1
    Wd=(eval_legendre(ell-1,c1)-eval_legendre(ell+1,c1))/((2*ell+1)*A1)
    Wb=(eval_legendre(ell-1,c2)-eval_legendre(ell+1,c2))/((2*ell+1)*(1-c2))
    return Wd-(Wb*(1-c2)-Wd*A1)/((1-c2)-A1)
_cache={}
def sig(th):
    k=round(th,2)
    if k not in _cache: W=W_cdisk(th); _cache[k]=float(np.sqrt(np.sum((2*ell+1)/(4*np.pi)*W**2*CTT)))
    return _cache[k]
out={}
for name,S in (('slab1',S1),('slab2',S2),('slab2_fresh',S2[S2[:,7]==0]),('slab2_overlap',S2[S2[:,7]==1])):
    th=S[:,2]; s=np.array([sig(t/2) for t in th]); N=len(S)
    q=th>=np.percentile(th,75)
    st=np.sqrt((s**2).sum())/N; stq=np.sqrt((s[q]**2).sum())/q.sum()
    w=1/s**2; siv=1/np.sqrt(w.sum())
    out[name]=dict(N=int(N),theta_med=float(np.median(th)),theta_q=[float(x) for x in np.percentile(th,[25,75])],
                   Delta_med=float(np.median(S[:,3])),sigma_void_med=float(np.median(s)),
                   sigma_stack=float(st),Z35=float(35/st),upper3=float(3*st),
                   Nq=int(q.sum()),sigma_stack_q=float(stq),Z35_q=float(35/stq),sigma_iv=float(siv))
for k,v in out.items(): print(k,{kk:(round(vv,2) if isinstance(vv,float) else vv) for kk,vv in v.items()})
# sigma per radius table for the figure
out['sig_curve']=[[t,round(sig(t),2)] for t in (1.5,2,2.5,3,4,5,6,7.5,9,10.5,12.5)]
out['cat1']=S1.tolist(); out['cat2']=S2.tolist()
json.dump(out,open('prereg.json','w'))

# ---- Statistic that keeps a profile falling to zero at the void edge: inner disk (theta/2) minus outer annulus (theta..1.5 theta) ----
def Wcap(th):
    c=np.cos(np.radians(th)); return (eval_legendre(ell-1,c)-eval_legendre(ell+1,c))/((2*ell+1)*(1-c)), 1-c
def W_io(th):
    Wd,_=Wcap(th/2); W1,A1=Wcap(th); W2,A2=Wcap(1.5*th)
    return Wd-(A2*W2-A1*W1)/(A2-A1)
_c2={}
def sig_io(th):
    k=round(th,2)
    if k not in _c2: W=W_io(th); _c2[k]=float(np.sqrt(np.sum((2*ell+1)/(4*np.pi)*W**2*CTT)))
    return _c2[k]
# what the compensated (equal-area ring) statistic would retain of a profile T0(1-(r/R)^2): disk mean 0.875 T0, ring mean 0.625 T0
retain_comp=(0.875-0.625)/0.875
io={}
for name,S in (('slab1',S1),('slab2_fresh',S2[S2[:,7]==0]),('slab2',S2)):
    th=S[:,2]; s=np.array([sig_io(t) for t in th]); N=len(S); q=th>=np.percentile(th,75)
    st=np.sqrt((s**2).sum())/N; stq=np.sqrt((s[q]**2).sum())/q.sum()
    io[name]=dict(N=int(N),sigma_void_med=float(np.median(s)),sigma_stack=float(st),Z35=float(35/st),Z30=float(30/st),upper3=float(3*st),Nq=int(q.sum()),sigma_stack_q=float(stq),Z35_q=float(35/stq),upper3_q=float(3*stq))
    print('io',name,{kk:(round(vv,2) if isinstance(vv,float) else vv) for kk,vv in io[name].items()})
print('compensated statistic keeps',round(retain_comp,3),'of the inner value for a profile falling to zero at R')
D=json.load(open('prereg.json')); D['io']=io; D['retain_comp']=retain_comp
D['sig_io_curve']=[[t,round(sig_io(t),2)] for t in (3,4,5,6,7.5,9,10.5,12.5,15,18,21,25)]
json.dump(D,open('prereg.json','w'))
