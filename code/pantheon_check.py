"""Workbench Entry 31, check: the static Hubble laws against the Pantheon supernova sample (Scolnic et al. 2018; 1,048 supernovae, z = 0.010 to 2.26).
Input: redshift and corrected peak magnitude of every supernova (pantheon_rows.py, read from VizieR J/ApJ/859/101/fullz).
Errors here are STATISTICAL ONLY: the scatter of the supernovae themselves, measured in redshift slices about a flexible smooth curve.
The survey-to-survey calibration errors (the 'systematic covariance') are not in this table, so the chi-square differences below overstate the significance;
the robustness section asks instead how large a calibration error would have to be to rescue each static law.
Every model gets a free additive constant (supernova luminosity and H0 are not known separately)."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy import integrate, optimize
from pantheon_rows import z, m, e_low
N=len(z); rng=np.random.default_rng(31)
# ---- 1. empirical scatter about a flexible curve (5 log z + quartic in log(1+z)) ----
x=np.log10(1+z); base=5*np.log10(z)
def flex(w):
    A=np.vander(x,5); c=np.linalg.lstsq(A*w[:,None],(m-base)*w,rcond=None)[0]; return m-base-A@c
edges=np.quantile(z,np.linspace(0,1,15)); mid=np.array([np.median(z[(z>=a)&(z<=b)]) for a,b in zip(edges[:-1],edges[1:])])
sig=np.full(N,0.15)
for it in range(4):
    r=flex(1/sig); s=np.array([1.4826*np.median(np.abs(r[(z>=a)&(z<=b)]-np.median(r[(z>=a)&(z<=b)]))) for a,b in zip(edges[:-1],edges[1:])])
    srms=np.array([r[(z>=a)&(z<=b)].std(ddof=1) for a,b in zip(edges[:-1],edges[1:])])
    sig=np.interp(np.log(z),np.log(mid),srms)
r=flex(1/sig); chi_flex=float(np.sum((r/sig)**2))
print(f"{N} supernovae; scatter about a smooth curve in 14 redshift slices (mag): "+" ".join(f"{v:.3f}" for v in srms))
print(f"   slice median redshifts: "+" ".join(f"{v:.3f}" for v in mid))
print(f"   catalogue's own error for the 319 supernovae below z = 0.16: median {np.median(e_low):.3f}; the measured scatter there is {np.sqrt(np.mean(srms[:4]**2)):.3f}")
print(f"   flexible curve (5 free numbers): chi2 = {chi_flex:.0f} for {N} supernovae\n")
w=1/sig**2
def chi2(mu_model,sel=None,step=None):
    """chi2 with a free constant; optional free step for supernovae with z below `step`."""
    s=np.ones(N,bool) if sel is None else sel; d=(m-mu_model)[s]; ww=w[s]
    if step is None:
        a=np.sum(ww*d)/np.sum(ww); return float(np.sum(ww*(d-a)**2)),float(a),0.0
    lo=z[s]<step; a1=np.sum(ww[lo]*d[lo])/np.sum(ww[lo]); a2=np.sum(ww[~lo]*d[~lo])/np.sum(ww[~lo])
    return float(np.sum(ww[lo]*(d[lo]-a1)**2)+np.sum(ww[~lo]*(d[~lo]-a2)**2)),float(a2),float(a1-a2)
mu=lambda dl: 5*np.log10(dl)
E=lambda zz,om: np.sqrt(om*(1+zz)**3+1-om)
zg=np.linspace(0,2.4,2401)
def dl_lcdm(om,zz=z):
    dc=integrate.cumulative_trapezoid(1/E(zg,om),zg,initial=0); return (1+zz)*np.interp(zz,zg,dc)
LAWS={"flat LCDM, matter density fitted":None,
      "static, Lerner: d = cz/H0, energy loss only":lambda t: t*np.sqrt(1+t),
      "static, d = cz/H0, photons also stretched in time":lambda t: t*(1+t),
      "static, tired light: d = (c/H0) ln(1+z), energy loss only":lambda t: np.log(1+t)*np.sqrt(1+t),
      "static, d = (c/H0) ln(1+z), photons also stretched in time":lambda t: np.log(1+t)*(1+t),
      "empty expanding universe (coasting)":lambda t: t*(1+t/2),
      "expanding, matter only (no dark energy)":lambda t: 2*(1+t)*(1-1/np.sqrt(1+t))}
def fit_om(sel=None,step=None):
    f=lambda o: chi2(mu(dl_lcdm(o)),sel,step)[0]; g=optimize.minimize_scalar(f,bounds=(0.05,0.9),method="bounded")
    lo=optimize.brentq(lambda o:f(o)-g.fun-1,0.05,g.x); hi=optimize.brentq(lambda o:f(o)-g.fun-1,g.x,0.9); return g.x,lo,hi,g.fun
om,om_lo,om_hi,c_l=fit_om()
print(f"flat LCDM, statistical errors only: matter density {om:.3f} (+{om_hi-om:.3f}, -{om-om_lo:.3f}); chi2 = {c_l:.0f}")
out={"N":N,"om":[om,om_lo,om_hi],"chi_flex":chi_flex,"slices":{"z":mid.tolist(),"rms":srms.tolist()},"models":{}}
def table(sel=None,step=None,label="all supernovae"):
    o,_,_,c0=fit_om(sel,step); n=N if sel is None else int(sel.sum()); res={}
    print(f"\n--- {label}: {n} supernovae"+(f", free step at z = {step}" if step else "")+f"; LCDM matter density {o:.3f}, chi2 = {c0:.0f} ---")
    for name,fn in LAWS.items():
        if fn is None: c,a,st=chi2(mu(dl_lcdm(o)),sel,step)
        else: c,a,st=chi2(mu(fn(z)),sel,step)
        print(f"   {name:60s} chi2 {c:8.0f}   delta {c-c0:+8.0f}"+(f"   step needed {st:+.3f} mag" if step else "")); res[name]=dict(chi2=c,dchi2=c-c0,step=st,offset=a)
    return res
out["models"]["all"]=table()
# ---- 2. the dimming exponent: d_L = d (1+z)^q ----
def fit_q(dfun,sel=None,step=None):
    f=lambda q: chi2(mu(dfun(z)*(1+z)**q),sel,step)[0]; g=optimize.minimize_scalar(f,bounds=(-1,3),method="bounded")
    lo=optimize.brentq(lambda q:f(q)-g.fun-1,-1,g.x); hi=optimize.brentq(lambda q:f(q)-g.fun-1,g.x,3); return g.x,g.x-lo,hi-g.x,g.fun
q1=fit_q(lambda t:t); q2=fit_q(lambda t:np.log(1+t))
print(f"\nd = cz/H0, brightness falling as an extra (1+z)^(2q): q = {q1[0]:.3f} (+{q1[2]:.3f}, -{q1[1]:.3f}), chi2 {q1[3]:.0f}; energy loss only is 0.5, with stretched time 1.0")
print(f"d = (c/H0) ln(1+z):                                  q = {q2[0]:.3f} (+{q2[2]:.3f}, -{q2[1]:.3f}), chi2 {q2[3]:.0f}")
out["q"]=list(q1); out["q_log"]=list(q2)
# ---- 3. robustness: what would have to be wrong with the data ----
out["models"]["step"]=table(step=0.1,label="free calibration step between supernovae below and above z = 0.1")
out["models"]["z>0.1"]=table(sel=z>0.1,label="only supernovae beyond z = 0.1 (no nearby anchor)")
out["models"]["z<0.7"]=table(sel=z<0.7,label="only supernovae nearer than z = 0.7")
out["models"]["z>0.3"]=table(sel=z>0.3,label="only supernovae beyond z = 0.3")
# a smooth drift: add b*z to the magnitudes (a calibration or evolution error growing with redshift); how big must b be for each static law?
def drift(fn):
    f=lambda b: chi2(mu(fn(z))+b*z)[0]; g=optimize.minimize_scalar(f,bounds=(-2,2),method="bounded"); return g.x,g.fun
print("\nA drift in supernova brightness proportional to redshift (magnitudes per unit z) that each law would need, and the chi2 left:")
out["drift"]={}
for name,fn in LAWS.items():
    if fn is None: continue
    b,c=drift(fn); print(f"   {name:60s} {b:+.3f} mag per unit z   chi2 {c:8.0f}  (delta {c-c_l:+.0f})"); out["drift"][name]=[b,c,c-c_l]
bl,cl2=optimize.minimize_scalar(lambda b: chi2(mu(dl_lcdm(om))+b*z)[0],bounds=(-2,2),method="bounded").x,None
print(f"   (the same drift fitted on top of flat LCDM: {bl:+.3f} mag per unit z)")
# ---- 4. scatter-only error model check: bootstrap of the chi2 difference Lerner minus LCDM, and of q ----
dL=[];qs=[]
mu_l=mu(dl_lcdm(om)); mu_s=mu(LAWS["static, Lerner: d = cz/H0, energy loss only"](z))
for b in range(2000):
    i=rng.integers(0,N,N); ww=w[i]
    def c(mm):
        d=(m-mm)[i]; a=np.sum(ww*d)/np.sum(ww); return np.sum(ww*(d-a)**2)
    dL.append(c(mu_s)-c(mu_l))
dL=np.array(dL); print(f"\nBootstrap (2,000 resamplings of the supernovae): chi2(Lerner) - chi2(LCDM at {om:.3f}) = {np.median(dL):.0f}, 68% range {np.percentile(dL,16):.0f} to {np.percentile(dL,84):.0f}, lowest {dL.min():.0f}")
out["boot"]=[float(np.median(dL)),float(np.percentile(dL,16)),float(np.percentile(dL,84)),float(dL.min())]
# ---- 5. binned residuals for the figure and the table ----
be=np.r_[np.quantile(z[z<1.0],np.linspace(0,1,21))[:-1],1.0,1.25,1.45,2.3]; rows=[]
a_l=chi2(mu_l)[1]
for a,b in zip(be[:-1],be[1:]):
    s=(z>=a)&(z<b) if b<2.3 else (z>=a); ww=w[s]; zc=float(np.exp(np.sum(ww*np.log(z[s]))/ww.sum())); d=float(np.sum(ww*(m-mu_l-a_l)[s])/ww.sum()); e=float(1/np.sqrt(ww.sum())); rows.append([zc,d,e,int(s.sum())])
out["bins"]=rows
print("\nBinned residuals from the best LCDM curve: z, mean residual (mag), error, count")
for rrow in rows: print(f"   {rrow[0]:.3f}  {rrow[1]:+.3f} +- {rrow[2]:.3f}  ({rrow[3]})")
zz=np.geomspace(0.01,2.3,80); ref=mu(dl_lcdm(om,zz))+a_l; out["curves"]={}
for key,name in (("lerner","static, Lerner: d = cz/H0, energy loss only"),("lerner_td","static, d = cz/H0, photons also stretched in time"),("log_td","static, d = (c/H0) ln(1+z), photons also stretched in time"),("coast","empty expanding universe (coasting)")):
    out["curves"][key]=[[float(a),float(b)] for a,b in zip(zz,mu(LAWS[name](zz))+out["models"]["all"][name]["offset"]-ref)]
    pick=np.interp([0.03,0.5,1.0,1.5,2.26],zz,np.array(out["curves"][key])[:,1]); print(f"{name:60s} minus LCDM at z = 0.03, 0.5, 1.0, 1.5, 2.26: "+" ".join(f"{v:+.3f}" for v in pick))
# the 12 supernovae beyond z = 1.3
s=z>1.3; print(f"\n{int(s.sum())} supernovae beyond z = 1.3: mean residual from LCDM {np.sum(w[s]*(m-mu_l-a_l)[s])/w[s].sum():+.3f} +- {1/np.sqrt(w[s].sum()):.3f}; from Lerner's law {np.sum(w[s]*(m-mu_s-out['models']['all']['static, Lerner: d = cz/H0, energy loss only']['offset'])[s])/w[s].sum():+.3f}")
# ---- 6. the same 'no nearby anchor' question asked of the JLA binned distances, which carry their full covariance (statistical and systematic) ----
import io, contextlib
_cwd=os.getcwd(); os.chdir(os.path.dirname(os.path.abspath(__file__)))
with contextlib.redirect_stdout(io.StringIO()): import static_hubble as sh
os.chdir(_cwd)
def jla(sel):
    zb=sh.zb[sel]; C=sh.C[np.ix_(sel,sel)]; Ci=np.linalg.inv(C); one=np.ones(len(zb)); y=sh.mub[sel]
    def c2(mm):
        d=y-mm; a=(one@Ci@d)/(one@Ci@one); r=d-a; return float(r@Ci@r)
    dl=lambda o: np.array([(1+t)*integrate.quad(lambda u:1/E(u,o),0,t)[0] for t in zb])
    g=optimize.minimize_scalar(lambda o:c2(mu(dl(o))),bounds=(0.02,0.95),method="bounded")
    fq=lambda q:c2(mu(zb*(1+zb)**q)); gq=optimize.minimize_scalar(fq,bounds=(-1,3),method="bounded")
    qlo=optimize.brentq(lambda q:fq(q)-gq.fun-1,-1,gq.x); qhi=optimize.brentq(lambda q:fq(q)-gq.fun-1,gq.x,3)
    return dict(n=int(len(zb)),om=float(g.x),chi_lcdm=float(g.fun),chi_lerner=c2(mu(zb*np.sqrt(1+zb))),chi_log_td=c2(mu(np.log(1+zb)*(1+zb))),chi_coast=c2(mu(zb*(1+zb/2))),q=[float(gq.x),float(gq.x-qlo),float(qhi-gq.x)])
print("\nJLA binned distances with full covariance (statistical + systematic), restricted in redshift:")
out["jla"]={}
for lab,sel in (("all 31 points",sh.zb>0),("z > 0.1 (16 points)",sh.zb>0.1),("z > 0.3 (10 points)",sh.zb>0.3),("z < 0.7 (27 points)",sh.zb<0.7),("z < 0.3 (21 points)",sh.zb<0.3)):
    r=jla(sel); out["jla"][lab]=r
    print(f"   {lab:22s} n={r['n']:2d}  LCDM matter {r['om']:.3f} chi2 {r['chi_lcdm']:.1f} | Lerner {r['chi_lerner']:.1f} (delta {r['chi_lerner']-r['chi_lcdm']:+.1f}) | log+stretched {r['chi_log_td']:.1f} ({r['chi_log_td']-r['chi_lcdm']:+.1f}) | coasting {r['chi_coast']:.1f} ({r['chi_coast']-r['chi_lcdm']:+.1f}) | q = {r['q'][0]:.3f} +{r['q'][2]:.3f} -{r['q'][1]:.3f}")
# Pantheon q in the same cuts (statistical only)
out["q_cuts"]={}
for lab,sel in (("all",None),("z > 0.1",z>0.1),("z > 0.3",z>0.3),("z < 0.7",z<0.7)):
    qq=fit_q(lambda t:t,sel); out["q_cuts"][lab]=list(qq); print(f"   Pantheon, {lab:8s}: q = {qq[0]:.3f} +{qq[2]:.3f} -{qq[1]:.3f}")
# ---- 7. the two halves separately: does one dimming power fit both? (each half has its own free constant, so no calibration seam enters) ----
print("\nPantheon in two disjoint halves, each with its own constant:")
out["halves"]={}
for lab,sel in (("z < 0.3",z<0.3),("z > 0.3",z>0.3)):
    qq=fit_q(lambda t:t,sel); o=fit_om(sel); cs=chi2(mu(LAWS["static, Lerner: d = cz/H0, energy loss only"](z)),sel)[0]
    out["halves"][lab]=dict(n=int(sel.sum()),q=list(qq),om=[o[0],o[0]-o[1],o[2]-o[0]],chi_lcdm=o[3],chi_lerner=cs)
    print(f"   {lab}: {int(sel.sum())} supernovae; q = {qq[0]:.3f} +{qq[2]:.3f} -{qq[1]:.3f};  LCDM matter density {o[0]:.3f} +{o[2]-o[0]:.3f} -{o[0]-o[1]:.3f} (chi2 {o[3]:.0f});  Lerner chi2 {cs:.0f} (delta {cs-o[3]:+.0f})")
h=out["halves"]; dq=h["z < 0.3"]["q"][0]-h["z > 0.3"]["q"][0]; eq=np.hypot(h["z < 0.3"]["q"][1],h["z > 0.3"]["q"][2]); do=h["z > 0.3"]["om"][0]-h["z < 0.3"]["om"][0]; eo=np.hypot(h["z < 0.3"]["om"][2],h["z > 0.3"]["om"][1])
print(f"   the two halves disagree on q by {dq:.3f} +- {eq:.3f} ({dq/eq:.1f} sigma); on the matter density by {do:.3f} +- {eo:.3f} ({do/eo:.1f} sigma)")
out["halves"]["dq"]=[float(dq),float(eq)]; out["halves"]["dom"]=[float(do),float(eo)]
# Lerner's law levelled on the distant half only, shown against the nearby supernovae
sel=z>0.3; a_far=chi2(mu_s,sel)[1]; near=z<0.1; d_near=np.sum(w[near]*(m-mu_s-a_far)[near])/w[near].sum(); e_near=np.hypot(1/np.sqrt(w[near].sum()),1/np.sqrt(w[sel].sum()))
print(f"   Lerner's law levelled on z > 0.3: the {int(near.sum())} supernovae below z = 0.1 then lie {d_near:+.3f} +- {e_near:.3f} mag from it (negative = brighter than the law says)")
o_far=fit_om(sel); mu_f=mu(dl_lcdm(o_far[0])); a_f=chi2(mu_f,sel)[1]; d_near_l=np.sum(w[near]*(m-mu_f-a_f)[near])/w[near].sum()
print(f"   flat LCDM levelled and fitted on z > 0.3 (matter {o_far[0]:.3f}): the same nearby supernovae lie {d_near_l:+.3f} mag from it")
out["anchor"]=dict(lerner=[float(d_near),float(e_near)],lcdm=float(d_near_l),om_far=float(o_far[0]))
out["curves"]["lerner_far"]=[[float(a),float(b)] for a,b in zip(zz,mu(LAWS["static, Lerner: d = cz/H0, energy loss only"](zz))+a_far-ref)]
# the 12 most distant supernovae against LCDM with the matter density preferred beyond z = 0.1
o1=fit_om(z>0.1); mu1=mu(dl_lcdm(o1[0])); a1=chi2(mu1,z>0.1)[1]; s=z>1.3
print(f"   12 supernovae beyond z = 1.3 against LCDM fitted beyond z = 0.1 (matter {o1[0]:.3f}): {np.sum(w[s]*(m-mu1-a1)[s])/w[s].sum():+.3f} +- {1/np.sqrt(w[s].sum()):.3f}")
out["far12"]=dict(all=float(np.sum(w[s]*(m-mu_l-a_l)[s])/w[s].sum()),z01=float(np.sum(w[s]*(m-mu1-a1)[s])/w[s].sum()),err=float(1/np.sqrt(w[s].sum())))
json.dump(out,open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"pantheon_check.json"),"w"),indent=1)
