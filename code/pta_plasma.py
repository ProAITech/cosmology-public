"""Workbench Entry 117: could the nanohertz 'gravitational-wave background' be plasma?

Pulsar timing arrays see a common red signal with the Hellings-Downs (HD) correlation between pulsar pairs. The largest noise in pulsar timing is
plasma: free electrons between us and each pulsar delay the pulses by 4.149 ms x DM x (nu/GHz)^-2, and the interstellar column and the solar
wind both change with time. Three tests:
 1. Size: plasma delays against the background's own timing residuals.
 2. Colour: plasma delays depend on radio frequency; gravitational waves do not.
 3. Pattern: a foreground that changes the delay of each pulsar according to its own direction (the solar wind, or any shared screen) gives a
    correlation matrix of the separable form g(a) g(b). HD depends only on the angle between the two pulsars and is not separable. How much of HD
    can a separable pattern reproduce, at best and for the solar wind?"""
import os, json
import numpy as np
from scipy import integrate
HERE=os.path.dirname(os.path.abspath(__file__)); rng=np.random.default_rng(117)
YR=3.15576e7; FYR=1/YR; out={}
# ---- 1. sizes
def gwb_rms(A,gamma=13/3,T=15.0):
    P=lambda f: A**2/(12*np.pi**2)*FYR**-3*(f/FYR)**-gamma          # timing-residual power spectral density, s^2/Hz
    return np.sqrt(integrate.quad(P,1/(T*YR),100*FYR,limit=200)[0])
for A in (2.4e-15,2.1e-15):
    out[f"gwb_rms_ns_{A}"]=gwb_rms(A)*1e9; print(f"GWB A = {A:.1e}, gamma 13/3, 15 years: rms timing residual {gwb_rms(A)*1e9:.0f} ns")
PC_AU=4.8481e-6; K=4.149e-3                                          # pc per AU; s per (pc cm^-3) at 1 GHz
def sw_dm(eps,n0=7.0):                                               # spherical 1/r^2 wind, line of sight at elongation eps (rad)
    return n0*PC_AU*(np.pi-eps)/np.sin(eps)
for deg in (5,20,60,90,150):
    d=sw_dm(np.radians(deg)); print(f"   solar wind, n0 = 7 cm^-3, elongation {deg:3d} deg: DM = {d:.2e} pc cm^-3, delay {K*d*1e9:.0f} ns at 1 GHz, {K*d/1.4**2*1e9:.0f} ns at 1.4 GHz")
out["sw_delay_90deg_1400_ns"]=K*sw_dm(np.pi/2)/1.4**2*1e9
ism=1e-3; out["ism_delay_ns_1400"]=K*ism/1.4**2*1e9
print(f"   interstellar DM change of 1e-3 pc cm^-3 (typical over years): {K*ism/1.4**2*1e9:.0f} ns at 1.4 GHz")
# ---- 2. colour: residual of an unmodelled DM change when timing at two frequencies and fitting nu^-2
# (a statement of method, not a number: an achromatic signal survives a nu^-2 fit unchanged; a plasma signal is removed)
# ---- 3. pattern
def hd(c):
    x=(1-c)/2; x=np.clip(x,1e-12,1); return 0.5-0.25*x+1.5*x*np.log(x)
def iso(n):
    v=rng.normal(size=(n,3)); return v/np.linalg.norm(v,axis=1)[:,None]
eps_obl=np.radians(23.44)
def ecl_lat(v):                                                      # v in equatorial coordinates; ecliptic pole direction
    pole=np.array([0,-np.sin(eps_obl),np.cos(eps_obl)]); return np.arcsin(np.clip(v@pole,-1,1))
def g_sw(beta,n=720):                                                # annual mean solar-wind DM geometry factor at ecliptic latitude beta
    lam=np.linspace(0,2*np.pi,n,endpoint=False); ce=np.cos(beta)[:,None]*np.cos(lam)[None,:]; e=np.arccos(np.clip(ce,-1,1))
    return np.mean((np.pi-e)/np.sin(np.maximum(e,1e-3)),axis=1)
def offdiag(M): i,j=np.triu_indices(M.shape[0],1); return M[i,j]
def corr(a,b): a=a-a.mean(); b=b-b.mean(); return float(a@b/np.sqrt((a@a)*(b@b)))
def frac(P,H):                                                       # fraction of HD's pair-to-pair structure a pattern P captures (best scaling + offset)
    X=np.c_[np.ones(len(P)),P]; c=np.linalg.lstsq(X,H,rcond=None)[0]; r=H-X@c; h=H-H.mean(); return 1-float(r@r)/float(h@h)
N=67; res={"solar_wind":[],"ecliptic_quadrupole":[],"best_separable":[],"sw_negative_pairs":[]}
for it in range(300):
    v=iso(N); C=v@v.T; H=offdiag(hd(C)); b=ecl_lat(v)
    g=g_sw(b); Psw=offdiag(np.outer(g,g)); res["solar_wind"].append(frac(Psw,H)); res["sw_negative_pairs"].append(float(np.mean(Psw<0)))
    q=3*np.sin(b)**2-1; res["ecliptic_quadrupole"].append(frac(offdiag(np.outer(q,q)),H))
    # best separable (rank-1) pattern: leading eigenvector of the HD matrix with its diagonal set to the mean pair value
    Hm=hd(C); np.fill_diagonal(Hm,0); w,U=np.linalg.eigh(Hm); u=U[:,-1]; res["best_separable"].append(frac(offdiag(np.outer(u,u)),H))
    u2=U[:,0]; res["best_separable"][-1]=max(res["best_separable"][-1],frac(offdiag(np.outer(u2,u2)),H))
out["pattern"]={k:dict(median=float(np.median(v)),p95=float(np.percentile(v,95)),max=float(np.max(v))) for k,v in res.items()}
print(f"\nFraction of the Hellings-Downs pair pattern a separable plasma pattern can reproduce ({N} pulsars, 300 random arrays):")
for k,lab in (("solar_wind","solar wind (annual mean geometry)"),("ecliptic_quadrupole","best ecliptic pattern (3 sin^2 b - 1)"),("best_separable","best separable pattern of any shape")):
    s=out["pattern"][k]; print(f"   {lab:42s} median {s['median']:.3f}; 95% {s['p95']:.3f}; max {s['max']:.3f}")
print(f"   solar-wind pairs with negative correlation: {np.median(res['sw_negative_pairs']):.3f} (HD is negative for separations of 50 to 121 deg)")
# HD Legendre content: share of the quadrupole and of one quadrupole mode
ls=np.arange(2,60); cl=(2*ls+1)/((ls+2)*(ls+1)*ls*(ls-1)); lam=cl/(2*ls+1)     # per-mode eigenvalue, up to a constant
tot=np.sum((2*ls+1)*lam**2); out["hd_one_mode_share"]=float(lam[0]**2/tot); out["hd_quadrupole_share"]=float(5*lam[0]**2/tot)
print(f"   HD power in the quadrupole: {5*lam[0]**2/tot:.3f}; in any single mode: {lam[0]**2/tot:.3f} (full-sky limit of a separable pattern)")
# angular profile of the solar-wind pattern for the figure (binned over random arrays)
bins=np.linspace(0,180,16); prof=[[] for _ in range(15)]; hdp=[[] for _ in range(15)]
for it in range(60):
    v=iso(N); C=v@v.T; i,j=np.triu_indices(N,1); ang=np.degrees(np.arccos(np.clip(C[i,j],-1,1))); g=g_sw(ecl_lat(v)); P=(np.outer(g,g))[i,j]
    # scale the solar-wind pattern to the same mean-square as HD, after removing nothing: a common plasma signal with HD's total power
    Hh=hd(C)[i,j]; P=P*np.sqrt(np.mean(Hh**2)/np.mean(P**2))
    k=np.digitize(ang,bins)-1
    for q in range(15): prof[q]+=list(P[k==q]); hdp[q]+=list(Hh[k==q])
mid=0.5*(bins[1:]+bins[:-1]); out["profile"]=dict(angle=mid.tolist(),sw=[float(np.mean(p)) for p in prof],sw_sd=[float(np.std(p)) for p in prof],hd=[float(np.mean(p)) for p in hdp])
za=np.linspace(0.5,180,200); out["hd_curve"]=dict(angle=za.tolist(),hd=hd(np.cos(np.radians(za))).tolist())
json.dump(out,open(os.path.join(HERE,"pta_plasma.json"),"w"),indent=1)
