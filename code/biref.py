"""Workbench Entry 116: is the rotation of the CMB's polarization Faraday rotation by magnetized plasma?

The CMB's E and B polarization patterns correlate as if every polarization direction on the sky were rotated by the same small angle,
beta. Planck gave 0.30-0.35 deg; ACT DR6 0.21 deg; a joint 2026 analysis 0.277 +- 0.057 deg (4.8 sigma), 0.236 +- 0.067 with dust
contamination limited. Magnetized plasma rotates polarization by Faraday rotation, angle = RM lambda^2. Two tests:
 1. Frequency: Faraday rotation scales as nu^-2. Planck measured beta band by band from 30 to 353 GHz (Eskilt 2022, Table 3).
 2. Symmetry: a rotation by the same angle in every direction needs the line-of-sight field to have the same sign in every direction.
Data: Eskilt 2022 (A&A 662, A10), Table 3, f_sky = 0.92, foreground EB ignored; errors treated as independent."""
import os, json
import numpy as np
from scipy import optimize, integrate
HERE=os.path.dirname(os.path.abspath(__file__))
nu=np.array([30,44,70,100,143,217,353.]); b=np.array([0.75,0.27,0.82,0.42,0.20,0.28,-0.19]); e=np.array([0.68,0.68,0.45,0.16,0.12,0.12,0.38])
w=1/e**2; out={"bands":dict(nu=nu.tolist(),beta=b.tolist(),err=e.tolist())}
# constant
c=np.sum(w*b)/w.sum(); ec=1/np.sqrt(w.sum()); chic=float(np.sum(w*(b-c)**2)); out["const"]=dict(beta=c,err=ec,chi=chic)
print(f"Constant angle: beta = {c:.3f} ± {ec:.3f} deg, chi2 = {chic:.1f} for 6")
# Faraday: beta = beta150 (nu/150)^-2
def chi_n(n):
    t=(nu/150)**n; a=np.sum(w*b*t)/np.sum(w*t*t); return float(np.sum(w*(b-a*t)**2)),a,1/np.sqrt(np.sum(w*t*t))
cf,af,ef=chi_n(-2); out["faraday"]=dict(beta150=af,err=ef,chi=cf,pred30=af*(30/150)**-2,pred353=af*(353/150)**-2)
print(f"Faraday (nu^-2): beta(150 GHz) = {af:.3f} ± {ef:.3f} deg, chi2 = {cf:.1f} for 6; it predicts {af*25:.2f} deg at 30 GHz (measured 0.75 ± 0.68)")
# power law, profile in n
ns=np.linspace(-3,2,2001); ch=np.array([chi_n(n)[0] for n in ns]); i=int(np.argmin(ch)); nb=ns[i]; ok=ns[ch<ch[i]+1]
out["power"]=dict(n=float(nb),lo=float(ok.min()),hi=float(ok.max()),chi=float(ch[i]),dchi_faraday=float(cf-ch[i]))
print(f"Power law nu^n: n = {nb:.2f} (+{ok.max()-nb:.2f} −{nb-ok.min():.2f}); Faraday n = −2 worse by Δχ² = {cf-ch[i]:.1f} ({np.sqrt(max(cf-ch[i],0)):.1f}σ)")
print("   Eskilt's own fits: n = −0.35 (+0.48 −0.47) with 30 and 44 GHz; −0.37 (+0.49 −0.64) with the dust model")
# HFI only (100-353 GHz), where the signal sits
m=nu>=100; ns2=np.linspace(-4,3,1401)
def chi_n2(n):
    t=(nu[m]/150)**n; a=np.sum(w[m]*b[m]*t)/np.sum(w[m]*t*t); return float(np.sum(w[m]*(b[m]-a*t)**2))
ch2=np.array([chi_n2(n) for n in ns2]); j=int(np.argmin(ch2)); ok2=ns2[ch2<ch2[j]+1]
out["power_hfi"]=dict(n=float(ns2[j]),lo=float(ok2.min()),hi=float(ok2.max()),dchi_faraday=float(chi_n2(-2)-ch2[j]))
print(f"HFI bands only: n = {ns2[j]:.2f} (+{ok2.max()-ns2[j]:.2f} −{ns2[j]-ok2.min():.2f}); Faraday worse by Δχ² = {chi_n2(-2)-ch2[j]:.1f}")
# what Faraday rotation would need: RM for beta = 0.277 deg at 150 GHz
lam2=(2.99792458e8/150e9)**2; beta=0.277; RM=np.radians(beta)/lam2; out["RM_needed"]=float(RM)
print(f"\nFor 0.277 deg at 150 GHz, RM = {RM:.0f} rad/m2, the same sign in every direction")
# a field threading the last-scattering plasma: RM = 0.812 n_e B L (cm^-3, microgauss, pc) with n_e L = tau / sigma_T and tau ~ 1 behind
# the visibility peak; B redshifts with the plasma so the comoving field enters directly.
sigT=6.6524e-25; PC=3.0857e18; neL=1/sigT/PC; RM_per_nG=0.812*neL*1e-3
out["RM_per_nG_comoving"]=float(RM_per_nG); out["beta_per_nG_150"]=float(np.degrees(RM_per_nG*lam2)); out["B_needed_nG"]=float(RM/RM_per_nG)
print(f"A comoving field B0 along the line of sight at last scattering: RM = {RM_per_nG:.0f} rad/m2 per nG, {np.degrees(RM_per_nG*lam2):.3f} deg at 150 GHz; "
      f"0.277 deg needs {RM/RM_per_nG:.1f} nG pointing toward us (or away) in every direction")
# the symmetry argument, made numerical: a field with no monopole has zero all-sky mean of B . n on any sphere around us.
# Check with random divergence-free fields: the sky-mean RM over the sky-rms RM, for many realisations
rng=np.random.default_rng(116)
def sky_mean_ratio(nmodes=300, kmax=6.0, R=1.0, npix=20000):
    k=rng.normal(size=(nmodes,3))*kmax; A=rng.normal(size=(nmodes,3)); A-= (np.sum(A*k,1)/np.sum(k*k,1))[:,None]*k   # transverse amplitudes: div B = 0
    ph=rng.uniform(0,2*np.pi,nmodes); u=rng.normal(size=(npix,3)); u/=np.linalg.norm(u,axis=1)[:,None]
    rs=np.linspace(0.02,R,60); rm=np.zeros(npix)
    for r in rs:
        x=u*r; Bv=np.cos(x@k.T+ph)@A; rm+=np.sum(Bv*u,1)*(rs[1]-rs[0])
    return abs(rm.mean())/rm.std()
rat=[sky_mean_ratio() for _ in range(40)]; out["sky_mean_over_rms"]=dict(median=float(np.median(rat)),p95=float(np.percentile(rat,95)))
print(f"Random divergence-free fields, line-of-sight integral to the shell: |sky mean| / sky rms = {np.median(rat):.3f} (95%: {np.percentile(rat,95):.3f});"
      f" pixel noise alone gives about {1/np.sqrt(20000):.3f}")
# Galactic Faraday screen: an upper bound on its monopole at CMB frequencies
for RMg in (10,100):
    out[f"galactic_{RMg}"]={str(f):float(np.degrees(RMg*(2.99792458e8/(f*1e9))**2)) for f in (30,70,150,353)}
    print(f"   a Galactic RM of {RMg} rad/m2 rotates by "+", ".join(f"{f} GHz: {np.degrees(RMg*(2.99792458e8/(f*1e9))**2):.4f} deg" for f in (30,70,150,353)))
fg=np.geomspace(25,400,60); out["curve"]=dict(nu=fg.tolist(),faraday=(af*(fg/150)**-2).tolist(),const=c)
json.dump(out,open(os.path.join(HERE,"biref.json"),"w"),indent=1)
