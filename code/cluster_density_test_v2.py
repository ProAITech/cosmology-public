"""Workbench Entry 24: were galaxy clusters denser in the past?  A test that uses no distances.
Two things can be measured for a cluster without knowing how far away it is:
  T  : the temperature of its gas, from the shape of its X-ray spectrum;
  y0 : the depth of the shadow it casts on the CMB (the central Compton parameter), which is the pressure of the gas summed along the line of sight.
At a fixed temperature a cluster has a fixed depth of gravitational well. If the universe was denser when the light left, clusters of that
temperature were denser too, and y0 is larger: for self-similar clusters y0 ~ E(z)^0.92 T^1.62 with ACT's mass scaling (y0 ~ E^2 M^1.08, M ~ T^1.5/E).
In a universe that does not evolve, y0 at fixed T is the same at every redshift: exponent m = 0.
Data: clusters_y_T.txt = 125 ACT DR5 clusters with an X-ray temperature (columns: z, kT in keV, log10 y0 in arbitrary units, ACT signal-to-noise,
source of the temperature: A = ACCEPT/Chandra, XC = XMM Cluster Survey, XS = 2XMMi/SDSS, B = BAX compilation, M = Maughan et al. 2008, Chandra),
matched in the browser (sz_shadows.js). maughan2008_matched.txt holds all 38 ACT clusters with a Maughan et al. temperature, 27 of them also in the first list.
y0 is rebuilt from the catalogue as E(z)^2 x (uncorrected mass)^1.08, which undoes the catalogue's own conversion and returns the measured signal."""
import numpy as np, json
rows=[r.split() for r in open("clusters_y_T_v2.txt").read().strip().split(";")]
z=np.array([float(r[0]) for r in rows]); T=np.array([float(r[1]) for r in rows]); ly=np.array([float(r[2]) for r in rows]); snr=np.array([float(r[3]) for r in rows]); src=np.array([r[4] for r in rows])
E=np.sqrt(0.3*(1+z)**3+0.7); lE=np.log10(E); lT=np.log10(T/5)
print(len(z),"clusters; z",z.min(),"-",z.max(),"; kT median",np.median(T),"keV; sources",{s:int((src==s).sum()) for s in sorted(set(src))})
def ols(X,y):
    X=np.asarray(X); b,res,rk,sv=np.linalg.lstsq(X,y,rcond=None); r=y-X@b; s2=r@r/(len(y)-X.shape[1]); cov=s2*np.linalg.inv(X.T@X); return b,np.sqrt(np.diag(cov)),np.sqrt(s2)
def fit(mask,label):
    b,e,rms=ols(np.c_[np.ones(mask.sum()),lT[mask],lE[mask]],ly[mask]); b2,e2,_=ols(np.c_[np.ones(mask.sum()),lE[mask]],ly[mask]-1.62*lT[mask])
    print(f"  {label:34s} n={mask.sum():3d}  temperature exponent {b[1]:.2f}+-{e[1]:.2f}  m = {b[2]:+.2f}+-{e[2]:.2f}   (exponent fixed at 1.62: m = {b2[1]:+.2f}+-{e2[1]:.2f})")
    return {"n":int(mask.sum()),"alpha":b[1],"alpha_se":e[1],"m":b[2],"m_se":e[2],"m_fixed":b2[1],"m_fixed_se":e2[1],"rms":rms}
out={}
print("fit: log y0 = a + alpha log T + m log E(z)   [expanding, self-similar: m = 0.92;  no evolution: m = 0]")
all_=np.ones(len(z),bool); out["all"]=fit(all_,"all clusters")
out["snr6"]=fit(snr>=6,"ACT signal-to-noise >= 6"); out["T4"]=fit(T>=4,"kT >= 4 keV"); out["chandra"]=fit((src=="A")|(src=="M"),"Chandra only (ACCEPT, Maughan)"); out["notbax"]=fit(src!="B","without the BAX compilation"); out["zlt06"]=fit(z<0.6,"z < 0.6"); out["zlt045"]=fit(z<0.45,"z < 0.45")
m,se=out["all"]["m"],out["all"]["m_se"]; print(f"\nall clusters: m = {m:.2f} +- {se:.2f}: {m/se:.1f} sigma from no evolution, {abs(m-0.92)/se:.1f} sigma from the self-similar value")
rng=np.random.default_rng(3); bs=[]
for _ in range(2000):
    i=rng.integers(0,len(z),len(z)); b,_,_=ols(np.c_[np.ones(len(i)),lT[i],lE[i]],ly[i]); bs.append(b[2])
print(f"bootstrap: m = {np.mean(bs):.2f} +- {np.std(bs):.2f}; share of resamples with m <= 0: {np.mean(np.array(bs)<=0):.4f}")
out["boot"]=[float(np.mean(bs)),float(np.std(bs)),float(np.mean(np.array(bs)<=0))]
# the 38 clusters with temperatures from one uniform Chandra analysis (Maughan et al. 2008), as an independent set of temperatures
mr=[r.split() for r in open("maughan2008_matched.txt").read().strip().split(";")]
mz=np.array([float(r[0]) for r in mr]); mT=np.array([float(r[1]) for r in mr]); mly=np.array([float(r[2]) for r in mr]); mlE=np.log10(np.sqrt(0.3*(1+mz)**3+0.7))
b,e,_=ols(np.c_[np.ones(len(mz)),np.log10(mT/5),mlE],mly); b2,e2,_=ols(np.c_[np.ones(len(mz)),mlE],mly-1.62*np.log10(mT/5))
print(f"\nMaughan et al. temperatures only: n={len(mz)}, z {mz.min()}-{mz.max()}: temperature exponent {b[1]:.2f}+-{e[1]:.2f}, m = {b[2]:+.2f}+-{e[2]:.2f} (fixed 1.62: {b2[1]:+.2f}+-{e2[1]:.2f})")
out["maughan"]={"n":len(mz),"alpha":b[1],"alpha_se":e[1],"m":b[2],"m_se":e[2],"m_fixed":b2[1],"m_fixed_se":e2[1]}
# attempted scale-up with eROSITA survey temperatures (eRASS1, Bulbul et al. 2024), numbers from the browser run
print("eRASS1 attempt: 3,389 clusters have a temperature; 375 are ACT clusters; median temperature error 64%; fitted temperature exponent 0.42+-0.25 (expected 1.6);")
print("                107 lie beyond z = 0.45, of which 2 have a temperature error below 30%. Not usable for this test.")
out["erass1"]={"with_kT":3389,"matched":375,"median_rel_err":0.637,"alpha":0.422,"alpha_se":0.248,"z_gt_045":107,"z_gt_045_good":2}
# binned residuals at the self-similar temperature exponent, for the figure
b2,_,_=ols(np.c_[np.ones(len(z)),lE],ly-1.62*lT); a0=b2[0]; res=ly-1.62*lT-a0
edges=[0.03,0.1,0.2,0.3,0.45,0.95]; bins=[]
for lo,hi in zip(edges,edges[1:]):
    k=(z>=lo)&(z<hi); bins.append([float(np.median(z[k])),int(k.sum()),float(res[k].mean()),float(res[k].std(ddof=1)/np.sqrt(k.sum())),float(lE[k].mean())])
    print(f"  z {lo}-{hi}: n={k.sum():2d}  mean residual {res[k].mean():+.3f} +- {res[k].std(ddof=1)/np.sqrt(k.sum()):.3f}   (E = {10**lE[k].mean():.2f})")
out["bins"]=bins
# how many clusters a decisive version needs
s=out["all"]["m_se"]*np.sqrt(len(z)); print(f"\nerror scales as {s:.1f}/sqrt(N) for this redshift spread: +-0.10 needs about {int((s/0.10)**2)} clusters, or fewer if they reach higher redshift")
out["n_for_0.1"]=int((s/0.10)**2)
# forecast: clusters added between z = 0.5 and 1.2 with pointed-observation temperatures
rms=out["all"]["rms"]
for n2 in (40,80,160):
    z2=np.linspace(0.5,1.2,n2); l2=np.log10(np.sqrt(0.3*(1+z2)**3+0.7)); la=np.r_[lE,l2]; print(f"adding {n2} clusters at z = 0.5-1.2: error on m about {rms/np.sqrt(((la-la.mean())**2).sum()):.2f}")
    out.setdefault("forecast",{})[str(n2)]=float(rms/np.sqrt(((la-la.mean())**2).sum()))
zz=np.linspace(0,0.95,39); out["curve"]={"z":zz.tolist(),"logE":np.log10(np.sqrt(0.3*(1+zz)**3+0.7)).tolist()}
json.dump(out,open("cluster_density_test_v2.json","w"),indent=1,default=float)
