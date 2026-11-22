"""Workbench Entry 29, second part: variants for a pure-hydrogen target, other spectra, and what the yield means for the model's budget.
Run after deuterium_yield.py (reads deuterium_yield.json)."""
import numpy as np, json, sys, os
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import deuterium_yield as dy
J=json.load(open("deuterium_yield.json")); E=np.array(J["E"]); out={}
Tf=np.geomspace(10,1e5,400)
def avg(curve,spec):
    w=spec(Tf)*np.gradient(Tf); w/=w.sum(); c=np.array(curve); return [float((w*np.interp(np.log(Tf),np.log(E),c[:,j],left=0)).sum()) for j in range(4)], float((w*Tf).sum())
base=lambda T:(T+980.0)**-2.7
print("pure hydrogen target, variants (deuterons per cosmic-ray proton = slow neutrons + direct deuterons):")
out["H_variants"]={}
for name,kw in (("ionised gas, stops twice as fast",dict(plasma=2.0)),("inelastic p-p x0.8",dict(f_inel=0.8)),("inelastic p-p x1.2",dict(f_inel=1.2)),("neutron in 70% of p-p pion events",dict(pn_frac=0.70)),("leading nucleon keeps more",dict(lead=(0.7,0.95)))):
    M=dy.Medium(Y=0.0,**kw); yc=dy.yield_curve(M,"p",E,1200); a,_=avg(yc,base); out["H_variants"][name]=a[0]+a[1]; print(f"  {name:36s} {a[0]+a[1]:.2f}")
print("\nother spectra (mean energy; pure hydrogen; hydrogen and helium; per alpha):")
out["spectra"]={}
for name,sp in (("(E+980)^-2.7, the model's mean",base),("softer, (E+300)^-2.7",lambda T:(T+300.0)**-2.7),("harder, nothing below 300 MeV",lambda T:np.where(T>=300,(T+980.0)**-2.7,0.0)),("Lerner's mean of 1.5 GeV: (E+1050)^-2.7",lambda T:(T+1050.0)**-2.7)):
    h,m=avg(J["cases"]["hydrogen only"]["curve"],sp); he,_=avg(J["cases"]["hydrogen and helium (Y = 0.25)"]["curve"],sp); al,_=avg(J["alpha"]["curve"],sp)
    out["spectra"][name]=dict(mean=m,H=h[0]+h[1],HHe=he[0]+he[1],alpha=al[0]+al[1]); print(f"  {name:40s} mean {m/1000:.2f} GeV   {h[0]+h[1]:.2f}   {he[0]+he[1]:.2f}   {al[0]+al[1]:.2f}")
c=J["cases"]; H=c["hydrogen only"]["avg"]; HHe=c["hydrogen and helium (Y = 0.25)"]["avg"]; A=J["alpha"]["avg"]; apr=0.085
tot=[HHe[j]+apr*A[j] for j in range(4)]
allv=[v["avg"][0]+v["avg"][1] for k,v in c.items() if k!="hydrogen only"]
print(f"\nsummary, deuterons per cosmic-ray proton if every slow neutron is captured by hydrogen:")
print(f"  hydrogen target, protons only:            {H[0]+H[1]:.2f}   (variants {min(out['H_variants'].values()):.2f} to {max(out['H_variants'].values()):.2f})")
print(f"  hydrogen and helium, protons only:        {HHe[0]+HHe[1]:.2f}   (variants {min(allv):.2f} to {max(allv):.2f})")
print(f"  hydrogen and helium, with {apr} alphas:   {tot[0]+tot[1]:.2f}")
print(f"  the model's figure:                       0.9")
out["summary"]=dict(H=H[0]+H[1],HHe=HHe[0]+HHe[1],HHe_range=[min(allv),max(allv)],H_range=[min(out['H_variants'].values()),max(out['H_variants'].values())],with_alpha=tot[0]+tot[1],alpha=A[0]+A[1])
print("\nneutrons lost to 3He (capture cross section 5330 b against 0.332 b for hydrogen):")
for r in (0,1e-5,3e-5):
    f=1/(1+5330/0.332*r); print(f"  3He/H = {r:.0e}: {100*f:.0f}% captured by hydrogen -> {tot[0]*f+tot[1]:.2f} deuterons per proton (with helium and alphas), {H[0]*f+H[1]:.2f} (hydrogen only)"); out.setdefault("he3",[]).append([r,f,tot[0]*f+tot[1],H[0]*f+H[1]])
he3=(tot[2]+tot[3]); print(f"\nby-product: 3He + 3H (which decays to 3He) per deuteron: {he3/(tot[0]+tot[1]):.2f}; observed 3He/D: about {1.1e-5/2.5e-5:.2f} (Galactic 3He/H 1.1e-5, D/H 2.5e-5)"); out["he3_per_D"]=he3/(tot[0]+tot[1])
# lithium of Entry 25 with these yields
li=json.load(open("lithium_thick_target.json" if os.path.exists("lithium_thick_target.json") else "../gole/lithium_thick_target.json")); print("\nlithium per deuteron of Entry 25, recomputed with these deuterium yields (it used 0.9):")
lo,hi=li["LiD_range"]; mid=li["cases"]["(E+980)^-2.7, mean 1.4 GeV/n | cold gas"]["LiD"]; y=tot[0]+tot[1]
print(f"  model spectrum, cold gas: {mid:.2e} -> {mid*0.9/y:.2e};  model-spectrum cases 1.0-2.0e-5 -> {1.0e-5*0.9/y:.1e} to {2.0e-5*0.9/y:.1e};  observed 0.64-0.90e-5"); out["LiD"]=dict(old=mid,new=mid*0.9/y,lo=1.0e-5*0.9/y,hi=2.0e-5*0.9/y)
# energy
Emean=1400.0; Etot=Emean*(1+4*apr); perD=Etot/(tot[0]+tot[1]); perD_H=Emean/(H[0]+H[1]); print(f"\ncosmic-ray energy spent per deuteron: {perD/1000:.2f} GeV (helium present, alphas included), {perD_H/1000:.1f} GeV (hydrogen only); the model: 1.5 GeV / 0.9 = {1.5/0.9:.2f} GeV")
DH=1.4e-5; print(f"cosmic-ray energy needed per hydrogen atom for the model's pregalactic D/H of {DH}: {1000*DH*perD:.0f} keV (model: 20 keV = 0.5% of 4 MeV); in pure hydrogen {1000*DH*perD_H:.0f} keV")
fus=0.25*26.73/4*1000/0.75    # keV released per hydrogen atom of the final mix when a quarter of the mass has become helium (26.73 MeV per helium nucleus)
print(f"energy released in making the helium: {fus/1000:.2f} MeV per remaining hydrogen atom (Y = 0.25); cosmic-ray share needed: {100*1000*DH*perD/fus:.2f}% (helium present) to {100*1000*DH*perD_H/fus:.2f}% (hydrogen only); with the model's 4 MeV: {100*1000*DH*perD/4000:.2f}% to {100*1000*DH*perD_H/4000:.2f}%")
mw=(0.7e41/2.3e44,1.0e41/1.5e44)    # cosmic-ray power 0.7-1.0e41 erg/s (Strong et al. 2010), starlight 1.5-2.3e44 erg/s
sb=(3e40/2.2e43,6e40/2.2e43)          # one supernova per 50-100 solar masses formed, 1e51 erg each, a tenth into cosmic rays; starlight 2.2e43 erg/s per solar mass a year
print(f"cosmic-ray share of stellar power: Milky Way today {100*mw[0]:.3f}-{100*mw[1]:.3f}%; a young starburst (10% of supernova energy) {100*sb[0]:.2f}-{100*sb[1]:.2f}%")
need=(1000*DH*perD/4000,1000*DH*perD_H/fus); print(f"needed / Milky Way: {need[0]/mw[1]:.0f} to {need[1]/mw[0]:.0f} times; needed / starburst: {need[0]/sb[1]:.1f} to {need[1]/sb[0]:.1f} times")
out["energy"]=dict(perD=perD,perD_H=perD_H,keV_per_H=1000*DH*perD,keV_per_H_pureH=1000*DH*perD_H,fusion_keV_per_H=fus,share=[1000*DH*perD/4000,1000*DH*perD/fus,1000*DH*perD_H/4000,1000*DH*perD_H/fus],mw=list(mw),starburst=list(sb),need_over_mw=[need[0]/mw[1],need[1]/mw[0]],need_over_sb=[need[0]/sb[1],need[1]/sb[0]])
# capture gamma rays
nH=0.75*4.2e-31/1.6726e-24; u=nH*2.5e-5*2.223e6; mev=4*np.pi/2.998e10*3e-3*1e6*np.log(3/0.5)
print(f"\ncapture gamma rays: one 2.223 MeV photon per deuteron made by capture; if all escaped, {u:.1e} eV/cm3; the observed background from 0.5 to 3 MeV holds about {mev:.1e} eV/cm3: at most {100*mev/u:.0f}% may escape")
out["capture_gamma"]=dict(u=u,mev=mev,max_escape=mev/u)
# does the neutron live long enough to be captured?  capture rate = n_H x (0.332 b x 2200 m/s); free-neutron lifetime 880 s
print("\nneutron decay against capture, at the depth where the cascade develops (grey atmosphere, electron-scattering opacity 0.35 cm2/g, P = g x column):")
out["decay"]=[]
for name,Teff,logg in (("hot main-sequence star, 30,000 K, log g 4.0",3.0e4,4.0),("B star, 20,000 K, log g 4.0",2.0e4,4.0),("hot supergiant, 25,000 K, log g 3.0",2.5e4,3.0)):
    for m in (60.0,200.0):
        Tg=Teff*(0.75*0.35*m+0.5)**0.25; ntot=10**logg*m/(1.380649e-16*Tg); nHy=0.75*0.62*ntot; tc=1/(nHy*0.332e-24*2.2e5); f=880/(880+tc)
        print(f"  {name:44s} column {m:4.0f} g/cm2: T = {Tg:6.0f} K, hydrogen {nHy:.1e} per cm3, capture time {tc:5.0f} s, captured before decay {100*f:.0f}%"); out["decay"].append([name,m,Tg,nHy,tc,f])
Tg=1.2e4; nHy=0.9*10**4.44*60/(1.380649e-16*Tg)*0.75; tc=1/(nHy*0.332e-24*2.2e5); print(f"  the Sun, 60 g/cm2 down (about 12,000 K): capture time {tc:.0f} s, captured {100*880/(880+tc):.0f}%")
fs=[d[5] for d in out["decay"][:4]]; print(f"hot main-sequence stars: {100*min(fs):.0f}-{100*max(fs):.0f}% captured -> yields {tot[0]*min(fs)+tot[1]:.2f}-{tot[0]*max(fs)+tot[1]:.2f} (helium and alphas), {HHe[0]*min(fs)+HHe[1]:.2f}-{HHe[0]*max(fs)+HHe[1]:.2f} (helium, protons only), {H[0]*min(fs)+H[1]:.2f}-{H[0]*max(fs)+H[1]:.2f} (hydrogen only)")
out["decay_yields"]=dict(f=[min(fs),max(fs)],all=[tot[0]*min(fs)+tot[1],tot[0]*max(fs)+tot[1]],he=[HHe[0]*min(fs)+HHe[1],HHe[0]*max(fs)+HHe[1]],h=[H[0]*min(fs)+H[1],H[0]*max(fs)+H[1]])
dy_hot=out["decay_yields"]; eD=[Etot/dy_hot["all"][1],Etot/dy_hot["all"][0]]; eDH=[Emean/dy_hot["h"][1],Emean/dy_hot["h"][0]]
print(f"\non hot main-sequence stars: cosmic-ray energy per deuteron {eD[0]/1000:.2f}-{eD[1]/1000:.2f} GeV (helium present, alphas), {eDH[0]/1000:.1f}-{eDH[1]/1000:.1f} GeV (hydrogen only); per hydrogen atom {1000*DH*eD[0]:.1f}-{1000*DH*eD[1]:.1f} keV and {1000*DH*eDH[0]:.0f}-{1000*DH*eDH[1]:.0f} keV")
sh=[1000*DH*eD[0]/4000,1000*DH*eD[1]/fus]; shH=[1000*DH*eDH[0]/4000,1000*DH*eDH[1]/fus]
print(f"share of fusion energy: {100*sh[0]:.2f}-{100*sh[1]:.2f}% (helium present), {100*shH[0]:.2f}-{100*shH[1]:.2f}% (hydrogen only); over a starburst's share: {sh[0]/sb[1]:.1f}-{sh[1]/sb[0]:.1f} times; over the Milky Way's: {sh[0]/mw[1]:.0f}-{sh[1]/mw[0]:.0f} times")
out["hot"]=dict(eD=eD,eDH=eDH,keV=[1000*DH*eD[0],1000*DH*eD[1]],keV_H=[1000*DH*eDH[0],1000*DH*eDH[1]],share=sh,share_H=shH,over_sb=[sh[0]/sb[1],sh[1]/sb[0]],over_mw=[sh[0]/mw[1],sh[1]/mw[0]])
print(f"lithium per deuteron on hot stars: range of Entry 25's cases {1.0e-5*0.9/dy_hot['all'][1]:.1e}-{2.0e-5*0.9/dy_hot['all'][0]:.1e}; model spectrum, cold gas {mid*0.9/dy_hot['all'][1]:.2e}-{mid*0.9/dy_hot['all'][0]:.2e}")
out["LiD_hot"]=[1.0e-5*0.9/dy_hot["all"][1],2.0e-5*0.9/dy_hot["all"][0]]; out["he3_per_D_hot"]=[he3/dy_hot["all"][1],he3/dy_hot["all"][0]]; print("3He per deuteron on hot stars:",[round(x,2) for x in out["he3_per_D_hot"]])
json.dump(out,open("deuterium_budget.json","w"),indent=1)
