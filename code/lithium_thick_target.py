"""Workbench Entry 25: lithium and deuterium from cosmic rays, thin target against thick target.
The classic objection to making deuterium with cosmic rays: the same cosmic rays make far too much lithium (alpha + alpha -> 6Li, 7Li, 7Be).
That objection is worked out for cosmic rays in thin interstellar gas. In Lerner's GOLE model the cosmic rays are trapped near stars
and stop in dense stellar material. This script works out the lithium yield of alpha particles that slow down and stop in a thick
target of hydrogen and helium, and the resulting Li/D ratio when deuterium is made at the model's own stated rate.

Inputs
  alpha + alpha cross sections: fits of Mercer et al. (2001, nucl-ex/0104005), E = alpha lab energy in MeV:
      mass 7 (7Li + 7Be): 510 exp(-0.0354 E) mb above the 34.7 MeV threshold; mass 6 (6Li + 6He): 0.014 + 75 exp(-0.0159 E) mb above 44.7 MeV.
      The fits rest on data from 60 MeV up; below that they are extrapolations.
  stopping: Bethe formula for a cold gas of hydrogen and helium (X = 0.75, Y = 0.25); a fully ionised plasma stops roughly twice as fast (variant).
  nuclear losses on the way down: sigma(alpha + p) = 100 mb, sigma(alpha + alpha) = 270 mb (assumed; varied by +-50%).
  cosmic-ray spectrum: dN/dE ~ (E + 980 MeV/n)^-2.7, which has the mean energy the model quotes (1.4 GeV per nucleon); two other shapes as variants.
  deuterium: 0.9 deuterons per cosmic-ray proton (the model's own figure, not checked here); alphas per proton 0.085 (helium/hydrogen by number for Y = 0.25).
Observed: D/H = 2.5e-5; Li/H = 1.6e-10 (Spite plateau, A(Li) = 2.2) to 2.3e-10 (A(Li) = 2.36, the value the model quotes)."""
import numpy as np, json
me=0.510999; u=931.494; K=0.307075
X,Y=0.75,0.25
comp=[(X,1/1.008,19.2e-6),(Y,2/4.0026,41.8e-6)]              # mass fraction, Z/A, mean excitation energy (MeV)
nH=X/(1.008*1.66054e-24); nHe=Y/(4.0026*1.66054e-24)          # atoms per gram of the mixture
def S(En,plasma=1.0):                                         # energy loss of an alpha, MeV per g/cm^2, at En MeV per nucleon
    g=1+En/u; b2=1-1/g**2
    return plasma*sum(w*K*4*za/b2*(np.log(2*me*b2*g*g/I)-b2) for w,za,I in comp)
sig7=lambda E: np.where(E>=34.7,510*np.exp(-0.0354*E),0.0)    # mb, E = total alpha energy (MeV)
sig6=lambda E: np.where(E>=44.7,0.014+75*np.exp(-0.0159*E),0.0)
def yields(E0n,plasma=1.0,fnuc=1.0):
    """6Li and 7Li nuclei made per alpha entering a thick target with E0n MeV per nucleon"""
    En=np.linspace(8.0,E0n,4000) if E0n>8.5 else np.array([8.0,8.5]); E=4*En
    s=S(En,plasma); dx=np.gradient(E)/s                         # g/cm^2 per step
    path=np.cumsum(dx[::-1])[::-1]-dx                           # grammage already crossed when the alpha has slowed to E
    surv=np.exp(-fnuc*(nH*100e-27+nHe*270e-27)*path)
    return float(np.sum(nHe*sig6(E)*1e-27*surv*dx)), float(np.sum(nHe*sig7(E)*1e-27*surv*dx)), float(path[0])
print("one alpha stopping in a thick target:")
print("  start (MeV/n)   range (g/cm^2)   7Li+7Be per alpha   6Li+6He per alpha   6/7")
tab=[]
for e0 in (20,40,70,100,200,300,500,1000,3000):
    y6,y7,rng=yields(e0); tab.append([e0,rng,y7,y6]); print(f"  {e0:6d}        {rng:9.2f}        {y7:.2e}           {y6:.2e}        {y6/y7:.2f}")
def spectrum_average(shape,plasma=1.0,fnuc=1.0):
    E0=np.logspace(np.log10(9),np.log10(2e4),260); w=shape(E0)*np.gradient(E0); w/=w.sum()
    ys=np.array([yields(e,plasma,fnuc)[:2] for e in E0]); return float((w*ys[:,0]).sum()), float((w*ys[:,1]).sum()), float((w*E0).sum())
shapes={"(E+980)^-2.7, mean 1.4 GeV/n":lambda E:(E+980.0)**-2.7,"harder: nothing below 300 MeV/n":lambda E:np.where(E>=300,(E+980.0)**-2.7,0.0),"softer: (E+300)^-2.7":lambda E:(E+300.0)**-2.7}
D_per_p=0.9; a_per_p=0.085; LiD_obs=(1.6e-10/2.5e-5,2.3e-10/2.55e-5)
print(f"\nobserved Li/D: {LiD_obs[0]:.1e} to {LiD_obs[1]:.1e};  the model's own prediction 1.13e-5;  thin interstellar gas (Wright): 1/80 = {1/80:.1e}")
out={"per_alpha":tab,"cases":{}}
for name,sh in shapes.items():
    for lab,pl,fn in (("cold gas",1.0,1.0),("ionised plasma, stopping x2",2.0,1.0),("nuclear losses x0.5",1.0,0.5),("nuclear losses x1.5",1.0,1.5)):
        y6,y7,mean=spectrum_average(sh,pl,fn); r=a_per_p*y7/D_per_p
        if lab=="cold gas": print(f"\nspectrum {name} (mean {mean:.0f} MeV/n):")
        print(f"   {lab:28s} 7Li per alpha {y7:.2e}   Li/D = {r:.1e}   ({r/LiD_obs[0]:.1f} to {r/LiD_obs[1]:.1f} x observed)   6Li/7Li made = {y6/y7:.2f}")
        out["cases"][name+" | "+lab]={"y7":y7,"y6":y6,"LiD":r,"ratio67":y6/y7}
vals=[c["LiD"] for c in out["cases"].values()]; r67=[c["ratio67"] for c in out["cases"].values()]
print(f"\nrange over all cases: Li/D = {min(vals):.1e} to {max(vals):.1e};  6Li/7Li at production {min(r67):.2f} to {max(r67):.2f}")
for lim,star in ((0.007,"HD 84937"),(0.006,"HD 140283"),(0.017,"LP 815-43")):
    print(f"  {star}: 6Li/7Li < {lim} (2 sigma) -> at least {100*(1-lim/min(r67)):.1f}% to {100*(1-lim/max(r67)):.1f}% of the 6Li destroyed, a factor {min(r67)/lim:.0f} to {max(r67)/lim:.0f}")
out["LiD_range"]=[min(vals),max(vals)]; out["r67_range"]=[min(r67),max(r67)]; out["LiD_obs"]=list(LiD_obs)
json.dump(out,open("lithium_thick_target.json","w"),indent=1)
