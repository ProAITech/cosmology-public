"""Workbench Entry 113: the neutron-lifetime puzzle in the early plasma. Does primordial helium choose between the bottle and the beam?

Bottle experiments count neutrons that survive in a trap: UCNtau 2025 gives 877.83 s (+-0.3). Beam experiments count protons from decays in a
cold beam: NIST 2013 gives 887.7 +- 2.2 s. The 9.9 s gap (about 4 sigma) is the neutron-lifetime puzzle. A J-PARC beam that counts electrons
instead of protons gives 877.2 +- 4.2 s (2025), siding with the bottles.
One proposed resolution: about 1% of neutrons decay into dark particles. Then the beam measures the beta-decay lifetime, the bottle the total.
In the early plasma the neutron-proton conversion rates are fixed by the beta-decay matrix element, so they follow the BEAM lifetime, and the
dark channel removes a few neutrons on top.

The network of Entry 53 (bbn_net.py) is used differentially: slopes from the network, anchored to the full-code prediction of Entry 102,
Yp = 0.2467 +- 0.0002 at 878.4 s, eta from the CMB, N_eff = 3.044."""
import os, sys, json, inspect, textwrap
HERE=os.path.dirname(os.path.abspath(__file__))
CODE=HERE if os.path.exists(os.path.join(HERE,"bbn_net.py")) else os.path.abspath(os.path.join(HERE,"..","artifact-files","4379c2ac-1195-4a41-978e-eb44b945a265","code"))
sys.path.insert(0,CODE)
import numpy as np
import bbn_net as bn
# a copy of the network's run() with one extra term: neutrons lost to a dark channel at rate ldark (per second)
src=textwrap.dedent(inspect.getsource(bn.run))
src=src.replace("def run(eta10=6.10,","def run_dark(ldark=0.0, eta10=6.10,",1)
src=src.replace("d[1+N]=-w-r1+r3-r5+r6-r10","d[1+N]=-w-r1+r3-r5+r6-r10-ldark*Y[N]",1)
assert "ldark*Y[N]" in src
exec(src,bn.__dict__); run_dark=bn.run_dark
ETA=6.124                                               # eta10 from omega_b = 0.02237 (273.78 x omega_b, rounded)
def checked(r):
    assert r['t_end']>1e5, "integration stopped early"   # the stiff solver can end early on a glitch; every run is checked
    return r
out={"eta10":ETA}
BOTTLE=(877.83,0.30); BEAM=(887.7,2.2); JPARC=(877.2,4.2); PDG=(878.4,0.5)
print(f"eta10 = {ETA:.3f}.  Network runs (Yp, D/H, Li/H):")
runs={}
for tau in (868.4,878.4,888.4):
    r=checked(bn.run(eta10=ETA,nnu=3.0,tau=tau)); runs[tau]=r; print(f"   tau = {tau:6.1f} s: Yp {r['Yp']:.5f}  D/H {r['D']*1e5:.4f}e-5  Li/H {r['Li7']*1e10:.3f}e-10")
sY=(runs[888.4]['Yp']-runs[868.4]['Yp'])/20; sD=(runs[888.4]['D']/runs[868.4]['D']-1)/20; sL=(runs[888.4]['Li7']/runs[868.4]['Li7']-1)/20
curv=(runs[888.4]['Yp']+runs[868.4]['Yp']-2*runs[878.4]['Yp'])/100
print(f"   dYp/dtau = {sY*1e4:.3f}e-4 per s (curvature {curv:.1e});  dlnD/dtau = {sD*100:.3f}% per s;  dlnLi/dtau = {sL*100:.3f}% per s")
r3=checked(bn.run(eta10=ETA,nnu=3.3,tau=878.4)); sN=(r3['Yp']-runs[878.4]['Yp'])/0.3; sDN=(r3['D']/runs[878.4]['D']-1)/0.3
print(f"   dYp/dN_eff = {sN:.4f};  dlnD/dN_eff = {sDN*100:.2f}%")
out.update(slope_Y=sY,slope_lnD=sD,slope_lnLi=sL,slope_YN=sN,slope_lnDN=sDN)
# dark decay: weak rates follow the beam lifetime; the dark channel is 1/tau_bottle - 1/tau_beam
ld=1/BOTTLE[0]-1/BEAM[0]; rd=checked(run_dark(ldark=ld,eta10=ETA,nnu=3.0,tau=BEAM[0])); rb=checked(bn.run(eta10=ETA,nnu=3.0,tau=BEAM[0]))
dark_shift=rd['Yp']-rb['Yp']
print(f"   dark channel: branching {ld*BOTTLE[0]*100:.2f}%, lifetime {1/ld:,.0f} s; it lowers Yp by {-dark_shift:.5f} relative to beam-only weak rates")
out["dark"]=dict(branch=ld*BOTTLE[0],tau_dark=1/ld,dYp=dark_shift,Yp_rel=rd['Yp']-runs[878.4]['Yp'],D_rel=rd['D']/runs[878.4]['D']-1)
# predictions anchored to the full codes
Y0,eY0=0.2467,0.0002
def pred(tau): return Y0+sY*(tau-PDG[0])
P={"bottle":pred(BOTTLE[0]),"beam":pred(BEAM[0]),"jparc":pred(JPARC[0]),"dark":Y0+out["dark"]["Yp_rel"]}
eP={"bottle":np.hypot(eY0,sY*BOTTLE[1]),"beam":np.hypot(eY0,sY*BEAM[1]),"jparc":np.hypot(eY0,sY*JPARC[1]),"dark":np.hypot(eY0,sY*BEAM[1])}
DATA={"four routes (Entry 102)":(0.24572,0.00126),"LBT 2026 alone":(0.2458,0.0013)}
out["pred"]=P; out["epred"]=eP; out["vs"]={}
print("\nPredicted primordial helium, N_eff = 3.044, and the distance from the measurements:")
for k in P:
    row={}
    for dk,(y,e) in DATA.items(): row[dk]=(y-P[k])/np.hypot(e,eP[k])
    out["vs"][k]=row; print(f"   {k:7s} Yp = {P[k]:.5f} ± {eP[k]:.5f}   "+"; ".join(f"{dk}: {v:+.1f}σ" for dk,v in row.items()))
# helium as a neutron clock
out["clock"]={}
for dk,(y,e) in DATA.items():
    t=PDG[0]+(y-Y0)/sY; et=np.hypot(e,eY0)/abs(sY)
    out["clock"][dk]=[t,et]; print(f"\nHelium as a clock ({dk}): tau = {t:.1f} ± {et:.1f} s;  bottle {(BOTTLE[0]-t)/np.hypot(et,BOTTLE[1]):+.1f}σ, beam {(BEAM[0]-t)/np.hypot(et,BEAM[1]):+.1f}σ, J-PARC {(JPARC[0]-t)/np.hypot(et,JPARC[1]):+.1f}σ")
# with N_eff free: add the CMB's own N_eff (Planck + ACT DR6 + DESI: 2.86 +- 0.13) and deuterium
# deuterium: measured 2.527 +- 0.030 (e-5); prediction at standard inputs 2.51 +- 0.07 (spread of full codes, nuclear rates); fit lnD only for N_eff
NE=(2.86,0.13); Dm=(2.527,0.030); Dp=(2.51,0.07)
taus=np.linspace(840,920,801); nes=np.linspace(2.3,3.6,651); best=None; chi=np.zeros((len(taus),len(nes)))
for i,t in enumerate(taus):
    yv=Y0+sY*(t-PDG[0])+sN*(nes-3.044); dv=Dp[0]*(1+sD*(t-PDG[0])+sDN*(nes-3.044))
    chi[i]=((DATA["four routes (Entry 102)"][0]-yv)/np.hypot(0.00126,eY0))**2+((nes-NE[0])/NE[1])**2+((Dm[0]-dv)/np.hypot(Dm[1],Dp[1]))**2
i,j=np.unravel_index(np.argmin(chi),chi.shape); pt=chi.min(axis=1)-chi.min()
ok=taus[pt<1]; out["free_neff"]=dict(tau=float(taus[i]),neff=float(nes[j]),lo=float(ok.min()),hi=float(ok.max()),
    chi_bottle=float(np.interp(BOTTLE[0],taus,pt)),chi_beam=float(np.interp(BEAM[0],taus,pt)))
print(f"\nN_eff free (CMB 2.86 ± 0.13, deuterium, helium): tau = {taus[i]:.0f} (+{ok.max()-taus[i]:.0f} −{taus[i]-ok.min():.0f}) s, N_eff = {nes[j]:.2f};"
      f" Δχ² beam {out['free_neff']['chi_beam']:.2f}, bottle {out['free_neff']['chi_bottle']:.2f}")
# how precise must helium be to tell the two lifetimes apart at 3 sigma (N_eff = 3.044)
need=abs(sY*(BEAM[0]-BOTTLE[0]))/3; out["need_3sigma"]=need
print(f"\nThe two lifetimes differ by {sY*(BEAM[0]-BOTTLE[0]):.4f} in Yp; telling them apart at 3σ needs Yp to ±{need:.4f} with the theory at ±0.0002")
# ---- heat after nucleosynthesis: the photodisintegration window for lithium (Poulin & Serpico 2015)
# Photons between 1.59 MeV (7Be -> 3He + 4He) and 2.22 MeV (deuterium) can destroy 7Be and spare deuterium. Below the pair threshold
# m_e^2/(22 T) they are lost only by Compton scattering on electrons. Order of magnitude here: the chance that a photon breaks a 7Be nucleus
# before Compton scattering drops it below 1.59 MeV, and the energy that must be injected to destroy two thirds of the 7Be.
RE2=7.9407877e-26; ME=0.51099895; MU=3.016029*4.002603/(3.016029+4.002603)*931.494; Q=1.5866; ALPHA=1/137.036
def sig_kn(E):                                          # Klein-Nishina total cross-section, cm^2
    k=E/ME; return 2*np.pi*RE2*((1+k)/k**2*(2*(1+k)/(1+2*k)-np.log(1+2*k)/k)+np.log(1+2*k)/(2*k)-(1+3*k)/(1+2*k)**2)
def sig_be(Eg, S_keVb=0.34):                            # 7Be(g,a)3He from 3He(a,g)7Be ground-state capture by detailed balance
    E=Eg-Q
    if E<=0: return 0.0
    tpe=2*np.pi*4*ALPHA*np.sqrt(MU/(2*E)); scap=S_keVb*1e-3/E*np.exp(-tpe)*1e-24    # cm^2 (S in MeV b / E in MeV)
    return 0.25*(2*MU*E/Eg**2)*scap
XBE=3.9e-10; NE_NH=1.16; NH_NG=6.12e-10*0.753           # 7Be/H after nucleosynthesis; electrons and hydrogen per photon
out["loophole"]={}
for E0 in (1.8,2.0,2.2):
    P=XBE*sig_be(E0)/(NE_NH*sig_kn(E0))*1.5               # about 1.5 scatterings before falling below threshold
    Nph=0.65*XBE/P; zeta=Nph*NH_NG*E0                       # injected photons per H, and MeV per CMB photon
    out["loophole"][str(E0)]=dict(sig_be_mb=sig_be(E0)*1e27,P=P,photons_per_H=Nph,zeta_MeV=zeta)
    print(f"   E0 = {E0} MeV: sigma(7Be) = {sig_be(E0)*1e27*1e3:.1f} microbarn, Klein-Nishina {sig_kn(E0)*1e24:.3f} b; chance per photon {P:.1e};"
          f" photons per H {Nph:.1e}; energy {zeta:.1e} MeV per CMB photon")
print("   Poulin & Serpico's window: 1.3e-6 to 1e-3 MeV per CMB photon, lifetimes 1e4 to 1e6 s; the 2-sigma entropy bound 0.022 is not reached")
# can the plasma itself supply such photons? Its hard photons after 1e4 s come from captures still running (d+p -> 3He + 5.5 MeV gamma,
# 3He + 4He -> 7Be + gamma, t + 4He -> 7Li + gamma). Count every nucleus made after 1e4 s, an upper bound on the hard photons.
rs=checked(bn.run(eta10=ETA,nnu=3.0,tau=878.4,dense=True)); sol=rs['sol']
def at(t):
    y=sol.sol(t); Y=np.maximum(y[1:],0.0); return Y
t1=1e4; tend=rs['t_end']; Y1=at(t1); Y2=at(tend)
made=max(Y2[4]-Y1[4],0)+max(Y2[6]+Y2[7]-Y1[6]-Y1[7],0)+max(Y2[3]-Y1[3],0)       # 3He, A = 7 and 3H gained after 1e4 s, per baryon
T9_1=float(sol.sol(t1)[0])
out["plasma_photons"]=dict(T_keV=T9_1*86.17,per_baryon=float(made),per_H_needed=float(out["loophole"]["2.0"]["photons_per_H"]))
print(f"   the plasma at 1e4 s (T = {T9_1*86.17:.0f} keV) makes at most {made:.1e} capture photons per baryon from then on; the window needs {out['loophole']['2.0']['photons_per_H']:.0e} per H")
# curve for the figure: Yp against tau
tg=np.linspace(860,900,41); out["curve"]=dict(tau=tg.tolist(),Y=[float(pred(t)) for t in tg])
json.dump(out,open(os.path.join(HERE,"neutron_bbn.json"),"w"),indent=1)
