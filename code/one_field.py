"""Workbench Entry 115: does a nearly scale-invariant field of ~45 pG clump recombination enough to ease the Hubble tension?

Entry 14 found one primordial field that could do several jobs: a nearly scale-invariant spectrum normalised to about 45 pG on 1 Mpc reaches
the seed LOFAR's filaments prefer, supplies the Milky Way's vertical flux (Entry 54) and fills the voids (Entry 109). Its role in the Hubble
tension rests on clumping. Clumping is made on scales below the photon mean free path at recombination, and Jedamzik & Saveliev (2019) found
the scale inside that range does not matter. So what counts is the field's power between the dissipation scale and the photon mean free path.

For a spectrum P_B ~ k^n, the energy per ln k goes as k^(n+3). Normalised the Planck way (Gaussian smoothing on 1 Mpc), the power in the
clumping band depends steeply on how close n is to -3. The clumping strength b follows from Jedamzik & Saveliev's simulations,
b = 0.15 at 52.5 pG and b < 0.119 at 47 pG: locally b ~ B^2.1, steepening to B^4 for weak fields. H0 follows from the notebook's track
(Entry 13, clump_track.json). Limits: the 2026 three-telescope bound read in the notebook's model, b < 0.08 to 0.21 (Entry 105)."""
import os, json
import numpy as np
from scipy import integrate
HERE=os.path.dirname(os.path.abspath(__file__))
TRACK=json.load(open(os.path.join(HERE,"..","clump","clump_track.json"))) if os.path.exists(os.path.join(HERE,"..","clump","clump_track.json")) else json.load(open(os.path.join(HERE,"clump_track.json")))
tb=np.array([t['b'] for t in TRACK]); tH=np.array([t['H0'] for t in TRACK])
H0_of=lambda b: float(np.interp(b,tb,tH))
# photon mean free path at recombination, comoving (x_e = 1 just before, z = 1100); dissipation scale ~0.5 kpc (Jedamzik & Saveliev)
nb0=2.51e-7; Y=0.245; z=1100; ne=nb0*(1-Y/2)*(1+z)**3; sigT=6.6524e-25; MPC=3.0857e24
l_gamma=(1+z)/(ne*sigT)/MPC                                # comoving Mpc
out={"l_gamma_Mpc":l_gamma}
print(f"photon mean free path before recombination: {l_gamma*1e3/(1+z):.2f} kpc physical, {l_gamma:.2f} Mpc comoving")
k_IR=2*np.pi/14000.0                                       # largest scale: today's horizon
def band_ratio(n, l_g=l_gamma, l_d=0.5e-3, lam=1.0):
    """(field in the clumping band / field smoothed on lam)^2 for P_B ~ k^n."""
    p=n+3; kg=2*np.pi/l_g; kd=2*np.pi/l_d
    num=integrate.quad(lambda lk: np.exp(p*lk),np.log(kg),np.log(kd))[0]
    den=integrate.quad(lambda lk: np.exp(p*lk-np.exp(2*lk)*lam*lam),np.log(k_IR),np.log(50/lam),limit=200)[0]
    return num/den
def b_of(Bc,pexp):
    return 0.119*(Bc/47.0)**pexp
out["rows"]=[]
print("\nField 45 pG on 1 Mpc. Clumping-band field, clumping b (local exponent 2.1 / weak-field 4) and H0 along the track:")
for n in (-3.2,-3.1,-3.05,-3.0,-2.95,-2.9,-2.8):
    rlo,rmid,rhi=[np.sqrt(band_ratio(n,l_g=lg,l_d=ld)) for lg,ld in ((3.0,2e-3),(l_gamma,0.5e-3),(1.0,0.2e-3))]
    Bc=45*rmid; b2=b_of(Bc,2.1); b4=b_of(Bc,4.0)
    row=dict(n=n,ratio=[rlo,rmid,rhi],Bc=Bc,b=[b4,b2],H0=[H0_of(b4),H0_of(b2)],Bc_range=[45*rlo,45*rhi])
    out["rows"].append(row)
    print(f"   n = {n:5.2f}: band field {Bc:5.1f} pG ({45*rlo:4.1f}-{45*rhi:5.1f} for the band edges), b = {b4:.3f}-{b2:.3f}, H0 = {H0_of(b4):.2f}-{H0_of(b2):.2f}")
# what the 2026 limit allows for a 45 pG field: the steepest tilt with b below 0.08 or 0.21
for blim in (0.08,0.21):
    for pexp in (2.1,4.0):
        Bmax=47*(blim/0.119)**(1/pexp)
        ns=np.linspace(-3.3,-2.6,701); r=np.array([np.sqrt(band_ratio(n)) for n in ns]); ok=ns[45*r<Bmax]
        out[f"nmax_b{blim}_p{pexp}"]=float(ok.max()); out[f"Bmax_b{blim}_p{pexp}"]=Bmax
        print(f"   b < {blim} (exponent {pexp}): band field below {Bmax:.1f} pG, so a 45 pG field needs n below {ok.max():.3f}")
# the needed H0: what b and band field would reach 71 (the direct one-step average 71.4 +- 1.2, Entry 56)
b71=float(np.interp(71.0,tH,tb)); out["b_for_71"]=b71
print(f"\nReaching H0 = 71 needs b = {b71:.2f}, beyond both readings of the 2026 limit; a 45 pG field reaches it only for n above "
      f"{[r['n'] for r in out['rows'] if r['b'][1]>=b71][0] if any(r['b'][1]>=b71 for r in out['rows']) else 'none listed'} (exponent 2.1)")
# curve for the figure
ng=np.round(np.arange(-3.20,-2.849,0.01),3); cv=[]
for n in ng:
    Bc=45*np.sqrt(band_ratio(n)); b2=b_of(Bc,2.1); b4=b_of(Bc,4.0); bl=min(b2,b4); bh=max(b2,b4) if Bc<47 else b2
    cv.append([float(n),float(Bc),float(min(H0_of(bl),72.3)),float(min(H0_of(bh),72.3)),float(bl),float(bh)])
out["curve"]=cv
json.dump(out,open(os.path.join(HERE,"one_field.json"),"w"),indent=1)
