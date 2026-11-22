"""Workbench Entry 23: could stars have made the helium?  An energy account.
Lerner's Galactic Origin of Light Elements (GOLE) has all helium (mass fraction 0.248) made by stars in young galaxies, within about
200 Myr, and has the released energy, thermalised by an intergalactic fog, become the CMB. Entry 22 found no fog.
So where would the light be, and how much of it is there?"""
import math, json
eV=1.602176634e-19; c=2.99792458e8; Mpc_cm=3.0857e24; Gyr=3.15576e16; Lsun=3.828e26; Msun=1.98847e30
out={}
# --- energy released per nucleon turned into helium
Q_He=26.731        # MeV per 4He nucleus formed from 4 protons
nu_loss=0.02       # share carried off by neutrinos (pp chain)
q=Q_He/4*(1-nu_loss)                       # MeV of light per nucleon
print(f"Light released per nucleon fused to helium: {q:.2f} MeV ({100*q/931.494:.2f}% of its rest energy)")
Y=0.248            # GOLE's helium mass fraction (0.248 +- 0.009); measured values are 0.245-0.250
# --- baryon density: Lerner's own figure and the standard one
n_lerner=1.2e-7/(1-Y)                      # nucleons per cm^3, from his 'observed H density about 1.2e-7 per cm^3'
n_std=0.0224*1.87834e-29/1.6726e-24        # Omega_b h^2 = 0.0224
u_cmb=7.5657e-15*2.7255**4/1e7*1e7/1.602176634e-12   # erg/cm^3 -> eV/cm^3
u_cmb=7.5657e-15*2.7255**4/1.602176634e-12
I_ebl,dI=50.0,6.4                          # nW m^-2 sr^-1: 24+-4 optical plus 26+-5 infrared, from galaxy counts (Driver et al. 2016)
u_ebl=4*math.pi*I_ebl*1e-9/c/eV/1e6; du_ebl=u_ebl*dI/I_ebl
print(f"CMB energy density: {u_cmb:.3f} eV/cm^3")
print(f"All starlight and dust glow, 0.1-1000 micron: {I_ebl:.0f} +- {dI:.0f} nW/m^2/sr = {u_ebl:.4f} +- {du_ebl:.4f} eV/cm^3  ({100*u_ebl/u_cmb:.1f}% of the CMB)")
rows=[]
for name,n in (("Lerner's hydrogen density",n_lerner),("standard baryon density",n_std)):
    u=n*Y*q*1e6
    print(f"\n{name}: {n:.2e} nucleons/cm^3 -> energy released in making the helium {u:.3f} eV/cm^3")
    print(f"   = {u/u_cmb:.2f} x the CMB;  {u/u_ebl:.0f} x the starlight background")
    # at today's total stellar output
    j=1.6e35                                # W/Mpc^3, 0.1-1000 micron (Driver et al. 2008)
    t=u*eV*Mpc_cm**3/j/Gyr
    print(f"   at today's output of all stars ({j:.1e} W/Mpc^3) it would take {t:.0f} Gyr, {t/13.97:.0f} Hubble times")
    # light emitted when the oldest stars formed, 13 Gyr ago, in a static universe where light loses energy with H0 = 70
    for label,f in (("exponential loss",math.exp(13/13.97)),("linear distance law",1+13/13.97)):
        print(f"   emitted 13 Gyr ago and degraded ({label}, factor {f:.2f}): {u/f:.3f} eV/cm^3 = {u/f/u_ebl:.0f} x the starlight background")
    rows.append({"name":name,"n":n,"u":u,"x_cmb":u/u_cmb,"x_ebl":u/u_ebl,"t_Gyr":t,"u_13Gyr":[u/math.exp(13/13.97),u/(1+13/13.97)]})
    print(f"   helium that the starlight background can pay for: Y = {Y*u_ebl/u:.3f} (no energy loss), {Y*2*u_ebl/u:.3f} if the light has lost half its energy")
out["rows"]=rows; out["u_cmb"]=u_cmb; out["u_ebl"]=[u_ebl,du_ebl]
# --- the bright phase GOLE requires
LM=Y*q/931.494*c**2/(0.2*Gyr)/(Lsun/Msun)
print(f"\nMaking Y = {Y} in 200 Myr needs {LM:.0f} Lsun per Msun of gas; a galaxy of 6e10 Msun would shine at {LM*6e10:.1e} Lsun for 200 Myr")
out["L_over_M"]=LM; out["L_MW"]=LM*6e10
json.dump(out,open("helium_budget.json","w"),indent=1)
