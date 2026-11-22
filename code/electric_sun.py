"""Workbench Entry 71: could the Sun be powered from outside by an electric current (Juergens 1972, 1982; the 'electric Sun')?
Whatever the circuit, current I and voltage V must satisfy I V = L.  Two measurements bound V and the particle flux; a third,
solar neutrinos, asks where the energy is made."""
import json, math, os
HERE = os.path.dirname(os.path.abspath(__file__))
L = 3.828e26; AU = 1.495978707e11; e = 1.602176634e-19; c = 2.998e8; MeV = 1.602176634e-13; RE = 6.371e6
A1 = 4 * math.pi * AU**2; S = L / A1
out = {"solar_constant": S}
print(f"Energy that must cross every sphere around the Sun, per m2 at 1 au: {S:.0f} W/m2")
# --- Juergens's own numbers
I_J, V_J = 4e16, 1e10
j = I_J / A1; flux = j / e
print(f"\n1. Juergens's circuit: I = {I_J:.0e} A at V = {V_J:.0e} V gives I V = {I_J*V_J:.1e} W (L = {L:.2e})")
print(f"   current density at 1 au {j:.2e} A/m2 = {flux:.1e} electrons m^-2 s^-1 = {flux/1e4:.1e} cm^-2 s^-1")
print(f"   if the volt drop is in the outer heliosphere (as he places it), each arrives with ~{V_J/1e9:.0f} GeV: energy flux {flux*V_J*e:.0f} W/m2")
P_earth = math.pi * RE**2 * flux * V_J * e
print(f"   Earth would intercept {P_earth:.1e} W of 10-GeV electrons; all auroral precipitation is ~1e10-1e11 W")
u_cr = 1.0 * e * 1e6           # J/m3: about 1 eV/cm3, all Galactic cosmic rays, Voyager 1 beyond the heliopause
F_cr = c * u_cr / 4            # one-sided energy flux of an isotropic population
print(f"   all Galactic cosmic rays together carry ~{F_cr:.1e} W/m2 through a surface; the circuit needs {S/F_cr:.0e} times that, in electrons alone")
out["juergens"] = dict(j=j, flux=flux, P_earth=P_earth, ratio_to_all_CR=S / F_cr)
# --- if the drop is near the Sun instead, the solar wind protons fall through it
print("\n2. If the voltage drop is near the Sun, protons leaving it fall through it.")
mp = 1.67262192e-27
for v in (300e3, 450e3, 800e3):
    print(f"   solar wind at {v/1e3:.0f} km/s: {0.5*mp*v*v/e/1e3:.2f} keV per proton")
Vmax = 0.5 * mp * (800e3)**2 / e
Imin = L / Vmax; jmin = Imin / A1
print(f"   so the potential the wind falls through is at most ~{Vmax/1e3:.1f} kV (Parker Solar Probe: 50-300 V beyond 0.1-0.2 au).")
print(f"   with V <= {Vmax/1e3:.1f} kV the current must be >= {Imin:.1e} A: {jmin:.2f} A/m2 at 1 au = {jmin/e/1e4:.1e} charges cm^-2 s^-1,")
nv = 5e6 * 450e3
print(f"   against a whole solar-wind particle flux of {nv/1e4:.1e} cm^-2 s^-1: {jmin/e/nv:.0e} times larger")
out["near"] = dict(Vmax=Vmax, Imin=Imin, jmin=jmin, ratio_to_wind=jmin / e / nv)
# --- neutrinos
print("\n3. Where the energy is made: pp neutrinos")
Q_eff = 26.73 - 0.53          # MeV per helium nucleus left as heat and light, pp-I (two neutrinos carry ~0.53 MeV)
n_he = L / (Q_eff * MeV); n_pp = 2 * n_he
F_pred = n_pp / (A1 * 1e4)
print(f"   L / {Q_eff:.1f} MeV = {n_he:.2e} helium per s, so {n_pp:.2e} pp neutrinos per s: {F_pred:.2e} cm^-2 s^-1 at Earth")
print("   Borexino measures 6.1 +- 0.5 e10 cm^-2 s^-1 (pp, below 0.42 MeV); neutrino-inferred luminosity 3.89e33 erg/s vs 3.846e33 in light")
out["pp_pred"] = F_pred
# Could beam fusion at the surface make them?  Weak p+p cross-section, largest possible value
S0 = 4.01e-25                  # MeV barn (Adelberger et al. 2011)
EG = 0.4935                    # MeV, Gamow energy for p + p
Ebest = EG / 4; sig_max = S0 / Ebest * math.exp(-2) * 1e-24      # cm2, at the peak of S/E exp(-sqrt(EG/E))
col = 1.4e21                   # H atoms per cm2: about the stopping column of a 1 MeV proton in hydrogen (2.3e-3 g/cm2)
p_max = sig_max * col
need = n_pp / p_max
print(f"   largest weak p+p cross-section {sig_max:.1e} cm2 (at {Ebest*1e3:.0f} keV); a beam proton stopping in hydrogen fuses with probability <= {p_max:.0e}")
print(f"   to make {n_pp:.1e} reactions per s a surface beam needs >= {need:.1e} protons per s, carrying >= {need*Ebest*MeV/L:.0e} times the Sun's luminosity")
out["beam"] = dict(sig_max=sig_max, p_max=p_max, protons_per_s=need, power_over_L=need * Ebest * MeV / L)
json.dump(out, open(os.path.join(HERE, "electric_sun.json"), "w"), indent=1)
