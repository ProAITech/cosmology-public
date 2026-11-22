"""Workbench Entry 72: the 'electric comet' (a rocky body discharging in the Sun's field, its water made from solar-wind protons
and oxygen stripped from surface minerals) against Rosetta's measurements of comet 67P at perihelion (13 August 2015, 1.24 au)."""
import json, math, os
HERE = os.path.dirname(os.path.abspath(__file__))
au = 1.24; mp = 1.67262192e-27; e = 1.602176634e-19
V = 18.7e9; R = (3 * V / (4 * math.pi)) ** (1 / 3); A = math.pi * R**2          # nucleus volume (OSIRIS), equal-volume radius
Q_kg = 1e8 / 86400; m_w = 18.015 * 1.66054e-27; Q = Q_kg / m_w                    # ~1e8 kg/day near perihelion (Hansen et al. 2016)
n = 5e6 / au**2; v = 4e5                                                         # solar wind at 1 au: 5 per cm3, 400 km/s
flux = n * v; Nsw = flux * A
print(f"Nucleus: equal-volume radius {R/1e3:.2f} km, cross-section {A/1e6:.1f} km2")
print(f"Water production near perihelion: {Q_kg:.0f} kg/s = {Q:.1e} molecules/s")
print(f"\n1. Hydrogen supply. Solar-wind protons striking the nucleus: {Nsw:.1e}/s, enough for at most {Nsw/2:.1e} water molecules/s")
print(f"   shortfall: {Q/(Nsw/2):.0e}.  Even every proton crossing a disc 1,000 km across: {flux*math.pi*(5e5)**2/2:.1e}/s, still {Q/(flux*math.pi*(5e5)**2/2):.0e} short")
Lsub = 2.83e6; P_sub = Q_kg * Lsub
S = 1361 / au**2; P_sun = S * A * (1 - 0.06)
P_sw = 0.5 * n * mp * v**3 * A
print(f"\n2. Energy. Sublimating the water takes {P_sub:.1e} W.")
print(f"   sunlight absorbed by the nucleus: {P_sun:.1e} W ({100*P_sub/P_sun:.0f}% of it would do)")
print(f"   solar-wind kinetic power on the nucleus: {P_sw:.1e} W, {P_sub/P_sw:.0e} times too little; a 3 kV drop (Entry 71) on those protons: {Nsw*e*3e3:.1e} W")
DH_c = (5.3e-4, 0.7e-4); DH_ps = (2.1e-5, 0.5e-5)
print(f"\n3. Isotopes. Water D/H measured {DH_c[0]:.1e} +- {DH_c[1]:.0e}; hydrogen from the Sun carries at most the protosolar {DH_ps[0]:.1e}")
print(f"   (the Sun has burned its deuterium, so the wind carries less). Water built from it would be >= {DH_c[0]/DH_ps[0]:.0f} times too poor in deuterium")
rho = 533.0
for rr in (2500, 3000, 3500):
    print(f"4. Density {rho:.0f} kg/m3: a body of rock at {rr} kg/m3 would need {100*(1-rho/rr):.0f}% empty space")
out = dict(R=R, A=A, Q=Q, Q_kg=Q_kg, shortfall_H=Q / (Nsw / 2), P_sub=P_sub, P_sun=P_sun, frac_sun=P_sub / P_sun, P_sw=P_sw, shortfall_E=P_sub / P_sw,
           DH_ratio=DH_c[0] / DH_ps[0], porosity_rock=[1 - rho / r for r in (2500, 3000, 3500)])
json.dump(out, open(os.path.join(HERE, "electric_comet.json"), "w"), indent=1)
