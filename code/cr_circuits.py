"""cr_circuits.py - Workbench Entry 134
Can cosmic circuits and their double layers make the cosmic rays? Voltage, power and reach against the measured spectrum.
Inputs (all from the cited papers or earlier entries):
  Galaxy's rotation EMF 2.7e15 V and circuit power 5.3e32 W (Entries 18 and 107).
  Carlqvist (1986): filamentary double layers reach ~1e11 eV (solar circuit) and ~1e14 eV (galactic circuit) for protons.
  Galactic cosmic-ray luminosity 6.6-8.1e40 erg/s (Strong et al. 2010).
  LHAASO proton knee near 3 PeV (Science Bulletin 2025). Auger dipole above 8 EeV, about 115 deg from the Galactic centre (ApJ 976, 48).
  Highest energy, 2.4e20 eV (Telescope Array 2023). M87 horizon field 1-30 G (EHT 2021).
Unipolar-inductor estimate for a spinning black hole: V ~ c B r_g (order of magnitude, maximal spin).
"""
import json, math
c, G, Msun = 2.998e8, 6.674e-11, 1.989e30
V_gal, P_gal = 2.7e15, 5.3e32
L_cr = (6.6e40*1e-7, 8.1e40*1e-7)                        # W
E_dl_sun, E_dl_gal = 1e11, 1e14
knee_p, E_max, E_dipole = 3e15, 2.4e20, 8e18
Z = {'p': 1, 'He': 2, 'O': 8, 'Fe': 26}
out = {'galaxy_reach_eV': {k: V_gal*z for k, z in Z.items()},
       'power_shortfall': [L/P_gal for L in L_cr],
       'dl_gal_below_knee': knee_p/E_dl_gal,
       'needed_V': {'p': E_max, 'Fe': E_max/26},
       'galaxy_short_by': {'p': E_max/V_gal, 'Fe': E_max/26/V_gal}}
rg = G*6.5e9*Msun/c**2
out['M87'] = dict(r_g_m=rg, V=[c*B*1e-4*rg for B in (1, 30)], Fe_eV=[26*c*B*1e-4*rg for B in (1, 30)])
out['BH_1e9_1e4G_V'] = c*1.0*G*1e9*Msun/c**2
out['ladder'] = [
    ['Solar double layer', 9, 11.0, 'acc'],
    ['Galactic double layer', 9, 14.0, 'acc'],
    ['Galaxy’s whole voltage', 9, round(math.log10(V_gal), 3), 'acc'],
    ['… iron nuclei (×26)', 9, round(math.log10(26*V_gal), 3), 'acc'],
    ['M87 black hole', 9, round(math.log10(out['M87']['V'][1]), 3), 'bh', round(math.log10(out['M87']['V'][0]), 3)],
    ['… iron nuclei (×26)', 9, round(math.log10(out['M87']['Fe_eV'][1]), 3), 'bh', round(math.log10(out['M87']['Fe_eV'][0]), 3)]]
out['marks'] = [['proton knee', 15.48], ['Auger dipole from', 18.9], ['highest measured', 20.38]]
json.dump(out, open('cr_circuits.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in out.items() if k not in ('ladder',)}, indent=1))
