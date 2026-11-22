"""massflux.py - Workbench Entry 97
Where does magnetism beat gravity? The mass-to-flux ratio from interstellar clouds to the cosmic web.
A field B threading a layer of surface density Sigma can hold it up against its own gravity only if
   Sigma < Sigma_crit = B / (2 pi sqrt(G))        (Nakano & Nakamura 1978)
   lambda = Sigma / Sigma_crit = Sigma [Msun/pc^2] / (2.95 B [uG])  (= 3.8e-21 N_H / B, Crutcher's form with 1.4 m_H per H)
lambda < 1: the field can control the gas (subcritical); lambda > 1: gravity wins (supercritical).
For each class of object: measured fields and measured surface densities, for the gas alone and for all the gravitating mass.
Plasma cosmology scales laboratory plasma up to galaxies (Peratt), which needs lambda well below 1 on galactic scales.
Values and sources are in ROWS; ranges are carried through as low/high.
"""
import json, math, os
HERE = os.path.dirname(os.path.abspath(__file__))
G = 6.674e-8; MSUN_PC2 = 1.989e33/(3.0857e18)**2           # g/cm^2 per Msun/pc^2
K = 1e-6/(2*math.pi*math.sqrt(G))/MSUN_PC2                 # Sigma_crit for 1 uG, in Msun/pc^2
MH = 1.6735e-24*1.4/MSUN_PC2                               # Msun/pc^2 per hydrogen atom per cm^2 (with helium)
h = 0.70
# name, B (uG) central, low, high; Sigma_gas central, low, high; Sigma_total central, low, high (Msun/pc^2; None if not separable); source note
ROWS = [
 ('Cold atomic clouds', 6.0, 3.0, 10.0, 1e20.__mul__(3)*MH, 1e20*MH, 1e21*MH, None, None, None,
  'Heiles & Troland 2005: median field 6 uG in cold HI; N_H 1e20-1e21 cm^-2'),
 ('Molecular clouds', 50.0, 15.0, 300.0, 3e22*MH, 1e22*MH, 1e23*MH, None, None, None,
  'Crutcher & Kemball 2019: N_H >~ 1e22 cm^-2, B ~ 10 uG (n/300)^0.65; on average 2-3 times supercritical'),
 ('Galactic disc near the Sun', 6.0, 5.0, 7.0, 13.7, 12.0, 15.5, 68.0, 64.0, 72.0,
  'Field ~6 uG (Entry 16); gas 13.7 Msun/pc^2 (McKee et al. 2015); all mass within 1.1 kpc of the plane 68 +- 4 (Bovy & Rix 2013)'),
 ('Starburst centres', 75.0, 50.0, 100.0, 3e3, 1e3, 1e4, None, None, None,
  'Beck 2015: 50-100 uG in central starbursts; Kennicutt 1998: starburst gas 1e3-1e4 Msun/pc^2 (stars add more)'),
 ('Milky Way halo, out to 100 kpc', 0.73, 0.3, 1.0, None, None, None, 4.1e11/(math.pi*1e10), 2.75e11/(math.pi*1e10), 5.4e11/(math.pi*1e10),
  'Halo field 0.73 uG (Entry 16); mass within ~100 kpc 2.75e11 (Jiao et al. 2023) to 5.4e11 Msun (Wang et al. 2022)'),
 ('Rich cluster (Coma)', 2.0, 1.0, 4.7, 0.13*1.88e15/h/(math.pi*(1.99e6/h)**2), 0.10*1.88e15/h/(math.pi*(1.99e6/h)**2), 0.16*1.88e15/h/(math.pi*(1.99e6/h)**2),
  1.88e15/h/(math.pi*(1.99e6/h)**2), 1.32e15/h/(math.pi*(1.99e6/h)**2), 2.53e15/h/(math.pi*(1.99e6/h)**2),
  'Bonafede et al. 2010: 4.7 uG at the centre, falling outward; Kubo et al. 2007: M200 = 1.88e15/h Msun inside 1.99/h Mpc; gas 10-16% of the mass'),
 ('Cosmic-web filaments', 0.03, 0.011, 0.06, 0.157*1.6e13/(7.1*2.5/h**2*1e12), 0.157*1.3e13/(7.1*2.5/h**2*1e12), 0.157*1.9e13/(7.1*2.5/h**2*1e12),
  1.6e13/(7.1*2.5/h**2*1e12), 1.3e13/(7.1*2.5/h**2*1e12), 1.9e13/(7.1*2.5/h**2*1e12),
  'Fields 11-15 nG (Carretti et al. 2025) to 30-60 nG (Vernstrom et al. 2021); 1.6e13 Msun in 7.1 x 2.5 /h Mpc (Epps & Hudson 2017); baryons 15.7% of it'),
]
PERATT = ('Peratt-type galaxy filament', 100.0, 30.0, 300.0, 1e11/35e3**2, 3e10/35e3**2, 3e11/35e3**2,
          'Peratt, Peter & Snell 1990: filaments ~35 kpc wide with ~100 uG fields; a galaxy of 1e11 Msun (range 3e10-3e11) across that width')
lam = lambda S, B: S/(K*B)
out = {'K_Msun_pc2_per_uG': round(K, 3), 'rows': []}
for r in ROWS:
    name, B, Bl, Bh, Sg, Sgl, Sgh, St, Stl, Sth, note = r
    row = dict(name=name, B=[B, Bl, Bh], note=note)
    if Sg is not None: row.update(Sgas=[round(Sg, 4), round(Sgl, 4), round(Sgh, 4)], lam_gas=[round(lam(Sg, B), 3), round(lam(Sgl, Bh), 3), round(lam(Sgh, Bl), 3)])
    if St is not None: row.update(Stot=[round(St, 4), round(Stl, 4), round(Sth, 4)], lam_tot=[round(lam(St, B), 3), round(lam(Stl, Bh), 3), round(lam(Sth, Bl), 3)])
    out['rows'].append(row)
    print('%-32s B %-6g gas lambda %-24s total lambda %s' % (name, B, row.get('lam_gas'), row.get('lam_tot')))
n, B, Bl, Bh, S, Sl, Sh, note = PERATT
out['peratt'] = dict(name=n, B=[B, Bl, Bh], Stot=[round(S, 2), round(Sl, 2), round(Sh, 2)], lam_tot=[round(lam(S, B), 3), round(lam(Sl, Bh), 3), round(lam(Sh, Bl), 3)], note=note,
                     field_needed_for_lambda1_at_measured_halo=round(4.1e11/(math.pi*1e10)/K, 1))
print('Peratt', out['peratt'])
# plasma that is magnetically dominated, for scale: Earth's plasma sheet and the solar corona
out['controls'] = {'plasma sheet (n 0.3 cm^-3 over 10 Earth radii, 20 nT)': lam(0.3*6.4e9*MH, 200.0),
                   'solar corona (n 1e8 cm^-3 over 0.1 solar radius, 5 G)': lam(1e8*7e9*MH, 5e6)}
print({k: '%.1e' % v for k, v in out['controls'].items()})
json.dump(out, open(os.path.join(HERE, 'massflux.json'), 'w'), indent=1)
