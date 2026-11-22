"""void_fields.py - Workbench Entry 109: what could fill the cosmic voids with the field the blazars require?

Lower bound (Fermi-LAT & H.E.S.S. 2023, missing GeV cascades of TeV blazars): B > 7.1e-16 G for a 1 Mpc coherence length with
blazar duty cycles as short as 10 yr; > 1.8e-14 G for 1e4 yr; > 3.9e-14 G for 1e7 yr.
Upper bounds on a primordial field (Jedamzik & Saveliev 2019, recombination clumping): 47 pG scale-invariant, 8.9 pG Batchelor.
Candidates: the Biermann battery (Kulsrud et al. 1997, ~1e-21 G by z ~ 3 in collapsing regions); a primordial field; galactic
outflows (volume filling <~20%, Seller & Sigl 2025 summarising Garcia et al. 2021 and others); currents along cosmic filaments
(toroidal field B = mu0 I/(2 pi r)); the notebook's filament-current limit, 3e17 A (Entries 68, 75).
"""
import json, os, math
HERE = os.path.dirname(os.path.abspath(__file__))
MU0 = 4e-7*math.pi; MPC = 3.0857e22; G_T = 1e-4                  # 1 G = 1e-4 T
bounds = {'10 yr': 7.1e-16, '1e4 yr': 1.8e-14, '1e7 yr': 3.9e-14}
def I_needed(BG, r_mpc): return 2*math.pi*r_mpc*MPC*BG*G_T/MU0      # A
def B_at(I, r_mpc): return MU0*I/(2*math.pi*r_mpc*MPC)/G_T          # G
out = {'bounds': bounds, 'cmb_primordial_upper': {'scale-invariant': 47e-12, 'Batchelor': 8.9e-12}, 'biermann': 1e-21}
out['current_needed'] = {k: {f'{r} Mpc': float('%.2g' % I_needed(b, r)) for r in (5, 10, 20)} for k, b in bounds.items()}
Ilim = 3e17
out['filament_limit'] = dict(I=Ilim, B_at_10Mpc=float('%.2g' % B_at(Ilim, 10)), B_at_5Mpc=float('%.2g' % B_at(Ilim, 5)),
                             factor_above_needed=float('%.2g' % (Ilim/I_needed(bounds['10 yr'], 10))))
out['biermann_shortfall'] = float('%.2g' % (bounds['10 yr']/1e-21))
out['primordial_window_decades'] = dict(scale_invariant=round(math.log10(47e-12/bounds['10 yr']), 2), batchelor=round(math.log10(8.9e-12/bounds['10 yr']), 2))
# what a current of the needed size involves: drift speed and driving voltage
I = I_needed(bounds['10 yr'], 10); a = 1*MPC; ne = 10.0                 # filament radius 1 Mpc, n_e = 1e-5 cm^-3
J = I/(math.pi*a**2); vd = J/(ne*1.602e-19)
ell = 30*MPC; L = MU0*ell/(2*math.pi)*math.log(10)                     # inductance of a 30 Mpc line out to ten times its radius
V = L*I/(1e9*3.156e7)
out['needed_current_physics'] = dict(I=float('%.2g' % I), J=float('%.2g' % J), drift_m_s=float('%.2g' % vd), inductance_H=float('%.2g' % L), emf_V_1Gyr=float('%.2g' % V),
                                     milky_way_emf_V=2.7e15)
# magnetic energy density at the bound, against the CMB
uB = (bounds['10 yr'])**2/(8*math.pi); ucmb = 4.17e-13
out['energy'] = dict(uB_erg_cm3=float('%.2g' % uB), over_cmb=float('%.2g' % (uB/ucmb)))
# field contrast between a sight line near a filament and the middle of a void, for a current-made field
out['contrast_2_vs_15_Mpc'] = round(15/2, 1)
json.dump(out, open(os.path.join(HERE, 'void_fields.json'), 'w'), indent=1)
print(json.dumps(out, indent=1))
