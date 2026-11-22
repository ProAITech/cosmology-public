"""Workbench Entry 85: black holes or plasmoids? Ring sizes of Sgr A* and M87* against masses weighed by stars.
GR fixes the scale of a black hole's image with no free parameter: the shadow's diameter is 2*sqrt(27) = 10.39
gravitational radii GM/c^2 for a non-spinning hole and between about 9.6 and 10.4 for any spin and viewing angle;
the bright ring the EHT measures sits a few per cent outside it (EHT calibrates ring/theta_g ~ 10-11.5 with simulations).
theta_g = GM/(c^2 D). Masses here come from stellar orbits (Sgr A*) and stellar dynamics (M87), not from the EHT.
Run: python3 -I shadow_scale.py > shadow_scale.json
"""
import json, math
G, c, Msun, pc = 6.67430e-11, 2.99792458e8, 1.98847e30, 3.0856776e16
uas = math.pi/180/3600*1e-6
def theta_g(M, D_pc): return G*M*Msun/c**2/(D_pc*pc)/uas
out = {}
# Sgr A*: GRAVITY 2022 (A&A 657, L12): M = 4.297 +- 0.012 e6, R0 = 8277 +- 9 pc (stat); EHT 2022 ring 51.8 +- 2.3 uas
tg = theta_g(4.297e6, 8277)
out['SgrA'] = {'mass': 4.297e6, 'D_pc': 8277, 'theta_g_uas': round(tg, 3), 'ring_uas': [51.8, 2.3],
               'ring_over_theta_g': [round(51.8/tg, 2), round(2.3/tg, 2)], 'shadow_pred_a0_uas': round(2*math.sqrt(27)*tg, 1)}
# M87*: ring 43.9 +- 0.6 uas (2017-2021, EHT 2025); stellar-dynamical masses
ring, ering = 43.9, 0.6
m87 = []
for name, M, lo, hi, D in [('Gebhardt et al. 2011 (stars)', 6.6e9, 0.4e9, 0.4e9, 17.9e6), ('Liepold, Ma & Walsh 2023 (stars)', 5.37e9, 0.25e9, 0.37e9, 16.8e6),
                           ('Simon, Cappellari & Hartke 2024 (stars)', 8.7e9, 1.77e9, 1.77e9, 16.8e6), ('Walsh et al. 2013 (gas)', 3.5e9, 0.3e9, 0.9e9, 17.9e6)]:
    t = theta_g(M, D); r = ring/t
    # error of the ratio: mass error (asymmetric -> symmetrised) and ring error
    er = r*math.hypot(ering/ring, 0.5*(lo+hi)/M)
    m87.append({'source': name, 'mass': M, 'D_pc': D, 'theta_g_uas': round(t, 3), 'ring_over_theta_g': [round(r, 2), round(er, 2)], 'shadow_pred_a0_uas': round(2*math.sqrt(27)*t, 1)})
out['M87'] = {'ring_uas': [ring, ering], 'years': '2017, 2018, 2021', 'masses': m87}
out['mass_ratio_M87_to_SgrA'] = round(6.0e9/4.297e6)
out['GR_shadow_range_rg'] = [9.6, 10.39]
print(json.dumps(out))
