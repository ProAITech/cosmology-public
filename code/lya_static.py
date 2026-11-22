"""Workbench Entry 82: does the Lyman-alpha forest thicken with distance?
In a static universe that is the same everywhere, absorbers have the same column densities and line widths (in km/s)
at every distance, and the number per unit proper length is constant. Seen at observed wavelength, an absorber at
distance l absorbs at velocity coordinate v = c ln(lambda_obs/lambda_alpha) = c ln(1+z), so the number of absorbers per
unit velocity is  nu(z) = n (1+z) (dl/dz) / c.  For lines scattered at random (Poisson), the mean transmitted flux is
exp(-nu * W_v), so the effective optical depth tau_eff is proportional to nu exactly, saturation included:
    tired light / logarithmic law, d = (c/H0) ln(1+z):  dl/dz = (c/H0)/(1+z)  ->  tau_eff constant
    Lerner's law, d = cz/H0:                            dl/dz = c/H0          ->  tau_eff proportional to (1+z)
The measured tau_eff(z) = -ln <F> is compared with both.
Data: Kirkman et al. 2007 (HST FOS, low-density IGM, 0 < z < 1.6: 0.037 +- 0.004 absorbed at z = 1, power-law index 1.01);
Becker et al. 2013 (SDSS composites, 2.15 <= z <= 4.85, Table 2); Bosman et al. 2022 (XQR-30, 4.9 <= z <= 6.1, Table 4).
Run: python3 -I lya_static.py > lya_static.json
"""
import json, math
import numpy as np
kirk = [(z, 0.037*((1+z)/2.0)**1.01) for z in (0.2, 0.6, 1.0, 1.4)]          # D_A of the low-density IGM
kirk = [(z, -math.log(1-d), 0.21*(-math.log(1-d))) for z, d in kirk]          # tau_eff, ~21% error per interval
becker = [(2.15, 0.8806, 0.0103), (2.45, 0.7968, 0.0089), (2.95, 0.6966, 0.0084), (3.45, 0.5548, 0.0071),
          (3.95, 0.4255, 0.0071), (4.45, 0.3216, 0.0094), (4.85, 0.2225, 0.0151)]
bosman = [(4.9, 0.171, 0.014), (5.0, 0.1581, 0.0086), (5.2, 0.1222, 0.0050), (5.4, 0.0801, 0.0055), (5.6, 0.0447, 0.0035),
          (5.8, 0.0172, 0.0025), (6.0, 0.0089, 0.0031)]
def tau(rows): return [(z, -math.log(F), e/F) for z, F, e in rows]
B, X = tau(becker), tau(bosman)
out = {'kirkman_lowdensity_tau': kirk, 'becker_tau': B, 'bosman_tau': X}
def slope(rows):
    x = np.log([1+r[0] for r in rows]); y = np.log([r[1] for r in rows]); w = 1/np.array([r[2]/r[1] for r in rows])**2
    W = np.diag(w); A = np.vstack([np.ones_like(x), x]).T; cov = np.linalg.inv(A.T@W@A); p = cov@A.T@W@y
    return round(float(p[1]), 2), round(float(math.sqrt(cov[1, 1])), 2)
out['power_law_index'] = {'0<z<1.6 (Kirkman, low-density gas)': [1.01, None], '2.15-4.85 (Becker)': slope(B), '4.9-6.0 (Bosman)': slope(X)}
# what a static universe needs: the absorbing gas per unit length, relative to the anchor, so that tau matches
anchors = {'z=1 (HST)': kirk[2], 'z=2.15 (SDSS)': B[0]}
need = {}
for an, (za, ta, _) in anchors.items():
    rows = []
    for z, t, e in [r for r in kirk + B + X if r[0] > za]:
        rows.append([z, round(t, 4), round(t/ta, 2), round((t/ta)/((1+z)/(1+za)), 2)])
    need[an] = rows
out['static_needs_gas_factor'] = {'columns': ['z', 'tau_eff', 'needed by the logarithmic law (tau ratio)', 'needed by Lerner law (tau ratio over (1+z) ratio)'], 'by_anchor': need}
# distances in Lerner's law at the far end, for the text
out['lerner_distance_Gpc_at'] = {str(z): round(z*299792.458/70/1000, 1) for z in (2, 4, 6)}
print(json.dumps(out))
