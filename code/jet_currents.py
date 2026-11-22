"""jet_currents.py - Workbench Entry 142
How much current do black-hole jets carry? From one jet (3C 303) to 102: core-shift fields turned into currents.
Data, read in a browser from VizieR and copied below exactly as served:
  MOJAVE IX (Pushkarev et al. 2012, J/A+A/545/A113/table4): source, optical class (Q quasar, B BL Lac), z, fastest apparent speed
  beta_app, upper-limit flag on B1, B1 = magnetic field at 1 pc (G) from the frequency-dependent core shift under equipartition,
  B1 ~ 0.042 Omega_rnu^(3/4) (1+z)^(1/2) (1+beta_app^2)^(1/8), with half-opening angle phi = 0.6 deg and theta ~ 1/Gamma assumed.
  MOJAVE XIV (Pushkarev et al. 2017, J/MNRAS/468/4992/table4), matched by B1950 name: intrinsic opening angle (deg), Doppler factor,
  Lorentz factor, viewing angle (deg); blank where not measured.
Current: if the field at 1 pc is mostly toroidal (rotation-measure gradients across many jets say it has a large toroidal part),
  Ampere's law gives the current enclosed within the jet radius r = 1 pc x tan(phi):
     I' = 2 pi r B1 / mu0           (B1 is the field in the jet's own frame)
     I  = Gamma I'                  (a toroidal field is larger by Gamma in the galaxy's frame, B_lab = Gamma B')
  Main reading: the catalogue's own assumptions, phi = 0.6 deg and Gamma = (1 + beta_app^2)^(1/2).
  Second reading (55 jets): measured intrinsic opening angle phi = alpha_int/2, B1 rescaled as phi^(-1/4), and the fitted Gamma.
Power: the Poynting flux of a relativistic current-carrying jet, P = (mu0 c / 4 pi) I^2 = 30 ohm x I^2.
"""
import json, math, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
RAW = """0003-066,B,0.347,2.89,<,0.20,3.2,5.1,6.4,11.1;0106+013,Q,2.099,26.50,<,0.79,1.2,18.2,12.4,2.8;0119+115,Q,0.570,17.10,,1.63,,,,;0133+476,Q,0.859,12.98,,0.77,1.4,20.5,16.9,2.7;0149+218,Q,1.320,18.55,,1.69,,,,;0202+149,Q,0.405,6.41,,0.48,1.4,15.0,15.9,3.8;0212+735,Q,2.367,7.64,,1.27,1.4,8.4,6.8,6.6;0215+015,Q,1.715,34.16,,1.41,1.2,27.1,18.2,1.8;0224+671,Q,0.523,11.63,,0.75,,,,;0234+285,Q,1.207,12.26,,1.71,1.3,16.0,26.9,3.3;0333+321,Q,1.259,12.76,,1.95,0.4,22.0,14.9,2.3;0336-019,Q,0.852,22.36,,0.92,1.7,17.2,26.0,3.1;0403-132,Q,0.571,19.69,,1.54,,,,;0420-014,Q,0.914,7.36,,1.41,0.3,19.7,10.0,0.7;0454+844,B,0.112,0.14,,0.27,,,,;0528+134,Q,2.070,19.20,,1.60,0.4,30.9,17.4,1.2;0552+398,Q,2.363,0.36,<,0.36,,,,;0605-085,Q,0.872,19.79,,0.84,1.2,7.5,72.0,3.4;0607-157,Q,0.324,3.93,,0.68,,,,;0716+714,B,0.310,10.06,,0.49,1.6,10.8,21.5,4.6;0730+504,Q,0.720,14.06,,1.47,,,,;0736+017,Q,0.191,14.32,<,0.20,1.9,8.5,12.6,6.4;0738+313,Q,0.631,10.76,,0.81,,,,;0748+126,Q,0.889,18.37,,0.84,,,,;0754+100,B,0.266,14.40,,0.88,1.1,5.5,2.9,4.0;0804+499,Q,1.436,1.83,,0.48,0.1,35.2,17.6,0.1;0805-077,Q,1.837,50.60,,2.71,1.1,14.9,60.4,2.5;0823+033,B,0.506,17.80,,0.83,,,,;0827+243,Q,0.940,22.01,,1.18,1.6,13.0,21.6,4.0;0829+046,B,0.174,10.09,,0.34,,,,;0836+710,Q,2.218,25.38,,1.93,0.5,16.1,21.9,3.4;0851+202,B,0.306,15.17,<,0.28,1.3,16.8,9.7,2.3;0859-140,Q,1.339,16.47,,1.70,,,,;0906+015,Q,1.024,20.68,,1.61,,,,;0917+624,Q,1.446,15.57,,1.09,,,,;0923+392,Q,0.695,4.29,<,0.33,2.8,4.3,3.1,12.4;0945+408,Q,1.249,18.60,,1.10,1.4,6.3,35.6,5.2;0953+254,Q,0.712,11.52,<,0.42,,,,;1015+359,Q,1.226,12.46,,0.92,,,,;1036+054,Q,0.473,6.15,,0.74,,,,;1038+064,Q,1.265,11.87,,1.19,,,,;1045-188,Q,0.595,8.57,,0.86,,,,;1101+384,B,0.031,0.82,,0.10,,,,;1127-145,Q,1.184,14.18,,0.84,0.7,21.9,19.1,2.6;1150+812,Q,1.250,7.09,,0.68,,,,;1156+295,Q,0.729,24.73,,1.17,0.8,28.2,24.8,2.0;1219+044,Q,0.965,2.35,,0.81,,,,;1222+216,Q,0.432,21.10,,0.90,1.4,5.2,42.9,5.3;1253-055,Q,0.536,20.57,<,0.42,0.5,23.8,20.8,2.4;1302-102,Q,0.278,5.41,,0.70,,,,;1308+326,Q,0.997,27.17,,0.96,2.4,15.3,22.5,3.5;1334-127,Q,0.539,10.26,,1.23,,,,;1413+135,B,0.247,1.80,,0.44,0.2,12.1,6.2,0.9;1458+718,Q,0.904,7.04,,0.51,,,,;1502+106,Q,1.839,14.77,,0.69,3.0,11.9,20.0,4.4;1504-166,Q,0.876,4.31,,0.66,,,,;1508-055,Q,1.191,18.64,,1.52,,,,;1510-089,Q,0.360,20.14,,0.73,1.1,16.5,22.4,3.4;1532+016,Q,1.420,14.11,,1.14,,,,;1538+149,B,0.605,8.73,,0.49,4.1,4.3,11.2,10.5;1546+027,Q,0.414,12.08,<,0.32,,,,;1606+106,Q,1.226,18.91,,0.79,0.9,24.8,17.4,2.1;1611+343,Q,1.397,14.11,,0.66,1.6,13.6,42.4,3.1;1633+382,Q,1.814,29.45,,1.62,0.9,21.3,30.7,2.6;1637+574,Q,0.751,10.61,,0.71,1.1,13.9,13.6,4.1;1641+399,Q,0.593,19.27,,1.20,1.4,7.7,28.0,5.1;1642+690,Q,0.751,16.65,<,0.48,,,,;1652+398,B,0.033,0.21,,0.10,,,,;1655+077,Q,0.621,14.45,,0.47,,,,;1726+455,Q,0.717,1.82,<,0.28,,,,;1730-130,Q,0.902,35.69,,1.69,1.1,10.6,14.2,5.2;1749+096,B,0.322,6.84,,0.33,2.7,11.9,8.6,4.5;1749+701,B,0.770,6.03,,1.04,,,,;1803+784,B,0.680,8.97,<,0.39,0.3,12.1,6.2,0.9;1807+698,B,0.051,0.10,,0.12,1.4,1.1,1.0,7.3;1823+568,B,0.664,20.85,,0.74,0.7,6.3,41.1,4.8;1828+487,Q,0.692,13.65,,0.69,2.0,5.6,18.1,7.4;1849+670,Q,0.657,30.63,<,0.52,1.6,8.1,32.5,4.7;1901+319,Q,0.635,2.67,,0.87,,,,;1928+738,Q,0.302,8.43,,0.54,2.2,1.9,16.2,14.2;1936-155,Q,1.657,2.60,,1.31,,,,;2005+403,Q,1.736,12.21,,2.34,,,,;2037+511,Q,1.686,3.30,<,0.45,,,,;2113+293,Q,1.514,1.40,,1.16,,,,;2121+053,Q,1.941,13.29,,1.43,1.4,15.2,9.4,2.9;2128-123,Q,0.501,6.94,,0.98,,,,;2131-021,B,1.285,20.02,,1.02,,,,;2134+004,Q,1.932,5.94,,1.31,0.3,16.0,8.2,0.9;2136+141,Q,2.427,5.43,<,0.55,1.8,8.2,4.6,4.0;2145+067,Q,0.990,2.52,<,0.34,0.4,15.5,8.1,1.4;2155-152,Q,0.672,18.11,,1.89,,,,;2200+420,B,0.069,10.57,,0.09,3.5,7.2,7.3,8.0;2201+171,Q,1.076,2.55,,1.55,1.6,10.0,20.5,4.9;2201+315,Q,0.295,7.87,,0.95,2.0,6.6,8.2,8.5;2209+236,Q,1.125,3.43,<,0.39,,,,;2223-052,Q,1.404,17.34,,1.47,0.8,15.9,17.9,3.6;2227-088,Q,1.560,8.14,,1.44,0.5,15.8,8.2,1.4;2230+114,Q,1.037,15.41,,2.12,0.8,15.5,17.9,3.7;2243-123,Q,0.632,5.49,,0.78,,,,;2251+158,Q,0.859,14.19,,1.13,0.1,32.9,16.6,0.4;2345-167,Q,0.576,13.45,,0.90,,,,;2351+456,Q,1.986,18.01,,1.72,,,,"""
pc, mu0 = 3.0857e16, 4e-7*math.pi
rows = []
for item in RAW.split(';'):
    p = item.split(',')
    o = dict(src=p[0], cls=p[1], z=float(p[2]), bapp=float(p[3]), lim=p[4] == '<', B1=float(p[5]),
             aint=float(p[6]) if p[6] else None, delta=float(p[7]) if p[7] else None, Gam=float(p[8]) if p[8] else None, theta=float(p[9]) if p[9] else None)
    phi = 0.6; Gmin = math.sqrt(1 + o['bapp']**2)
    r = pc*math.tan(math.radians(phi)); Ip = 2*math.pi*r*o['B1']*1e-4/mu0
    o.update(I_rest=Ip, I_lab=Gmin*Ip, Gmin=Gmin)
    if o['aint'] and o['Gam'] and o['aint'] > 0:
        phi2 = o['aint']/2; B2 = o['B1']*(phi2/0.6)**(-0.25); r2 = pc*math.tan(math.radians(phi2))
        o['I_rest2'] = 2*math.pi*r2*B2*1e-4/mu0; o['I_lab2'] = o['Gam']*o['I_rest2']
    rows.append(o)
P = lambda I: 30.0*I*I*1e7                                   # erg/s
def st(a): a = np.array(a); return [float(np.percentile(a, 5)), float(np.median(a)), float(np.percentile(a, 95)), float(a.min()), float(a.max())]
det = [o for o in rows if not o['lim']]
Q = [o for o in det if o['cls'] == 'Q']; Bl = [o for o in det if o['cls'] == 'B']
m2 = [o for o in det if 'I_lab2' in o]
res = dict(n=len(rows), n_det=len(det), n_lim=len(rows) - len(det), nQ=len(Q), nB=len(Bl), n_second=len(m2),
           B1_Q=st([o['B1'] for o in Q]), B1_B=st([o['B1'] for o in Bl]),
           I_rest_all=st([o['I_rest'] for o in det]), I_lab_all=st([o['I_lab'] for o in det]), I_lab_Q=st([o['I_lab'] for o in Q]), I_lab_B=st([o['I_lab'] for o in Bl]),
           Gmin=st([o['Gmin'] for o in det]),
           I_lab2=st([o['I_lab2'] for o in m2]), I_rest2=st([o['I_rest2'] for o in m2]),
           P_lab_all=st([P(o['I_lab']) for o in det]), P_lab2=st([P(o['I_lab2']) for o in m2]),
           frac_lab_1e17_1e19=sum(1e17 <= o['I_lab'] <= 1e19 for o in det)/len(det),
           frac_lab_above_1e18=sum(o['I_lab'] >= 1e18 for o in det)/len(det),
           ref_3C303=3e18, ref_MW_halo=9e16, ref_peratt=[1e17, 1e19])
res['hist'] = {c: [round(math.log10(o['I_lab']), 3) for o in det if o['cls'] == c] for c in ('Q', 'B')}
res['hist_lim'] = [round(math.log10(o['I_lab']), 3) for o in rows if o['lim']]
res['hist_rest'] = [round(math.log10(o['I_rest']), 3) for o in det]
res['hist2'] = [round(math.log10(o['I_lab2']), 3) for o in m2]
json.dump(res, open(os.path.join(HERE, 'jet_currents.json'), 'w'), indent=1)
for k, v in res.items():
    if k not in ('hist', 'hist_lim', 'hist2', 'hist_rest'): print(k, v)
