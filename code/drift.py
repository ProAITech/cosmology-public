"""Workbench Entry 120: watching redshifts change - the direct test of expansion.

In an expanding universe the redshift of a distant source changes slowly as we watch (Sandage 1962, Loeb 1998):
    dz/dt0 = H0 (1 + z) - H(z),   or, as a velocity,   dv/dt0 = c H0 [1 - E(z)/(1 + z)].
In any static universe, tired light, plasma redshift or Lerner's, a source at fixed distance keeps its redshift: the drift is zero.
So is it in a coasting universe (H = H0 (1 + z)). This script:
 1. computes the signal for Planck LCDM and DESI DR2's evolving dark energy, and where it changes sign;
 2. puts the published measurements on the same scale;
 3. adds our own Galaxy's pull: the Sun accelerates toward the Galactic centre (Gaia EDR3), which shifts every redshift by up to 0.73 cm/s/yr,
    the size of the cosmological signal, as a dipole on the sky; projects it on the ESPRESSO target;
 4. scales the ESPRESSO precision to ask when a static universe and LCDM can be told apart at 5 sigma, with the VLT and with the ELT;
 5. estimates what a plasma-redshift static universe adds: drifts from sources moving through the redshifting plasma.
"""
import json, math, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
C = 299792.458; YR = 3.15576e7; MPC = 3.0857e19
def H0s(H0): return H0/MPC*YR          # km/s/Mpc -> per year
def E_lcdm(z, Om=0.315): return np.sqrt(Om*(1+z)**3 + 1 - Om)
def E_cpl(z, Om, w0, wa):
    de = (1+z)**(3*(1+w0+wa))*np.exp(-3*wa*z/(1+z)); return np.sqrt(Om*(1+z)**3 + (1-Om)*de)
def vdot(z, E, H0): return C*1e5*H0s(H0)*(1 - E/(1+z))        # cm/s per year
z = np.linspace(0, 5, 501)
models = {
  'Planck LCDM': (lambda zz: E_lcdm(zz, 0.315), 67.4),
  'DESI DR2 w0waCDM (CMB + DES 5-yr SNe)': (lambda zz: E_cpl(zz, 0.319, -0.752, -0.86), 66.7),
  'DESI DR2 w0waCDM (CMB + Pantheon+)': (lambda zz: E_cpl(zz, 0.311, -0.838, -0.62), 67.5),
  'Einstein-de Sitter (matter only)': (lambda zz: E_lcdm(zz, 0.9999), 67.4),
}
out = {'z': z[::10].round(2).tolist(), 'curves': {}}
for k, (E, H0) in models.items():
    v = vdot(z, E(z), H0); out['curves'][k] = v[::10].round(4).tolist()
    m = z > 0.1; zc = z[m][1:][np.diff(np.sign(v[m])) != 0]
    print(f"{k:42s} v-dot at z = 0.5, 1, 3.962, 5: " + ", ".join(f"{vdot(q, E(q), H0):+.3f}" for q in (0.5, 1, 3.962, 5)) + f" cm/s/yr; peak {v.max():+.3f} at z = {z[v.argmax()]:.2f}; zero at z = {zc[0] if len(zc) else float('nan'):.2f}")
    out[k] = dict(v05=vdot(0.5, E(0.5), H0), v1=vdot(1, E(1), H0), v396=vdot(3.962, E(3.962), H0), v5=vdot(5, E(5), H0), zero=float(zc[0]) if len(zc) else None, peak=float(v.max()), zpeak=float(z[v.argmax()]))
print("static universe (any tired light, plasma redshift, Lerner), coasting universe: 0 at every redshift")
# 2. measurements, in cm/s/yr
L = out['Planck LCDM']
meas = [('Darling 2012, HI 21 cm, 10 absorbers, z 0.09-0.69, 13.5 yr', 0.4, -550, 220),
        ('Trost et al. 2025, ESPRESSO, z 3.962, 0.9 yr', 3.962, -143, 509),
        ('Trost et al. 2026, ESPRESSO, z 3.962, 2 yr', 3.962, -350, 360)]
out['meas'] = []
for n, zz, v, s in meas:
    sig = abs(vdot(zz, E_lcdm(zz), 67.4))
    print(f"   {n}: {v:+.0f} ± {s:.0f} cm/s/yr; LCDM {vdot(zz, E_lcdm(zz), 67.4):+.3f}; error is {s/sig:,.0f} times the signal")
    out['meas'].append([n, zz, v, s, s/sig])
print("   Cooke 2020, Lya cells: below 6,500 cm/s/yr (2 sigma)")
# 3. the Sun's acceleration (Gaia EDR3): (2.32 ± 0.16) x 1e-10 m/s^2 toward RA 269.1, Dec -31.6
a = 2.32e-10*YR*100; ea = 0.16e-10*YR*100                     # cm/s per year
def unit(ra, de): ra, de = math.radians(ra), math.radians(de); return np.array([math.cos(de)*math.cos(ra), math.cos(de)*math.sin(ra), math.sin(de)])
na = unit(269.1, -31.6); ns = unit(82.316, -43.864)          # J052915.80-435152.0
cth = float(na @ ns); proj = -a*cth                           # observer moving toward the source lowers its redshift
out.update(a_obs=a, a_err=ea, cos_target=cth, proj_target=proj)
print(f"   Sun's acceleration: {a:.3f} ± {ea:.3f} cm/s/yr; toward J0529-4351 the angle is {math.degrees(math.acos(cth)):.0f}°, so it adds {proj:+.3f} cm/s/yr (LCDM signal there {L['v396']:+.3f})")
# 4. when can static and LCDM be told apart? ESPRESSO anchor: sigma(dv) = 4.45 m/s from 12 h (Trost et al. 2025)
s12 = 445.0                                                   # cm/s for 12 h, VLT
kELT = math.sqrt(978/49.3)                                    # collecting areas, ELT vs one VLT unit telescope
out['kELT'] = kELT
def sig_vdot(T, hours_per_yr, k, uniform=False):
    H = hours_per_yr*T; s = s12*math.sqrt(12/H)/k/T          # two epochs at the ends of T years, hours split equally
    return s*(math.sqrt(12)/2 if uniform else 1)
target = abs(L['v396'])
out['clock'] = {}
for lab, k in (('VLT/ESPRESSO', 1.0), ('ELT/ANDES', kELT)):
    for hpy in (50, 100):
        for uni in (False, True):
            T = next(t for t in np.arange(1, 400, 0.5) if target/sig_vdot(t, hpy, k, uni) >= 5)
            out['clock'][f"{lab}|{hpy}|{uni}"] = float(T)
            print(f"   {lab}, {hpy} h a year, {'spread evenly' if uni else 'at the two ends'}: 5 sigma between static and LCDM at z = 3.96 after {T:.0f} years")
# check against the paper's own estimate: ELT, 5,400 h over 54 years
for uni in (False, True):
    print(f"   ELT, 5,400 h over 54 yr ({'even' if uni else 'ends'}): {target/sig_vdot(54, 100, kELT, uni):.1f} sigma (paper: 99% detection)")
out['curve_sig'] = {lab: [[t, sig_vdot(t, 100, k, True)] for t in (2, 5, 10, 20, 30, 40, 50, 60, 80, 100, 150, 200)] for lab, k in (('VLT/ESPRESSO', 1.0), ('ELT/ANDES', kELT))}
# differential: a low-z absorber against a high-z forest on the same sight line doubles the contrast and cancels the Sun's pull
dd = L['v05'] - L['v396']; out['diff'] = dd
print(f"   difference z = 0.5 minus z = 3.96: {dd:.3f} cm/s/yr in LCDM, 0 in a static universe; the Sun's pull cancels")
# 5. plasma redshift: source moving through the plasma at v_r changes the column, so z drifts
ne = (67.4/MPC)/(C*1e5)/3.326e-25; out['ne_plasma'] = ne     # electrons per cm^3 that give H0 = 67.4
for vr in (300, 1000):
    dzdt = (vr/C)*(67.4/MPC)*YR                               # d ln(1+z)/dt = (H0/c) v_r
    print(f"   plasma redshift (n_e = {ne:.1e} cm^-3): a source moving at {vr} km/s through it drifts by ±{C*1e5*dzdt:.4f} cm/s/yr, random in sign")
    out[f'plasma_{vr}'] = C*1e5*dzdt
# CMB cooling
out['Tdot'] = 2.7255*H0s(67.4)
print(f"   the CMB cools by {out['Tdot']:.2e} K per year in an expanding universe; FIRAS's absolute scale is good to about 1 mK")
json.dump(out, open(os.path.join(HERE, 'drift.json'), 'w'), indent=1)
