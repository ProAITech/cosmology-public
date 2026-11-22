"""Workbench Entry 126: too big to form? How far matter travelled to build the cosmic web.

Lerner's argument (The Big Bang Never Happened, 1991): walls and voids spaced hundreds of millions of light-years apart would need matter to
travel ~270 million light-years at ~1000 km/s, taking several times the age of the universe. In the expanding model, matter does not
assemble structures by travelling across them: the pattern is present at the start as small density differences (seen in the CMB at 1e-5)
and grows in place; each parcel moves only far enough to change the local density. This script computes, from the Planck 2018 linear
spectrum of Entry 101:
 1. the rms displacement of matter since the beginning (Zel'dovich), 3D and per axis, and the share that moved more than given distances;
 2. how far matter must move to empty a void or build a wall of given size and contrast;
 3. the rms peculiar velocity today and the expected bulk flow in spheres of 50-200 Mpc/h, beside the Cosmicflows-4 measurements;
 4. Lerner's own timescale for comparison."""
import os, sys, json, math, contextlib, io
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
CODE = HERE if os.path.exists(os.path.join(HERE, 'progenitors.py')) else os.path.abspath(os.path.join(HERE, '..', 'artifact-files', '4379c2ac-1195-4a41-978e-eb44b945a265', 'code'))
sys.path.insert(0, CODE)
with contextlib.redirect_stdout(io.StringIO()): import progenitors as P
k, Pk, h, Om = P.kk, P.Pk, P.h, P.Om
f0 = Om**0.55; H0 = 100*h
def W_th(x): return 3*(np.sin(x) - x*np.cos(x))/x**3
def sig_psi(R=None):                     # 3D rms displacement (Mpc), optionally smoothed with a top hat of radius R (Mpc)
    W2 = 1.0 if R is None else W_th(k*R)**2
    return math.sqrt(np.trapezoid(Pk*W2, k)/(2*math.pi**2))
out = {}
s3 = sig_psi(); s1 = s3/math.sqrt(3)
out.update(sig3=s3, sig1=s1)
print(f"rms displacement since the beginning (linear, z = 0): {s3:.1f} Mpc in 3D, {s1:.1f} Mpc per axis")
# Maxwell distribution of |psi|: share moving more than d
from scipy.special import erfc
def frac_gt(d, s1=s1):
    x = d/s1; return float(erfc(x/math.sqrt(2)) + math.sqrt(2/math.pi)*x*math.exp(-x*x/2))
for d in (10, 20, 41, 83):
    out[f'frac_{d}'] = frac_gt(d)
    print(f"   share of matter that moved more than {d} Mpc: {frac_gt(d):.1e}")
# displacement smoothed on large scales: how much of it is coherent over R
for Rh in (8, 25, 50, 100, 150):
    R = Rh/h; sp = sig_psi(R); sd = P.sigma(R)
    out[f'psi_{Rh}'] = sp; out[f'delta_{Rh}'] = sd
    print(f"   top hat {Rh} Mpc/h ({R:.0f} Mpc): coherent displacement {sp:.1f} Mpc (3D rms), linear density contrast {sd:.3f}; building that contrast moves matter about {sd*R/3:.1f} Mpc")
# voids and walls
for R in (15, 30, 50):
    for dv in (-0.8, -0.9):
        dr = R*(1 - (1+dv)**(1/3)); out[f'void_{R}_{dv}'] = dr
        print(f"   void of radius {R} Mpc emptied to δ = {dv}: matter at its edge moved {dr:.1f} Mpc")
# peculiar velocities today and bulk flows
sv1 = f0*H0*s1/h*h                                      # km/s: f H0 x (displacement in Mpc) with H0 in km/s/Mpc
print(f"rms peculiar velocity per axis (linear): {f0*H0*s1:.0f} km/s; 3D {f0*H0*s3:.0f} km/s")
out['v1'] = f0*H0*s1; out['v3'] = f0*H0*s3
out['bulk'] = {}
for Rh in (50, 100, 150, 200):
    R = Rh/h; vb = f0*H0*sig_psi(R)
    # Gaussian window for comparison
    Wg2 = np.exp(-(k*R)**2); vg = f0*H0*math.sqrt(np.trapezoid(Pk*Wg2, k)/(2*math.pi**2))
    out['bulk'][Rh] = (vb, vg)
    print(f"   expected bulk flow (3D rms) in a sphere of {Rh} Mpc/h: {vb:.0f} km/s (top hat), {vg:.0f} km/s (Gaussian of that radius)")
# measured
meas = [('Watkins et al. 2023 (CF4, 150 Mpc/h)', None, None, 'probability below 0.03%'),
        ('Whitford, Howlett & Davis 2023 (CF4, minimum variance, 173 Mpc/h)', 428, 108, '0.11%'),
        ('Nusser & Tully 2026 (CF4 dipole, shell 120-160 Mpc/h)', 628, 82, '3.4σ in that shell; 2.3σ over all shells')]
out['meas'] = meas
# Lerner's timescale and the actual history
for d_mly, v in ((270, 1000), (75, 1000)):
    d = d_mly*0.3066; t = d*3.0857e19/v/3.15576e16
    print(f"Lerner: {d_mly} million light-years ({d:.0f} Mpc) at {v} km/s takes {t:.0f} Gyr")
    out[f'lerner_{d_mly}'] = t
# growth: linear growth since recombination
Dr = P.D(1100); print(f"linear growth since z = 1100: D(0)/D(1100) = {1/Dr:.0f}")
out['growth'] = 1/Dr
out['cum'] = [[d, frac_gt(d)] for d in range(0, 91)]
bc = []
for Rh in range(40, 211, 10):
    R = Rh/h; Wg2 = np.exp(-(k*R)**2)
    bc.append([Rh, f0*H0*sig_psi(R), f0*H0*math.sqrt(np.trapezoid(Pk*Wg2, k)/(2*math.pi**2))])
out['bulk_curve'] = bc
json.dump(out, open(os.path.join(HERE, 'flows.json'), 'w'), indent=1)
