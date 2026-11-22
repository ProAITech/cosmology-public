"""Workbench Entry 124: can laboratory pinches be scaled up to galaxies? The dimensionless numbers.

Plasma cosmology rests on scaling: what a z-pinch does in the laboratory, filaments of current do at the size of galaxies (Alfven; Peratt 1986).
Ideal magnetohydrodynamics has no built-in length, so a laboratory plasma and a galaxy behave alike if both are ideal and share their
dimensionless numbers (Ryutov et al. 2001): dissipation numbers (Lundquist S, magnetic Reynolds Rm) much larger than 1 in both; the
collision mean free path short or long in both alike; the same plasma beta and Mach numbers; and the system much larger than the ion
inertial length and gyroradius in both, so that two-fluid and kinetic effects do not enter. This script computes those numbers for
laboratory pinches, a tokamak, two space plasmas where Birkeland currents are measured, and three galactic and intergalactic plasmas.
Gravity is added as one more number: magnetic pressure against the pressure needed to hold the matter in orbit."""
import json, math, os
HERE = os.path.dirname(os.path.abspath(__file__))
E = 1.602176634e-19; ME = 9.1093837e-31; MP = 1.67262e-27; EPS0 = 8.8541878e-12; MU0 = 4e-7*math.pi; C = 2.99792458e8; KPC = 3.0857e19
# name: n_e (m^-3), T (eV), B (T), L (m), flow speed v (m/s), ion mass number A, orbital speed for gravity (m/s) and total density (kg/m^3) or None
SYS = [
 ('Dense plasma focus pinch',            1e25, 1000, 300,    1e-3,   1e5, 2, None),
 ('Wire-array z-pinch (Z machine)',      1e27, 1000, 3000,   1e-3,   3e5, 1, None),
 ('Tokamak core',                        1e20, 1e4,  5,      2,      1e3, 2, None),
 ('Auroral current region',              1e6,  1000, 1e-5,   1e7,    1e5, 1, None),
 ('Solar coronal loop',                  1e15, 86,   1e-2,   1e7,    1e5, 1, None),
 ('Warm ionized gas, Milky Way disc',    1e5,  0.7,  5e-10,  1*KPC,  1e4, 1, (2.2e5, 0.1*1.989e30/(3.0857e16)**3)),
 ('Hot halo, Milky Way',                 1e3,  172,  3e-10,  10*KPC, 2e5, 1, (2.2e5, 0.01*1.989e30/(3.0857e16)**3)),
 ('Intergalactic filament',              1e2,  86,   3e-12,  1000*KPC, 3e5, 1, (3e5, 10*0.27*8.5e-27)),
]
def lnL(n, T):
    nc = n/1e6
    return 24 - math.log(math.sqrt(nc)/T) if T > 10 else 23 - math.log(math.sqrt(nc)*T**-1.5)
rows = []
for name, n, T, B, L, v, A, grav in SYS:
    mi = A*MP; ll = lnL(n, T)
    lD = math.sqrt(EPS0*T*E/(n*E*E)); de = C/math.sqrt(n*E*E/(EPS0*ME)); di = de*math.sqrt(mi/ME)
    rho_i = math.sqrt(T*E/mi)/(E*B/mi)
    mfp = 1.44e17*T*T/(n*ll)
    eta = 5.2e-5*ll/T**1.5; D = eta/MU0
    vA = B/math.sqrt(MU0*n*mi)
    S = L*vA/D; Rm = L*v/D
    beta = 2*n*T*E/(B*B/(2*MU0))
    Mg = None
    if grav:
        vc, rho_tot = grav; Mg = (B*B/(2*MU0))/(rho_tot*vc*vc)
    r = dict(name=name, n=n, T=T, B=B, L=L, lnL=ll, L_lD=L/lD, L_de=L/de, L_di=L/di, L_rhoi=L/rho_i, Kn=mfp/L, S=S, Rm=Rm, beta=beta, vA_c=vA/C, Mg=Mg)
    rows.append(r)
    print(f"{name:36s} L/d_i {L/di:9.1e}  L/rho_i {L/rho_i:9.1e}  mfp/L {mfp/L:9.1e}  S {S:9.1e}  Rm {Rm:9.1e}  beta {beta:8.1e}" + (f"  magnetic/orbital {Mg:8.1e}" if Mg else ""))
# what a particle simulation can span: a 3D grid of N^3 cells resolving the electron skin depth with ~2 cells
for N in (512, 4096):
    print(f"   a {N}^3 particle simulation spans about {N/2:.0f} electron skin depths; with the true mass ratio that is {N/2/math.sqrt(1836):.0f} ion inertial lengths")
lab = rows[0]; gal = rows[6]
print(f"   gap, dense plasma focus to Milky Way halo: L/d_i x{gal['L_di']/lab['L_di']:.0e}, S x{gal['S']/lab['S']:.0e}")
json.dump(rows, open(os.path.join(HERE, 'scaling.json'), 'w'), indent=1)
