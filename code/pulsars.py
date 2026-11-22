"""Workbench Entry 125: are pulsars rotating neutron stars, or electric discharges on ordinary remnants?

Healy and Peratt (1995) modelled pulsar radio emission as plasma processes in a magnetosphere. Electric Universe writers extend this:
'a normal stellar remnant undergoing periodic discharges will suffice', with no neutron star. Whatever makes the pulses, each pulsar is a
clock carried by a body, and the clock's timing weighs and sizes that body. Tests, none of which needs to know how the pulses are made:
 1. Mass: Shapiro delay gives PSR J0740+6620 2.08 +- 0.07 Msun; a white dwarf cannot exceed the Chandrasekhar mass, about 1.4 Msun.
 2. Orbit size: the double pulsar's two bodies (1.34 and 1.25 Msun) orbit every 2.45 hours, 1.26 solar radii apart.
 3. GW170817: two bodies of 1.17-1.60 Msun entered LIGO's band at 24 Hz, 400 km apart, and spiralled for about 3,000 cycles.
 4. Rotation: if a 716 Hz pulsar's period is its spin, the body must be denser than 3 pi/(G P^2); a white dwarf cannot spin faster than seconds.
 5. Rotation shown directly: the double pulsar's eclipse is modulated at half of pulsar B's 2.8 s period, as B's magnetosphere turns.
 6. Radii: NICER's pulse-profile fits give 12-13 km.
"""
import json, math, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
G = 6.674e-11; MSUN = 1.98847e30; RSUN = 6.957e8; C = 2.99792458e8
out = {}
# 4. spin limits
for P in (1.396e-3, 2.947e-3, 22.7e-3):
    rho = 3*math.pi/(G*P*P)
    Rmax = (G*1.4*MSUN*P*P/(4*math.pi**2))**(1/3)
    out[f'spin_{P}'] = (rho, Rmax)
    print(f"spin period {P*1e3:.3f} ms: mean density above {rho/1e3:.1e} g/cm^3; a 1.4 Msun body must be smaller than {Rmax/1e3:.0f} km (Newtonian mass shedding)")
def R_wd(M, Mch=1.44):                                      # Nauenberg 1972
    x = M/Mch; return 0.0112*RSUN*math.sqrt(max(x**(-2/3) - x**(2/3), 0))
for M in (0.6, 1.0, 1.25, 1.34, 1.4):
    R = R_wd(M); rho = M*MSUN/(4/3*math.pi*R**3); Pmin = math.sqrt(3*math.pi/(G*rho))
    print(f"   white dwarf of {M} Msun: radius {R/1e3:.0f} km, mean density {rho/1e3:.1e} g/cm^3, fastest spin {Pmin:.2f} s")
    out[f'wd_{M}'] = (R, rho, Pmin)
# 2. double pulsar
Mt = (1.338185 + 1.248868)*MSUN; Pb = 2.4535*3600
a = (G*Mt*Pb**2/(4*math.pi**2))**(1/3); out['dp_a'] = a
print(f"double pulsar: separation {a/1e3:.2e} km = {a/RSUN:.2f} solar radii; two 1.3 Msun main-sequence stars would each be about {1.3**0.8:.2f} solar radii")
# 3. GW170817 separation at 24 Hz (gravitational-wave frequency = twice the orbital)
Mgw = 2.74*MSUN
for f in (24, 100, 400):
    r = (G*Mgw/(math.pi*f)**2)**(1/3); out[f'gw_r_{f}'] = r
    print(f"   GW170817 at {f} Hz: separation {r/1e3:.0f} km")
for R in (R_wd(1.37), 12e3):
    fc = (1/math.pi)*math.sqrt(G*Mgw/(2*R)**3); out[f'fcontact_{R:.0f}'] = fc
    print(f"   two {R/1e3:.0f} km bodies of 2.74 Msun in total touch at {fc:.2f} Hz")
# chirp cycles from 24 Hz to contact (leading order)
Mc = 1.188*MSUN; tc = G*Mc/C**3
K = (1/(32*math.pi**(8/3)))*tc**(-5/3)
N = K*(24**(-5/3) - 1600**(-5/3)); out['gw_cycles'] = N
print(f"   leading-order cycles from 24 Hz to contact for point masses: {N:.0f} (observed: about 3,000 in the analysed band)")
# tidal deformability -> radius: Lambda = (2/3) k2 (R c^2/G M)^5
for k2 in (0.08, 0.12):
    R = (800*1.5/k2)**(1/5)*G*1.4*MSUN/C**2; out[f'R_lambda_{k2}'] = R
    print(f"   Lambda(1.4 Msun) <= 800 with k2 = {k2}: R <= {R/1e3:.1f} km")
lam_wd = (2/3)*0.1*(R_wd(1.37)*C**2/(G*1.37*MSUN))**5; out['lam_wd'] = lam_wd
print(f"   a 1.37 Msun white dwarf would have Lambda about {lam_wd:.0e}")
# 1. J0740 mass against Chandrasekhar
print(f"PSR J0740+6620: 2.08 ± 0.07 Msun, {(2.08-1.44)/0.07:.1f} sigma above 1.44 Msun")
out['j0740_sig'] = (2.08-1.44)/0.07
# figure curves: white-dwarf branch, main sequence, spin limit
Ms = np.linspace(0.2, 1.43, 60)
out['wd_curve'] = [[float(m), float(R_wd(m)/1e3)] for m in Ms]
out['ms_curve'] = [[float(m), float(m**0.8*RSUN/1e3)] for m in np.linspace(0.3, 2.5, 30)]
out['spin_curve'] = [[float(m), float((G*m*MSUN*1.396e-3**2/(4*math.pi**2))**(1/3)/1e3)] for m in np.linspace(0.2, 2.5, 30)]
json.dump(out, open(os.path.join(HERE, 'pulsars.json'), 'w'), indent=1)
