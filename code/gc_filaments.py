"""gc_filaments.py - Workbench Entry 110: the Galactic-centre radio filaments as twisted flux tubes. What current, and how to see it.

Inputs (published): equipartition fields along the ~1000 MeerKAT filaments 100-400 uG (Yusef-Zadeh et al. 2022);
ALMA 100 GHz polarimetry of the Radio Arc, free of Faraday rotation: the field is turned from the filament by a constant
20-60 deg in one filament (PF1) and -20 to -80 deg in its neighbour (PF2) (Salem et al. 2026);
VLA rotation measures along the Arc -500 to -5800 rad/m^2, no sign reversal reported across a filament (Pare et al. 2019).
Reading tested here: the turn is the pitch psi of a twisted field, B_phi = B sin(psi), B_z = B cos(psi), carried by a uniform
axial current I = 2 pi r B_phi / mu0 in a filament of radius r (assumed 0.1-0.5 pc).
Prediction: a rotation-measure profile across the filament, RM(x) = 0.812 n_e B_phi (x/r) 2 sqrt(r^2 - x^2) (uG, cm^-3, pc),
odd in x, peak 0.812 n_e B_phi r at x = r/sqrt(2), opposite in sign for opposite handedness.
"""
import json, os, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
MU0 = 4e-7*math.pi; PC = 3.0857e16
def current(B_uG, psi_deg, r_pc): return 2*math.pi*r_pc*PC*B_uG*1e-10*math.sin(math.radians(psi_deg))/MU0
grid = [current(B, p, r) for B in (100, 200, 400) for p in (20, 40, 60) for r in (0.1, 0.25, 0.5)]
base = dict(B=200, psi=40, r=0.25)
out = {'I_range': [float('%.2g' % min(grid)), float('%.2g' % max(grid))], 'I_base': float('%.2g' % current(**{'B_uG': 200, 'psi_deg': 40, 'r_pc': 0.25}))}
# forces: magnetic pressure vs thermal; Bennett pinch per unit length vs thermal content
B = 200e-10; pB = B**2/(2*MU0); n = 10e6; T = 1e4; pth = 2*n*1.380649e-23*T
r = 0.25*PC; I = out['I_base']; pinch = MU0*I**2/(8*math.pi); NkT = math.pi*r**2*pth
out['pressure'] = dict(magnetic_Pa=float('%.2g' % pB), thermal_Pa_n10_T1e4=float('%.2g' % pth), ratio=float('%.2g' % (pB/pth)),
                       bennett_J_per_m=float('%.2g' % pinch), thermal_J_per_m=float('%.2g' % NkT), pinch_over_thermal=float('%.2g' % (pinch/NkT)))
# rotation-measure step across the filament
def rm_peak(ne, B_uG, psi, r_pc): return 0.812*ne*B_uG*math.sin(math.radians(psi))*r_pc
out['rm_peak'] = dict(low=round(rm_peak(1, 100, 20, 0.1), 1), base_ne1=round(rm_peak(1, 200, 40, 0.25), 1), base_ne3=round(rm_peak(3, 200, 40, 0.25), 1),
                      base_ne10=round(rm_peak(10, 200, 40, 0.25), 1), high=round(rm_peak(10, 400, 60, 0.5), 1))
x = np.linspace(-1.6, 1.6, 161)
prof = {}
for ne in (1, 3, 10):
    Bphi = 200*math.sin(math.radians(40))
    rm = np.where(np.abs(x) < 1, 0.812*ne*Bphi*x*2*np.sqrt(np.clip(1 - x**2, 0, None))*0.25, 0.0)
    prof[str(ne)] = [round(float(v), 1) for v in rm]
out['profile'] = dict(x=[round(float(v), 2) for v in x], rm=prof, B=200, psi=40, r=0.25)
# electron drift needed, if thermal electrons carry it
J = I/(math.pi*r**2)
out['drift'] = dict(J=float('%.2g' % J), v_ne1_m_s=float('%.2g' % (J/(1e6*1.602e-19))))
json.dump(out, open(os.path.join(HERE, 'gc_filaments.json'), 'w'), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != 'profile'}, indent=1))
