"""halo_loop.py - Workbench Entry 107: where the Milky Way's halo current flows, and where its energy goes.

The halo toroids of Xu & Han (2024), fitted to pulsar and extragalactic rotation measures:
  B_phi(R, z) = sign(z) B0 (|z|/z0) exp(-(|z| - z0)/z0) exp(-((R - R0)/R_T)^2),
  B0 = 0.73 +- 0.12 uG, z0 = 3.0 kpc (fit capped; -0.65), R0 = 7.97 (+1.49, -0.93) kpc, R_T = 5.31 (+3.91, -1.75) kpc.
Ampere's law with an azimuthal field only: the poloidal current runs along contours of R B_phi (its stream function),
  J_R = -(1/mu0) dB_phi/dz,  J_z = (1/(mu0 R)) d(R B_phi)/dR,
  current through a disc of radius R at height z:  I_z = 2 pi R B_phi / mu0.
Sense (Entry 18): inward through the disc region, away from the plane at small R, in each hemisphere.
Energy: the field's energy B^2/2mu0 over the halo, against the e.m.f. of the disc turning through its own 0.3 uG vertical field
(2.7e15 V, Entry 18), and against Ohmic loss in the halo's warm ionized gas.
"""
import json, os, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
MU0 = 4e-7*math.pi; KPC = 3.0857e19; UG = 1e-10
P0 = dict(B0=0.73, z0=3.0, R0=7.97, RT=5.31)
def Bphi(R, z, p=P0):
    az = np.abs(z)
    return np.sign(z)*p['B0']*(az/p['z0'])*np.exp(-(az - p['z0'])/p['z0'])*np.exp(-((R - p['R0'])/p['RT'])**2)   # uG
def Rpeak(p): return (p['R0'] + math.sqrt(p['R0']**2 + 2*p['RT']**2))/2
def Iz(R, z, p=P0): return 2*math.pi*R*KPC*abs(Bphi(np.array(R), np.array(z), p))*UG/MU0        # A, through a disc of radius R at height z
def summary(p):
    rp = Rpeak(p); Iv = float(Iz(rp, p['z0'], p)); Ir = 2*Iv
    cap1 = float(Iz(1.0, p['z0'], p)); cap2 = float(Iz(2.0, p['z0'], p))
    return dict(R_peak=round(rp, 2), I_loop=float('%.3g' % Iv), I_radial_disc=float('%.3g' % Ir), cap_1kpc=float('%.3g' % cap1), cap_2kpc=float('%.3g' % cap2),
                cap_1kpc_frac=round(cap1/Iv, 3), cap_2kpc_frac=round(cap2/Iv, 3))
out = {'model': P0, 'base': summary(P0)}
# parameter ranges within the published errors
var = []
for B0 in (0.61, 0.73, 0.85):
    for z0 in (2.35, 3.0):
        for R0 in (7.04, 7.97, 9.46):
            for RT in (3.56, 5.31, 9.22):
                var.append(summary(dict(B0=B0, z0=z0, R0=R0, RT=RT)))
out['ranges'] = {k: [min(v[k] for v in var), max(v[k] for v in var)] for k in ('R_peak', 'I_loop', 'cap_1kpc_frac', 'cap_2kpc_frac')}
# grid, current density, drift speed
R = np.linspace(0.05, 25, 500); Z = np.linspace(-12, 12, 481); RR, ZZ = np.meshgrid(R, Z)
B = Bphi(RR, ZZ)*UG
dBdz = np.gradient(B, Z*KPC, axis=0); dRBdR = np.gradient(RR*KPC*B, R*KPC, axis=1)
JR = -dBdz/MU0; JZ = dRBdR/(MU0*RR*KPC); J = np.hypot(JR, JZ)
ne = 0.016*np.exp(-np.abs(ZZ)/1.7)*1e6        # warm ionized gas, m^-3 (scale height 1.7 kpc as used with the fit)
vd = J/(ne*1.602e-19)
m = (RR > 2) & (RR < 20) & (np.abs(ZZ) > 0.2) & (np.abs(ZZ) < 10)
out['J'] = dict(max=float('%.3g' % J[m].max()), typical=float('%.3g' % np.median(J[m])), drift_max_m_s=float('%.3g' % vd[m].max()), drift_typ_m_s=float('%.3g' % np.median(vd[m])))
# energy and power
dV = 2*np.pi*RR*KPC*(R[1] - R[0])*KPC*(Z[1] - Z[0])*KPC
W = float(np.sum(B**2/(2*MU0)*dV))
V_emf = 2.7e15; P = out['base']['I_radial_disc']*V_emf
eta = 5.2e-5*20/(1e4*8.617e-5)**1.5            # Spitzer resistivity at 1e4 K, ohm m
ohm = float(np.sum(J**2*eta*dV))
out['energy'] = dict(W_J=float('%.3g' % W), W_erg=float('%.3g' % (W*1e7)), emf_V=V_emf, power_W=float('%.3g' % P), build_Myr=round(W/P/3.156e13, 1),
                     resistivity_ohm_m=float('%.3g' % eta), ohmic_W=float('%.3g' % ohm), ohmic_over_power=float('%.3g' % (ohm/P)),
                     ohmic_decay_Gyr=float('%.3g' % (W/ohm/3.156e16)))
# current lines: contours of R*B_phi in each hemisphere, with the current between neighbouring lines
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
psi = RR*Bphi(RR, ZZ)                                   # kpc uG
pmax = float(np.abs(psi).max()); levels = [f*pmax for f in (0.2, 0.35, 0.5, 0.65, 0.8, 0.92)]
cs = plt.contour(RR, ZZ, psi, levels=sorted([-l for l in levels] + levels))
lines = []
for lev, segs in zip(cs.levels, cs.allsegs):
    for sgm in segs:
        if len(sgm) < 5: continue
        step = max(1, len(sgm)//80); pts = sgm[::step]
        lines.append(dict(level=round(float(lev)/pmax, 2), pts=[[round(float(a), 2), round(float(b), 2)] for a, b in pts]))
out['lines'] = lines
out['psi_max_kpc_uG'] = round(pmax, 3); out['I_per_line'] = float('%.3g' % (0.15*out['base']['I_loop']))   # between lines 0.2-0.35-0.5-0.65-0.8
json.dump(out, open(os.path.join(HERE, 'halo_loop.json'), 'w'), indent=1)
for k in ('base', 'ranges', 'J', 'energy'): print(k, json.dumps(out[k]))
print(len(lines), 'lines')

# orient each current line in the sense of Entry 18: inward (toward smaller R) where the line comes closest to the plane
for ln_ in out['lines']:
    p = ln_['pts']; i = min(range(len(p)), key=lambda k: abs(p[k][1]))
    j = (i + 1) % len(p)
    if p[j][0] > p[i][0]: ln_['pts'] = p[::-1]
json.dump(out, open(os.path.join(HERE, 'halo_loop.json'), 'w'), indent=1)
print('oriented', [round(l['level'], 2) for l in out['lines']])
