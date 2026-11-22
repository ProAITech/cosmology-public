"""Workbench Entry 123: are the Fermi and eROSITA bubbles a Galactic pinch?

Above and below the Galactic centre sit two gamma-ray bubbles about 8-10 kpc tall (Fermi), inside two X-ray bubbles about 14 kpc tall
(eROSITA), with polarized radio lobes whose ridges wind around them (S-PASS). Their hourglass shape, centred on the Galaxy's axis, is read
in plasma cosmology as a z-pinch: a current along the axis whose own magnetic field squeezes the plasma. This script asks four questions:
 1. What current would hold the bubbles' plasma in a Bennett pinch, and what field would it make? (Measured: up to 15 microgauss.)
 2. What current do the measured fields allow, and does the winding of the field imply one?
 3. Is the edge a confining boundary or a shock? Rankine-Hugoniot against the published temperature and speed; the UV-absorption kinematics.
 4. Power: the bubbles' energy over their age, against Sgr A*, nuclear star formation and the Galaxy's own circuit (Entry 107).
"""
import json, math, os
HERE = os.path.dirname(os.path.abspath(__file__))
MU0 = 4e-7*math.pi; KB = 1.380649e-23; MP = 1.67262e-27; KPC = 3.0857e19; MYR = 3.15576e13
out = {}
# 1. Bennett pinch: mu0 I^2/(8 pi) = N_l k (Te + Ti), N_l = pi a^2 n_e
n_e = 1e-3*1e6; T = 4.5e6                                   # Miller & Bregman 2016: n ~ 1e-3 cm^-3, log T 6.60-6.70
for a_kpc in (2.0, 3.0, 4.0):
    a = a_kpc*KPC; Nl = math.pi*a*a*n_e
    I = math.sqrt(8*math.pi*Nl*KB*2*T/MU0); B = MU0*I/(2*math.pi*a)
    out[f'bennett_{a_kpc}'] = (I, B*1e10)
    print(f"Bennett pinch, radius {a_kpc} kpc, n_e 1e-3 cm^-3, T 4.5e6 K: I = {I:.2e} A, edge field {B*1e10:.1f} microgauss")
# plasma beta with the measured field
for Bmu in (6, 15):
    pB = (Bmu*1e-10)**2/(2*MU0); pth = 2*n_e*KB*T
    print(f"   {Bmu} microgauss: magnetic pressure {pB:.2e} Pa against thermal {pth:.2e} Pa (beta {pth/pB:.2f})")
    out[f'beta_{Bmu}'] = pth/pB
# 2. current allowed by the field, and current implied by a wound (toroidal) field
for a_kpc in (2.0, 3.0, 4.0):
    Imax = 2*math.pi*a_kpc*KPC*15e-10/MU0
    out[f'Imax_{a_kpc}'] = Imax
    print(f"   at {a_kpc} kpc a 15 microgauss toroidal field carries at most {Imax:.1e} A")
for Bphi in (2, 5):
    I = 2*math.pi*3*KPC*Bphi*1e-10/MU0; out[f'Iwound_{Bphi}'] = I
    print(f"   a wound component of {Bphi} microgauss around 3 kpc means an axial current of {I:.1e} A")
B_peratt = MU0*1e19/(2*math.pi*3*KPC)*1e10
out['B_peratt'] = B_peratt
print(f"   a Peratt-scale current of 1e19 A through a 3 kpc bubble makes {B_peratt:.0f} microgauss, {B_peratt/15:.0f} times the measured maximum")
# 3. shock: Rankine-Hugoniot, gamma 5/3, halo at 2e6 K (Miller & Bregman 2016)
g = 5/3; mu = 0.6; T1 = 2e6
cs = math.sqrt(g*KB*T1/(mu*MP))/1e3
def T_ratio(M): return (2*g*M*M - (g-1))*((g-1)*M*M + 2)/((g+1)**2*M*M)
def comp(M): return (g+1)*M*M/((g-1)*M*M + 2)
out['cs'] = cs; out['shock'] = []
for v in (413, 490, 720):
    M = v/cs; out['shock'].append((v, M, T_ratio(M)*T1, comp(M)))
    print(f"   shell speed {v} km/s: Mach {M:.2f} in a 2e6 K halo (sound speed {cs:.0f} km/s); shocked gas at {T_ratio(M)*T1:.2e} K, compressed {comp(M):.2f} times")
# inverse: Mach from the measured temperature 4.0-5.0e6 K
def M_of(Tr):
    lo, hi = 1.0, 10.0
    for _ in range(60):
        m = (lo+hi)/2; lo, hi = (m, hi) if T_ratio(m) < Tr else (lo, m)
    return (lo+hi)/2
out['M_from_T'] = [M_of(4.0e6/T1), M_of(5.0e6/T1)]
print(f"   measured 4.0-5.0e6 K implies Mach {out['M_from_T'][0]:.2f}-{out['M_from_T'][1]:.2f}, a speed of {out['M_from_T'][0]*cs:.0f}-{out['M_from_T'][1]*cs:.0f} km/s")
for Bh in (3, 6):
    print(f"   a halo field of {Bh} microgauss, compressed {comp(490/cs):.2f} times at the edge: {Bh*comp(490/cs):.0f} microgauss")
# 4. power
E_er = 1e56*1e-7                                            # J (eROSITA, Predehl et al. 2020)
need = {}
for age in (4.3, 6.0, 9.0):
    need[f'eROSITA 1e56 erg over {age} Myr'] = E_er/(age*MYR)
need['Miller & Bregman injection rate'] = 2.3e42*1e-7
for k, v in need.items(): print(f"   needed: {k}: {v:.1e} W")
LEdd = 1.26e31*4.3e6                                        # W, for Sgr A* (4.3e6 Msun)
src = {'Sgr A* at 1% of Eddington': 0.01*LEdd, 'nuclear star formation, 0.1 Msun/yr': 0.1/100*1e44/3.15576e7,
       'the Galaxy’s circuit (Entry 107)': 5.3e32}
for k, v in src.items(): print(f"   source: {k}: {v:.1e} W")
lo, hi = min(need.values()), max(need.values())
print(f"   the Galaxy's circuit falls {lo/5.3e32:.0f}-{hi/5.3e32:.0f} times short; the voltage it would need at the Bennett current: {lo/out['bennett_3.0'][0]:.1e}-{hi/out['bennett_3.0'][0]:.1e} V, against its 2.7e15 V")
print(f"   Sgr A* needs {lo/LEdd*100:.2f}-{hi/LEdd*100:.2f}% of its Eddington power; nuclear star formation at 0.1 Msun/yr gives {src['nuclear star formation, 0.1 Msun/yr']/lo*100:.0f}-{src['nuclear star formation, 0.1 Msun/yr']/hi*100:.0f}% of it")
out.update(need=need, src=src, LEdd=LEdd, V_need=(lo/out['bennett_3.0'][0], hi/out['bennett_3.0'][0]))
json.dump(out, open(os.path.join(HERE, 'bubbles.json'), 'w'), indent=1)
