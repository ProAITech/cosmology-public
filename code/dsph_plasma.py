"""dsph_plasma.py - Workbench Entry 141
Stars with nothing to push: can electromagnetic forces hold the dwarf spheroidal galaxies together?
Data: McConnachie (2012) catalogue of Local Group dwarfs, VizieR J/AJ/144/4/catalog, read in a browser; the 46 galaxies with a
measured stellar velocity dispersion, half-light radius and stellar mass are copied below exactly as served
(name, subgroup, type, distance kpc, half-light radius R_e pc, stellar mass 1e6 Msun, upper-limit flag, sigma km/s, +err, -err,
HI mass 1e6 Msun, HI flag, ellipticity). The catalogue's stellar masses assume M/L_V = 1.
For each dwarf spheroidal:
  mass within the 3D half-light radius  M_1/2 = 4 sigma^2 R_e / G = 930 sigma^2 R_e   (Wolf et al. 2010; r_1/2 = 4R_e/3)
  missing mass                         M_miss = M_1/2 - M_*/2 - M_HI/2
  (a) could charges on the stars do it? A star in hydrostatic balance carries the Pannekoek-Rosseland charge,
      q/m = 2 pi eps0 G m_p / e = 3.9e-29 C/kg. Field needed to hold a star on an orbit of radius r_1/2 at speed sigma:
      E = sigma^2 / (r q/m), or B = sigma / (r q/m) with the Lorentz force.
  (b) could the missing mass be plasma, which the fields could then act on? Spread uniformly inside r_1/2 at T ~ 10^4 K
      (the temperature that matches sigma), ionized: n_e = 1.2 n_H, n_H = M_miss X / (m_p V), X = 0.75.
      Emission measure through the centre EM = n_e n_p 2 r_1/2; WHAM's 3-sigma limit is 0.15 R, about 0.3 pc cm^-6.
      If neutral instead, it is HI: GBT limits are M_HI/M_dyn ~ 5e-5 (Spekkens et al. 2014).
  (c) could the dwarfs be unbound and flying apart? crossing time r_1/2/sigma against stellar ages of ~10 Gyr.
"""
import json, math, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
RAW = """Sagittarius dSph,MW,dSph,26,2587,21.00000,,11.4,0.7,0.7,,b,0.64
Segue I,MW,dSph,23,29,0.00034,,3.9,0.8,0.8,,,0.48
Ursa Major II,MW,dSph,32,149,0.00410,,6.7,1.4,1.4,,,0.63
Bootes II,MW,dSph,42,51,0.00100,,10.5,7.4,7.4,,,0.21
Segue II,MW,dSph,35,35,0.00086,,3.4,2.5,1.2,,,0.15
Willman 1,MW,dSph,38,25,0.00100,,4.3,2.3,1.3,,,0.47
Coma Berenices,MW,dSph,44,77,0.00370,,4.6,0.8,0.8,,,0.38
Bootes I,MW,dSph,66,242,0.02900,,2.4,0.9,0.5,,,0.39
Draco,MW,dSph,76,221,0.29000,,9.1,1.2,1.2,,,0.31
Ursa Minor,MW,dSph,76,181,0.29000,,9.5,1.2,1.2,,,0.56
Sculptor,MW,dSph,86,283,2.30000,,9.2,1.4,1.4,2.2e-01,e,0.32
Sextans I,MW,dSph,86,695,0.44000,,7.9,1.3,1.3,,,0.35
Ursa Major I,MW,dSph,97,319,0.01400,,7.6,1.0,1.0,,,0.80
Carina,MW,dSph,105,250,0.38000,,6.6,1.2,1.2,,,0.33
Hercules,MW,dSph,132,330,0.03700,,3.7,0.9,0.9,,,0.68
Fornax,MW,dSph,147,710,20.00000,,11.7,0.9,0.9,1.7e-01,f,0.30
Leo IV,MW,dSph,154,206,0.01900,,3.3,1.7,1.7,,,0.49
Canes Venatici II,MW,dSph,160,74,0.00790,,4.6,1.0,1.0,,,0.52
Leo V,MW,dSph,178,135,0.01100,,3.7,2.3,1.4,,,0.50
Canes Venatici I,MW,dSph,218,564,0.23000,,7.6,0.4,0.4,,,0.39
Leo II,MW,dSph,233,176,0.74000,,6.6,0.7,0.7,,,0.13
Leo I,MW,dSph,254,251,5.50000,,9.2,1.4,1.4,,,0.21
M32,M31,cE,805,110,320.00000,,92.0,5.0,5.0,,,0.25
Andromeda IX,M31,dSph,766,557,0.15000,,4.5,3.6,3.6,,,
NGC 205,M31,dEdSph,824,590,330.00000,,35.0,5.0,5.0,4.0e-01,,0.43
Andromeda I,M31,dSph,745,672,3.90000,,10.6,1.1,1.1,,,0.22
Andromeda III,M31,dSph,748,479,0.83000,,4.7,1.8,1.8,,,0.52
Andromeda XI,M31,dSph,759,157,0.04900,L,4.6,,,,,
Andromeda V,M31,dSph,773,315,0.39000,,11.5,5.3,4.4,,,0.18
Andromeda X,M31,dSph,701,265,0.09600,,3.9,1.2,1.2,,,0.44
Andromeda XII,M31,dSph,871,304,0.03100,,2.6,5.1,2.6,,,
NGC 147,M31,dEdSph,676,623,62.00000,,16.0,1.0,1.0,,,0.41
Andromeda XIV,M31,dSph,735,363,0.20000,,5.4,1.3,1.3,,,0.31
Andromeda XV,M31,dSph,631,222,0.49000,,11.0,7.0,5.0,,,0.00
Andromeda XIII,M31,dSph,912,207,0.04100,,9.7,8.9,4.5,,,
Andromeda II,M31,dSph,652,1176,7.60000,,7.3,0.8,0.8,,,0.20
NGC 185,M31,dEdSph,617,458,68.00000,,24.0,1.0,1.0,1.1e-01,,0.15
Andromeda VII,M31,dSph,762,776,9.50000,,9.7,1.6,1.6,,,0.13
LGS 3 Local Group Suspect 3,M31,dIrrdSph,769,470,0.96000,,7.9,5.3,2.9,3.8e-01,,0.20
Andromeda VI,M31,dSph,783,524,2.80000,,9.4,3.2,2.4,,,0.41
Andromeda XVI,M31,dSph,525,136,0.41000,L,10.0,,,,,0.00
Cetus,Rest,dSph,755,703,2.60000,,17.0,2.0,2.0,,,0.33
Leo T,Rest,dIrrdSph,417,120,0.14000,,7.5,1.6,1.6,2.8e-01,,0.00
WLM Wolf-Lundmark-Melotte,Rest,dIrr,933,2111,43.00000,,17.5,2.0,2.0,6.1e01,,0.65
Leo A,Rest,dIrr,798,499,6.00000,,9.3,1.3,1.3,1.1e01,,0.40
Tucana,Rest,dSph,887,284,0.56000,,15.8,4.1,3.1,,,0.48"""
G = 4.30091e-3                                 # pc Msun^-1 (km/s)^2
pc, Msun, mp, kB = 3.0857e16, 1.989e30, 1.6726e-27, 1.380649e-23
eps0, Gsi, e = 8.8541878e-12, 6.674e-11, 1.602176634e-19
qm = 2*math.pi*eps0*Gsi*mp/e                   # C/kg, Pannekoek-Rosseland
EM_lim = 0.3                                   # pc cm^-6, WHAM 3-sigma (0.15 R)
rows = []
for line in RAW.split('\n'):
    p = line.split(',')
    name, sub, typ = p[0], p[1], p[2]
    D, Re, Mst, lim, sig = float(p[3]), float(p[4]), float(p[5])*1e6, p[6] == 'L', float(p[7])
    MHI = float(p[10])*1e6 if p[10] else 0.0
    M12 = 4*sig**2*Re/G
    r12 = 4*Re/3
    miss = M12 - 0.5*Mst - 0.5*MHI
    o = dict(name=name, sub=sub, type=typ, D_kpc=D, Re_pc=Re, Mstar=Mst, sigma=sig, sigma_upper_limit=lim, MHI=MHI, M12=M12,
             ML12=M12/(0.5*Mst), miss=miss)
    if miss > 0:
        V = 4/3*math.pi*(r12*pc)**3
        nH = 0.75*miss*Msun/(mp*V)/1e6            # cm^-3
        ne = 1.2*nH
        o['nH_cm3'] = nH; o['EM'] = ne*nH*2*r12; o['EM_over_lim'] = o['EM']/EM_lim
        o['T_gas_K'] = 0.6*mp*(sig*1e3)**2/kB
    r = r12*pc; v = sig*1e3
    o['E_needed_Vm'] = v**2/(r*qm); o['B_needed_T'] = v/(r*qm)
    o['t_cross_Myr'] = r/v/3.156e13
    rows.append(o)
ds = [o for o in rows if o['type'] == 'dSph' and not o['sigma_upper_limit']]
def q(a, p): return float(np.percentile(a, p))
def st(k, sel): a = [o[k] for o in sel if k in o]; return [q(a, 5), q(a, 50), q(a, 95), min(a), max(a)]
res = dict(n_rows=len(rows), n_dsph=len(ds), qm_C_per_kg=qm,
           ML12=st('ML12', ds), M12=st('M12', ds), EM_over_lim=st('EM_over_lim', ds), nH=st('nH_cm3', ds), T_gas=st('T_gas_K', ds),
           E_needed=st('E_needed_Vm', ds), B_needed=st('B_needed_T', ds), t_cross=st('t_cross_Myr', ds),
           n_ML_gt_10=sum(o['ML12'] > 10 for o in ds), n_ML_gt_100=sum(o['ML12'] > 100 for o in ds), n_ML_lt_3=sum(o['ML12'] < 3 for o in ds),
           HI_detected=[(o['name'], o['MHI'], o['M12']) for o in rows if o['MHI'] > 0],
           others=[(o['name'], o['type'], round(o['ML12'], 1)) for o in rows if o['type'] != 'dSph'],
           table=[dict(name=o['name'], sub=o['sub'], Re=o['Re_pc'], sigma=o['sigma'], Mstar=o['Mstar'], M12=round(o['M12'], -3), ML12=round(o['ML12'], 1),
                       EM=round(o.get('EM', float('nan')), 0), B_T=o['B_needed_T'], tcross=round(o['t_cross_Myr'], 1)) for o in ds])
# Fig data: M_1/2 against stellar mass inside r_1/2, all 46
res['fig'] = [[o['name'], round(0.5*o['Mstar'], 1), round(o['M12'], 0), o['type'], o['sigma_upper_limit']] for o in rows]
json.dump(res, open(os.path.join(HERE, 'dsph_plasma.json'), 'w'), indent=1)
for k in ('n_rows', 'n_dsph', 'qm_C_per_kg', 'ML12', 'M12', 'EM_over_lim', 'nH', 'T_gas', 'E_needed', 'B_needed', 't_cross', 'n_ML_gt_10', 'n_ML_gt_100', 'n_ML_lt_3', 'HI_detected', 'others'):
    print(k, res[k])
for o in ds: print(f"{o['name']:<22} Re={o['Re_pc']:5.0f} s={o['sigma']:5.1f} M*={o['Mstar']:.2e} M12={o['M12']:.2e} ML={o['ML12']:8.1f} EM={o.get('EM',0):9.1f} B={o['B_needed_T']:.1e} tc={o['t_cross_Myr']:.0f}")
