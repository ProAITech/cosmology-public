"""Entry 86: is the CMB lensed from behind redshift 4?

Published lensing cross-correlations (CMB lensing x tracer at known redshift)
are compared with (a) the standard prediction and (b) the transmission that
Lerner's stated intergalactic fog would allow at CMB frequencies.
A photon absorbed and re-emitted by a fog loses its direction, so only the
fraction exp(-tau) of CMB photons emitted behind a lens carries its imprint.
Also: static-law distances to the lens redshifts, and what radio quasars at
z ~ 6 say about absorption that falls with frequency.
"""
import json, math
import numpy as np

C = 299792.458  # km/s
# ---------- distances ----------
H0_static = 70.0                 # Lerner's linear law, d = cz/H0
DH_static = C / H0_static / 1000  # Gpc, 4.283
H0p, Om = 67.66, 0.3111          # Planck 2018 (TT,TE,EE+lowE+lensing)
def E(z): return math.sqrt(Om*(1+z)**3 + 1-Om)
def chi_std(z, n=20000):
    zz = np.linspace(0, z, n+1); f = 1/np.sqrt(Om*(1+zz)**3 + 1-Om)
    return C/H0p/1000 * np.trapezoid(f, zz)
def d_linear(z): return DH_static*z
def d_tired(z): return DH_static*math.log(1+z)

# ---------- linear growth in flat LCDM ----------
def growth(z):
    # D(a) ∝ E(a) ∫_0^a da' / (a' E(a'))^3
    a = 1/(1+z)
    def D(a):
        aa = np.linspace(1e-6, a, 40001)
        Ea = np.sqrt(Om/aa**3 + 1-Om)
        return math.sqrt(Om/a**3 + 1-Om) * np.trapezoid(1/(aa*Ea)**3, aa)
    return D(a)/D(1.0)
s8_fid = 0.8102

# ---------- lensing measurements ----------
# Quaia x Planck PR4 lensing (GMV), Piccirilli et al. 2024, sigma8(z) per bin
quaia = [(0.69, 0.65, 0.10, 13.5, 8.7), (1.59, 0.39, 0.06, 20.0, 13.5), (2.72, 0.22, 0.06, 5.7, 5.1)]
rows = []
for z, s, e, sn, sn_pol in quaia:
    pred = s8_fid*growth(z)
    rows.append(dict(sample='Quaia quasars', z=z, sn=sn, sn_noTT=sn_pol, s8_meas=s, s8_err=e,
                     s8_pred=round(pred, 3), A=round(s/pred, 3), A_err=round(e/pred, 3)))
# HSC Lyman-break galaxies x Planck PR3 lensing, Miyatake et al. 2022: sigma8 (z=0 parameter)
mi = dict(sample='HSC Lyman-break galaxies', z=3.8, sn=5.1, s8_meas=0.46, s8_up=0.25, s8_lo=0.20, s8_pred=0.811)
mi['A'] = round(mi['s8_meas']/mi['s8_pred'], 3); mi['A_up'] = round(0.25/0.811, 3); mi['A_lo'] = round(0.20/0.811, 3)
# Combined Quaia amplitude relative to Planck
quaia_all = dict(A=0.97, err=0.09, A_noTT=1.12, err_noTT=0.09)

# ---------- Lerner's stated fog (Entry 22 table) ----------
# tau = 2.3 within 300 Mpc (radio lost tenfold); 7.6 by z = 2.2 +/- 1.0; 10 by z = 4.3 +/- 2; band 150 MHz to >200 GHz.
def tau_fog(z, zA=2.2, zB=4.3):
    d = d_linear(z); pts = [(0.3, 2.3), (d_linear(zA), 7.6), (d_linear(zB), 10.0)]
    if d <= pts[0][0]: return 2.3*d/0.3
    for (d1, t1), (d2, t2) in zip(pts, pts[1:]):
        if d <= d2: return t1 + (t2-t1)*math.log(d/d1)/math.log(d2/d1)
    (d1, t1), (d2, t2) = pts[-2], pts[-1]
    return t2 + (t2-t1)*math.log(d/d2)/math.log(d2/d1)
def T_band(z):  # transmission: central, most transparent (zA=3.2,zB=6.3), most opaque (zA=1.2,zB=2.3)
    return (math.exp(-tau_fog(z)), math.exp(-tau_fog(z, 3.2, 6.3)), math.exp(-tau_fog(z, 1.2, 2.3)))

for r in rows + [mi]:
    c, hi, lo = T_band(r['z']); r['fog_T'] = [float('%.3g' % c), float('%.3g' % hi), float('%.3g' % lo)]
    r['d_linear_Gpc'] = round(d_linear(r['z']), 2); r['d_tired_Gpc'] = round(d_tired(r['z']), 2); r['chi_std_Gpc'] = round(chi_std(r['z']), 2)
    # how many times larger the measured amplitude is than the most transparent fog allows
    r['ratio_to_fog_max'] = round(r['A']/r['fog_T'][1], 1)
# upper limits on CMB-band optical depth in front of each lens (2 sigma, measured amplitude vs no fog)
for r in rows:
    lo2 = r['A'] - 2*r['A_err']; r['tau_2sig_max'] = round(-math.log(min(lo2, 1.0)), 2) if lo2 > 0 else None
mi_lo2 = mi['A'] - 2*mi['A_lo']; mi['tau_2sig_max'] = round(-math.log(mi_lo2), 2) if mi_lo2 > 0 else None
mi['sig_vs_zero'] = 5.1

curve = []
for z in np.linspace(0.05, 5.0, 100):
    c, hi, lo = T_band(z); curve.append([round(float(z), 3), float('%.3g' % c), float('%.3g' % hi), float('%.3g' % lo)])

# ---------- radio quasars at z ~ 6 ----------
radio = {}
# grey fog in Lerner's band: radio loudness R = f(5 GHz)/f(4400 A) suppressed by exp(-tau) at z~6 (optical is out of band)
for z in (5.7, 6.0, 6.82):
    c, hi, lo = T_band(z); radio['z%.2f' % z] = dict(suppression_central=float('%.3g' % c), suppression_least=float('%.3g' % hi))
# absorption falling with frequency, tau ∝ nu^-2.1 (free-free): no turnover at 144 MHz -> tau(144 MHz) <~ 1
radio['tau150GHz_if_tau144MHz_is_1'] = float('%.3g' % (144e6/150e9)**2.1)
radio['flattening_144_1400_per_unit_tau144'] = round((1-(144/1400)**2.1)/math.log(1400/144), 3)
radio['RLF_z6'] = '8.1 (+5.0/-3.2)%'


# ---------- lensing ratio: CMB lensing / galaxy lensing around the same lenses (Prat et al. 2019) ----------
# DES Y1 redMaGiC lenses z 0.15-0.30, 0.30-0.45, 0.45-0.60; sources 0.63-0.90 and 0.90-1.30; SPT 150 GHz + Planck 143 GHz.
# The ratio depends on distances only, not on how much mass the lenses hold. Measured A = 1.1 +/- 0.1 relative to Planck LCDM.
# Galaxy shear is measured in visible light, outside the fog's band, so a fog divides the ratio by exp(-tau(z_lens)).
lenses = [0.225, 0.375, 0.525]; sources = [0.765, 1.10]
chi_star = chi_std(1090)
def beta(dl, ds): return (ds-dl)/ds
geo = {}
for name, dist, dstar in (('standard', chi_std, chi_star), ('tired light, CMB at infinity', d_tired, None), ('Lerner linear law, CMB at infinity', d_linear, None)):
    rs = []; rf = []
    for zl in lenses:
        for zs in sources:
            bc = 1.0 if dstar is None else beta(dist(zl), dstar)
            r = bc/beta(dist(zl), dist(zs)); rs.append(r); rf.append(r*math.exp(-tau_fog(zl, 3.2, 6.3)))
    geo[name] = dict(ratio=float(np.mean(rs)), ratio_with_least_fog=float(np.mean(rf)))
ref = geo['standard']['ratio']
for k in geo:
    geo[k]['A_no_fog'] = round(geo[k]['ratio']/ref, 3); geo[k]['A_with_fog'] = float('%.3g' % (geo[k]['ratio_with_least_fog']/ref))
    del geo[k]['ratio'], geo[k]['ratio_with_least_fog']
lensing_ratio = dict(A_meas=1.1, err=0.1, models=geo, tau_fog_least=[round(tau_fog(z, 3.2, 6.3), 2) for z in lenses],
                     tau_2sig_max=round(-math.log(1.1-2*0.1), 3))
print(lensing_ratio)

out = dict(lensing_ratio=lensing_ratio, DH_static_Gpc=round(DH_static, 3), rows=rows, miyatake=mi, quaia_all=quaia_all, fog_curve=curve, radio=radio,
           distances={str(z): dict(linear=round(d_linear(z), 2), tired=round(d_tired(z), 2), std=round(chi_std(z), 2)) for z in (0.69, 1.59, 2.72, 3.8, 4.0, 6.0, 6.82, 1100)})
json.dump(out, open(__file__.replace('.py', '.json'), 'w'), indent=1)
for r in rows + [mi]: print(r)
print(radio); print(out['distances'])
