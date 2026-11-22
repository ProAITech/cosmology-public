"""Entry 90: does the corona redden sunlight? Brynjolfsson's plasma redshift at the Sun.

Brynjolfsson (2004, astro-ph/0401420): photons crossing hot, sparse plasma lose energy
d ln E = -3.326e-25 N_e dx (half the Thomson cross-section per electron), for wavelengths
shorter than a cutoff lambda_0.5 = 3.18e-6 T / sqrt(N_e) cm (his Eq. 16). He holds that this,
not gravity, makes the Sun's redshift: photons "lose their gravitational redshift" in transit.
Here: the plasma redshift of sunlight through a standard corona (Baumbach-Allen), at disk centre
and for the whole disk; its sensitivity to changes of the corona; and the cutoff check.
"""
import json, math, os, numpy as np
os.chdir(os.path.dirname(os.path.abspath(__file__)))
C = 299792.458e3                 # m/s
RSUN = 6.957e10                  # cm
SIG = 3.326e-25                  # cm^2, Brynjolfsson's coefficient
def ne_ba(r):                    # Baumbach-Allen quiet corona, r in solar radii, cm^-3
    return 1e8*(2.99*r**-16 + 1.55*r**-6 + 0.036*r**-1.5)
def column(mu, scale=1.0, smax=215.0, n=200000):
    s = np.concatenate([np.linspace(0, 2, n//2, endpoint=False), np.linspace(2, smax, n//2)])
    r = np.sqrt(1 + 2*s*mu + s*s)
    y = ne_ba(r)*scale
    return float(np.trapezoid(y, s)*RSUN)
def z_of_mu(mu, scale=1.0): return SIG*column(mu, scale)
res = {}
res['column_centre_cm2'] = float('%.3g' % column(1.0))
res['z_centre'] = float('%.3g' % z_of_mu(1.0)); res['v_centre_ms'] = round(z_of_mu(1.0)*C, 1)
mus = np.linspace(0.02, 1.0, 120)
zs = np.array([z_of_mu(m) for m in mus])
u = 0.6                                                   # visible limb darkening
w = (1-u*(1-mus))*mus
res['z_disk'] = float('%.3g' % (np.trapezoid(zs*w, mus)/np.trapezoid(w, mus)))
res['v_disk_ms'] = round(res['z_disk']*C, 1)
res['v_by_mu_ms'] = {f'{m:.2f}': round(z_of_mu(m)*C, 1) for m in (1.0, 0.7, 0.5, 0.3, 0.1, 0.05)}
# Brynjolfsson's own Table 3 estimate at disk centre: 1.00e-6
res['brynjolfsson_centre_z'] = 1.00e-6; res['brynjolfsson_centre_v_ms'] = round(1.00e-6*C, 1)
# sensitivity: fractional change of the corona's density -> change in disk-integrated shift
res['dv_per_10pct_ms'] = {'Baumbach-Allen': round(0.1*res['v_disk_ms'], 1), 'Brynjolfsson centre': round(0.1*1.00e-6*C, 1)}
res['dv_for_modulations_ms'] = {f: [round((f-1)*res['v_disk_ms']), round((f-1)*1.00e-6*C)] for f in (1.3, 1.6, 4.3)}
# cutoff wavelength along the line of sight from disk centre, coronal temperature 1-2 MK
cut = {}
for T in (1.0e6, 2.0e6):
    for r in (1.01, 1.1, 1.3, 2.0):
        cut[f'T={T:.0e},r={r}'] = round(3.18e-6*T/math.sqrt(ne_ba(r))*1e4, 2)   # micrometres
res['cutoff_um'] = cut
# share of the centre column inside r < 1.3 (where the cutoff is several micrometres)
s = np.linspace(0, 0.3, 30001); res['share_column_below_1p3'] = round(float(np.trapezoid(ne_ba(1+s), s)*RSUN/column(1.0)), 3)
# observed: solar redshift 638 +/- 6 m/s (weak Fe lines, HARPS-LFC Moon spectra), GR 633.1 m/s
gr = 633.1; obs, err = 638.0, 6.0
res['gr_ms'] = gr; res['observed_ms'] = [obs, err]
res['plasma_allowed_on_top_of_GR_2sigma_ms'] = round(obs - gr + 2*err, 1)
res['plasma_allowed_fraction_of_BA'] = round((obs-gr+2*err)/res['v_disk_ms'], 3)
# Brynjolfsson: no net gravitational shift, so the whole observed shift would be plasma
res['plasma_only_needs_column_factor'] = round(obs/res['v_disk_ms'], 2)
# HARPS-N Sun-as-a-star 2015-2018: day-to-day MAD 1.07 m/s; range within about +/-6 m/s
res['harpsn_daily_mad_ms'] = 1.07
res['allowed_column_change_daily'] = round(1.07/res['v_disk_ms'], 4)
json.dump(res, open('solar_plasma_z.json', 'w'), indent=1)
for k, v in res.items(): print(k, v)
