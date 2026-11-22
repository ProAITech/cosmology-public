"""oldest.py - Workbench Entry 99
Is anything older than the universe? The oldest dated objects nearby, and the oldest stellar populations found at each redshift.
Lerner: the universe holds objects "10 times older" than the Big Bang allows. An unchanging, infinitely old universe has no
age ceiling, and galaxies seen far away should include old ones like nearby ellipticals. An expanding one caps the age of any
population seen at redshift z by the time since the Big Bang, t(z).
Ages (all as published): Milky Way globular clusters, best-fitting ages 11.5-13.5 Gyr, typical error 0.5-0.75 (Ying et al. 2025);
M92 13.80 +- 0.75 (Ying et al. 2023); NGC 6362 13.6 +- 0.5 (Oliveira et al. 2020); a halo white dwarf 13.9 +- 0.8 (Fouesneau et al. 2019);
ultrafaint dwarf galaxies 13.7-14.1 (Brown et al. 2014); the last four as collected by Cimatti & Moresco (2023).
Distant: z = 2.675, half its stars formed when the universe was 108 (+37/-39) Myr old (McConachie et al. 2025);
z = 3.205, stars formed ~1.5 Gyr before (Glazebrook et al. 2024); z = 4.658, ~200 Myr of star formation ending at z = 6.5 (+0.2/-0.5)
(Carnall et al. 2023); z = 7.29, a burst at z ~ 8-9, quenched 50-100 Myr before (Weibel et al. 2025).
"""
import json, math, os
import numpy as np
from scipy import integrate, optimize
HERE = os.path.dirname(os.path.abspath(__file__))
def age_of(z, H0=67.36, Om=0.3153):
    OL = 1 - Om; tH = 977.79/H0
    return 2/(3*math.sqrt(OL))*tH*math.asinh(math.sqrt(OL/Om)*(1+z)**-1.5)
def z_of(t, **kw): return optimize.brentq(lambda z: age_of(z, **kw) - t, 0.0, 200.0)
t0 = age_of(0)
LOCAL = [('Milky Way globular clusters, oldest (Ying et al. 2025)', 13.5, 0.6), ('M92 (Ying et al. 2023)', 13.80, 0.75), ('NGC 6362 (Oliveira et al. 2020)', 13.6, 0.5),
         ('halo white dwarf (Fouesneau et al. 2019)', 13.9, 0.8), ('ultrafaint dwarfs (Brown et al. 2014)', 13.9, 0.5)]
w = np.array([1/e**2 for _, a, e in LOCAL]); a = np.array([a for _, a, e in LOCAL])
mean = float(np.sum(w*a)/w.sum()); err = float(w.sum()**-0.5)
out = {'t0_planck': round(t0, 3), 'local': LOCAL, 'local_mean': [round(mean, 2), round(err, 2)]}
# H0 from the local ages: t0 = age + formation delay (0.1-0.4 Gyr, z_f 11-30), flat LCDM with Om = 0.30 (+- 0.02)
for Om in (0.28, 0.30, 0.32):
    for delay in (0.1, 0.4):
        tU = mean + delay
        H = optimize.brentq(lambda H: age_of(0, H0=H, Om=Om) - tU, 40, 120)
        out.setdefault('H0_from_ages', {})['Om=%.2f delay=%.1f' % (Om, delay)] = round(H, 1)
# what t0 SH0ES's H0 gives
out['t0_at_H0_73'] = {('Om=%.2f' % Om): round(age_of(0, H0=73.0, Om=Om), 2) for Om in (0.28, 0.30, 0.32)}
out['sigma_oldest_vs_H0_73'] = round((mean + 0.1 - age_of(0, H0=73.0, Om=0.30))/err, 2)
# distant: the oldest populations found, as age at the time seen and formation time after the Big Bang (Planck cosmology)
D = []
for name, z, tf_lo, tf_hi, note in (('McConachie et al. 2025', 2.675, 0.069, 0.145, 'half the stars formed at 108 (+37, -39) Myr'),
                                    ('Glazebrook et al. 2024 (ZF-UDS-7329)', 3.205, None, None, 'formed ~1.5 Gyr before, z ~ 11'),
                                    ('Carnall et al. 2023 (GS-9209)', 4.658, None, None, '~200 Myr of star formation ending at z = 6.5 (+0.2, -0.5)'),
                                    ('Weibel et al. 2025 (RUBIES-UDS-QG-z7)', 7.29, None, None, 'burst at z ~ 8-9, quenched 50-100 Myr before')):
    tz = age_of(z)
    if name.startswith('McConachie'): tf, lo, hi = 0.108, tf_lo, tf_hi
    elif name.startswith('Glazebrook'): tf = tz - 1.5; lo, hi = tf - 0.2, tf + 0.2
    elif name.startswith('Carnall'): tq = age_of(6.5); tf = tq - 0.1; lo, hi = age_of(7.0) - 0.2, age_of(6.3) - 0.0
    else: tf = 0.5*(age_of(9.0) + age_of(8.0)); lo, hi = age_of(9.0), age_of(8.0)
    D.append(dict(name=name, z=z, t_z=round(tz, 3), t_form=round(tf, 3), t_form_range=[round(lo, 3), round(hi, 3)], age=round(tz - tf, 3),
                  age_range=[round(tz - hi, 3), round(tz - lo, 3)], z_form=round(z_of(tf), 1), note=note))
    print('%-40s z %.3f  t(z) %.2f Gyr  formed at %.2f Gyr (z_f %.1f)  age %.2f Gyr' % (name, z, tz, tf, z_of(tf), tz - tf))
out['distant'] = D
print('local oldest mean %.2f +- %.2f Gyr; Planck t0 %.2f' % (mean, err, t0)); print(out['H0_from_ages']); print('t0 at H0 = 73:', out['t0_at_H0_73'], 'sigma', out['sigma_oldest_vs_H0_73'])
zz = np.round(np.concatenate([np.linspace(0, 1, 21), np.linspace(1.1, 9, 80)]), 3)
out['curve'] = [[float(z), round(age_of(z), 4), round(age_of(z, H0=73.0, Om=0.30), 4)] for z in zz]
json.dump(out, open(os.path.join(HERE, 'oldest.json'), 'w'), indent=1)
