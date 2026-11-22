"""Entry 157: are stars hot inside? The large frequency separation of solar-like oscillations against a cool interior.
The p-mode spacing is dnu = [2 * integral_0^R dr/c]^-1. If no layer were hotter than the photosphere, c <= c(Teff) everywhere,
so dnu <= c(Teff)/(2R): the cool-interior ceiling. Conversely the hottest layer must reach at least T_min = Teff (dnu/dnu_cool)^2
(ideal gas, mu = 0.6, gamma = 5/3; no magnetic support). Per-star work was done in the browser (astero.js); this script
re-derives the solar case, the magnetic escape, and summarises."""
import json, math, sys
d = json.load(open(sys.argv[1] if len(sys.argv) > 1 else 'astero_raw.json'))
k, mp, G = 1.380649e-23, 1.67262e-27, 6.674e-11
Rsun, Msun, Teff_sun, dnu_sun = 6.957e8, 1.989e30, 5772.0, 134.9e-6   # solar dnu from Sun-as-a-star helioseismology
def c_gas(T, mu=0.6, gam=5/3): return math.sqrt(gam*k*T/(mu*mp))
c_ph = c_gas(Teff_sun); dcool = c_ph/(2*Rsun)
ratio_sun = dnu_sun/dcool; Tmin_sun = Teff_sun*ratio_sun**2
tau = 1/(2*dnu_sun); cbar = Rsun/tau
rho_bar = Msun/(4/3*math.pi*Rsun**3)                     # kg/m^3
B_needed = cbar*math.sqrt(4e-7*math.pi*rho_bar)            # tesla, if waves were Alfvenic at the mean density
P_B = B_needed**2/(2*4e-7*math.pi)                         # Pa
T_vir = G*Msun*0.6*mp/(k*Rsun)
out = dict(sun=dict(c_photosphere_kms=round(c_ph/1e3, 2), dnu_cool_uHz=round(dcool*1e6, 2), ratio=round(ratio_sun, 1), Tmin_MK=round(Tmin_sun/1e6, 2),
                    sound_crossing_s=round(tau), mean_speed_kms=round(cbar/1e3), B_alfvenic_MG=round(B_needed*1e4/1e6, 1), magnetic_pressure_Pa=float(f'{P_B:.2g}'),
                    T_virial_MK=round(T_vir/1e6, 1)))
# dnu predicted for neutral photospheric gas (mu = 1.3) as a further-limiting case
out['sun']['dnu_cool_neutral_uHz'] = round(c_gas(Teff_sun, mu=1.3)/(2*Rsun)*1e6, 2)
out['fit'] = dict(slope_R=d['fit_logdnu_logR'][1], slope_sd=d['fit_slope_bootstrap_sd'], with_mass=d['fit_logdnu_logR_logM'][1:],
                  hydrostatic_expected=[-1.5, 0.5], cool_interior_expected=[-1.0, 0.0],
                  sigma_from_cool=round((d['fit_logdnu_logR'][1]-(-1.0))/d['fit_slope_bootstrap_sd'], 0))
out['dwarfs'] = d['dwarfs']; out['giants'] = d['giants']; out['n'] = d['n']; out['n_ratio_below_3'] = d['n_ratio_below_3']
out['share_below_3'] = round(d['n_ratio_below_3']/d['n']*100, 2)
print(json.dumps(out, indent=1))
json.dump(dict(summary=out, bins=d['bins'], hist=dict(start=d['hist_logT_start'], step=d['hist_step'], dwarfs=d['hist_dwarfs'], giants=d['hist_giants'])), open('astero.json', 'w'))
