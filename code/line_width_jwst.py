"""line_width_jwst.py - Workbench Entry 59
Two 2025-26 papers (Raikov, Tsymbal & Lovyagin, arXiv:2507.19651; Tsymbal, Raikov & Lovyagin, arXiv:2604.27867) find
that JWST galaxies keep one linear size at all redshifts if distance follows static tired light, and ask for a test
of sizes and masses that does not depend on the adopted cosmology. A line width is such a test: it does not change
with the distance assumed, while size and luminosity do.

 1. The catalogue result, as re-derived in a browser tab from VizieR J/A+A/691/A240 (astrodeep_sizes.js), typed in.
 2. Seventeen galaxies at redshift 4.0 to 7.6 with resolved sizes and line widths from JWST's high-resolution
    spectrograph (de Graaff et al. 2024; Saldana-Lopez et al. 2025), re-read under each distance law:
    size, stellar mass, gas mass, dynamical mass.
"""
import json, math
import numpy as np
from scipy import integrate

C, H0, OM = 299792.458, 70.0, 0.3
CH = C / H0                                      # Mpc
zg = np.linspace(0, 12, 12001); _dc = integrate.cumulative_trapezoid(1 / np.sqrt(OM * (1 + zg) ** 3 + 1 - OM), zg, initial=0)
DM = lambda z: CH * float(np.interp(z, zg, _dc))
LAWS = {   # angular-size distance and luminosity distance, Mpc
 "standard":                      (lambda z: DM(z) / (1 + z),        lambda z: DM(z) * (1 + z)),
 "tired light":                   (lambda z: CH * math.log(1 + z),   lambda z: CH * math.log(1 + z) * math.sqrt(1 + z)),
 "tired light, stretched time":   (lambda z: CH * math.log(1 + z),   lambda z: CH * math.log(1 + z) * (1 + z)),
 "Lerner's straight line":        (lambda z: CH * z,                 lambda z: CH * z * math.sqrt(1 + z)),
}
G = 4.30091e-3                                   # pc (km/s)^2 / Msun
ARC = math.pi / 180 / 3600
def shen(m):                                     # Shen et al. 2003, late types: half-light radius in kpc
    return 0.1 * m ** 0.14 * (1 + m / 3.98e10) ** 0.25
# name, z, r_e [arcsec], sigma [km/s], log M_dyn, log M*, SFR [Msun/yr], source
GAL = [
 ("JADES-NS-00016745", 5.566, 0.32, 55.0, 10.23, 8.80, 12.6, "dG24"), ("JADES-NS-00019606", 5.890, 0.117, 39.1, 9.17, 7.54, 1.46, "dG24"),
 ("JADES-NS-00022251", 5.799, 0.122, 39.0, 9.23, 8.21, 6.2, "dG24"), ("JADES-NS-00047100", 7.432, 0.130, 71.0, 9.85, 8.53, 5.1, "dG24"),
 ("JADES-NS-10016374", 5.504, 0.084, 62.0, 9.45, 7.86, 2.19, "dG24"), ("JADES-NS-20086025", 7.263, 0.36, 25.0, 10.31, 8.85, 24.0, "dG24"),
 ("1871-10", 7.5990, 0.14, 65.2, 9.54, 7.96, 6.0, "SL25"), ("1871-11", 7.6100, 0.17, 60.9, 9.77, 9.07, 9.7, "SL25"), ("1871-12", 7.5020, 0.07, 67.3, 9.41, 8.31, 19.9, "SL25"),
 ("1871-29", 5.0191, 0.13, 72.0, 9.96, 9.24, 10.7, "SL25"), ("1871-40", 7.0948, 0.14, 60.6, 9.73, 8.47, 9.2, "SL25"), ("1871-70", 7.0302, 0.12, 41.3, 9.42, 8.85, 1.7, "SL25"),
 ("1871-105", 6.5669, 0.07, 49.9, 9.34, 7.81, 2.8, "SL25"), ("1871-451", 4.4615, 0.08, 40.7, 9.37, 7.99, 1.7, "SL25"), ("1871-545", 4.0416, 0.12, 96.4, 10.11, 9.41, 66.7, "SL25"),
 ("1871-912", 6.7161, 0.18, 44.4, 9.68, 9.69, 27.1, "SL25"), ("1871-1032", 7.0889, 0.11, 37.3, 9.41, 9.07, 6.3, "SL25"),
]   # SFR for the SL25 rows: the larger of the 10 Myr and 100 Myr values, which gives the gas its best case
out = {"galaxies": [], "laws": {}}
print("How each law rescales a galaxy read with standard distances")
print("  z     law                           size    stars   gas     dynamical   dynamical/stars")
fac = {}
for z in (1.0, 2.0, 4.0, 6.0, 7.5):
    for name, (da, dl) in LAWS.items():
        fA, fL = da(z) / LAWS["standard"][0](z), dl(z) / LAWS["standard"][1](z)
        fgas = fA ** (2 - 2 / 1.4) * fL ** (2 / 1.4)
        fac.setdefault(name, []).append([z, fA, fL ** 2, fgas, fA, fA / fL ** 2])
        if name != "standard": print("  %.1f   %-28s  x%5.1f  x%5.2f  x%5.2f  x%5.1f      x%6.1f" % (z, name, fA, fL ** 2, fgas, fA, fA / fL ** 2))
out["factors"] = fac
rows = {k: [] for k in LAWS}
for nm, z, re_as, sig, lmd, lms, sfr, src in GAL:
    g = dict(name=nm, z=z, re_arcsec=re_as, sigma=sig, src=src)
    for law, (da, dl) in LAWS.items():
        fA, fL = da(z) / LAWS["standard"][0](z), dl(z) / LAWS["standard"][1](z)
        re = re_as * ARC * da(z) * 1e3                       # kpc
        ms = 10 ** lms * fL ** 2; md = 10 ** lmd * fA; s = sfr * fL ** 2
        sig_sfr = s / (2 * math.pi * re ** 2)                  # Msun / yr / kpc^2
        sig_gas = (sig_sfr / 2.5e-4) ** (1 / 1.4)              # Msun / pc^2 (Kennicutt 1998)
        mg = sig_gas * 2 * math.pi * (re * 1e3) ** 2
        vir = 5 * sig ** 2 * re * 1e3 / G                       # a plain virial mass, 5 sigma^2 r_e / G, for reference
        r = dict(re=re, mstar=ms, mgas=mg, mdyn=md, bar_over_dyn=(ms + mg) / md, dyn_over_star=md / ms, size_over_local=re / shen(ms), virial=vir,
                 surf=ms / (2 * math.pi * (re * 1e3) ** 2))
        g[law] = r; rows[law].append(r)
    out["galaxies"].append(g)
print("\nSeventeen galaxies, redshift %.1f to %.1f, medians [16th-84th percentile]" % (min(g[1] for g in GAL), max(g[1] for g in GAL)))
def pct(v): v = np.array(v); return np.median(v), np.percentile(v, 16), np.percentile(v, 84)
for law in LAWS:
    R = rows[law]; a = pct([r["re"] for r in R]); b = pct([np.log10(r["mstar"]) for r in R]); c_ = pct([r["dyn_over_star"] for r in R]); d = pct([r["bar_over_dyn"] for r in R]); e = pct([r["size_over_local"] for r in R]); f = pct([r["mgas"] / r["mstar"] for r in R])
    print("  %-28s half-light radius %5.2f kpc [%.2f-%.2f]; log stellar mass %.2f; size over local galaxies of that mass x%.1f [%.1f-%.1f]" % (law, a[0], a[1], a[2], b[0], e[0], e[1], e[2]))
    print("  %-28s dynamical/stars x%.0f [%.0f-%.0f]; gas/stars %.1f; (stars+gas)/dynamical %.3f [%.3f-%.3f]; galaxies with baryons below a tenth of dynamical: %d of 17"
          % ("", c_[0], c_[1], c_[2], f[0], d[0], d[1], d[2], sum(r["bar_over_dyn"] < 0.1 for r in R)))
    out["laws"][law] = dict(re=a, logmstar=b, dyn_over_star=c_, bar_over_dyn=d, size_over_local=e, gas_over_star=f, n_below_tenth=int(sum(r["bar_over_dyn"] < 0.1 for r in R)),
                            n_above_half=int(sum(r["bar_over_dyn"] > 0.5 for r in R)))
# check of the published dynamical masses against a plain virial estimate
v = [math.log10(r["virial"] / r["mdyn"]) for r in rows["standard"]]; print("\n  plain virial mass 5 sigma^2 r_e / G over the published dynamical mass, standard distances: median x%.2f" % 10 ** np.median(v))
# line width a local galaxy of the same size and stellar mass would need if baryons were all the mass: sigma^2 = G M_bar / (5 r_e)
print("\nLine width if stars and gas were all the mass, sigma = sqrt(G M / 5 r_e), against the measured one")
for law in LAWS:
    s = [math.sqrt(G * (r["mstar"] + r["mgas"]) / (5 * r["re"] * 1e3)) for r in rows[law]]; m = [g[3] for g in GAL]
    print("  %-28s predicted median %.0f km/s, measured median %.0f; measured/predicted median x%.1f" % (law, np.median(s), np.median(m), np.median(np.array(m) / np.array(s))))
    out["laws"][law]["sigma_pred"] = float(np.median(s)); out["laws"][law]["sigma_ratio"] = float(np.median(np.array(m) / np.array(s)))
out["sigma_measured_median"] = float(np.median([g[3] for g in GAL]))
# a quantity no distance law can change: star formation rate over stellar mass (both scale as luminosity distance squared)
ssfr = np.array([g[6] / 10 ** g[5] for g in GAL]) * 1e9          # per Gyr
dbl = 1e3 / ssfr                                                   # Myr to double the stellar mass at the present rate
local = 10 ** (-0.35 * (8.53 - 10) - 9.83) * 1e9                   # Salim et al. 2007 star-forming sequence at the sample's median mass, per Gyr
print("\nStar formation rate over stellar mass: median %.0f per Gyr [%.0f-%.0f]; doubling time median %.0f Myr [%.0f-%.0f]; local star-forming galaxies of 10^8.5 Msun: %.2f per Gyr, doubling in %.1f Gyr"
      % (np.median(ssfr), np.percentile(ssfr, 16), np.percentile(ssfr, 84), np.median(dbl), np.percentile(dbl, 16), np.percentile(dbl, 84), local, 1 / local))
out["ssfr"] = dict(median_per_gyr=float(np.median(ssfr)), p16=float(np.percentile(ssfr, 16)), p84=float(np.percentile(ssfr, 84)), doubling_myr=float(np.median(dbl)), local_per_gyr=float(local))
# stellar surface density inside the half-light radius
for law in LAWS:
    v = [r["surf"] for r in rows[law]]; print("  %-28s stellar surface density, median %.2g Msun/pc^2" % (law, np.median(v))); out["laws"][law]["surf"] = float(np.median(v))
# mass-size plane for the figure
out["plane"] = {law: [[round(math.log10(r["mstar"]), 3), round(math.log10(r["re"]), 3)] for r in rows[law]] for law in LAWS}
out["shen"] = [[m, round(math.log10(shen(10 ** m)), 3)] for m in np.arange(6.5, 11.01, 0.25)]

# ---------------- the catalogue, re-derived in the browser run (astrodeep_sizes.js) ------------------------------
out["catalogue"] = dict(
 source="VizieR J/A+A/691/A240 (Merlin et al. 2024), seven fields; spectroscopic redshift, flag < 80, CLASS_STAR < 0.9", n=6316, pixel_arcsec=0.03, star_r50_pix=3.05,
 cols=["median z", "n", "median diameter 2 r50 [arcsec]", "after removing the point-source size in quadrature [arcsec]", "that, standard distances [kpc]", "that, tired light [kpc]", "share within 1.5 x point-source size"],
 rows=[[0.35,635,0.591,0.562,2.59,3.27,0.02],[0.74,1686,0.51,0.476,3.45,5.53,0.02],[1.23,1119,0.462,0.424,3.52,7.06,0.02],[1.75,763,0.425,0.383,3.24,8.0,0.08],[2.24,622,0.406,0.363,2.98,8.77,0.11],
       [2.74,317,0.343,0.29,2.3,7.99,0.2],[3.41,539,0.326,0.27,1.99,8.23,0.26],[4.42,256,0.293,0.229,1.51,8.04,0.41],[5.43,211,0.292,0.228,1.38,8.81,0.42],[6.55,91,0.274,0.204,1.13,8.56,0.52],
       [7.29,47,0.266,0.193,0.97,8.55,0.53],[8.78,22,0.23,0.138,0.63,6.73,0.73]],
 raw_median_diameter_kpc=dict(z=[0.35,0.74,1.23,1.75,2.24,2.74,3.41,4.42,5.43,6.55,7.29,8.78], standard=[2.74,3.71,3.84,3.59,3.35,2.72,2.39,1.93,1.77,1.5,1.34,1.06],
                              tired=[3.44,5.94,7.68,8.88,9.83,9.43,9.97,10.29,11.27,11.53,11.74,11.24], lerner=[4.03,7.96,11.74,15.41,18.7,19.63,22.94,27.6,32.99,37.73,41.07,45.22]),
 exponents=dict(note="diameter proportional to (1+z)^a, least squares on single galaxies; the papers: standard -0.48 (all), -0.87 (z>1), -1.00 (z>2)",
                all=dict(standard=-0.38, tired=0.67, lerner=1.26, angle=-0.46), z1=dict(standard=-0.77, tired=0.33, lerner=0.94, angle=-0.47),
                z2=dict(standard=-0.91, tired=0.22, lerner=0.85, angle=-0.42), z3=dict(standard=-0.83, tired=0.32, lerner=0.96, angle=-0.25),
                z2_deconvolved=dict(standard=-1.09, tired=0.04, lerner=0.67)),
 surface_brightness=dict(note="detection-band flux over pi r50^2; 5th, 50th, 95th percentile of log10 by redshift bin",
                         rows=[["0.1-0.5",622,-2.96,-2.27,-1.17],["0.5-1",1686,-2.83,-1.99,-1.04],["1-1.5",1119,-2.88,-1.9,-0.99],["1.5-2",763,-3.22,-1.98,-1.07],["2-3",939,-3.41,-2.25,-1.43],
                               ["3-4",539,-3.83,-2.7,-1.87],["4-5",256,-3.87,-2.97,-2.18],["5-6",211,-3.94,-3.17,-2.4],["6-7",91,-3.78,-3.13,-2.33],["7-8.5",53,-3.86,-3.25,-2.65]],
                         exponent_all=-1.58, exponent_z1=-2.64))
json.dump(out, open("line_width_jwst.json", "w"), indent=1, default=float)
