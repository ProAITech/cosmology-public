"""Workbench Entry 69: does the filament gas density decide the filament-current test of Entry 68?
(1) The field it is compared with was itself inferred from rotation measures with an assumed density, so both scale as 1/n_e:
    compare in rotation-measure units instead.  (2) What the measurements say the density is."""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
S17 = 0.1652            # rest-frame on-axis step, rad/m^2, per 1e17 A at n_e = 1e-5 cm^-3 (Entry 13)
# current estimates toward the denser end, 1e17 A at n_e = 1e-5 (Entry 68, combine_filaments.py)
est = {"LOFAR, nothing lost": (0.284, 0.400), "LOFAR, half lost": (0.568, 0.800), "SPICE-RACS DR1": (-1.454, 1.294)}
def comb(*xs):
    w = [1 / s**2 for _, s in xs]; m = sum(a * x for a, (x, _) in zip(w, xs)) / sum(w); return m, sum(w) ** -0.5
est["both, LOFAR half lost"] = comb(est["LOFAR, half lost"], est["SPICE-RACS DR1"])
est["both, LOFAR nothing lost"] = comb(est["LOFAR, nothing lost"], est["SPICE-RACS DR1"])
RRM_F = (0.71, 0.07)    # Carretti et al. 2022: rest-frame RRM per filament, from which B_f = 32 nG at n_e = 1e-5
out = {"step": {}}
print("1. The same test in rotation-measure units (rest frame), which no assumed density enters")
for k, (m, s) in est.items():
    sm, ss = m * S17, s * S17; lim = max(sm, 0) + 2 * ss
    print(f"   {k:26s} coherent step {sm:+.3f} +- {ss:.3f} rad/m2; 95% limit {lim:.2f} = {100*lim/RRM_F[0]:.0f}% of the 0.71 rad/m2 per filament")
    out["step"][k] = dict(step=sm, sd=ss, lim=lim, frac=lim / RRM_F[0])
# 2. mean electron density, fully ionised, Omega_b h^2 = 0.0224, Y = 0.245
rho_b0 = 0.0224 * 1.8785e-29; ne0 = rho_b0 * (1 - 0.245 + 0.245 / 2) / 1.67262e-24
z = 0.035; ne = ne0 * (1 + z) ** 3
print(f"\n2. Mean electron density today {ne0:.2e} cm^-3; at z = {z} (half the crossings are nearer) {ne:.2e}")
meas = [("SZ, Tanimura+2020, beta model, centre", 19.0, 27.3, 12.1, "Planck y; 24,544 DisPerSE filaments 30-100 Mpc long, z 0.2-0.6; core radius 1.5 Mpc"),
        ("SZ, Tanimura+2020, mean within 1 Mpc", 19.0 * 0.828 + 0.828 - 1, None, None, "beta = 2/3 profile averaged over 1 Mpc with r_c = 1.5 Mpc"),
        ("SZ, Tanimura+2020, constant density to 5 Mpc", 6.3, 0.9, 0.8, ""),
        ("X-ray, Zhang+2024, mean of emitting gas", 10 ** 1.88 - 1, None, None, "eRASS1; 7,817 DisPerSE filaments 20-100 Mpc, median z 0.47; emission-weighted"),
        ("X-ray, Zhang+2024, centre", 10 ** 2.28 - 1, None, None, "")]
out["density"] = {}
for name, d, up, dn, note in meas:
    n = (1 + d) * ne; rng = "" if up is None else f" ({(1+d-dn)*ne:.1e} to {(1+d+up)*ne:.1e})"
    print(f"   {name:46s} overdensity {d:6.1f} -> n_e {n:.1e} cm^-3{rng}   {note}")
    out["density"][name] = dict(delta=d, ne=n)
print("   Faraday rotation weighs n_e itself; X-ray weighs n_e^2, so clumped gas reads high in X-rays and the SZ values are the relevant ones.")
m, s = est["both, LOFAR half lost"]
for n in (4.1e-6, 1e-5, 1.9e-5):
    print(f"   absolute 95% limit on the current at n_e = {n:.1e}: {(max(m,0)+2*s)*1e-5/n:.1f} x 1e17 A")
json.dump(out, open(os.path.join(HERE, "density_free.json"), "w"), indent=1)
