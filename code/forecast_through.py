"""forecast_through.py - Workbench Entry 57
How big is the rotation-measure signal of a current running through a galaxy along its spin axis, and what
would it take to see it?

1. The signal. Uniform gas (the proposal's assumption) against halo gas that thins outward.
2. The same-sense ring field of the Milky Way's own disc, read as a current through the plane.
3. Forecast: error on the current per galaxy and for stacks, for rotation-measure grids of different density.

Results of the catalogue runs (through_current_stack.js, in a browser tab) are typed in at the end and written
to through_results.json with everything computed here.
"""
import json, math
import numpy as np

MU0 = 4e-7 * math.pi
PC = 3.0857e16                      # m
K_UGPC = MU0 * 1e18 / (2 * math.pi) * 1e10 / PC   # mu0 I / 2 pi for 1e18 A, in microgauss parsec
KRM = 0.812 * K_UGPC                # rad/m2 per cm^-3 per unit of the dimensionless line integral
N0_THIN = 1.35e-2                   # Miller & Bregman 2015, as fitted with solar metallicity
N0_BASE = 1.1e-4 * 48.2 ** 1.5      # same slope, pinned to the density at the Large Magellanic Cloud (Salem+ 2015)
N0_Z03 = 1.35e-2 / 0.3              # Miller & Bregman with their preferred metallicity, 0.3 solar

def F(theta_deg, i_deg, n=20000):
    """Dimensionless line-of-sight integral for a line current along the spin axis, gas density ~ r^-1.5."""
    th, i = math.radians(theta_deg), math.radians(i_deg)
    ph = -math.pi / 2 + (np.arange(n) + 0.5) * math.pi / n
    cp, sp = np.cos(ph), np.sin(ph)
    q = math.sin(th) * math.cos(i) * cp - sp * math.sin(i)
    return float(np.sum(cp ** 1.5 * math.cos(th) * math.sin(i) / (math.cos(th) ** 2 * cp ** 2 + q ** 2)) * math.pi / n)

def F2mean(i_deg):
    th = (np.arange(90) + 0.5)
    return float(np.mean([F(t, i_deg, 40000 if t > 80 else 6000) ** 2 for t in th]))

out = {}

# ---- 1. the signal -------------------------------------------------------------------------------------
print("mu0 I / 2 pi for 1e18 A: %.4g microgauss pc ; ring field at 10 kpc: %.2f microgauss" % (K_UGPC, K_UGPC / 1e4))
cap = 0.812 * 1e-5 * K_UGPC * math.pi
print("uniform gas, 1e-5 cm^-3, endless path: RM = 0.812 n mu0 I / 2 = %.3f rad/m2 per 1e18 A, at any distance" % cap)
rows = []
for x, Bprop, RMprop in [(10, 6.5, 53), (30, 2.2, 18), (100, 0.65, 5.3), (300, 0.22, 1.8)]:
    L = 1000.0
    rm_uni = 0.812 * 1e-5 * K_UGPC * 2 * math.atan(L / (2 * x))
    rm_thin = KRM * N0_THIN * x ** -1.5 * F(0, 90)
    rm_base = KRM * N0_BASE * x ** -1.5 * F(0, 90)
    rows.append(dict(x=x, proposal=RMprop, uniform=round(rm_uni, 2), over=round(RMprop / rm_uni, 1),
                     halo_thin_major=round(rm_thin, 2), halo_base_major=round(rm_base, 2),
                     halo_base_minor=round(KRM * N0_BASE * x ** -1.5 * math.pi, 2)))
    print("x %4d kpc: proposal %5.1f | uniform gas, 1 Mpc path %.2f (proposal is %.1fx too big) | halo gas, edge-on: major axis %.2f (thin) %.2f (LMC-pinned), beside the minor axis %.2f"
          % (x, RMprop, rm_uni, RMprop / rm_uni, rm_thin, rm_base, KRM * N0_BASE * x ** -1.5 * math.pi))
out["signal"] = dict(cap_uniform=round(cap, 3), rows=rows, N0_thin=N0_THIN, N0_base=round(N0_BASE, 4), N0_Z03=round(N0_Z03, 4),
                     n_at_48kpc=dict(thin=N0_THIN * 48.2 ** -1.5, base=N0_BASE * 48.2 ** -1.5, z03=N0_Z03 * 48.2 ** -1.5))
print("gas density at 48.2 kpc: thin %.2e, LMC-pinned %.2e, MB15 at 0.3 solar %.2e cm^-3" % (N0_THIN * 48.2 ** -1.5, 1.1e-4, N0_Z03 * 48.2 ** -1.5))
print("F: edge-on major axis %.4f (analytic 1.7480), beside the minor axis -> pi sin^1.5 i; i=77: major %.3f" % (F(0, 90), F(0, 77)))
prof = []
for b in [10, 15, 20, 30, 40, 60, 80, 100, 150, 200, 300]:
    prof.append([b, round(KRM * N0_THIN * b ** -1.5 * F(0, 90), 3), round(KRM * N0_BASE * b ** -1.5 * F(0, 90), 3),
                 round(KRM * N0_BASE * b ** -1.5 * math.pi, 3), round(0.812 * 1e-5 * K_UGPC * 2 * math.atan(1000 / (2 * b)), 3)])
out["profile_cols"] = ["b_kpc", "halo_thin_major", "halo_base_major", "halo_base_minor", "uniform_1e-5_1Mpc"]
out["profile"] = prof
# the proposal's decision rule: "no antisymmetric term above about 2 rad/m2 at 30 kpc bounds I below roughly 1e16 A"
r30 = [KRM * N0_THIN * 30 ** -1.5 * F(0, 90), KRM * N0_BASE * 30 ** -1.5 * math.pi]
print("2 rad/m2 at 30 kpc corresponds to %.2g - %.2g x 1e18 A in halo gas, %.2g x 1e18 A in uniform 1e-5 gas"
      % (2 / r30[1], 2 / r30[0], 2 / (0.812 * 1e-5 * K_UGPC * 2 * math.atan(1000 / 60))))
out["rule_2radm2_at_30kpc_A"] = dict(halo_low=2 / r30[1] * 1e18, halo_high=2 / r30[0] * 1e18, uniform=2 / (0.812 * 1e-5 * K_UGPC * 2 * math.atan(1000 / 60)) * 1e18)

# ---- 2. the Milky Way's disc field as a current through the plane ---------------------------------------
kpc = 1e3 * PC
for B, R in [(2.0, 8.2), (1.5, 8.2)]:
    I = 2 * math.pi * R * kpc * B * 1e-10 / MU0
    print("same-sense disc field %.1f microgauss at %.1f kpc: %.2g A crossing the plane inside that circle" % (B, R, I))
out["mw_disc_current_A"] = [2 * math.pi * 8.2 * kpc * B * 1e-10 / MU0 for B in (2.0, 1.5)]
out["A_per_microgauss_kpc"] = 2 * math.pi * kpc * 1e-10 / MU0

# ---- 3. forecast ---------------------------------------------------------------------------------------
def sumT2(D_mpc, i_deg, rho_deg, N0, bin_=30.0, bout=200.0):
    """Sum of template^2 over the sight lines a grid of rho_deg sources per square degree puts in the ring."""
    kpc_per_deg = D_mpc * 1e3 * math.pi / 180
    return rho_deg / kpc_per_deg ** 2 * (KRM * N0) ** 2 * 2 * math.pi * F2mean(i_deg) * (1 / bin_ - 1 / bout)

f2 = {i: F2mean(i) for i in (55, 60, 65, 70, 77, 80, 85, 90)}
print("mean of F^2 over angle:", {k: round(v, 2) for k, v in f2.items()})
def f2of(i):
    ks = sorted(f2); return float(np.interp(i, ks, [f2[k] for k in ks]))
def sig(D, i, rho, s, N0, bin_=30.0, bout=200.0):
    kpd = D * 1e3 * math.pi / 180
    return s / math.sqrt(rho / kpd ** 2 * (KRM * N0) ** 2 * 2 * math.pi * f2of(i) * (1 / bin_ - 1 / bout))

surveys = {"NVSS (1.1 per deg2, 14 rad/m2)": (1.1, 14.0), "LoTSS (0.43, 7)": (0.43, 7.0),
           "POSSUM (40, 8)": (40.0, 8.0), "SKA-like (100, 6)": (100.0, 6.0)}
fc = {}
for name, (rho, s) in surveys.items():
    fc[name] = {str(D): round(sig(D, 80, rho, s, N0_BASE), 3) for D in (3.5, 10, 20, 30)}
    print("%-32s error on I per galaxy [1e18 A], LMC-pinned halo, 30-200 kpc: " % name, fc[name])
out["per_galaxy_sigma"] = fc

# the 118 large spirals of the sample (blue magnitude below -19.5): distance [Mpc], inclination, south of the equator?
D = [float(x) for x in "19 48.8 27.8 32.6 47.5 17.8 17.1 16.9 16.9 17.5 17.4 25.5 13 14.3 13.6 20 17.2 33.9 22.2 20.4 33.7 12 21.9 24.4 20.4 22 22.4 26.7 23.3 27.4 26.7 43.2 14.1 65.3 28.4 12.8 24.6 17 23.5 26.4 27.8 21.9 17 17 17 25 16.8 16.8 6.8 17 16.8 16.8 9.7 16.8 20 16.8 25.9 15.1 38.9 21.5 21.1 44.3 21.3 18.7 29 7.2 27.3 36.6 35.1 33.5 24 37.8 31 52.6 32.6 26.4 35.8 34.5 28.8 42.2 29.4 31.8 32.6 26.7 30.6 28.5 31.6 39.8 14.9 39.2 17.5 31.2 25.4 36.7 30.2 61.3 25.5 28.8 60.9 59.1 25.3 48.8 60.9 36.6 38.8 52 34.1 20.5 45 18.3 14.3 20.1 20.9 30.3 33.4 17.6 28.9 26.1".split()]
inc = [float(x) for x in "82 65 69 56 67 67 58 59 59 56 66 82 80 57 62 71 58 59 61 74 73 67 59 83 90 63 73 61 78 56 58 78 81 57 64 60 90 63 62 63 69 56 57 62 70 57 78 85 70 65 71 59 90 65 69 57 63 64 59 73 79 61 64 64 74 57 57 57 63 70 90 84 58 59 66 74 90 64 77 57 90 60 66 82 81 61 68 71 90 59 69 60 67 66 67 63 62 78 80 90 78 62 68 79 62 64 82 63 64 65 73 76 75 65 73 68 69 69".split()]
south = [c == "1" for c in "110011111111111100001000000000010010100010000000000000101101" + "1100101111100000000101101010000010010111111111111101100111"]
assert len(D) == len(inc) == len(south) == 118
print("large spirals: %d, of them south of the equator: %d, median distance %.1f Mpc" % (len(D), sum(south), float(np.median(D))))

def stack(sel, rho, s, N0, bin_=30.0, bout=200.0):
    sg = np.array([sig(d, i, rho, s, N0, bin_, bout) for d, i, k in zip(D, inc, sel) if k])
    signed = 1 / math.sqrt(np.sum(1 / sg ** 2))                 # spin sense known, current tied to it
    s_I2 = 1 / math.sqrt(np.sum(1 / (2 * sg ** 4)))             # sense unknown: only I^2 adds
    return dict(n=int(len(sg)), median_per_galaxy=float(np.median(sg)), signed_sigma=signed, signed_3sigma=3 * signed,
                unsigned_95=math.sqrt(1.645 * s_I2))
st = {}
allsel = [True] * 118
for tag, sel, rho, s in [("POSSUM, 56 southern large spirals", south, 40, 8), ("SKA-like grid, all 118", allsel, 100, 6),
                         ("NVSS, all 118", allsel, 1.1, 14), ("LoTSS density, all 118", allsel, 0.43, 7)]:
    for dens, N0 in [("LMC-pinned", N0_BASE), ("thin", N0_THIN)]:
        r = stack(sel, rho, s, N0); st[tag + " | " + dens] = {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}
        print("%-34s %-10s n=%3d  per galaxy (median) %.2f  signed stack: sigma %.3f, 3 sigma %.3f  | sense unknown, 95%%: %.2f   [1e18 A]"
              % (tag, dens, r["n"], r["median_per_galaxy"], r["signed_sigma"], r["signed_3sigma"], r["unsigned_95"]))
r15 = stack(south, 40, 8, N0_BASE, 15.0, 200.0)
print("POSSUM, southern large spirals, sight lines from 15 kpc: signed 3 sigma %.3f, unsigned 95%% %.2f" % (r15["signed_3sigma"], r15["unsigned_95"]))
st["POSSUM, 56 southern, from 15 kpc | LMC-pinned"] = {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r15.items()}
out["stacks"] = st
# how many galaxies at 10 Mpc for 1e17 A (3 sigma, signed) and 1e16 A
s10 = sig(10, 80, 40, 8, N0_BASE)
print("POSSUM-like grid, one edge-on spiral at 10 Mpc: +-%.2f x 1e18 A. For 3 sigma at 1e17 A: %d such galaxies; at 1e16 A: %d"
      % (s10, math.ceil((3 * s10 / 0.1) ** 2), math.ceil((3 * s10 / 0.01) ** 2)))
out["n_needed_10Mpc"] = dict(per_galaxy=s10, for_1e17=math.ceil((3 * s10 / 0.1) ** 2), for_1e16=math.ceil((3 * s10 / 0.01) ** 2))
# the proposal's five targets with a POSSUM-like grid (none of the five is in POSSUM's sky: all but one lie well north)
targets = {"NGC 891": (9.6, 89, 42.3), "NGC 4631": (7.4, 86, 32.5), "NGC 5775": (28.9, 86, 3.5), "NGC 4565": (11.9, 87, 26.0), "NGC 3628": (11.3, 87, 13.6)}
tg = {k: round(sig(d, i, 40, 8, N0_BASE), 3) for k, (d, i, dec) in targets.items()}
comb = 1 / math.sqrt(sum(1 / v ** 2 for v in tg.values()))
print("five named targets, POSSUM-like grid:", tg, " together +-%.3f (signed)" % comb)
out["five_targets"] = dict(per_galaxy=tg, combined_signed=round(comb, 3), declinations={k: v[2] for k, v in targets.items()})

# single nearby galaxies
n253 = sig(3.5, 78, 40, 8, N0_BASE)
print("NGC 253 alone (3.5 Mpc, inclination 78) in POSSUM: +-%.3f, 3 sigma %.2f x 1e18 A (LMC-pinned); thin: 3 sigma %.2f" % (n253, 3 * n253, 3 * sig(3.5, 78, 40, 8, N0_THIN)))
# M31: scale the error actually found with LoTSS (62 sight lines in the 30-150 kpc ring, +-1.0 for the thin halo) to a dense grid
ring_deg2 = math.pi * ((150 / 13.44) ** 2 - (30 / 13.44) ** 2)
n_dense = 40 * ring_deg2
m31_thin = 1.0 * math.sqrt(62 / n_dense)
print("M31 ring 30-150 kpc = %.0f deg2; at 40 per deg2: %d sources; error %.3f (thin), %.3f (LMC-pinned); 3 sigma %.2f and %.3f x 1e18 A"
      % (ring_deg2, n_dense, m31_thin, m31_thin * N0_THIN / N0_BASE, 3 * m31_thin, 3 * m31_thin * N0_THIN / N0_BASE))
out["single"] = dict(ngc253_sigma=n253, m31_ring_deg2=ring_deg2, m31_sources=n_dense, m31_sigma_thin=m31_thin, m31_sigma_base=m31_thin * N0_THIN / N0_BASE)
# the published excess near minor axes (Heesen+ 2023: 3.7 rad/m2 within 100 kpc), read as a through-current's step
for nm, N0 in [("LMC-pinned", N0_BASE), ("thin", N0_THIN)]:
    lo, hi = 3.7 / (KRM * N0 * 30 ** -1.5 * math.pi), 3.7 / (KRM * N0 * 100 ** -1.5 * math.pi)
    print("3.7 rad/m2 beside the minor axis at 30-100 kpc = %.2g to %.2g x 1e18 A (%s)" % (lo, hi, nm))
out["heesen_equiv_base"] = [3.7 / (KRM * N0_BASE * 30 ** -1.5 * math.pi), 3.7 / (KRM * N0_BASE * 100 ** -1.5 * math.pi)]

# ---- catalogue runs (through_current_stack.js), currents in 1e18 A for the thin halo ----------------------
runs = {
  "sample": dict(galaxies=228, median_distance_Mpc=17.1, nvss_residuals=28127, lotss=2461),
  "stack_all": dict(
     nvss=dict(sight_lines=1311, galaxies=142, pairs=30621, I2=0.289, null_mean=0.862, null_sd=0.73, p=0.77, limit95=1.1, recovery=0.85),
     lotss=dict(sight_lines=241, galaxies=42, pairs=2344, I2=0.585, null_mean=1.055, null_sd=0.69, p=0.77, limit95=1.35, recovery=1.0),
     weight_share_top8_nvss="NGC 2976 0.34, IC 2574 0.20, NGC 4244 0.07, NGC 247 0.07, NGC 4236 0.06, NGC 5253 0.05, NGC 2366 0.04, NGC 253 0.03: 85% in eight galaxies, seven of them dwarfs",
     weight_share_top5_lotss="NGC 4244 0.47, NGC 2976 0.22, NGC 4605 0.12, IC 2574 0.12, NGC 4144 0.06: 99% in five dwarfs"),
  "stack_large_spirals": dict(n_galaxies=118,
     nvss=dict(sight_lines_by_ring_15_30_60_100_200=[2, 8, 17, 127], I2=6.43, null_median=9.0, limit95=">5"),
     lotss=dict(sight_lines_by_ring_15_30_60_100_200=[0, 1, 3, 27], I2=3.65, null_median=3.1, limit95=4.7)),
  "M31": dict(distance_Mpc=0.77, inclination=77, receding_PA=38,
     nvss_median_filter=dict(sight_lines=472, ring_kpc=[30, 150], rms=22.0, A=-0.71, recovery=0.76, A_corrected=-0.93,
                             null_sd=1.0, null_sd_corrected=1.3, null="180 placements at the same Galactic latitude"),
     lotss_raw_quadratic=dict(sight_lines=62, rms=13.2, A=-0.50, formal=0.85, null_sd_elsewhere=1.0, rotation_sd=1.7, recovery=1.0),
     nvss_raw_quadratic=dict(A=-2.35, null_sd=2.65), nvss_raw_cubic=dict(A=-0.85, null_sd=1.57),
     bins_cols=["x_from", "x_to", "n", "mean_residual", "error", "mean_filtered_template_per_1e18A_thin"],
     bins=[[-150,-100,63,-0.6,2.22,0.008],[-100,-60,72,-6.54,2.72,0.129],[-60,-30,63,-2.71,2.91,0.672],[-30,0,52,-1.35,2.95,-1.933],
           [0,30,63,2.4,2.86,1.802],[30,60,45,-3.16,3.49,-0.304],[60,100,60,1.59,2.77,-0.099],[100,150,54,0.3,2.91,-0.033]]),
}
scale = N0_THIN / N0_BASE
print("thin -> LMC-pinned halo: multiply currents by %.3f" % scale)
print("stack of all 228: below %.2f (NVSS), %.2f (LoTSS) x 1e18 A for the thin halo; %.2f, %.2f for the LMC-pinned one"
      % (1.1, 1.35, 1.1 * scale, 1.35 * scale))
print("M31: NVSS %.2f +- %.2f, LoTSS %.2f +- %.1f to %.1f (thin); x%.3f: %.2f +- %.2f and %.2f +- %.2f to %.2f"
      % (-0.93, 1.3, -0.5, 1.0, 1.7, scale, -0.93 * scale, 1.3 * scale, -0.5 * scale, 1.0 * scale, 1.7 * scale))
out["runs"] = runs
out["thin_to_base"] = scale
json.dump(out, open("through_results.json", "w"), indent=1)
