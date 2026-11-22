"""Workbench Entry 17: what ejection of quasars from nearby galaxies predicts for Gaia, and what the
measured outward drift allows. Measured numbers come from arp_proper_motion.js (run in a browser tab on
Milliquas v8 x Gaia DR3 and on SDSS DR16Q x Gaia DR2)."""
import json, numpy as np
rng = np.random.default_rng(17)
KPC_PER_KMS_YR = 1.0227e-9          # 1 km/s for 1 yr = 1.0227e-9 kpc
KV = 4.740                          # km/s per (mas/yr x kpc)
C = 299792.458

# ---------- measured (5-sigma-clipped primary sample; error = largest of formal, jackknife, shifted-sky scatter)
meas = {
  "bins_kpc": [[0,50],[50,100],[100,200],[200,500],[500,1000]],
  "pairs":    [1621, 4713, 18264, 132138, 478984],
  "galaxies": [421, 542, 624, 761, 834],
  "mu_r":     [0.24, -0.29, 2.21, -0.38, -0.31],      # micro-arcsec / yr, outward
  "mu_r_err": [6.95, 5.23, 2.26, 0.87, 0.59],
  "mu_t":     [9.05, -6.61, -0.99, 0.97, 0.37],       # sideways (null)
  "mu_t_err": [7.84, 5.76, 2.28, 1.09, 1.34],
  "v":        [40, -183, -5, 9, -27],                 # km/s outward if at the galaxy's distance
  "v_err":    [353, 221, 141, 23, 30],
  "invt":     [-2.338, -2.271, -0.138, 0.027, -0.030],# 1e-9 / yr
  "invt_err": [7.892, 2.847, 0.878, 0.063, 0.043],
  "in200": {"pairs":24598,"galaxies":648,"mu_r":1.59,"mu_r_err":1.99,"mu_t":-1.42,"mu_t_err":2.10,"v":-39,"v_err":105,"invt":-0.285,"invt_err":0.880},
  "in500": {"pairs":156736,"galaxies":777,"mu_r":-0.07,"mu_r_err":0.81,"mu_t":0.60,"mu_t_err":1.01,"v":3,"v_err":23,"invt":0.021,"invt_err":0.064},
  "in50":  {"pairs":1621,"galaxies":421,"v":40,"v_err":353},
}
variants = {  # within 200 kpc: [pairs, galaxies, mu_r, err, v, err]
  "galaxies nearer than 11 Mpc": [8523, 39, 2.92, 3.67, -46, 112],
  "galaxies at 11-40 Mpc":       [16075, 609, 0.68, 2.37, 42, 217],
  "without NGC 55":              [24370, 647, 2.32, 1.92, 53, 89],
  "brightest galaxies (M_B < -20.5)": [5453, 148, 3.54, 3.87, 174, 175],
  "quasars with z > 1.5":        [14628, 612, -0.53, 2.61, -60, 100],
  "SDSS spectra only":           [20719, 510, 1.84, 2.12, 102, 91],
  "SDSS DR16Q with Gaia DR2 motions": [19705, 510, -3.28, 4.36, -14, 185],
}
out = {"measured": meas, "variants": variants}
up = lambda m, e: m + 2*e           # ~97.7% one-sided; quoted as a 95% limit, conservatively

# ---------- 1. ejection at constant speed, random direction, age uniform 0..T (Bell's model as tested by Tang & Zhang)
def model(vej, T=1e8, R=(50, 200, 500), n=4_000_000):
    t = rng.uniform(0, T, n); ci = rng.uniform(-1, 1, n); si = np.sqrt(1-ci*ci)
    rp = vej*t*si*KPC_PER_KMS_YR; res = {}
    for r in R:
        m = rp < r
        res[r] = {"share_inside": float(m.mean()), "mean_vperp": float((vej*si[m]).mean()), "mean_invt_per_Gyr": float((1e9/t[m]).mean())}
    return res
out["model"] = {}
print("Ejection model: mean sideways speed of ejected quasars seen inside R (km/s)")
for vej in (11000, 40000, 80000):
    r = model(vej); out["model"][str(vej)] = r
    print(f"  v_ej {vej:6d} km/s:", {k: round(v['mean_vperp']) for k, v in r.items()}, " share of all ejected inside R:", {k: round(v['share_inside'],3) for k, v in r.items()})

# ---------- 2. largest share f of the quasars near these galaxies that can have been ejected
print("\n95%% upper limits on outward speed: <50 kpc %.0f, <200 kpc %.0f, <500 kpc %.0f km/s" % (up(40,353), up(-39,105), up(3,23)))
lim = {50: up(40,353), 200: up(-39,105), 500: up(3,23)}
out["limits_v"] = lim; out["f_max"] = {}
for vej in (11000, 40000, 80000):
    fm = {r: lim[r]/out["model"][str(vej)][r]["mean_vperp"] for r in (50, 200, 500)}
    out["f_max"][str(vej)] = fm
    print(f"  v_ej {vej:6d}: ejected share below", {k: f"{100*v:.1f}%" for k, v in fm.items()})

# ---------- 3. travel time, independent of the galaxy's distance: mu = theta / t
for k in ("in200", "in500"):
    lim_i = up(meas[k]["invt"], meas[k]["invt_err"])*1e-9
    out[k+"_travel_time_Gyr_if_all_ejected"] = 1/lim_i/1e9
    print(f"{k}: f x <1/t> < {lim_i:.2e} /yr  ->  travel time > {1/lim_i/1e9:.2f} Gyr if all were ejected; share < {100*lim_i*1e8:.1f}% if every trip took 1e8 yr")
# speed implied by that travel time at 200 kpc
print("  200 kpc in %.2f Gyr = %.0f km/s" % (out["in200_travel_time_Gyr_if_all_ejected"], 200/(out["in200_travel_time_Gyr_if_all_ejected"]*1e9)/KPC_PER_KMS_YR))

# ---------- 4. the Doppler version (redshift = recession speed): motion across the sky
for D in (3.6, 10, 25):
    mu = 0.8*C/(KV*D*1000)
    print(f"Doppler version, v = 0.8c sideways at {D} Mpc: {mu:.1f} mas/yr (typical Gaia error 0.51 mas/yr -> {mu/0.51:.0f} sigma)")
out["doppler_mu_mas_yr"] = {str(D): 0.8*C/(KV*D*1000) for D in (3.6, 10, 25)}
print("Arp-type speeds: 11,000 km/s at 10 Mpc = %.0f micro-arcsec/yr; 25,000 km/s at 7.6 Mpc = %.0f" % (11000/(KV*1e4)*1000, 25000/(KV*7600)*1000))
# fast movers
print("Quasars moving at >5 sigma: all sky 157 of 446,841 = %.3f%%; within 200 kpc of the galaxies 5 of 20,922 = %.3f%% (expected at the all-sky rate %.1f)" % (157/446841*100, 5/20922*100, 20922*157/446841))
print("  SDSS DR16Q + Gaia DR2: all sky 186 of 354,426 = %.3f%%; near galaxies 9 of 17,470 = %.3f%%" % (186/354426*100, 9/17470*100))

# ---------- 5. Arp's own showcase: the pair across NGC 4258 (z = 0.398 and 0.653/0.654)
z1, z2 = 0.398, 0.654
ratio = (1+z2)/(1+z1); beta = (ratio-1)/(ratio+1)
print(f"\nNGC 4258 pair: redshift ratio {ratio:.4f} -> +-{beta:.4f} c = +-{beta*C:.0f} km/s along the line of sight if the difference is Doppler")
D = 7.6e3; mu, e = -0.102, 0.087     # mas/yr, error-weighted mean outward motion of the two (Gaia DR3)
v, ev = KV*D*mu, KV*D*e
print(f"  measured outward motion {mu*1000:.0f} +- {e*1000:.0f} micro-arcsec/yr = {v:.0f} +- {ev:.0f} km/s at 7.6 Mpc; 95% limit {up(v,ev):.0f} km/s")
ang = np.degrees(np.arctan(up(v,ev)/(beta*C)))
print(f"  ejection axis within {ang:.1f} deg of the line of sight; chance for a random direction {100*(1-np.cos(np.radians(ang))):.2f}%")
print(f"  for a typical direction (60 deg from the line of sight) the same speed would show as {np.tan(np.radians(60))*beta*C/(KV*D)*1000:.0f} micro-arcsec/yr of sideways drift each")
out["ngc4258"] = {"beta": beta, "v_los": beta*C, "mu": mu, "mu_err": e, "v": v, "v_err": ev, "v_lim": up(v,ev), "angle_deg": ang, "chance": 1-np.cos(np.radians(ang))}

# ---------- 6. M82: Burbidge, Burbidge & Arp (2003) give the expected drift themselves
Dm = 3.53e3
print("\nM82 at 3.53 Mpc: 0.1c sideways = %.0f micro-arcsec/yr (they quote about 2,000 at 3.25 Mpc: %.0f)" % (0.1*C/(KV*Dm)*1000, 0.1*C/(KV*3250)*1000))
mu12 = 12000/(KV*Dm)*1000
print("  their average line-of-sight speed, 12,000 km/s, as a sideways drift: %.0f; for random directions (pi/2 times larger): %.0f" % (mu12, mu12*np.pi/2))
m, e = 38, 146   # 7 confirmed quasars within 15 arcmin with Gaia motions, micro-arcsec/yr
print("  measured, 7 quasars within 15': %+d +- %d  -> %.1f sigma below %.0f, %.1f sigma below 0.1c; 95%% limit %.0f = %.0f km/s" % (m, e, (mu12-m)/e, mu12, (0.1*C/(KV*Dm)*1000-m)/e, up(m,e), up(m,e)/1000*KV*Dm))
out["m82"] = {"mu_0p1c": 0.1*C/(KV*Dm)*1000, "mu_12000": mu12, "mu": m, "err": e, "v_lim": up(m,e)/1000*KV*Dm}

# showcase fields: [n quasars, outward motion mas/yr, error, distance Mpc]
cases = {"M82, quasars within 10' (5 of 12 measured)": [5,-0.029,0.179,3.53], "M82, quasars within 15' (7 of 20 measured)": [7,0.038,0.146,3.53], "M82, all within 1 deg":[76,0.033,0.022,3.53],
         "NGC 4258 pair (z 0.398, 0.654)": [2,-0.102,0.087,7.6], "NGC 4258, all within 1 deg":[180,-0.007,0.018,7.6],
         "NGC 3628, two X-ray quasars (z 0.98, 0.99)":[2,0.148,0.252,11.0], "NGC 3628, all within 30'":[19,0.074,0.104,11.0],
         "NGC 3516, the five quasars":[5,0.040,0.073,38.0], "Arp 220, all within 45'":[42,0.045,0.038,78.0]}
out["cases"] = {}
print("\nShowcase fields (outward speed if at the galaxy's distance):")
for k,(n,m,e,Dq) in cases.items():
    v, ev = KV*Dq*1000*m, KV*Dq*1000*e
    out["cases"][k] = {"n":n,"mu":m,"err":e,"D":Dq,"v":v,"v_err":ev}
    print(f"  {k:44s} n={n:3d}  {m*1000:6.0f} +- {e*1000:4.0f} uas/yr  {v:7.0f} +- {ev:6.0f} km/s  ({m/e:+.1f} sigma)")
json.dump(out, open("arp_results.json","w"), indent=1)
