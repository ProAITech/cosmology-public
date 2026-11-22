# Entry 160: Gaia DR3 radial velocities in the escape-speed tail, checked against APOGEE DR17.
# Counts below were measured in the browser (apogee_check.js); this script records them and
# computes how unusual APOGEE's fastest reliable star is under the two Leonard-Tremaine fits of Entry 151.
import json, math
counts = {
  "apogee_rows_plx_gt_0p333": 496454, "apogee_unique_gaia_ids": 433064,
  "apogee_A": {"desc": "plx > 0.333 mas, plx/e_plx > 20, e_HRV < 0.5, SNR >= 20, not (Nvis > 1 and s_HRV > 1)",
               "n": 334439, "gt300": 433, "gt350": 86, "gt400": 30, "gt450": 14, "gt475": 10, "gt500": 8},
  "apogee_B": {"desc": "as A, and Nvis >= 2 with visit scatter s_HRV < 1 km/s", "n": 229000,
               "gt300": 240, "gt350": 41, "gt400": 13, "gt450": 0, "fastest": 428},
  "single_visit_gt350": {"n": 45, "with_gaia_rv": 21, "agree_10": 21},
  "fastest_single_visit_examples": [[-719.1, "near Eta Carinae, RA 161.29 Dec -59.97", 949],
                                     [669.2, "parallax 22.9 mas (44 pc)", 690], [578.2, "parallax 39.1 mas (26 pc)", 646]],
  "gaia_tail_query": {"n": 36308, "matched_apogee": 419, "bright_grvs_lt12": 232, "faint": 187},
  "gaia_vs_apogee_B": {"bright": {"n": 124, "agree_10": 123, "gt350_gaia_rv": 26, "gt350_apogee_rv": 26},
                       "faint": {"n": 100, "agree_10": 88, "gt350_gaia_rv": 8, "gt350_apogee_rv": 8}},
  "entry151_n_gt350": 709}
def p_max_below(x, N, vesc, k, vmin=350.0):
    # Leonard-Tremaine: P(v > x) = ((vesc - x)/(vesc - vmin))^(k+1); probability all N stars are below x
    q = ((vesc - x) / (vesc - vmin)) ** (k + 1)
    return (1 - q) ** N
def median_max(N, vesc, k, vmin=350.0):
    q = 1 - 0.5 ** (1 / N)
    return vesc - (vesc - vmin) * q ** (1 / (k + 1))
N = counts["apogee_B"]["gt350"]; xmax = counts["apogee_B"]["fastest"]
fits = {"halo_482": (482, 2.0), "main_515": (515, 3.2), "mond_519": (519, 3.4)}
res = {key: {"p_max_below_428": round(p_max_below(xmax, N, v, k), 3), "median_fastest": round(median_max(N, v, k), 1)} for key, (v, k) in fits.items()}
frac_bright = counts["gaia_vs_apogee_B"]["bright"]["agree_10"] / counts["gaia_vs_apogee_B"]["bright"]["n"]
frac_faint = counts["gaia_vs_apogee_B"]["faint"]["agree_10"] / counts["gaia_vs_apogee_B"]["faint"]["n"]
added = counts["gaia_vs_apogee_B"]["faint"]["gt350_apogee_rv"]
out = {"counts": counts, "fastest_test": res, "frac_agree": {"bright": round(frac_bright, 3), "faint": round(frac_faint, 3)},
       "rescue_fraction": round(added / counts["entry151_n_gt350"], 4),
       "error_scaling": round(math.sqrt(counts["entry151_n_gt350"] / N), 1)}
json.dump(out, open("vesc_apogee.json", "w"), indent=1)
print(json.dumps({k: out[k] for k in ("fastest_test", "frac_agree", "rescue_fraction", "error_scaling")}))
