# Entry 162: is the rotation-measure flip of Entry 161 pinned to the clouds? Fits ran in the browser (cloud_pin.js).
import json
prof = {"ori": "5.3 5.4 4.7 2.0 3.7 3.7 3.7 1.5 1.5 0.7 0.3 0.1 0.1 0.1 0.7 1.3 0.2 8.7 9.3 9.2 9.2 8.7 8.7 10.8 10.8 13.1 13.1 13.1 9.8 6.5 6.5 7.6 6.6 6.5 10.3 6.0 3.6 6.9 7.7 6.4 3.6 2.7 2.2 2.2 1.4 1.4 0.5 0.3 0.1",
        "per": "0.8 0.9 0.1 0.1 0.1 0.2 0.2 0.2 0.0 0.0 0.1 0.0 0.0 0.0 0.4 4.3 9.5 13.1 13.3 24.9 28.2 26.7 27.8 25.1 14.0 8.7 9.4 5.9 5.9 1.9 0.7 0.1 0.1 0.1 0.0 0.1 0.0 0.0 0.3 0.2 0.2 0.4 0.9 0.7 1.0 1.5 3.4 6.0 7.7"}
d0 = [round(-6 + 0.25 * i, 2) for i in range(49)]
res = {"d0": d0, "prof": {k: [float(x) for x in v.split()] for k, v in prof.items()},
 "best": {"ori": [0.25, 13.05, 21.5], "per": [-1.0, 28.18, 37.3]},
 "null_best_within_1deg": {"ori": [6, 51], "per": [7, 66]},
 "null_stronger_anywhere": {"ori": 0.647, "per": 0.333},
 "extent_w": [1, 2, 3, 5, 10], "extent": {"ori": [0.1, 3.88, 3.27, 7.5, 10.76], "per": [0.07, 4.92, 11.07, 11.57, 13.99]},
 "twozone_lt2_vs_2to10": {"ori": [17.7, 22.0], "per": [22.6, 35.3]},
 "orientation_best_rot": {"ori": [-30, 22.52, 13.05], "per": [40, 43.0, 28.18]},
 "sky_inflation_z": {"ori": 2.47, "per": 3.09}}
p1 = 6 / 51; p2 = 7 / 66; res["p_both_on_spine"] = round(p1 * p2, 4)
res["ext_minus_conf2"] = {k: round(res["extent"][k][4] - res["extent"][k][1], 2) for k in ("ori", "per")}
res["ext_minus_conf2_sky"] = {k: round(res["ext_minus_conf2"][k] / res["sky_inflation_z"][k] ** 2, 2) for k in ("ori", "per")}
json.dump(res, open("cloud_pin.json", "w"), indent=1)
print(res["p_both_on_spine"], res["ext_minus_conf2"], res["ext_minus_conf2_sky"])
