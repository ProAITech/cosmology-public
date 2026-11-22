"""Entry 187: collect result_<sample>_<map>.json into results.json and print a table.  Run from the repository root."""
import json, glob, os
HERE = os.path.dirname(os.path.abspath(__file__))
R = {}
for f in sorted(glob.glob(os.path.join(HERE, 'result_*.json'))):
    r = json.load(open(f)); R[f"{r['sample']}|{r['map']}"] = r
json.dump({'rule': 'Workbench section 187, frozen 21 November 2026 19:49 UTC (commit 7fcf22d, Workbench version 199)',
           'primary': R.get('P|smica_l5'), 'all': R}, open(os.path.join(HERE, 'results.json'), 'w'), indent=1)
for k, r in R.items():
    print(f"{k:20s} used {r['n_used']:5.0f}/{r['n_in']:5.0f}  C {r['C']:+6.2f} null {r['null_C'][0]:+5.2f}+-{r['null_C'][1]:4.2f} t_C {r['t_C']:+5.2f} p {r['p_C']:.3f} | "
          f"D {r['D']:+6.2f} null {r['null_D'][0]:+5.2f}+-{r['null_D'][1]:4.2f} t_D {r['t_D']:+5.2f} p {r['p_D']:.3f} | f {r['f_top']:.3f} "
          f"delta_top {r['delta_top']:+6.2f}+-{r['sd_delta_top']:4.2f} up95 {r['upper95_cooling_top']:5.2f} | {r['decision_minus15']}")
