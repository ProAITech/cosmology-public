"""Entry 179, after the result: does Entry 176's mass scaling repeat on the other half of the sky?

Entry 176 found the excess inside R500 rising with cluster mass as M500^0.86 (eRASS1, mass quarters: 176 rad^2/m^4 at 2.54e14 Msun).
Entry 179 replicated the excess around 1,101 MCXC and PSZ2 clusters at l < 179.5 deg, at 0.93 times the relation's value for their
median mass.  Here those clusters are split at their median M500, and each half is tested on its own as in other_half.py
(500 longitude rotations): the excess inside R500 is set against the eRASS1 relation at the median mass of the clusters holding
the sources there, and the two halves against each other.  Nothing here was pre-registered.
Run from the repository root:  .venv/bin/python -I entries/179-cluster-faraday-other-half/mass_halves.py     (about 5 minutes)
"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cluster_rm as c
import other_half as oh

log = lambda s: print(s, flush=True)
rows, src, _ = c.load_spice()
S, _ = c.sample(src, c.CFG['fg_radius'], c.CFG['min_neigh'])
cl, info = oh.load_clusters()
med = float(np.median(cl['m500']))
ref = oh.CFG['ref_176']
out = dict(split_m500_1e13=med, halves={})
for k, (name, sel) in enumerate((('lower', cl['m500'] < med), ('upper', cl['m500'] >= med))):
    sub = {kk: v[sel] for kk, v in cl.items()}
    rr, x, which = c.profile_test(S, sub, [[0.0, 1.0], [1.0, 2.0]], 500, oh.CFG['seed'] + 100 + k, log=log, label=name)
    inA = x < 1
    mA = float(np.median(sub['m500'][which[inA]]))
    pred = ref['excess_at_2p54e14'] * (mA / 25.4) ** ref['exponent']
    out['halves'][name] = dict(clusters=int(sel.sum()), A=rr[0], B=rr[1], median_m500_1e13_holding_A=mA, entry176_predicts=float(pred),
                               ratio=float(rr[0]['excess_variance'] / pred))
    log('%s half (%d clusters): excess inside R500 %+.1f +- %.1f from %d sources (t %+.2f); relation at %.1fe13 predicts %.1f (ratio %.2f)' % (
        name, sel.sum(), rr[0]['excess_variance'], rr[0]['null_sd'], rr[0]['n'], rr[0]['t'], mA, pred, rr[0]['excess_variance'] / pred))
lo, up = out['halves']['lower'], out['halves']['upper']
diff = up['A']['excess_variance'] - lo['A']['excess_variance']
sd = float(np.hypot(up['A']['null_sd'], lo['A']['null_sd']))
slope = float(np.log(max(up['A']['excess_variance'], 1e-9) / max(lo['A']['excess_variance'], 1e-9)) / np.log(up['median_m500_1e13_holding_A'] / lo['median_m500_1e13_holding_A'])) \
    if lo['A']['excess_variance'] > 0 and up['A']['excess_variance'] > 0 else None
out['difference'] = dict(upper_minus_lower=float(diff), sd=sd, t=float(diff / sd), predicted_difference=float(up['entry176_predicts'] - lo['entry176_predicts']),
                         implied_exponent=slope)
log('upper - lower: %+.1f +- %.1f (t %+.2f); the relation predicts %+.1f; implied exponent %s' % (
    diff, sd, diff / sd, up['entry176_predicts'] - lo['entry176_predicts'], slope and round(slope, 2)))
json.dump(out, open(os.path.join(HERE, 'mass_halves.json'), 'w'), indent=1, default=float)
