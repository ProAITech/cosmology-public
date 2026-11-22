"""Entry 155: is lithium hidden below the surface of old stars? NGC 6397, Lind et al. (2009), VizieR J/A+A/503/545 table4.
n6397.csv: per star 'Teff,logL*100,A(Li)*100,e*100,upper-limit flag,A(Na)*100[u = Na upper limit]' in table order;
rows 0-304 are the first sample, 305-467 the second (stars with a second identifier); 14 stars appear in both.
Teff-scale sensitivity of the rise was computed in the browser from table2 (alternative scales) with dA(Li)/dTeff = 0.07 dex per 100 K."""
import json, sys, numpy as np
r = [x.split(',') for x in open(sys.argv[1] if len(sys.argv) > 1 else 'n6397.csv').read().strip().split(';')]
T = np.array([int(x[0]) for x in r]); lL = np.array([int(x[1]) for x in r])/100; Li = np.array([int(x[2]) for x in r])/100
eLi = np.array([int(x[3]) for x in r])/100; lim = np.array([int(x[4]) for x in r]) == 1
Na = np.array([int(x[5].rstrip('u')) for x in r])/100; Nau = np.array([x[5].endswith('u') for x in r])
blk = np.array([0]*305+[1]*163); assert len(r) == 468
dup = [(137,337),(151,384),(131,334),(139,340),(118,316),(130,329),(144,350),(153,387),(136,336),(125,325),(117,314),(124,318),(122,317),(145,359)]
dd = np.array([Li[b]-Li[a] for a, b in dup])
out = {'n': 468, 'dup_offset_B_minus_A': [round(float(dd.mean()), 3), round(float(dd.std(ddof=1)/np.sqrt(len(dd))), 3), len(dd)]}
keep = np.ones(len(r), bool)
LiA = Li.copy()
for a, b in dup: LiA[a] = (Li[a]+Li[b])/2; keep[b] = False
ok = keep & ~lim
def stat(s, y=LiA):
    v = y[s]; return dict(n=int(s.sum()), T=round(float(T[s].mean())), logL=round(float(lL[s].mean()), 2), Li=round(float(v.mean()), 3), se=round(float(v.std(ddof=1)/np.sqrt(len(v))), 3), sd=round(float(v.std(ddof=1)), 3))
# profile along the sequence (both samples)
edges = [-0.1, 0.1, 0.25, 0.33, 0.40, 0.47, 0.54, 0.60, 0.64, 0.68, 0.72, 0.80, 0.90, 1.0, 1.3, 1.6, 1.9, 2.6]
out['profile'] = [dict(lo=lo, hi=hi, **stat(ok & (lL >= lo) & (lL < hi))) for lo, hi in zip(edges[:-1], edges[1:]) if (ok & (lL >= lo) & (lL < hi)).sum() > 1]
# the rise, within each homogeneous sample (no duplicate averaging, so samples stay independent)
okr = ~lim
def rise(b, lo0, hi0):
    a = okr & (blk == b) & (lL >= lo0) & (lL < hi0); c = okr & (blk == b) & (lL >= 0.54) & (lL < 0.64)
    return dict(ref=stat(a, Li), sgb=stat(c, Li), rise=round(float(Li[c].mean()-Li[a].mean()), 3),
                err=round(float(np.hypot(Li[c].std(ddof=1)/np.sqrt(c.sum()), Li[a].std(ddof=1)/np.sqrt(a.sum()))), 3))
out['rise_sample1_TO'] = rise(0, 0.25, 0.47)
out['rise_sample2_MS'] = rise(1, -0.1, 0.25)
w = [1/out[k]['err']**2 for k in ('rise_sample1_TO', 'rise_sample2_MS')]
cm = sum(wi*out[k]['rise'] for wi, k in zip(w, ('rise_sample1_TO', 'rise_sample2_MS')))/sum(w)
out['rise_combined'] = [round(cm, 3), round(1/np.sqrt(sum(w)), 3)]
# first and second generation (sodium), sample 1
gen = {}
for lab, lo0, hi0 in (('TO', 0.25, 0.47), ('SGB', 0.47, 0.64)):
    s = okr & (blk == 0) & (lL >= lo0) & (lL < hi0) & (Na > 0) & ~Nau
    med = float(np.median(Na[s])); gen[lab] = dict(Na_median=med, Na_poor=stat(s & (Na < med), Li), Na_rich=stat(s & (Na >= med), Li))
s1 = okr & (blk == 0) & (Na > 0) & ~Nau
a1 = s1 & (lL >= 0.25) & (lL < 0.47); c1 = s1 & (lL >= 0.54) & (lL < 0.64)
med = float(np.median(Na[s1 & (lL >= 0.25) & (lL < 0.64)]))
pa, pc = a1 & (Na < med), c1 & (Na < med)
gen['first_gen_rise'] = [round(float(Li[pc].mean()-Li[pa].mean()), 3), round(float(np.hypot(Li[pc].std(ddof=1)/np.sqrt(pc.sum()), Li[pa].std(ddof=1)/np.sqrt(pa.sum()))), 3), int(pa.sum()), int(pc.sum())]
out['generations'] = gen
# Teff-scale robustness (browser, table2 alternative scales): rise, error, mean dT at TO/ref and SGB
out['teff_scales'] = {'sample1': {'O09 (adopted)': [0.052, 0.014, 0, 0], 'RM05': [0.063, 0.014, -95, -81], 'A96/99': [0.086, 0.014, -137, -88]},
                      'sample2': {'O09 (adopted)': [0.039, 0.013, 0, 0], 'RM05': [0.091, 0.013, -155, -80], 'A96/99': [0.122, 0.013, -207, -88]}}
# post-dilution plateau and the gap
plat = ok & (lL >= 1.0) & (lL < 1.9); out['rgb_plateau'] = stat(plat)
sg = ok & (lL >= 0.54) & (lL < 0.64); out['sgb_peak'] = stat(sg)
to = ok & (lL >= 0.25) & (lL < 0.47); out['turnoff'] = stat(to)
out['bbn'] = 2.74
out['gap_from_TO'] = round(2.74-out['turnoff']['Li'], 3)
out['visible_share_of_gap'] = [round(out['rise_combined'][0]/out['gap_from_TO'], 3), round(0.122/out['gap_from_TO'], 3)]
out['dilution_dex'] = round(out['sgb_peak']['Li']-out['rgb_plateau']['Li'], 3)
print(json.dumps(out, indent=1))
pts = [[int(T[i]), round(float(lL[i]), 2), round(float(LiA[i]), 2), int(lim[i]), int(blk[i])] for i in range(len(r)) if keep[i]]
json.dump(dict(summary=out, pts=pts), open('li6397.json', 'w'))
