"""Entry 156: the lithium rise from turnoff to subgiants in a second cluster, M30, and the subgiant plateau in NGC 6752.
m30.csv  (Gruyters et al. 2016, VizieR J/A+A/589/A61 table5): 'Teff,logg*100,A(Li)*100,e*100' per star, adopted (isochrone) temperatures.
n6752.csv (Gruyters et al. 2014, VizieR J/A+A/567/A72): 'Teff,logg*100,A(Li)*100,e*100,upper-limit flag,population (P/I/E)'.
Colour and H-alpha temperature scales for M30 (table1) were applied in the browser with dA(Li)/dTeff = 0.07 dex per 100 K.
NGC 6397 numbers come from Entry 155 (li6397.json)."""
import json, sys, numpy as np
def load(f):
    r = [x.split(',') for x in open(f).read().strip().split(';')]
    T = np.array([int(x[0]) for x in r]); g = np.array([int(x[1]) for x in r])/100; Li = np.array([int(x[2]) for x in r])/100
    lim = np.array([len(x) > 4 and x[4] == '1' for x in r]); pop = np.array([x[5] if len(x) > 5 else '-' for x in r])
    return T, g, Li, lim, pop
def st(v): return dict(n=int(len(v)), mean=round(float(v.mean()), 3), se=round(float(v.std(ddof=1)/np.sqrt(len(v))), 3), sd=round(float(v.std(ddof=1)), 3))
out = {}
T, g, Li, lim, pop = load(sys.argv[1] if len(sys.argv) > 1 else 'm30.csv')
assert len(T) == 119
outl = (g >= 3.7) & (Li <= 1.5)
TO = (g >= 4.0) & ~outl; SG = (g >= 3.70) & (g < 4.0) & ~outl; RG = (g >= 2.0) & (g < 3.2)
out['m30'] = dict(n=119, excluded_li_poor=int(outl.sum()), TO=dict(**st(Li[TO]), T=round(float(T[TO].mean()))), SGB=dict(**st(Li[SG]), T=round(float(T[SG].mean()))),
                  rise=round(float(Li[SG].mean()-Li[TO].mean()), 3), err=round(float(np.hypot(Li[SG].std(ddof=1)/np.sqrt(SG.sum()), Li[TO].std(ddof=1)/np.sqrt(TO.sum()))), 3),
                  giants=dict(**st(Li[RG]), T=round(float(T[RG].mean()))))
out['m30']['dilution'] = round(out['m30']['SGB']['mean']-out['m30']['giants']['mean'], 3)
out['m30']['teff_scales'] = {'isochrone (adopted)': [0.070, 0.040, 0, 0], 'Ramírez–Meléndez colours': [0.032, 0.041, 69, 14], 'Alonso colours': [0.069, 0.041, -155, -158], 'H-alpha profiles': [0.153, 0.058, -266, -148]}
edges = [0.5, 2.0, 2.8, 3.2, 3.4, 3.5, 3.6, 3.65, 3.7, 3.8, 3.9, 3.95, 4.0, 4.05, 4.25]
out['m30']['profile'] = [dict(lo=lo, hi=hi, T=round(float(T[s].mean())), **st(Li[s])) for lo, hi in zip(edges[:-1], edges[1:]) for s in [(g >= lo) & (g < hi) & ~outl] if s.sum() > 1]
pts30 = [[int(a), round(float(b), 2), round(float(c), 2)] for a, b, c in zip(T, g, Li)]
T6, g6, Li6, lim6, pop6 = load(sys.argv[2] if len(sys.argv) > 2 else 'n6752.csv')
assert len(T6) == 124
ok6 = ~lim6; pl6 = ok6 & (g6 >= 3.85)
out['n6752'] = dict(n=124, upper_limits=int(lim6.sum()), sgb_plateau=dict(**st(Li6[pl6]), T=round(float(T6[pl6].mean()))),
                    by_population={p: st(Li6[pl6 & (pop6 == p)]) for p in 'PIE'})
edges6 = [3.4, 3.5, 3.6, 3.65, 3.7, 3.75, 3.8, 3.85, 3.9, 3.95, 4.0, 4.05]
out['n6752']['profile'] = [dict(lo=lo, hi=hi, T=round(float(T6[s].mean())), **st(Li6[s])) for lo, hi in zip(edges6[:-1], edges6[1:]) for s in [(g6 >= lo) & (g6 < hi) & ok6] if s.sum() > 1]
pts6752 = [[int(a), round(float(b), 2), round(float(c), 2), int(d), str(e)] for a, b, c, d, e in zip(T6, g6, Li6, lim6, pop6)]
# combine with NGC 6397 (Entry 155)
E = json.load(open(sys.argv[3] if len(sys.argv) > 3 else 'li6397.json'))['summary']
r97, e97 = E['rise_combined']; r30, e30 = out['m30']['rise'], out['m30']['err']
w = np.array([1/e97**2, 1/e30**2]); comb = float((w*np.array([r97, r30])).sum()/w.sum()); ec = float(1/np.sqrt(w.sum()))
chi = float(((np.array([r97, r30])-comb)**2*w).sum())
out['combined_rise'] = dict(value=round(comb, 3), err=round(ec, 3), chi2=round(chi, 2), sig=round(comb/ec, 1))
out['harris_feh'] = {'NGC 6397': -1.95, 'NGC 6752': -1.55, 'M30': -2.12}
out['clusters'] = {'NGC 6397': dict(TO=E['turnoff']['Li'], SGB=E['sgb_peak']['Li'], giants=E['rgb_plateau']['Li']),
                   'M30': dict(TO=out['m30']['TO']['mean'], SGB=out['m30']['SGB']['mean'], giants=out['m30']['giants']['mean']),
                   'NGC 6752': dict(TO=None, SGB=out['n6752']['sgb_plateau']['mean'], giants=None)}
bbn = 2.74
out['gap_from_TO'] = {k: round(bbn-v['TO'], 3) for k, v in out['clusters'].items() if v['TO']}
print(json.dumps({k: v for k, v in out.items() if k not in ('m30', 'n6752')}, indent=1))
print(json.dumps({k: v for k, v in out['m30'].items() if k != 'profile'}))
print(json.dumps({k: v for k, v in out['n6752'].items() if k != 'profile'}))
json.dump(dict(summary=out, pts30=pts30, pts6752=pts6752), open('li3gc.json', 'w'))
