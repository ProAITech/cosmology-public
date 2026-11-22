"""Entry 153: where Entry 150's brightness effect comes from.
Re-derives the regression of stacked satellite offset on wing asymmetry from the
browser table (diagonal weights). Jackknife intercept errors come from the browser,
where the regression was repeated in each of the 20 sky zones."""
import json, math, sys
d = json.load(open(sys.argv[1] if len(sys.argv) > 1 else 'grz_why_raw.json'))
cols = d['sub_cols']; S = [dict(zip(cols, r)) for r in d['sub']]
def wls(x, y, w):
    Sw = sum(w); Sx = sum(a*b for a, b in zip(w, x)); Sy = sum(a*b for a, b in zip(w, y))
    Sxx = sum(a*b*b for a, b in zip(w, x)); Sxy = sum(a*b*c for a, b, c in zip(w, x, y))
    D = Sw*Sxx - Sx*Sx
    c0 = (Sxx*Sy - Sx*Sxy)/D; c1 = (Sw*Sxy - Sx*Sy)/D
    return c0, c1, math.sqrt(Sxx/D), math.sqrt(Sw/D)
out = {}
for stat, err in [('mean500','e500'),('mean1000','e1000'),('mean1500','e1500'),('fit_mu','e_fit')]:
    x = [s['A_2000_4000'] for s in S]; y = [s[stat] for s in S]; w = [1/s[err]**2 for s in S]
    c0, c1, e0, e1 = wls(x, y, w)
    chi_line = sum(wi*(yi-c0-c1*xi)**2 for wi, xi, yi in zip(w, x, y))
    ybar = sum(wi*yi for wi, yi in zip(w, y))/sum(w)
    chi_flat = sum(wi*(yi-ybar)**2 for wi, yi in zip(w, y))
    jk = d['reg_browser_jackknife'][stat]
    out[stat] = dict(c0=round(c0,1), e0_formal=round(e0,1), e0_jk=jk['e0'], slope=round(c1,4), es_formal=round(e1,4), es_jk=jk['es'],
                     chi2_line=round(chi_line,1), chi2_flat=round(chi_flat,1), dof_line=len(S)-2,
                     z_dm=round((c0-d['pred_entry150_K0']['dark_matter'])/jk['e0'],2),
                     z_pl=round((c0-d['pred_entry150_K0']['plasma_only'])/jk['e0'],2))
# Entry 150's margin series against the subsample line (fit statistic)
r = out['fit_mu']; m = d['margin']
out['margin_vs_line'] = [dict(mg=q['mg'], A=q['A'], mu=q['mu'], e=q['e'], line=round(r['c0']+r['slope']*q['A'],1),
                              pull=round((q['mu']-(r['c0']+r['slope']*q['A']))/q['e'],2)) for q in m]
xs = [q['A'] for q in m]; ys = [q['mu'] for q in m]
n = len(xs); mx = sum(xs)/n; my = sum(ys)/n
sl = sum((a-mx)*(b-my) for a, b in zip(xs, ys))/sum((a-mx)**2 for a in xs)
out['margin_series_slope'] = round(sl, 4); out['margin_series_intercept'] = round(my - sl*mx, 1)
# radius-combined apparent bright/faint difference (independent pair sets, shared clusters)
R = d['apparent_split_by_radius']; w = [1/q['ed']**2 for q in R]
out['apparent_diff_all_R'] = [round(sum(wi*q['diff'] for wi, q in zip(w, R))/sum(w),1), round(1/math.sqrt(sum(w)),1)]
# z-scaling of the luminous wing asymmetry
W = d['wings']; lo, hi = W['lumLoZ_2000_4000'], W['lumHiZ_2000_4000']
ratio = lo[0]/hi[0]; eratio = ratio*math.hypot(lo[1]/lo[0], hi[1]/hi[0])
out['lum_wing_lo_over_hi'] = [round(ratio,2), round(eratio,2)]
q0lo, q0hi = S[0], S[1]
out['q0_wing_lo_over_hi'] = [round(q0lo['A_2000_4000']/q0hi['A_2000_4000'],2), round(abs(q0lo['A_2000_4000']/q0hi['A_2000_4000'])*math.hypot(q0lo['eA']/q0lo['A_2000_4000'], q0hi['eA']/q0hi['A_2000_4000']),2)]
# luminous minus faint wing contrast
out['wing_contrast'] = [round(W['fnt_2000_4000'][0]-W['lum_2000_4000'][0],1), round(math.hypot(W['fnt_2000_4000'][1], W['lum_2000_4000'][1]),1)]
print(json.dumps(out, indent=1))
json.dump(out, open('grz_why.json','w'), indent=1)
