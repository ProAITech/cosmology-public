/* Workbench Entry 80: does the quasar Hubble diagram favour a static universe?
   Data: 2,421 quasars with rest-frame 2500 A and 2 keV fluxes (Lusso et al. 2020, VizieR J/A+A/642/A150).
   Run in a browser tab on a VizieR mirror:  await QH.load();  QH.run();
   The relation L_X = beta L_UV^gamma, written in fluxes with the catalogue's convention L = 4 pi dL^2 F, is
       log F_X = gamma log F_UV + beta' + 2 (gamma - 1) log dL(z).
   1. Cosmology-free: inside a narrow redshift bin every quasar has (nearly) the same dL, so the slope of log F_X on
      log F_UV is gamma itself. If gamma drifts with redshift, the relation cannot be used as a ruler.
   2. The diagram: with gamma, beta' and an intrinsic scatter delta fitted, compare the shape of dL(z) for:
      flat LCDM (Om 0.3, and Om free), Lerner's static law d = cz/H0 with energy loss (dL = z sqrt(1+z)), tired light
      in static space (dL = ln(1+z) sqrt(1+z)), the logarithm with stretched time (dL = ln(1+z)(1+z)), Milne
      (dL = z(1+z/2)) and Einstein-de Sitter. H0 only shifts beta'. For a static model the catalogue's 'dL' is the
      effective distance with L_nu = 4 pi d^2 (1+z) F_cat, which is d sqrt(1+z): energy loss without time dilation.
   Likelihood: Gaussian, variance e_FX^2 + gamma^2 e_FUV^2 + delta^2.                                            */
window.QH = (function () {
  const Q = {};
  Q.load = async function () {
    const ft = async u => { for (let k = 0; k < 8; k++) { try { const r = await fetch(u); if (r.ok) { const t = await r.text(); if (!t.includes('STDB_ERR')) return t; } } catch (e) {} await new Promise(r => setTimeout(r, 1500)); } return 'FAIL'; };
    const t = await ft('/viz-bin/asu-tsv?-source=J/A%2BA/642/A150&-out=z,logFUV,e_logFUV,logFX,e_logFX,Group&-out.max=unlimited&-oc.form=dec');
    const rows = t.split('\n').filter(l => l.length && !l.startsWith('#')).slice(3).map(l => l.split('\t').map(Number)).filter(r => r.length >= 6 && r.every(Number.isFinite));
    Q.d = rows; return rows.length;
  };
  const E = (z, om) => Math.sqrt(om * Math.pow(1 + z, 3) + 1 - om);
  const dc = (z, om) => { const n = 400, h = z / n; let s = 0; for (let i = 0; i <= n; i++) { const w = (i === 0 || i === n) ? 1 : (i % 2 ? 4 : 2); s += w / E(i * h, om); } return s * h / 3; };
  Q.models = {
    'LCDM Om=0.3': z => (1 + z) * dc(z, 0.3),
    'Lerner static, d=cz/H0': z => z * Math.sqrt(1 + z),
    'tired light, d=ln(1+z)': z => Math.log(1 + z) * Math.sqrt(1 + z),
    'logarithm with stretched time': z => Math.log(1 + z) * (1 + z),
    'Milne (empty)': z => z * (1 + z / 2),
    'Einstein-de Sitter': z => 2 * (1 + z) * (1 - 1 / Math.sqrt(1 + z))
  };
  // weighted fit of y = g x + b with scatter delta, iterating weights; returns -2 ln L
  Q.fit = function (x, y, ex, ey) {
    const n = x.length;
    const solve = (dl, g0) => { let g = g0, b = 0; for (let it = 0; it < 6; it++) { let S = 0, Sx = 0, Sy = 0, Sxx = 0, Sxy = 0; for (let i = 0; i < n; i++) { const w = 1 / (ey[i] ** 2 + g * g * ex[i] ** 2 + dl * dl); S += w; Sx += w * x[i]; Sy += w * y[i]; Sxx += w * x[i] * x[i]; Sxy += w * x[i] * y[i]; } g = (S * Sxy - Sx * Sy) / (S * Sxx - Sx * Sx); b = (Sy - g * Sx) / S; }
      let m2 = 0; for (let i = 0; i < n; i++) { const v = ey[i] ** 2 + g * g * ex[i] ** 2 + dl * dl, r = y[i] - g * x[i] - b; m2 += r * r / v + Math.log(2 * Math.PI * v); } return { g, b, m2 }; };
    let lo = 0.01, hi = 1.0, best = null; // golden-section on delta
    const gr = (Math.sqrt(5) - 1) / 2; let a = lo, c = hi, x1 = c - gr * (c - a), x2 = a + gr * (c - a), f1 = solve(x1, 0.6), f2 = solve(x2, 0.6);
    for (let k = 0; k < 60; k++) { if (f1.m2 < f2.m2) { c = x2; x2 = x1; f2 = f1; x1 = c - gr * (c - a); f1 = solve(x1, f2.g); } else { a = x1; x1 = x2; f1 = f2; x2 = a + gr * (c - a); f2 = solve(x2, f1.g); } }
    best = f1.m2 < f2.m2 ? Object.assign(f1, { delta: x1 }) : Object.assign(f2, { delta: x2 });
    // slope error from the Fisher information at the best delta
    let S = 0, Sx = 0, Sxx = 0; for (let i = 0; i < n; i++) { const w = 1 / (ey[i] ** 2 + best.g * best.g * ex[i] ** 2 + best.delta ** 2); S += w; Sx += w * x[i]; Sxx += w * x[i] * x[i]; }
    best.eg = Math.sqrt(S / (S * Sxx - Sx * Sx)); best.n = n; return best;
  };
  Q.run = function (zmin = 0, zmax = 10, nb = 12, ref = 'LCDM Om=0.3') {
    const d = Q.d.filter(r => r[0] >= zmin && r[0] < zmax).sort((a, b) => a[0] - b[0]);
    const out = { n: d.length, bin_reference: ref };
    // 1. slope in redshift bins of equal count. Inside a bin the small spread in distance would pull the slope toward 1
    //    (both fluxes fall together), so each quasar is moved to the bin's median redshift with a reference distance
    //    law; the result is checked with a static reference law as well (the two differ only at second order in a bin).
    out.bins = [];
    const R = Q.models[ref];
    for (let k = 0; k < nb; k++) {
      const s = d.slice(Math.floor(k * d.length / nb), Math.floor((k + 1) * d.length / nb));
      const zm = s[s.length >> 1][0], Lm = 2 * Math.log10(R(zm));
      const sh = s.map(r => 2 * Math.log10(R(r[0])) - Lm);
      const xs = s.map((r, i) => r[1] + sh[i]), xm = xs.reduce((a, v) => a + v, 0) / xs.length;
      const f = Q.fit(xs.map(v => v - xm), s.map((r, i) => r[3] + sh[i]), s.map(r => r[2]), s.map(r => r[4]));
      out.bins.push({ z_lo: s[0][0], z_hi: s[s.length - 1][0], z_med: s[s.length >> 1][0], n: s.length, gamma: +f.g.toFixed(3), e_gamma: +f.eg.toFixed(3), delta: +f.delta.toFixed(3) });
    }
    // 2. the diagram, model by model
    out.models = {};
    const x = d.map(r => r[1]), y = d.map(r => r[3]), ex = d.map(r => r[2]), ey = d.map(r => r[4]);
    const fitModel = dl => { // y = g x + b + 2(g-1) log dl  ->  y + 2 log dl = g (x + 2 log dl) + b
      const L = d.map(r => 2 * Math.log10(dl(r[0])));
      return Q.fit(x.map((v, i) => v + L[i]), y.map((v, i) => v + L[i]), ex, ey);
    };
    for (const [nm, f] of Object.entries(Q.models)) { const r = fitModel(f); out.models[nm] = { m2lnL: +r.m2.toFixed(2), gamma: +r.g.toFixed(4), e_gamma: +r.eg.toFixed(4), delta: +r.delta.toFixed(4), beta: +r.b.toFixed(3) }; }
    const scan = []; for (let om = 0; om <= 1.0001; om += 0.05) { const r = fitModel(z => (1 + z) * dc(z, om)); scan.push([+om.toFixed(2), +r.m2.toFixed(2)]); }
    out.LCDM_Om_scan = scan;
    return out;
  };
  // binned diagram for the figure: intercepts y - g x in equal-count bins, with gamma from one model
  Q.binned = function (g, nb = 20) {
    const d = Q.d.slice().sort((a, b) => a[0] - b[0]); const out = [];
    for (let k = 0; k < nb; k++) {
      const s = d.slice(Math.floor(k * d.length / nb), Math.floor((k + 1) * d.length / nb));
      const v = s.map(r => r[3] - g * r[1]); const m = v.reduce((a, b) => a + b, 0) / v.length;
      const sd = Math.sqrt(v.reduce((a, b) => a + (b - m) ** 2, 0) / (v.length - 1));
      out.push([+(s.reduce((a, r) => a + Math.log10(1 + r[0]), 0) / s.length).toFixed(4), +m.toFixed(4), +(sd / Math.sqrt(v.length)).toFixed(4), s.length]);
    }
    return out;
  };
  return Q;
})();
/* ---- Added during the round: the comparisons, subsets, selection checks and mocks reported in Entry 80 ---- */
(function (QH) {
// general weighted linear fit with intrinsic scatter (first column's coefficient is gamma)
QH.lin = function (cols, y, ex, ey) {
  const n = y.length, k = cols.length;
  const solve = (dl, g) => { let p; for (let it = 0; it < 5; it++) { const M = Array.from({length:k},()=>new Array(k).fill(0)), v = new Array(k).fill(0);
      for (let i = 0; i < n; i++) { const w = 1/(ey[i]**2 + g*g*ex[i]**2 + dl*dl); for (let a = 0; a < k; a++) { v[a] += w*cols[a][i]*y[i]; for (let b = 0; b < k; b++) M[a][b] += w*cols[a][i]*cols[b][i]; } }
      const A = M.map((r,i)=>r.concat([v[i]])); for (let c = 0; c < k; c++) { let piv = c; for (let r = c+1; r < k; r++) if (Math.abs(A[r][c]) > Math.abs(A[piv][c])) piv = r; [A[c],A[piv]]=[A[piv],A[c]]; for (let r = 0; r < k; r++) if (r!==c) { const f = A[r][c]/A[c][c]; for (let j = c; j <= k; j++) A[r][j] -= f*A[c][j]; } }
      p = A.map((r,i)=>r[k]/r[i]); g = p[0]; }
    let m2 = 0; for (let i = 0; i < n; i++) { let mu = 0; for (let a = 0; a < k; a++) mu += p[a]*cols[a][i]; const vv = ey[i]**2 + g*g*ex[i]**2 + dl*dl; m2 += (y[i]-mu)**2/vv + Math.log(2*Math.PI*vv); } return {p, m2}; };
  const gr = (Math.sqrt(5)-1)/2; let a = 0.01, c = 1.0, x1 = c-gr*(c-a), x2 = a+gr*(c-a), f1 = solve(x1,0.6), f2 = solve(x2,0.6);
  for (let it = 0; it < 50; it++) { if (f1.m2 < f2.m2) { c = x2; x2 = x1; f2 = f1; x1 = c-gr*(c-a); f1 = solve(x1, f2.p[0]); } else { a = x1; x1 = x2; f1 = f2; x2 = a+gr*(c-a); f2 = solve(x2, f1.p[0]); } }
  return f1.m2 < f2.m2 ? Object.assign(f1,{delta:x1}) : Object.assign(f2,{delta:x2});
};
// -2 ln L of each law relative to LCDM, for quasars with zmin <= z < zmax; evolve adds a slope change pivoting at the mean luminosity
QH.compare = function (zmin, zmax, evolve) {
  const d = QH.d.filter(r => r[0] >= zmin && r[0] < zmax), t0 = d.reduce((a,r)=>a+Math.log10(1+r[0]),0)/d.length, out = {};
  for (const [nm, f] of Object.entries(QH.models)) {
    const L = d.map(r => 2*Math.log10(f(r[0]))), X = d.map((r,i)=>r[1]+L[i]), Y = d.map((r,i)=>r[3]+L[i]), Xm = X.reduce((a,v)=>a+v,0)/X.length;
    const cols = [X, d.map(()=>1)]; if (evolve) cols.push(X.map((v,i)=>(v-Xm)*(Math.log10(1+d[i][0])-t0)));
    const r = QH.lin(cols, Y, d.map(r=>r[2]), d.map(r=>r[4])); out[nm] = {m2: r.m2, g: r.p[0], g1: evolve ? r.p[2] : 0};
  }
  const b = out['LCDM Om=0.3'].m2;
  return {n: d.length, rel: Object.fromEntries(Object.entries(out).map(([k,v])=>[k.slice(0,14), [+(v.m2-b).toFixed(1), +v.g.toFixed(3), +v.g1.toFixed(3)]]))};
};
// the ultraviolet-bright half: split at the median log F_UV in each of 40 equal-count redshift slices
QH.halves = function (rows) { const s = rows.slice().sort((a,b)=>a[0]-b[0]), br = [], fa = [];
  for (let k = 0; k < 40; k++) { const q = s.slice(Math.floor(k*s.length/40), Math.floor((k+1)*s.length/40)); const m = q.map(r=>r[1]).sort((a,b)=>a-b)[q.length>>1]; q.forEach(r => (r[1] >= m ? br : fa).push(r)); }
  return {bright: br, faint: fa}; };
// mocks: real z and F_UV; F_X from a fixed relation in LCDM (gamma 0.6, scatter delta, beta matched below z = 1);
// optional X-ray floor = 2nd percentile of the real log F_X in 24 equal-count redshift slices (redrawn until above it)
QH.makeMock = function (seed0) {
  const full = QH.d.slice(), srt = full.slice().sort((a,b)=>a[0]-b[0]), NS = 24, floors = [];
  for (let k = 0; k < NS; k++) { const s = srt.slice(Math.floor(k*srt.length/NS), Math.floor((k+1)*srt.length/NS)); const fx = s.map(r=>r[3]).sort((a,b)=>a-b); floors.push([s[0][0], s[s.length-1][0], fx[Math.floor(0.02*fx.length)]]); }
  const floorAt = z => { for (const f of floors) if (z <= f[1]) return f[2]; return floors[NS-1][2]; };
  let seed = seed0 || 7; const rand = () => { seed = (seed*16807) % 2147483647; return seed/2147483647; }; const gauss = () => { let u = 0; for (let i = 0; i < 12; i++) u += rand(); return u - 6; };
  const lcdm = QH.models['LCDM Om=0.3'];
  const res = full.filter(r=>r[0]<1).map(r=>{const L=2*Math.log10(lcdm(r[0])); return (r[3]+L) - 0.6*(r[1]+L);}).sort((a,b)=>a-b), beta = res[res.length>>1];
  QH.mock = function (useFloor, delta) { const out = []; for (const r of full) { const L = 2*Math.log10(lcdm(r[0])); for (let t = 0; t < 30; t++) { const fx = 0.6*(r[1]+L) + beta + delta*gauss() - L; if (!useFloor || fx >= floorAt(r[0])) { out.push([r[0], r[1], r[2], fx, r[4], r[5]]); break; } } } return out; };
  return {beta, floors};
};
/* Reported: QH.run(); QH.compare(0,10), (0,1.5), (1.5,10), (0,10,true); QH.compare on QH.halves(QH.d).bright / .faint and on Group 5
   (SDSS-4XMM); QH.run(0,10,12,'tired light, d=ln(1+z)') for the slope check; 20 mocks with and 20 without the floor, delta 0.24,
   one random stream started at seed 7 (four trial mocks drawn first). */
})(window.QH);
