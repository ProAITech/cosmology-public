/* Entry 186: the cluster Faraday signal at 1.4 GHz. NVSS rotation measures around eRASS1 clusters.
   Run in a browser tab on https://vizier.cfa.harvard.edu (all data stay in browser memory; only results leave).
   Data (VizieR, /viz-bin/asu-tsv):
     NVSS RMs, Taylor, Stil & Sunstrum (2009), J/ApJ/702/1230/catalog: RAJ2000, DEJ2000, RM, e_RM (37,543 sources).
     eRASS1 clusters, Bulbul et al. (2024), J/A+A/685/A106/emain: RAJ2000, DEJ2000, zBest, pcont, R500 (kpc).
   Method (Entry 176's, at 1.4 GHz): sources at |b| >= 20 deg with at least 10 other NVSS RM sources within 3 deg; dRM = RM - median RM
   of those others. Clusters with pcont < 0.5, zBest >= 0.05, R500 > 0; theta500 = R500 / D_A(z), flat LCDM H0 = 67.74, Om = 0.3089.
   x = smallest theta / theta500 over all clusters (searched to 8 theta500). V = (1.4826 MAD dRM)^2 - mean e_RM^2 in each range of x.
   Null: the cluster list rotated in Galactic longitude by 500 random angles in 10-350 deg; t = (V - mean) / sd.
   Use: E186.load() then E186.synthetic() or E186.run(); results in E186.state.                                                      */
window.E186 = (function () {
  'use strict';
  const D2R = Math.PI / 180;
  const CFG = { fg: 3.0, minNeigh: 10, bcut: 20, xmax: 8, nrot: 500, seed: 20261186, H0: 67.74, Om: 0.3089,
    ranges: [[0, 0.5], [0.5, 1], [0, 1], [1, 2], [2, 4]], spice: { '0-0.5': [174.6, 21.5], '0.5-1': [79.7, 11.9], '0-1': [105.5, 10.8] } };
  const state = { log: [], stage: 'idle' };
  const log = m => state.log.push(m);
  const rng = seed => { let a = seed >>> 0; return () => { a = (a + 0x6D2B79F5) >>> 0; let t = Math.imul(a ^ (a >>> 15), 1 | a); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; }; };
  const MG = [[-0.0548755604, -0.8734370902, -0.4838350155], [0.4941094279, -0.4448296300, 0.7469822445], [-0.8676661490, -0.1980763734, 0.4559837762]];
  function gal(ra, dec) { const a = ra * D2R, d = dec * D2R, v = [Math.cos(d) * Math.cos(a), Math.cos(d) * Math.sin(a), Math.sin(d)];
    const x = MG[0][0] * v[0] + MG[0][1] * v[1] + MG[0][2] * v[2], y = MG[1][0] * v[0] + MG[1][1] * v[1] + MG[1][2] * v[2], z = MG[2][0] * v[0] + MG[2][1] * v[1] + MG[2][2] * v[2];
    return [((Math.atan2(y, x) / D2R) + 360) % 360, Math.asin(Math.max(-1, Math.min(1, z))) / D2R]; }
  const unit = (l, b) => [Math.cos(b * D2R) * Math.cos(l * D2R), Math.cos(b * D2R) * Math.sin(l * D2R), Math.sin(b * D2R)];
  function median(a) { const s = Float64Array.from(a).sort(), n = s.length; return n % 2 ? s[(n - 1) / 2] : 0.5 * (s[n / 2 - 1] + s[n / 2]); }
  function Vof(d, e) { const m = median(d); const ad = d.map(v => Math.abs(v - m)); let e2 = 0; for (const v of e) e2 += v * v; return Math.pow(1.4826 * median(ad), 2) - e2 / e.length; }
  function DA(z) { let s = 0; const n = 400; for (let i = 0; i < n; i++) { const zz = (i + 0.5) * z / n; s += 1 / Math.sqrt(CFG.Om * Math.pow(1 + zz, 3) + 1 - CFG.Om); } return 299792.458 / CFG.H0 * s * z / n / (1 + z); }
  async function tsv(src, cols, extra) {
    const r = await fetch('/viz-bin/asu-tsv?-source=' + src + '&-out.max=unlimited&-out=' + cols.join(',') + (extra || ''));
    const lines = (await r.text()).split('\n').filter(l => l.trim() && !l.startsWith('#'));
    const head = lines[0].split('\t').map(h => h.trim()); const rows = [];
    for (const l of lines.slice(3)) { const f = l.split('\t'); if (f.length < head.length) continue; const o = {}; head.forEach((h, i) => o[h] = f[i].trim()); rows.push(o); }
    return rows;
  }
  const sexa = (s, h) => { if (!/[ :]/.test(s.trim())) return parseFloat(s); const sg = /^-/.test(s.trim()) ? -1 : 1; const p = s.trim().replace(/^[+-]/, '').split(/[ :]+/).map(Number); return sg * (p[0] + p[1] / 60 + p[2] / 3600) * (h ? 15 : 1); };
  // ---- a simple grid on the unit sphere: cells of 1 deg in Galactic latitude and longitude
  function grid(L, B) { const g = new Map(); for (let i = 0; i < L.length; i++) { const k = Math.floor(B[i] + 90) * 400 + Math.floor(L[i]); if (!g.has(k)) g.set(k, []); g.get(k).push(i); } return g; }
  function near(g, L, B, V, l0, b0, rdeg) {   // indices within rdeg of (l0, b0)
    const out = [], u0 = unit(l0, b0), c = Math.cos(rdeg * D2R), nb = Math.ceil(rdeg) + 1;
    for (let bb = Math.floor(b0 + 90) - nb; bb <= Math.floor(b0 + 90) + nb; bb++) {
      if (bb < 0 || bb > 179) continue;
      const cb = Math.cos(Math.min(89.9, Math.abs(bb - 90) + 1) * D2R), nl = Math.min(180, Math.ceil(nb / Math.max(cb, 0.02)));
      for (let dl = -nl; dl <= nl; dl++) { const ll = ((Math.floor(l0) + dl) % 360 + 360) % 360; const a = g.get(bb * 400 + ll); if (!a) continue;
        for (const i of a) { const v = V[i]; if (v[0] * u0[0] + v[1] * u0[1] + v[2] * u0[2] >= c) out.push(i); } }
    }
    return out;
  }
  let S = null, C = null;
  async function load() {
    state.stage = 'loading';
    const nv = await tsv('J/ApJ/702/1230/catalog', ['RAJ2000', 'DEJ2000', 'RM', 'e_RM'], '&-oc.form=d');
    const L = [], B = [], RM = [], E = [];
    for (const o of nv) { const ra = sexa(o.RAJ2000, true), de = sexa(o.DEJ2000, false), rm = parseFloat(o.RM), e = parseFloat(o.e_RM);
      if (!isFinite(ra) || !isFinite(de) || !isFinite(rm) || !(e > 0)) continue; const [l, b] = gal(ra, de); L.push(l); B.push(b); RM.push(rm); E.push(e); }
    const V = L.map((l, i) => unit(l, B[i])), g = grid(L, B);
    // foreground: median RM of the other sources within 3 deg; keep |b| >= 20 with >= 10 others
    const keep = [], d = [], e = [];
    for (let i = 0; i < L.length; i++) { if (Math.abs(B[i]) < CFG.bcut) continue; const nb = near(g, L, B, V, L[i], B[i], CFG.fg).filter(j => j !== i);
      if (nb.length < CFG.minNeigh) continue; keep.push(i); d.push(RM[i] - median(nb.map(j => RM[j]))); e.push(E[i]); }
    S = { l: keep.map(i => L[i]), b: keep.map(i => B[i]), d, e }; S.V = S.l.map((l, i) => unit(l, S.b[i])); S.g = grid(S.l, S.b);
    const cl = await tsv('J/A%2bA/685/A106/emain', ['RAJ2000', 'DEJ2000', 'zBest', 'pcont', 'R500'], '&-oc.form=d');
    const cL = [], cB = [], th = [];
    for (const o of cl) { const ra = parseFloat(o.RAJ2000), de = parseFloat(o.DEJ2000), z = parseFloat(o.zBest), pc = parseFloat(o.pcont), r5 = parseFloat(o.R500);
      if (!isFinite(ra) || !isFinite(de) || !isFinite(z) || !(pc < 0.5) || !(z >= 0.05) || !(r5 > 0)) continue; const [l, b] = gal(ra, de);
      cL.push(l); cB.push(b); th.push(r5 / 1000 / DA(z) / D2R); }
    C = { l: cL, b: cB, th };
    state.info = { nvss_rows: nv.length, nvss_usable: L.length, sources: d.length, cluster_rows: cl.length, clusters: th.length,
      median_e: median(e), median_theta500_arcmin: 60 * median(th) };
    state.stage = 'loaded'; log('loaded ' + JSON.stringify(state.info)); return state.info;
  }
  function xdist(shift) {
    const n = S.d.length, x = new Float64Array(n).fill(Infinity);
    for (let k = 0; k < C.th.length; k++) { const l0 = (C.l[k] + shift) % 360, b0 = C.b[k], r = Math.min(CFG.xmax * C.th[k], 3.0);
      const u0 = unit(l0, b0);
      for (const i of near(S.g, S.l, S.b, S.V, l0, b0, r)) { const v = S.V[i]; const sep = 2 * Math.asin(Math.min(1, Math.hypot(v[0] - u0[0], v[1] - u0[1], v[2] - u0[2]) / 2)) / D2R; const xs = sep / C.th[k]; if (xs < x[i]) x[i] = xs; } }
    return x;
  }
  function Vin(x, d, e) { return CFG.ranges.map(([lo, hi]) => { const dd = [], ee = []; for (let i = 0; i < x.length; i++) if (x[i] >= lo && x[i] < hi) { dd.push(d[i]); ee.push(e[i]); } return [dd.length >= 10 ? Vof(dd, ee) : NaN, dd.length]; }); }
  async function test(d, e, nrot, seed) {
    const x0 = xdist(0), real = Vin(x0, d, e), R = rng(seed), nulls = CFG.ranges.map(() => []), counts = CFG.ranges.map(() => []);
    for (let r = 0; r < nrot; r++) { const xr = xdist(10 + 340 * R()); const v = Vin(xr, d, e); v.forEach(([V, n], j) => { if (isFinite(V)) nulls[j].push(V); counts[j].push(n); });
      if (r % 50 === 0) await new Promise(res => setTimeout(res, 0)); }
    return CFG.ranges.map(([lo, hi], j) => { const a = nulls[j], mu = a.reduce((s, v) => s + v, 0) / a.length, sd = Math.sqrt(a.reduce((s, v) => s + (v - mu) ** 2, 0) / (a.length - 1));
      return { range: lo + '-' + hi, n: real[j][1], n_rot: counts[j].reduce((s, v) => s + v, 0) / counts[j].length, V: real[j][0], null_mean: mu, null_sd: sd,
        t: (real[j][0] - mu) / sd, excess: real[j][0] - mu, upper95: real[j][0] - mu + 1.645 * sd }; });
  }
  function decide(R) {
    const by = Object.fromEntries(R.map(r => [r.range, r])), A = by['0-1'], I = by['0-0.5'];
    const det = A.n < 30 ? 'not tested (fewer than 30 sources inside R500)' : A.t >= 3 ? 'cluster signal seen at 1.4 GHz' : A.t < 2 ? 'not seen at 1.4 GHz' : 'unresolved';
    const sp = CFG.spice['0-0.5'], tD = (I.excess - sp[0]) / Math.hypot(I.null_sd, sp[1]);
    const dep = I.n < 20 ? 'not tested (fewer than 20 sources inside 0.5 R500)' : tD >= 3 ? 'inner excess larger at 1.4 GHz than at 0.9 GHz: depolarization selection supported'
      : tD >= 2 ? 'inner excess leans larger at 1.4 GHz (2 <= t_D < 3)' : tD <= -2 ? 'inner excess smaller at 1.4 GHz' : 'no difference detected between 1.4 and 0.9 GHz inside 0.5 R500';
    return { detection: det, depolarization: dep, t_D: tD, D: I.excess - sp[0] };
  }
  async function synthetic() {   // simulated RMs at the real positions and errors: none, and a cluster term inside R500
    state.stage = 'synthetic'; const out = {};
    const x0 = xdist(0), g = rng(186);
    const gauss = () => { let u = 0, v = 0; while (u === 0) u = g(); v = g(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v); };
    for (const [name, add] of [['none', x => 0], ['cluster 400/100', x => x < 0.5 ? 400 : x < 1 ? 100 : 0]]) {
      const d = S.d.map((_, i) => gauss() * Math.sqrt(70 + add(x0[i])) + gauss() * S.e[i]);
      const R = await test(d, S.e, 200, CFG.seed + 7); out[name] = { ranges: R, decision: decide(R) };
      log('synthetic ' + name + ': ' + R.map(r => r.range + ' n ' + r.n + ' ex ' + r.excess.toFixed(1) + ' +- ' + r.null_sd.toFixed(1) + ' t ' + r.t.toFixed(2)).join('; ') + ' | ' + JSON.stringify(out[name].decision));
    }
    state.synthetic = out; state.stage = 'synthetic done'; return out;
  }
  async function run() {
    state.stage = 'running'; const R = await test(S.d, S.e, CFG.nrot, CFG.seed);
    state.result = { ranges: R, decision: decide(R), info: state.info, spice_entry176: CFG.spice };
    state.stage = 'done'; log('run: ' + JSON.stringify(state.result.decision)); return state.result;
  }
  return { load, synthetic, run, state, CFG };
})();
