// quasar_dilation.js - Workbench Entry 135. Runs in a browser on VizieR J/MNRAS/514/164/table1 (Stone et al. 2022):
// 190 SDSS Stripe 82 quasars with ~20-year light curves (SDSS + PS1 + DES), damped-random-walk timescales in g, r, i.
// Model: log tau_obs = a + n log(1+z) + b log(lambda_rest / 4000 A) + c log(M_BH / 1e9 Msun)  [+ d (log L_bol - 46)]
// Expansion predicts n = 1 (clocks at redshift z run slow by 1+z); a static universe predicts n = 0.
// Fit: weighted least squares with an intrinsic scatter found by maximum likelihood; uncertainty by bootstrap over quasars.
const QD = (() => {
  const URL = 'https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=J/MNRAS/514/164/table1&-out.all&-out.max=unlimited&-out.form=|';
  function parse(tsv) {
    const L = tsv.split('\n').filter(l => l && !l.startsWith('#')); const H = L[0].split('|').map(s => s.trim());
    const col = n => H.indexOf(n);
    return L.slice(3).map(l => { const r = l.split('|').map(s => s.trim()); const g = n => parseFloat(r[col(n)]);
      const q = { z: g('z'), m: g('log(Mbh)'), lb: g('log(Lbol)'), bands: [] };
      for (const b of ['g', 'r', 'i']) {
        const t = g('log(tauObs' + b + ')'), lo = g('e_log(tauObs' + b + ')'), hi = g('E_log(tauObs' + b + ')'), lr = g('lRest' + b);
        if (isFinite(t) && isFinite(lr)) q.bands.push({ b, t, e: (Math.abs(lo) + Math.abs(hi)) / 2 || 0.1, lr });
      }
      return q; }).filter(q => isFinite(q.z));
  }
  function rows(Q, opts) {                               // design rows: [y, err, x...]
    const R = [];
    for (const [qi, q] of Q.entries()) for (const s of q.bands) {
      if (opts.band && s.b !== opts.band) continue;
      if (opts.needM && !isFinite(q.m)) continue;
      if (opts.needL && !isFinite(q.lb)) continue;
      const x = [1, Math.log10(1 + q.z)];
      if (opts.lam) x.push(Math.log10(s.lr / 4000));
      if (opts.mass) x.push(q.m - 9);
      if (opts.lum) x.push(q.lb - 46);
      R.push({ y: s.t, e: s.e, x, qi });
    }
    return R;
  }
  function solve(A, b) {                                 // Gaussian elimination
    const n = b.length, M = A.map((r, i) => [...r, b[i]]);
    for (let c = 0; c < n; c++) { let p = c; for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
      [M[c], M[p]] = [M[p], M[c]]; for (let r = 0; r < n; r++) if (r !== c) { const f = M[r][c] / M[c][c]; for (let k = c; k <= n; k++) M[r][k] -= f * M[c][k]; } }
    return M.map((r, i) => r[n] / r[i]);
  }
  function wls(R, s2) {
    const k = R[0].x.length, A = Array.from({ length: k }, () => Array(k).fill(0)), b = Array(k).fill(0);
    for (const r of R) { const w = 1 / (r.e * r.e + s2); for (let i = 0; i < k; i++) { b[i] += w * r.x[i] * r.y; for (let j = 0; j < k; j++) A[i][j] += w * r.x[i] * r.x[j]; } }
    return solve(A, b);
  }
  function fit(R) {                                       // profile the intrinsic scatter
    let best = null;
    for (let s = 0.02; s <= 0.6; s += 0.01) { const s2 = s * s, beta = wls(R, s2); let ll = 0;
      for (const r of R) { const m = r.x.reduce((a, x, i) => a + x * beta[i], 0), v = r.e * r.e + s2; ll += -0.5 * ((r.y - m) ** 2 / v + Math.log(v)); }
      if (!best || ll > best.ll) best = { ll, beta, s }; }
    return best;
  }
  function boot(Q, opts, nb, seed) {
    let st = seed >>> 0; const rnd = () => { st = (st * 1664525 + 1013904223) >>> 0; return st / 4294967296; };
    const base = fit(rows(Q, opts)), ns = [];
    for (let k = 0; k < nb; k++) { const Qb = Q.map(() => Q[Math.floor(rnd() * Q.length)]); ns.push(fit(rows(Qb, opts)).beta); }
    const pick = j => { const v = ns.map(b => b[j]).sort((a, b) => a - b); return [v[Math.floor(0.16 * nb)], v[Math.floor(0.5 * nb)], v[Math.floor(0.84 * nb)], v[Math.floor(0.025 * nb)], v[Math.floor(0.975 * nb)]]; };
    return { beta: base.beta, scatter: base.s, n: rows(Q, opts).length, ci: base.beta.map((_, j) => pick(j)), frac_n_le0: ns.filter(b => b[1] <= 0).length / nb };
  }
  return { URL, parse, rows, fit, boot };
})();
// driver used on 11 Nov 2026 (seed 7, 2,000 bootstrap draws):
//   Q = QD.parse(await (await fetch(QD.URL)).text())
//   QD.boot(Q, {}, 2000, 7)                                   // n alone (no controls), all bands
//   QD.boot(Q, {lam:true}, 2000, 7); QD.boot(Q, {lam:true, mass:true, needM:true}, 2000, 7)
//   QD.boot(Q, {lam:true, lum:true, needL:true}, 2000, 7); QD.boot(Q, {lam:true, mass:true, lum:true, needM:true, needL:true}, 2000, 7)
//   per band: add band:'g' | 'r' | 'i'
// two-stage fit used for the entry: stage 1 puts 8 equal-count redshift bins in as fixed effects, so covariate slopes come only
// from quasars at similar redshift; stage 2 regresses the corrected log tau on log(1+z). Fixed-relation runs subtract published slopes.
QD.two = (Q, cov, nbin) => {
  const zs = Q.map(q => q.z).slice().sort((a, b) => a - b), edges = [];
  for (let k = 1; k < nbin; k++) edges.push(zs[Math.floor(k * zs.length / nbin)]);
  const bin = z => edges.filter(e => z >= e).length, R1 = [];
  for (const q of Q) for (const s of q.bands) { const x = []; for (let k = 0; k < nbin; k++) x.push(bin(q.z) === k ? 1 : 0);
    if (cov.includes('lam')) x.push(Math.log10(s.lr / 4000)); if (cov.includes('m')) x.push(q.m - 9); if (cov.includes('L')) x.push(q.lb - 46); R1.push({ y: s.t, e: s.e, x }); }
  const coef = QD.fit(R1).beta.slice(nbin), R2 = [];
  for (const q of Q) for (const s of q.bands) { let y = s.t, j = 0;
    if (cov.includes('lam')) y -= coef[j++] * Math.log10(s.lr / 4000); if (cov.includes('m')) y -= coef[j++] * (q.m - 9); if (cov.includes('L')) y -= coef[j++] * (q.lb - 46);
    R2.push({ y, e: s.e, x: [1, Math.log10(1 + q.z)] }); }
  const f2 = QD.fit(R2); return { coef, n: f2.beta[1] };
};
QD.fixed = (Q, b, c, d) => { const R = []; for (const q of Q) for (const s of q.bands) R.push({ y: s.t - b * Math.log10(s.lr / 4000) - c * (q.m - 9) - d * (q.lb - 46), e: s.e, x: [1, Math.log10(1 + q.z)] }); return QD.fit(R); };
// bootstraps: resample quasars with replacement (500-1000 draws), refit, take the 2.5/16/84/97.5 percentiles of n.
//   QD.two(Q, ['lam'], 8); QD.two(Q, ['lam','m'], 8); QD.two(Q, ['lam','L'], 8); QD.two(Q, ['lam','m','L'], 8); QD.fixed(Q, 0.17, 0.21, 0)
