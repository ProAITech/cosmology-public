// qso_align.js - Workbench Entry 132. Runs in a browser on the VizieR table J/A+A/441/915/qsos (Hutsemekers et al. 2005, 355 quasars).
// Alignment statistic S_D (Jain, Narain & Sarala 2004), as used by Hutsemekers et al.: for each quasar, the n_v nearest neighbours
// (itself included) have their polarization angles parallel-transported to its position; D_i = |mean exp(2i theta')|; S_D = mean D_i.
// Significance from shuffling the angles over the positions. Neighbours in 3D comoving space (flat LCDM, Om 0.3, H0 70) or on the sky.
const QA = (() => {
  const rad = Math.PI / 180;
  function parse(tsv) {                                   // '|' separated asu-tsv output
    const L = tsv.split('\n').filter(l => l && !l.startsWith('#')).slice(3).map(l => l.split('|').map(s => s.trim()));
    const hms = s => { const [h, m, x] = s.split(/\s+/).map(Number); return 15 * (h + m / 60 + x / 3600); };
    const dms = s => { const sg = s.startsWith('-') ? -1 : 1; const [d, m, x] = s.replace(/[+-]/, '').split(/\s+/).map(Number); return sg * (d + m / 60 + x / 3600); };
    return L.map(r => ({ name: r[1], z: +r[2], p: +r[3], ep: +r[4], pa: +r[5], epa: +r[6], ref: +r[7], ra: hms(r[9]), dec: dms(r[10]) }));
  }
  function Dc(z) { let s = 0, n = 400; for (let i = 0; i < n; i++) { const zz = (i + 0.5) * z / n; s += 1 / Math.sqrt(0.3 * (1 + zz) ** 3 + 0.7); } return 2997.92458 / 0.7 * s * z / n; }
  function vec(ra, dec) { return [Math.cos(dec * rad) * Math.cos(ra * rad), Math.cos(dec * rad) * Math.sin(ra * rad), Math.sin(dec * rad)]; }
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  function basis(ra, dec) {                              // north and east unit vectors at (ra, dec)
    const N = [-Math.sin(dec * rad) * Math.cos(ra * rad), -Math.sin(dec * rad) * Math.sin(ra * rad), Math.cos(dec * rad)];
    const E = [-Math.sin(ra * rad), Math.cos(ra * rad), 0]; return [N, E];
  }
  function rotate(v, k, ang) {                           // Rodrigues
    const c = Math.cos(ang), s = Math.sin(ang), kv = cross(k, v), kd = dot(k, v);
    return [0, 1, 2].map(j => v[j] * c + kv[j] * s + k[j] * kd * (1 - c));
  }
  function offset(a, b) {                                // PA rotation picked up by transporting a tangent vector from a to b
    const ra = vec(a.ra, a.dec), rb = vec(b.ra, b.dec), ax = cross(ra, rb), sn = Math.hypot(...ax);
    const [Na] = basis(a.ra, a.dec), [Nb, Eb] = basis(b.ra, b.dec);
    if (sn < 1e-12) return 0;
    const t = rotate(Na, ax.map(x => x / sn), Math.atan2(sn, dot(ra, rb)));
    return Math.atan2(dot(t, Eb), dot(t, Nb));          // radians
  }
  function prepare(Q, nvMax, mode) {
    const P = Q.map(q => { const r = vec(q.ra, q.dec), d = Dc(q.z); return { u: r, x: r.map(c => c * d) }; });
    return Q.map((q, i) => {
      const dist = Q.map((o, j) => [mode === '3d' ? Math.hypot(P[i].x[0] - P[j].x[0], P[i].x[1] - P[j].x[1], P[i].x[2] - P[j].x[2]) : Math.acos(Math.min(1, dot(P[i].u, P[j].u))), j]);
      dist.sort((a, b) => a[0] - b[0]);
      const nb = dist.slice(0, nvMax).map(d => d[1]);
      return { nb, off: nb.map(j => offset(Q[j], q)), r: dist.slice(0, nvMax).map(d => d[0]) };
    });
  }
  function SD(prep, ang, nv) {                           // ang: polarization angles in radians, indexed by position
    let s = 0;
    for (let i = 0; i < prep.length; i++) {
      let c = 0, sn = 0; const { nb, off } = prep[i];
      for (let k = 0; k < nv; k++) { const t = 2 * (ang[nb[k]] + off[k]); c += Math.cos(t); sn += Math.sin(t); }
      s += Math.hypot(c, sn) / nv;
    }
    return s / prep.length;
  }
  function shuffle(a, rnd) { const b = a.slice(); for (let i = b.length - 1; i > 0; i--) { const j = Math.floor(rnd() * (i + 1)); [b[i], b[j]] = [b[j], b[i]]; } return b; }
  function rng(seed) { let s = seed >>> 0; return () => { s = (s * 1664525 + 1013904223) >>> 0; return s / 4294967296; }; }
  // significance for several n_v at once, with the look-elsewhere correction from the same shuffles
  function test(Q, mode, nvs, nsh, seed) {
    const prep = prepare(Q, Math.max(...nvs), mode), ang = Q.map(q => q.pa * rad), rnd = rng(seed);
    const obs = nvs.map(nv => SD(prep, ang, nv));
    const sims = []; for (let s = 0; s < nsh; s++) { const a = shuffle(ang, rnd); sims.push(nvs.map(nv => SD(prep, a, nv))); }
    const p = nvs.map((nv, j) => (1 + sims.filter(r => r[j] >= obs[j]).length) / (nsh + 1));
    const pmin = Math.min(...p);
    // global p: how often a shuffle's best n_v is at least as extreme as the data's best
    const rankP = sims.map((r, s) => Math.min(...nvs.map((nv, j) => (1 + sims.filter(q => q[j] >= r[j]).length) / (nsh + 1))));
    const pglob = (1 + rankP.filter(x => x <= pmin).length) / (nsh + 1);
    return { mode, n: Q.length, nvs, obs, p, pmin, pglob, medR: nvs.map(nv => median(prep.map(x => x.r[nv - 1]))) };
  }
  const median = a => { const b = a.slice().sort((x, y) => x - y); return b[Math.floor(b.length / 2)]; };
  return { parse, test, prepare, SD, offset, Dc };
})();

// ---- shuffle counts reusing one neighbour table (what the runs below use) ----
QA.count = (Q, mode, nvs, nsh, seed) => {
  const rad = Math.PI / 180, prep = QA.prepare(Q, Math.max(...nvs), mode), ang = Q.map(q => q.pa * rad);
  let s = seed >>> 0; const rnd = () => { s = (s * 1664525 + 1013904223) >>> 0; return s / 4294967296; };
  const obs = nvs.map(nv => QA.SD(prep, ang, nv)), ex = nvs.map(() => 0), sims = nvs.map(() => []);
  for (let k = 0; k < nsh; k++) {
    const a = ang.slice(); for (let i = a.length - 1; i > 0; i--) { const j = Math.floor(rnd() * (i + 1)); [a[i], a[j]] = [a[j], a[i]]; }
    nvs.forEach((nv, j) => { const v = QA.SD(prep, a, nv); sims[j].push(v); if (v >= obs[j]) ex[j]++; });
  }
  const mu = sims.map(x => x.reduce((a, b) => a + b) / x.length), sd = sims.map((x, j) => Math.sqrt(x.reduce((a, b) => a + (b - mu[j]) ** 2, 0) / (x.length - 1)));
  return { obs, ex, nsh, sig: obs.map((o, j) => (o - mu[j]) / sd[j]) };   // sig: excess over the shuffles in their own standard deviations
};
// ---- Galactic interstellar polarization from Heiles (2000), VizieR II/226 ----
// per quasar: mean Stokes q,u of stars within 10 deg, beyond 200 pc, polarization detected at 2 sigma, angles transported to the quasar
QA.ism = (Q, S, R = 10, dmin = 200) => {
  const rad = Math.PI / 180, v = (ra, dec) => [Math.cos(dec * rad) * Math.cos(ra * rad), Math.cos(dec * rad) * Math.sin(ra * rad), Math.sin(dec * rad)];
  const St = S.filter(s => s.p >= 2 * s.ep && s.p > 0 && s.d >= dmin), Sv = St.map(s => v(s.ra, s.dec)), cr = Math.cos(R * rad);
  return Q.map(q => { const u = v(q.ra, q.dec); let qq = 0, uu = 0, n = 0;
    St.forEach((s, k) => { if (Sv[k][0] * u[0] + Sv[k][1] * u[1] + Sv[k][2] * u[2] >= cr) { const th = s.pa * rad + QA.offset(s, q); qq += s.p * Math.cos(2 * th); uu += s.p * Math.sin(2 * th); n++; } });
    return n >= 3 ? { q: qq / n, u: uu / n, n } : null; });
};
// subtract f times the ISM vector (rotated by phi[k] radians in Stokes space for the random-orientation control)
QA.correct = (Q, ism, f, phi) => Q.map((q, i) => i).filter(i => ism[i]).map((i, k) => {
  const rad = Math.PI / 180, m = ism[i], ph = phi ? phi[k] : 0, q = Q[i];
  const mq = m.q * Math.cos(ph) - m.u * Math.sin(ph), mu = m.q * Math.sin(ph) + m.u * Math.cos(ph);
  const Qs = q.p * Math.cos(2 * q.pa * rad) - f * mq, Us = q.p * Math.sin(2 * q.pa * rad) - f * mu;
  let pa = 0.5 * Math.atan2(Us, Qs) / rad; if (pa < 0) pa += 180; return { ...q, pa, thIsm: 0.5 * Math.atan2(m.u, m.q) };
});
// driver (as run on 11 Nov 2026 in the browser; seeds as listed)
//   Q = QA.parse(<asu-tsv of J/A+A/441/915/qsos, '|' separated>), S = Heiles table (RAJ2000, DEJ2000, Pol, e_Pol, PA, Dist)
//   QA.count(Q,'3d',[10,20,30,50,80],6000,23); QA.count(Q,'2d',[10,20,30,50,80],3000,29)
//   ism = QA.ism(Q,S); for f in [0,0.5,1,1.5,2,3]: C = QA.correct(Q,ism,f); mean cos 2(pa - thIsm); QA.count(C,'3d',[30,50],500,77)
//   controls: QA.correct(Q,ism,f,random phases) with f = 1 (10 runs, 600 shuffles) and f = 2 (8 runs, 400 shuffles)
