// Entry 182: region A1 on quasars measured after 2005 (Workbench section 182, rule frozen 21 Nov 2026 18:08 UTC, commit ff019f0).
// Run in a browser tab on https://vizier.cfa.harvard.edu (rows stay in browser memory).
// Sources: VizieR J/A+A/606/A101/table4 (Hutsemekers, Hall & Sluse 2017, EFOSC2 La Silla, 192 quasars) and
//          VizieR J/A+A/620/A68/table3 (Hutsemekers, Borguet, Sluse & Pelgrims 2018, FORS VLT, 86 quasars).
// Part 1 loads both tables and returns name, position and redshift only. Python (foreground.py) returns, for each row, the
// region flags, |b|, the 10" match to the 355 quasars of J/A+A/441/915 and the Planck foreground (q, u) in percent.
// Part 2 applies the frozen cuts, the correction and the binomial count.

async function part1() {
  const cols = 'Name,RAJ2000,DEJ2000,z,Filter,q,u,Pol,e_Pol,PA,e_PA';
  const rows = [];
  for (const [src, tag] of [['J/A%2bA/606/A101/table4', 'LS2017'], ['J/A%2bA/620/A68/table3', 'VLT2018']]) {
    const r = await fetch('/viz-bin/asu-tsv?-source=' + src + '&-out.max=unlimited&-out=' + cols);
    const lines = (await r.text()).split('\n').filter(l => l.trim() && !l.startsWith('#'));
    const head = lines[0].split('\t').map(h => h.trim());
    for (const l of lines.slice(3)) {
      const f = l.split('\t'); if (f.length < head.length - 1) continue;
      const o = { cat: tag }; head.forEach((h, i) => o[h] = (f[i] || '').trim()); rows.push(o);
    }
  }
  window.__a1 = rows;
  const sex = (s, h) => { const sg = s.startsWith('-') ? -1 : 1; const [a, b, c] = s.replace(/^[+-]/, '').split(/[ :]+/).map(Number); return sg * (a + b / 60 + c / 3600) * (h ? 15 : 1); };
  return rows.map((o, i) => [i, +sex(o.RAJ2000, true).toFixed(5), +sex(o.DEJ2000, false).toFixed(5), +o.z]);
}

function part2(fg) {   // fg: {i: [inA1, inA3, absb, match355, fq, fu, group]}
  const rows = window.__a1, D = Math.PI / 180;
  const inSec = (a, s) => s[0] > s[1] ? (a >= s[0] || a <= s[1]) : (a >= s[0] && a <= s[1]);
  const lgam = x => { const c = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5]; let y = x, t = x + 5.5; t -= (x + 0.5) * Math.log(t); let s = 1.000000000190015; for (const cc of c) s += cc / ++y; return -t + Math.log(2.5066282746310005 * s / x); };
  const binSF = (k, n, p) => { let s = 0; for (let j = k; j <= n; j++) s += Math.exp(lgam(n + 1) - lgam(j + 1) - lgam(n - j + 1) + j * Math.log(p) + (n - j) * Math.log(1 - p)); return Math.min(1, s); };
  const cpUpper = (k, n) => { if (k === n) return 1; let lo = k / n, hi = 1; for (let it = 0; it < 60; it++) { const m = (lo + hi) / 2; (1 - binSF(k + 1, n, m) >= 0.05 ? lo = m : hi = m); } return (lo + hi) / 2; };   // one-sided 95% upper bound: P(X <= k | m) = 0.05
  const corr = (p, pa, fq, fu, f) => { const q = p * Math.cos(2 * pa * D) - f * fq, u = p * Math.sin(2 * pa * D) - f * fu; let a = 0.5 * Math.atan2(u, q) / D; if (a < 0) a += 180; return [Math.hypot(q, u), a]; };
  const out = {};
  for (const [reg, idx, sec, mu] of [['A1', 0, [146, 46], 6], ['A3', 1, [84, 164], 124]]) {
    // candidates: in region, |b| >= 30, not within 10" of the 355, p >= 0.6, e_PA <= 14
    const cand = {};
    rows.forEach((o, i) => {
      const g = fg[i]; if (!g || !g[idx] || g[2] < 30 || g[3]) return;
      const p = +o.Pol, ePA = +o.e_PA, pa = +o.PA; if (o.Pol === '' || o.e_PA === '' || o.PA === '' || !(p >= 0.6) || !(ePA <= 14) || isNaN(pa)) return;
      const key = g[6]; if (!cand[key] || ePA < cand[key].ePA) cand[key] = { i, name: o.Name, cat: o.cat, p, ePA, pa, fq: g[4], fu: g[5], z: +o.z, filt: o.Filter };
    });
    const L = Object.values(cand), n = L.length;
    const res = { n, quasars: [] };
    for (const [tag, f] of [['measured', 0], ['f1', 1], ['f2', 2]]) {
      let k = 0, c = 0;
      for (const x of L) { const [pp, aa] = corr(x.p, x.pa, x.fq, x.fu, f); if (inSec(aa, sec)) k++; c += Math.cos(2 * (aa - mu) * D); if (tag === 'f1') x.pa1 = +aa.toFixed(1); }
      res[tag] = { k, P: n ? binSF(k, n, 80 / 180) : null, V: n ? Math.sqrt(2 * n) * c / n : null };
    }
    res.f1.upper95 = n ? cpUpper(res.f1.k, n) : null;
    if (reg === 'A1') {
      const r = res.f1;
      res.decision = n < 10 ? 'not tested' : r.P <= 0.00135 ? 'A1 confirmed on new quasars' : r.upper95 < 0.675 ? 'not repeated' : r.P <= 0.05 ? 'supported' : 'unresolved';
      res.combined_with_2005 = { k: r.k + 27, n: n + 40, P: binSF(r.k + 27, n + 40, 80 / 180) };
    }
    res.quasars = L.map(x => ({ name: x.name, cat: x.cat, z: x.z, filter: x.filt, p: x.p, e_PA: x.ePA, fg_p: +Math.hypot(x.fq, x.fu).toFixed(3) }));
    out[reg] = res;
  }
  return out;
}
