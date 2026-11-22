/* Entry 187: read the two catalogues in browser memory (VizieR, CfA mirror), select the galaxies Workbench section 187 names,
   and write them into the page as base64 records for decode.py.  Run in a tab on https://vizier.cfa.harvard.edu after
   pasting healpix.js.  Each galaxy is 6 bytes, big-endian: bits 47-20 the Galactic HEALPix RING index at Nside 4096 (28 bits),
   bits 19-4 round(z * 1e5) (16 bits), bits 3-0 flags: 1 large (Rp90 above the median of all N&A galaxies at 0.02 <= z < 0.05),
   2 dense (rho5 above the median of sample P), 4 the 2MRS control (else Nair & Abraham sample P), 8 unused.
   One line per 120 galaxies: "<line>:<count>:<checksum>:<base64>", checksum = sum over the line of (pix mod 999983) + 7 z5 + 13 flags, mod 1e9+7.
   Catalogues: Nair & Abraham (2010) J/ApJS/186/427/table2 (_RA, _DE, zs, TT, Rp90, rho5);
               2MRS, Huchra et al. (2012) J/ApJS/199/26/table3 (RAJ2000, DEJ2000, GLAT, cz, type, Riso).                     */
window.E187E = (function () {
  const D2R = Math.PI / 180, VZ = '/viz-bin/asu-tsv?-source=';
  const MG = [[-0.0548755604, -0.8734370902, -0.4838350155], [0.4941094279, -0.4448296300, 0.7469822445], [-0.8676661490, -0.1980763734, 0.4559837762]];
  const galvec = (ra, dec) => { const a = ra * D2R, d = dec * D2R, v = [Math.cos(d) * Math.cos(a), Math.cos(d) * Math.sin(a), Math.sin(d)]; return MG.map(r => r[0] * v[0] + r[1] * v[1] + r[2] * v[2]); };
  async function tsv(src, cols) {
    const r = await fetch(VZ + src + '&-out.max=unlimited&-out=' + cols.join(','));
    const L = (await r.text()).split('\n').filter(l => l.trim() && !l.startsWith('#'));
    const head = L[0].split('\t').map(h => h.trim()), rows = [];
    for (const l of L.slice(3)) { const f = l.split('\t'); if (f.length < head.length) continue; const o = {}; head.forEach((h, i) => o[h] = f[i].trim()); rows.push(o); }
    return rows;
  }
  const sexa = s => { s = s.trim(); if (!/[ :]/.test(s)) return parseFloat(s); const sg = /^-/.test(s) ? -1 : 1; const p = s.replace(/^[+-]/, '').split(/[ :]+/).map(Number); return sg * (p[0] + p[1] / 60 + (p[2] || 0) / 3600); };
  const median = a => { const s = Float64Array.from(a).sort(), n = s.length; return n % 2 ? s[(n - 1) / 2] : 0.5 * (s[n / 2 - 1] + s[n / 2]); };
  const st = { info: {} };
  function rec(v, z, flags) { return { pix: ang2pixRing(4096, Math.max(-1, Math.min(1, v[2])), Math.atan2(v[1], v[0])), z5: Math.round(z * 1e5), flags }; }
  async function build() {
    const na = await tsv('J/ApJS/186/427/table2', ['_RA', '_DE', 'zs', 'TT', 'Rp90', 'rho5']);
    const G = []; for (const o of na) { const ra = parseFloat(o._RA), de = parseFloat(o._DE), z = parseFloat(o.zs), T = o.TT === '' ? NaN : parseInt(o.TT, 10), R = parseFloat(o.Rp90), rho = parseFloat(o.rho5);
      if (!isFinite(ra) || !isFinite(de) || !isFinite(z)) continue; const v = galvec(ra, de); G.push({ v, z, T, R, rho, b: Math.asin(v[2]) / D2R }); }
    const inR = g => g.z >= 0.02 && g.z < 0.05;
    const Rmed = median(G.filter(g => inR(g) && isFinite(g.R)).map(g => g.R));
    const P = G.filter(g => inR(g) && Math.abs(g.b) >= 10 && g.T >= 3 && g.T <= 7);
    const rhoMed = median(P.filter(g => isFinite(g.rho)).map(g => g.rho));
    const out = P.map(g => rec(g.v, g.z, (isFinite(g.R) && g.R > Rmed ? 1 : 0) | (isFinite(g.rho) && g.rho > rhoMed ? 2 : 0)));
    const tm = await tsv('J/ApJS/199/26/table3', ['RAJ2000', 'DEJ2000', 'GLAT', 'cz', 'type', 'Riso']);
    let typed = 0, nA = 0;
    for (const o of tm) { const m = /^\s*(-?\d+)/.exec(o.type || ''); const ra = sexa(o.RAJ2000) * (/[ :]/.test(o.RAJ2000.trim()) ? 15 : 1), de = sexa(o.DEJ2000), cz = parseFloat(o.cz), gb = parseFloat(o.GLAT);
      if (!m || !isFinite(ra) || !isFinite(de) || !isFinite(cz)) continue; typed++; const T = parseInt(m[1], 10), z = cz / 299792.458;
      if (!(z >= 0.004 && z <= 0.02) || !(Math.abs(gb) >= 10) || T < 3 || T > 7) continue;
      const Rkpc = Math.pow(10, parseFloat(o.Riso)) * (cz / 100) * 1000 / 206264.8; if (!(Rkpc > 8.5)) continue;
      out.push(rec(galvec(ra, de), z, 4)); nA++; }
    st.info = { na_rows: na.length, na_usable: G.length, Rp90_median: Rmed, P: P.length, P_large: P.filter(g => isFinite(g.R) && g.R > Rmed).length,
      rho5_median_P: rhoMed, P_dense: P.filter(g => isFinite(g.rho) && g.rho > rhoMed).length, twomrs_rows: tm.length, twomrs_typed: typed, A: nA };
    // encode
    const A64 = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/', lines = [];
    for (let i = 0; i < out.length; i += 120) { const ch = out.slice(i, i + 120); let bin = '', cs = 0;
      for (const r of ch) { const V = r.pix * 1048576 + r.z5 * 16 + r.flags; let x = V; const by = []; for (let k = 0; k < 6; k++) { by.unshift(x % 256); x = Math.floor(x / 256); } bin += String.fromCharCode(...by);
        cs = (cs + (r.pix % 999983) + 7 * r.z5 + 13 * r.flags) % 1000000007; }
      lines.push((i / 120) + ':' + ch.length + ':' + cs + ':' + btoa(bin)); }
    st.lines = lines; st.total = out.length; return st.info;
  }
  function show(from, to) { document.body.innerHTML = '<article><pre>' + st.lines.slice(from, to).join('\n') + '</pre></article>'; return st.lines.slice(from, to).length; }
  return { build, show, st };
})();
