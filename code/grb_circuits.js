// grb_circuits.js - Workbench Entry 137. Runs in a browser on VizieR J/ApJ/893/77/table1 (Wang et al. 2020, 6,289 GRBs).
// For every burst with redshift, T90 and isotropic energy: the smallest magnetic field a circuit must hold if it stores the
// burst's energy inside the volume light can cross during the burst (a generous size: the whole rest-frame duration).
//   size l = c T90/(1+z);  W = E_iso (1 - cos theta_j), theta_j as published or 0.1 rad;  B_min = sqrt(2 mu0 W / (4/3 pi l^3))
const GRB = (() => {
  const URL = 'https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=J/ApJ/893/77/table1&-out=GRB,z,T90,Eiso,thetaj&-out.max=unlimited&-out.form=|';
  function parse(tsv) {
    const L = tsv.split('\n').filter(l => l && !l.startsWith('#')), H = L[0].split('|').map(s => s.trim());
    return L.slice(3).map(l => { const v = l.split('|').map(s => s.trim()), o = {}; H.forEach((h, i) => o[h] = v[i]); return o; });
  }
  function fields(G) {
    const c = 2.998e8, mu0 = 4e-7 * Math.PI, out = [];
    for (const g of G) { const z = parseFloat(g.z), T = parseFloat(g.T90), E = parseFloat(g.Eiso), th = parseFloat(g.thetaj);   // Eiso in 1e52 erg
      if (!(z > 0 && T > 0 && E > 0)) continue;
      const t = T / (1 + z), l = c * t, V = 4 / 3 * Math.PI * l ** 3, Wiso = E * 1e45, f = isFinite(th) && th > 0 ? 1 - Math.cos(th) : 1 - Math.cos(0.1);
      out.push({ grb: g.GRB, z, t, Wiso, Wcol: Wiso * f, Bcol_G: Math.sqrt(2 * mu0 * Wiso * f / V) * 1e4, Biso_G: Math.sqrt(2 * mu0 * Wiso / V) * 1e4, thMeasured: isFinite(th) }); }
    return out;
  }
  return { URL, parse, fields };
})();
// run: R = GRB.fields(GRB.parse(await (await fetch(GRB.URL)).text())); then quantiles and 0.5-dex histograms of log10 B.
