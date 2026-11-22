/* Workbench Entry 81: are distant galaxies younger? The 4000 A break at fixed velocity dispersion.
   Both quantities are free of any distance law: D4000n is a ratio of fluxes either side of 4000 A in the rest frame
   (Balogh et al. 1999 narrow definition), and the stellar velocity dispersion is a line width in km/s.
   In a static, stationary universe the galaxy population is the same at every distance, so D4000n at fixed sigma
   must not change with redshift.
   Near: SDSS DR18 (MPA-JHU galSpecIndx), run on https://skyserver.sdss.org/dr18/SearchTools/sql (same-origin SqlSearch).
   Far:  LEGA-C DR3 (van der Wel et al. 2021, VizieR J/ApJS/256/44), run on the VizieR CfA mirror.               */

// ---- near sample (SkyServer tab) ----
const SQL_NEAR = "SELECT TOP 120000 s.z, s.veldisp, s.veldispErr, i.d4000_n, i.d4000_n_err, i.lick_hd_a, i.lick_hd_a_err " +
  "FROM SpecObj s JOIN galSpecIndx i ON s.specObjID = i.specObjID WHERE s.class='GALAXY' AND s.zWarning=0 " +
  "AND s.z BETWEEN 0.07 AND 0.12 AND s.veldisp > 90 AND s.veldispErr < 40 AND i.d4000_n > 0.5 AND i.d4000_n < 3";
async function loadNear() {
  const t = await (await fetch('/dr18/SkyServerWS/SearchTools/SqlSearch?cmd=' + encodeURIComponent(SQL_NEAR) + '&format=csv')).text();
  return t.split('\n').slice(2).filter(l => l.trim().length).map(l => l.split(',').map(Number)).filter(r => r.length === 7 && r.every(Number.isFinite));
} // columns: z, sigma, e_sigma, D4000n, e, HdA, e

// ---- far sample (VizieR tab) ----
async function loadFar() {
  const u = '/viz-bin/asu-tsv?-source=J/ApJS/256/44&-out=Seq,Fprim,Fsp,zsp,sigma%2A,e_sigma%2A,D4000N,e_D4000N,HDA,e_HDA,Tcor,SN4000&-out.max=unlimited&-oc.form=dec';
  const t = await (await fetch(u)).text();
  const rows = t.split('\n').filter(l => l.length && !l.startsWith('#')).slice(3).map(l => l.split('\t').map(s => s.trim() === '' ? NaN : Number(s)));
  // primary targets, no spectroscopy warning, 0.55 < z < 1.05, sigma error < 40 km/s (as for SDSS), valid D4000n
  return rows.filter(r => r[1] === 1 && r[2] === 0 && r[3] > 0.55 && r[3] < 1.05 && r[4] > 0 && r[5] < 40 && r[6] > 0.5 && r[6] < 3);
} // columns: Seq, Fprim, Fsp, z, sigma, e_sigma, D4000n, e, HdA, e, Tcor (completeness weight), SN4000

// ---- statistics in bins of velocity dispersion: weighted quartiles, 90th percentile, bootstrap error of the median,
//      median H-delta_A, weighted fraction with D4000n > 1.85 ----
function STATS(rows, getS, getD, getH, getW) {
  const bins = [[100,150],[150,200],[200,250],[250,300],[300,450]], out = [];
  const wq = (arr, q) => { const a = arr.slice().sort((x,y)=>x[0]-y[0]); const tot = a.reduce((s,v)=>s+v[1],0); let c = 0; for (const [v,w] of a) { c += w; if (c >= q*tot) return v; } return a[a.length-1][0]; };
  for (const [lo,hi] of bins) { const s = rows.filter(r => getS(r) >= lo && getS(r) < hi); if (s.length < 10) { out.push([lo,hi,s.length]); continue; }
    const D = s.map(r => [getD(r), getW(r)]), H = s.filter(r=>Number.isFinite(getH(r))).map(r => [getH(r), getW(r)]);
    const fold = D.filter(v=>v[0] > 1.85).reduce((a,v)=>a+v[1],0)/D.reduce((a,v)=>a+v[1],0);
    const meds = []; for (let b = 0; b < 200; b++) { const bs = []; for (let i = 0; i < D.length; i++) bs.push(D[Math.floor(Math.random()*D.length)]); meds.push(wq(bs, 0.5)); } meds.sort((a,b)=>a-b);
    out.push([lo, hi, s.length, +wq(D,0.25).toFixed(3), +wq(D,0.5).toFixed(3), +wq(D,0.75).toFixed(3), +wq(D,0.9).toFixed(3), +((meds[168]-meds[32])/2).toFixed(3), +wq(H,0.5).toFixed(2), +fold.toFixed(3)]); }
  return out;
} // row: sigma_lo, sigma_hi, n, D4000n p25, median, p75, p90, error of median, median HdA, fraction > 1.85

/* Reported in Entry 81:
   near:  STATS(near, r=>r[1], r=>r[3], r=>r[5], r=>1) for all, z < 0.09 and z >= 0.10 (the larger-aperture comparison)
   far:   STATS(far.filter(r=>r[10]>0), r=>r[4], r=>r[6], r=>r[8], r=>r[10])   (completeness-weighted)
          STATS(far, ..., r=>1) unweighted, and for z < 0.8 and z >= 0.8 separately                                 */
