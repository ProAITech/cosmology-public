/* flare_circuits.js - Workbench Entry 139
   Solar flares as exploding double layers (Alfven & Carlqvist 1967): can one circuit current, crossing one double layer,
   supply the electrons a flare's hard X-rays require?
   Runs in a browser console on the CfA VizieR mirror (the two catalogues are read with asu-tsv queries, nothing is saved):
     J/ApJ/832/27  Aschwanden et al. 2016, Global energetics of solar flares III: nonthermal energies of 191 M/X flares (RHESSI, warm-target fits)
     J/ApJ/797/50  Aschwanden et al. 2014, Global energetics I: magnetic free and dissipated energies, length scale L, for 172 M/X flares
   The two papers share one event list, so flares are matched by sequence number (class and start date checked: 75 of 75 agree).
   Per flare:
     electrons      N_e >= E_nt / (2 E_c)              mean energy at most twice the cutoff for any spectrum steeper than delta = 3
     rate, current  I_e >= e N_e / t_SXR               soft X-ray duration, longer than the impulsive phase, so a lower bound
     power          P >= E_nt / t_SXR
     circuit        I = sqrt(2 E_mag / L_ind),  L_ind = mu0 R [ln(8R/a) - 7/4],  R = L (length scale), a = R/10 (R/3 and R/30 as checks),
                    E_mag = max(free, dissipated) magnetic energy: the largest current the stored energy allows for that size
     one layer      V = P / I  (the voltage one double layer must hold to pass the power with the circuit current)
   Compared with measured non-neutralized currents of flaring regions (Kontogiannis et al. 2017): 0.5-15 x 10^12 A.
*/
async function get(src, cols) {
  const r = await fetch('https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=' + src + '&-out=' + cols + '&-out.max=2000&-out.form=|');
  const L = (await r.text()).split('\n').filter(l => l && !l.startsWith('#'));
  const head = L[0].split('|').map(s => s.trim());
  return L.slice(3).map(l => l.split('|').map(s => s.trim())).filter(r => r.length === head.length && /^\d+$/.test(r[0]))
          .map(r => Object.fromEntries(head.map((h, i) => [h, r[i]])));
}
const e = 1.602176634e-19, keV = 1.602176634e-16, mu0 = 4e-7 * Math.PI;
function q(a, p) { const s = [...a].sort((x, y) => x - y); const i = (s.length - 1) * p, lo = Math.floor(i), hi = Math.ceil(i); return s[lo] + (s[hi] - s[lo]) * (i - lo); }
async function run() {
  const A = await get('J/ApJ/797/50', 'Seq,Fl,Cl,Pos,Ep,Efree,Ediss,Pdiss,L,T');
  const N = await get('J/ApJ/832/27', 'Seq,Obs,Cl,Pos,Dur,PkCts,LowE,UpE,Cutoff,Ewt,ERatio1,ERatio2');
  const mA = Object.fromEntries(A.map(r => [r.Seq, r]));
  const R = [];
  for (const r of N) {
    const Ewt = parseFloat(r.Ewt) * 1e-7, Ec = parseFloat(r.Cutoff), Dur = parseFloat(r.Dur);   // J, keV, s
    if (!(Ewt > 0 && Ec > 0 && Dur > 0)) continue;
    const Ne = Ewt / (2 * Ec * keV), Ie = e * Ne / Dur, P = Ewt / Dur;
    const o = { seq: r.Seq, cl: r.Cl, Ewt, Ec, Dur, Ne, Ie, P };
    const a = mA[r.Seq];
    if (a && a.Cl === r.Cl) {
      const E = Math.max(parseFloat(a.Efree), parseFloat(a.Ediss)) * 1e30 * 1e-7, Rr = parseFloat(a.L) * 1e6;
      if (E > 0 && Rr > 0) {
        for (const [k, asp] of [['10', 10], ['3', 3], ['30', 30]]) o['I' + k] = Math.sqrt(2 * E / (mu0 * Rr * (Math.log(8 * asp) - 1.75)));
        Object.assign(o, { Emag: E, R: Rr, V: P / o.I10, Vratio: P / o.I10 / (Ec * 1e3), nratio: Ie / o.I10 });
      }
    }
    R.push(o);
  }
  const lg = Math.log10, m = R.filter(o => o.I10), X = R.filter(o => o.cl.startsWith('X')), M = R.filter(o => o.cl.startsWith('M'));
  const st = a => [0.05, 0.5, 0.95].map(p => +lg(q(a, p)).toFixed(2));
  const hX = {}, hM = {};
  R.forEach(o => { const b = (Math.floor(lg(o.Ie) * 4) / 4).toFixed(2); const H = o.cl.startsWith('X') ? hX : hM; H[b] = (H[b] || 0) + 1; });
  return {
    n: R.length, nM: M.length, nX: X.length, n_matched: m.length,
    log_Ie_5_50_95: st(R.map(o => o.Ie)), log_Ie_med_M: +lg(q(M.map(o => o.Ie), 0.5)).toFixed(2), log_Ie_med_X: +lg(q(X.map(o => o.Ie), 0.5)).toFixed(2),
    log_Ie_min_max: [+lg(Math.min(...R.map(o => o.Ie))).toFixed(2), +lg(Math.max(...R.map(o => o.Ie))).toFixed(2)],
    n_above_1p5e13: R.filter(o => o.Ie > 1.5e13).length,
    log_Ne_med: +lg(q(R.map(o => o.Ne), 0.5)).toFixed(2), Ec_med_keV: +q(R.map(o => o.Ec), 0.5).toFixed(2), Ec_min_max: [Math.min(...R.map(o => o.Ec)), Math.max(...R.map(o => o.Ec))],
    Dur_med_s: Math.round(q(R.map(o => o.Dur), 0.5)), log_Ent_erg_med: +lg(q(R.map(o => o.Ewt * 1e7), 0.5)).toFixed(2), log_P_W_5_50_95: st(R.map(o => o.P)),
    I10_A_5_50_95: [0.05, 0.5, 0.95].map(p => +q(m.map(o => o.I10), p).toPrecision(3)), log_I3_5_50_95: st(m.map(o => o.I3)), log_I30_5_50_95: st(m.map(o => o.I30)),
    R_m_med: +q(m.map(o => o.R), 0.5).toPrecision(3), log_Emag_J_5_50_95: st(m.map(o => o.Emag)),
    V_V_5_50_95: [0.05, 0.5, 0.95].map(p => +q(m.map(o => o.V), p).toPrecision(3)), log_VoverEc_5_50_95: st(m.map(o => o.Vratio)), log_VoverEc_min: +lg(Math.min(...m.map(o => o.Vratio))).toFixed(2),
    log_IeoverI_5_50_95: st(m.map(o => o.nratio)), log_IeoverI_min: +lg(Math.min(...m.map(o => o.nratio))).toFixed(2),
    log_IeoverI_measured: { X_vs_1p5e13: +lg(q(X.map(o => o.Ie / 1.5e13), 0.5)).toFixed(2), M_vs_5e12: +lg(q(M.map(o => o.Ie / 5e12), 0.5)).toFixed(2), all_vs_1p5e13: +lg(q(R.map(o => o.Ie / 1.5e13), 0.5)).toFixed(2) },
    log_V_measured: { X_at_1p5e13: +lg(q(X.map(o => o.P / 1.5e13), 0.5)).toFixed(2), M_at_5e12: +lg(q(M.map(o => o.P / 5e12), 0.5)).toFixed(2) },
    Ent_over_Emag_med: +q(m.map(o => o.Ewt / o.Emag), 0.5).toFixed(3), n_Ent_gt_Emag: m.filter(o => o.Ewt > o.Emag).length,
    hist_log_Ie_quarter_dex: { X: hX, M: hM }
  };
}
run().then(r => console.log(JSON.stringify(r)));
