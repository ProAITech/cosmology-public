/* Workbench Entry 79: the event table behind gw_redshift.py.
   Run in a browser tab on https://gwosc.org/eventapi/html/GWTC/ (same-origin fetch); queried 11 November 2026.
   The cumulative GWTC list holds 391 candidates (GWTC-1, 2.1, 3, 4.1 and 5.0); 259 are binary black holes
   (m2_source > 3 Msun) with FAR < 1/yr, or with no FAR listed (the ten O1-O2 events, carried from GWTC-2.1 without one).
   GWOSC gives source-frame chirp masses and the redshift LVK computed from the distance; multiplying back by (1+z)
   recovers the detector-frame chirp mass, which is what the waveform measures.                                        */
const r = await fetch('/eventapi/json/GWTC/'); const GW = await r.json();
const ev = Object.values(GW.events);
const ok = e => e.chirp_mass_source != null && e.redshift != null && e.luminosity_distance != null;
const sel = ev.filter(e => ok(e) && ((e.far != null && e.far < 1) || (e.far == null && (e.p_astro == null || e.p_astro > 0.5)))
                       && e.mass_2_source != null && e.mass_2_source > 3).sort((a, b) => a.GPS - b.GPS);
// observing run from GPS time: 1 O1, 2 O2, 3 O3a, 4 O3b, 5 O4a, 6 O4b
const run = g => g < 1140000000 ? 1 : g < 1190000000 ? 2 : g < 1255000000 ? 3 : g < 1300000000 ? 4 : g < 1392000000 ? 5 : 6;
const ROWS = sel.map(e => [run(e.GPS), (e.chirp_mass_source*(1+e.redshift)).toFixed(2), Math.round(e.luminosity_distance),
  Math.round(-e.luminosity_distance_lower), Math.round(e.luminosity_distance_upper), e.network_matched_filter_snr.toFixed(1),
  (e.mass_2_source/e.mass_1_source).toFixed(2), e.redshift.toFixed(3)].join(',')).join(';');
// checksum against gw_rows.txt: 259 rows, sum of detector chirp masses 9617.72, sum of distance columns 1413685
ROWS;
