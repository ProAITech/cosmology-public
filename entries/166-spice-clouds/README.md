# Entry 166: currents along clouds in SPICE-RACS DR2

Does the Faraday-rotation flip across Orion A and Perseus stay pinned to each cloud's curved spine out to tens of parsecs, as a current along the cloud would make it, or does it die out beyond a few parsecs, as a field bent around the cloud would?
Rule frozen in the Workbench (section 166, "the rule, frozen first") in commit `1955eae` on 14 November 2026, before the file was opened. Two amendments, both before any statistic was computed, are in the same section (commits `c069c9f`, `9a58d95`).

## Result: not decided

By the rule as frozen and amended, **Orion A is unresolved and Perseus is not pinned**. Neither cloud is "current-like" or "bent-field-like".

| | Orion A (420 pc) | Perseus (300 pc) |
|---|---|---|
| sources in window | 556 | 476 |
| oriented amplitudes a_in, a_mid, a_out (rad/m²; zones 0–1.5°, 1.5–5°, 5–10°) | +13.5, +38.5, +53.2 | +13.8, +38.2, +35.0 |
| shifted-spine null, a_out (mean ± sd, 78 positions) | +18.9 ± 25.5 | +19.3 ± 27.9 |
| p_out, upper tail (current-like needs ≤ 0.05) | 0.10 | 0.15 |
| segments pinned (best step within ±1° of the spine, same sign) | 3 of 4 | 3 of 4 |
| p_pin (needs ≤ 0.05) | 0.013 | 0.063 |
| frozen decision | unresolved | not pinned |

- **Pinned, in Orion A.** The best step lies at +0.75°, −0.50° and −0.75° in three segments (the fourth, with the fewest sources, at +4.25°); no shifted position out of 78 reaches three same-sign pinned segments, so p = 1/79 = 0.013. In Perseus the three pinned segments are at −1.0°, −0.25° and −0.75°, the fourth at −1.5°, half a degree outside the rule's window; p = 0.063 misses the 0.05 line.
- **Reaches tens of parsecs, but not distinguishably.** The oriented amplitude stays positive and large out to 5–10° (37–73 pc for Orion A, 26–52 pc for Perseus), as a current predicts. But shifted copies of the same spine give +19 ± 26 and +19 ± 28 rad/m² there, because the sky carries antisymmetric structure of that size at that scale; the real values are 1.3σ and 0.6σ above the null mean. The formal errors (±8 and ±9 rad/m²) are 3 times too small. The profile (`profile_binned` in `results.json`) is a smooth antisymmetric curve crossing zero on the spine, about −40 rad/m² at −4.5° and +50 at +6.5° for Orion A, falling again toward ±10°.
- **Replication of Entries 161–162 (check a).** The model of Entry 161 on the straight spines gives the NVSS sign for both clouds (positive toward higher b) and the best step within 1.5° of NVSS's: Orion A A_in = +11.7 ± 3.7, A_out = +33.9 ± 3.7 (NVSS +9.7 ± 9.7, +23.8 ± 6.6), best step −0.50° (NVSS +0.25°); Perseus +14.1 ± 3.5, +39.5 ± 3.6 (NVSS +20.2 ± 8.6, +34.7 ± 9.5), best step −0.75° (NVSS −1.0°). Formal errors, so smaller than the sky allows.
- **Curved spine against straight chord (check b).** No preference. Orion A fits slightly better with the straight chord (χ² improvement 93.9 curved, 99.7 straight; its curved-minus-straight difference of −5.8 is exceeded at 90% of the shifted positions), Perseus slightly better with the curve (+5.4, exceeded at 14%). The data cannot say the flip follows the curve.
- **Orientation (check c).** The true orientation of the spine is the best of 13 rotations (−90° to +90°) for both clouds. *Exploratory, added after the results were seen and not part of the frozen rule:* against the same count made at the 78 shifted positions, rank 1 occurs with p = 0.051 (Orion A) and 0.063 (Perseus); the excess over the best rotated orientation has p = 0.038 and 0.051. Suggestive of a flip aligned with the cloud's long axis, not decisive, and chosen after the fact.
- **Background (check d).** Adding quadratic terms keeps a_out positive (+54.3 and +31.2 rad/m²; p_out 0.11 and 0.20).

What this means: the SPICE-RACS sample replicates the sign and position of the NVSS flips, and for Orion A the flip is pinned to the spine along its length (p = 0.013, about 2.2σ one-sided). It does not show whether the flip persists to tens of parsecs, because the large-scale RM structure of the sky at |b| ≈ 19° makes amplitudes of this size common. The current-versus-bent-field question is still open.

## Deviations, in the order they happened

1. **Table structure.** The file is one SPICE-RACS DR2 table (9,294,225 rows, `catalog_name` = SPICE-RACS DR2 for all), so no survey selection was needed. 333,463 rows have all four quality flags clear (`snr_flag`, `leakage_flag`, `channel_flag`, `stokesI_fit_flag` all False; True = flagged, as 96% of rows have `snr_flag` True); 333,174 of them have a finite RM and error. They are merged per `source_id` (weighted mean; 72,752 sources have more than one row) into 134,400 sources. The first read of the file after the header was one pass that kept columns and counted the values of the survey and flag columns (`meta.json`); no source was placed near a cloud.
2. **Amendment 1 (c069c9f).** The nearest-point coordinates of item 3 broke down far from a short, kinked polyline; replaced by chord-frame coordinates (s along the chord, d = n − h(t)). Found by testing on synthetic skies, before any real analysis.
3. **Amendment 2 (9a58d95).** The first run stopped at the gate: 556 and 476 sources, under the 1,000 required. The table has 3.7 sources per square degree, not the 6.7 I assumed (from earlier Workbench text) when setting the threshold. Threshold lowered to 400 after the counts were known and before any statistic was computed; `results_gate.json` is that first output.
4. **Implementation choices the rule left open:** σ² re-estimated as max(0, Σr²/(N − p) − mean e²) in 6 passes; MAD is the unscaled median absolute deviation; a merged source takes the position of its first row; the null shifts the frame centre, keeping the polyline's shape in the flat frame; a step fit needs at least 5 sources on each side. The rule's flat frame (x = Δl cos b₀) differs from true tangent-plane geometry by about 5–10% in the zone amplitudes at 10° from the centre (an independent astropy implementation gave Orion A [12.9, 36.4, 48.0] and Perseus [13.9, 37.3, 34.1] against [13.5, 38.5, 53.2] and [13.8, 38.2, 35.0]); the flat frame is what was frozen and is used.
5. **Exploratory orientation null** (above), added after the results were seen.

The near-identical near-zone amplitudes of the two clouds (a_in 13.5 and 13.8, a_mid 38.5 and 38.2) are in the data: the independent implementation reproduces them.

## Files and how to run

- `extract.py`: streams the table from the gzip once (about 50 s) and writes `work/compact.npz` (git-ignored, 632 MB) and `meta.json`. Run from the repository root: `.venv/bin/python -I entries/166-spice-clouds/extract.py`.
- `spice_clouds.py`: the test (about 15 s): `.venv/bin/python -I entries/166-spice-clouds/spice_clouds.py`. Deterministic; no random numbers.
- `results.json` (14 kB): the results above, with the null values of a_out, the segments, the checks and a binned profile for a figure. `results_gate.json`: the first run's output, which stopped at the gate.
- `meta.json`: counts of the table's survey, instrument and flag columns.

## Data

SPICE-RACS DR2, CSIRO Data Access Portal collection csiro:64891 (https://data.csiro.au/collection/csiro:64891), file `data/spice-racs-dr2/spice-racs.dr2.fits.gz`, 4,968,301,262 bytes, sha256 `a1a2f3c1d12a8a8b3aa8526d74b8e0115054f610e2e4186108f0f068005d0053`. One BINTABLE, 9,294,225 rows × 1030 bytes, 125 columns. Used: `l`, `b`, `rm`, `rm_err`, `source_id` and the four flags.

## What the cloud session should write up

Entry 166 with the result above stated plainly (not decided; Orion A pinned at p = 0.013 but the outer zone is not distinguishable from the sky's structure; Perseus misses the pinning line at p = 0.063; NVSS sign and position replicated; no preference for the curved spine), both amendments, the 3.7 rather than 6.7 sources per square degree (the Workbench's "seven times denser" in Entry 73 and in section 166's introduction should read about four times), a figure of the binned residual RM against signed distance for each cloud with the zone amplitudes and the null band, the Gap Map update (the cloud-current claim stays open, with a flip replicated and tied to the spine but not shown to reach tens of parsecs), a Ledger sentence, and the section 166 rule marked as run. What would decide it: a denser grid, or a model of the Galactic RM structure good enough to remove its large-scale antisymmetric part.
