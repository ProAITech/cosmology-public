# Entry 186: the cluster Faraday signal at 1.4 GHz (NVSS)

**Bottom line: no decision.** The test asks whether the small rotation-measure excess Entry 176 found near cluster centres at 0.9 GHz (+175 rad² m⁻⁴ inside 0.5 R500, against 500 or more predicted by the field fitted further out) is a depolarization selection. That reading predicts a larger excess at 1.4 GHz. NVSS rotation measures around the same 9,050 eRASS1 clusters give +224 ± 135 inside 0.5 R500 and +17 ± 68 inside R500. By the frozen rule the signal is *not seen at 1.4 GHz*, and the inner excess shows *no difference* from SPICE-RACS (t_D = +0.36). NVSS is about six times less sensitive here than SPICE-RACS (null sd 68 against 11 inside R500), so neither outcome contradicts Entry 176.

## Rule

Workbench section 186, "The cluster signal at 1.4 GHz: is the inner deficit depolarization? The rule, frozen first". It was fixed at 19:21 UTC on 21 November 2026 (commit `dc22cee`) and published as Workbench version 197. The synthetic run used simulated rotation measures at the real positions and errors and came first. No real rotation measure was set beside a cluster before the rule was fixed.

**Deviation.** The rule's commit was made on this branch at 19:21 UTC, but it reached GitHub only when the branch was pushed with the results. The record published before the run is the Workbench artifact, version 197 (about 19:22 UTC). The run began at 19:23. The commit is also on main as `c488c74`.

## Data

Read in browser memory from VizieR (vizier.cfa.harvard.edu, `/viz-bin/asu-tsv`). No file was downloaded.
- NVSS rotation measures, Taylor, Stil & Sunstrum (2009), `J/ApJ/702/1230/catalog`: 37,543 sources, of which 28,128 qualify (|b| ≥ 20°, at least 10 others within 3°). Median error 10.7 rad m⁻².
- eRASS1 clusters, Bulbul et al. (2024), `J/A+A/685/A106/emain`: 12,247 rows, 9,050 after Entry 176's cuts (pcont < 0.5, zBest ≥ 0.05, R500 > 0). Median θ500 is 3.0′, as in Entry 176.

## Result

| range (R500) | sources | around rotated clusters | V | null | t | excess (rad² m⁻⁴) | Entry 176, SPICE-RACS |
|---|---|---|---|---|---|---|---|
| < 0.5 | 71 | 33 | 376.8 | 152.5 ± 134.8 | +1.66 | +224 ± 135 | +175 ± 22 |
| 0.5–1 | 57 | 97 | 47.6 | 141.4 ± 77.6 | −1.21 | −94 ± 78 | +80 ± 12 |
| < 1 | 128 | 131 | 159.8 | 142.5 ± 67.7 | +0.25 | +17 ± 68 | +106 ± 11 |
| 1–2 | 249 | 376 | 94.2 | 143.9 ± 37.6 | −1.32 | −50 ± 38 | −1 ± 7 |
| 2–4 | 962 | 1,347 | 141.8 | 142.5 ± 19.2 | −0.04 | −1 ± 19 | −1 ± 5 |

- **Detection (x < 1): not seen at 1.4 GHz** (t = +0.25). The 95% upper bound is 129 rad² m⁻⁴, above SPICE-RACS's 106, so this does not contradict Entry 176.
- **Depolarization: no difference detected.** t_D = (224 − 175)/√(135² + 22²) = +0.36.
- **Reported only.**
  - At 0.5–1 R500: −94 ± 78, against SPICE-RACS's +80 ± 12; a 2.2σ difference.
  - Inside 0.5 R500 the 95% upper bound is 446, below the 500 or more that Entry 176's field model predicts there, by 2.0σ.
  - Counts: footprint-corrected, cluster radio galaxies are about two thirds of the 71 sources inside 0.5 R500.

## Reading

- The synthetic run assumed a field scatter of 70 rad² m⁻⁴. The real NVSS scatter around rotated clusters is 142–153, twice that. NVSS's two-band rotation measures carry wider tails, and the null spread grew accordingly: 135 inside 0.5 R500 against 81 in the synthetic run.
- The test can only detect an inner excess above about 400 at 1.4 GHz. It sees none above about 450 at 95%. That leans, at 2σ, against the field model's inner prediction being visible at 1.4 GHz. NVSS's 45″ beam (about 150 kpc at these clusters) still depolarizes sources behind rotation that varies across them, so the selection is only partly lifted.
- A decisive version needs a 1–2 GHz survey with RM errors near 1 rad m⁻² and a beam of a few arcseconds. POSSUM's band 2 or MeerKAT would do.

## Files

- `nvss_clusters.js`: the browser script; `E186.load()`, then `synthetic()` and `run()`.
- `results.json`.
- `figure.py`, `figure1.svg`: the excess against distance from the nearest cluster, NVSS beside SPICE-RACS.
