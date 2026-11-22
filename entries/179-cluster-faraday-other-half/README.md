# Entry 179: the cluster Faraday signal on the other half of the sky

**Result: replicated, by the frozen rule (Workbench section 179, commit `d467f31`, 17:13 UTC, 21 November 2026).**

Entry 176 found the rotation-measure scatter of SPICE-RACS sources raised behind eRASS1 clusters, which cover Galactic longitudes 180–360°. Here the same statistic is applied at longitudes below 179.5° to 1,101 clusters that eRASS1 did not select: 829 from MCXC (ROSAT X-ray) and 272 from Planck PSZ2 (SZ-selected, no MCXC counterpart). The sight lines and the clusters are new; only the radio survey and the statistic are shared.

| range | sources (rotated clusters) | V (rad²/m⁴) | null | t | excess |
|---|---|---|---|---|---|
| A: x < 1 | 80 (96) | 203.0 | 61.6 ± 20.0 | **+7.08** | +141.4 |
| B: 1–2 | 244 (279) | 70.6 | 60.9 ± 11.8 | +0.82 | +9.7 |
| C: 2–4 | 902 (1,074) | 54.5 | 60.1 ± 7.4 | −0.77 | −5.7 |

- Decision: A *replicated* (t ≥ +3); B and C *no excess* (|t| < 2).
- Shuffled null (ΔRM shuffled within latitude × error cells): t = +8.2, +1.2, −0.6.
- **Amplitude.** Entry 176's mass relation, 176 × (M500/2.54 × 10¹⁴ M☉)^0.86, predicts 152 rad²/m⁴ at the median mass of the clusters holding the sources inside R500 (2.15 × 10¹⁴ M☉). 141 ± 20 is measured, a ratio of 0.93.
- **Profile** (reported only):

| x | sources | excess (rad²/m⁴) | t | Entry 176 |
|---|---|---|---|---|
| 0–0.5 | 40 | +227.8 ± 41.4 | +5.5 | +174.6 ± 21.5 |
| 0.5–1 | 40 | +79.1 ± 19.7 | +4.0 | +79.7 ± 11.9 |
| 1–1.5 | 116 | +4.0 ± 18.1 | +0.2 | +5.6 ± 10.3 |
| 1.5–2 | 128 | +13.7 ± 15.5 | +0.9 | −3.2 ± 9.0 |
| 2–3 | 397 | −6.3 ± 9.9 | −0.6 | −2.3 ± 5.4 |
| 3–4 | 505 | −5.2 ± 9.7 | −0.5 | +0.1 ± 5.0 |

- MCXC clusters alone: t = +9.1 inside R500 (69 sources, excess +180). PSZ2 clusters without an MCXC counterpart hold only 11 sources inside R500, too few for the statistic.
- Clusters: median z 0.20, median M500 2.5 × 10¹⁴ M☉, median θ500 4.5′.

**Reading.** Entry 176's detection does not depend on the eRASS1 catalogue or on the half of the sky it covers. Clusters selected by ROSAT and by Planck, on sight lines Entry 176 never used, show the same raised scatter inside R500 at 7σ, nothing beyond R500, the same profile, and the amplitude Entry 176's mass relation predicts. Together the two entries detect cluster Faraday rotation at 0.5–1 R500 in about 10,000 clusters with no surplus of member sources (Entry 176), on both halves of the sky.

## Files

- `other_half.py`: `synthetic` and `run`.
- `cluster_rm.py`: Entry 176's, unchanged.
- `synthetic.json`: no effect gives t = −0.8; 12 rad/m² added inside R500 gives t = +9.8.
- `results.json`.
- Data:
  - `data/mcxc/table.tsv`: VizieR J/A+A/534/A109/mcxc, 551,285 bytes, sha256 9ff800289c01547bce757b48af21cffbf40585139ab3d9ce553962b70d2e4948.
  - `data/psz2/HFI_PCCS_SZ-union_R2.08.fits.gz`: https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/HFI_PCCS_SZ-union_R2.08.fits.gz, 155,324 bytes, sha256 afba0c3e2b7341579aa5b0212dd0ee4ae15ba19493329e2a4827d17db52053a5.
  - SPICE-RACS DR2 as Entry 176.

## The mass scaling on this half of the sky (after the result, not pre-registered)

`mass_halves.py`, `mass_halves.json`. The 1,101 clusters are split at their median M500 (2.55 × 10¹⁴ M☉), and each half is tested alone with 500 longitude rotations.

| half | clusters | sources inside R500 | excess (rad² m⁻⁴) | t | Entry 176's relation predicts | ratio |
|---|---|---|---|---|---|---|
| lower mass (median 1.4 × 10¹⁴ M☉ for the clusters holding sources) | 550 | 44 | +113 ± 27 | +4.2 | 106 | 1.07 |
| upper mass (4.0 × 10¹⁴ M☉) | 551 | 36 | +199 ± 30 | +6.7 | 260 | 0.77 |

- The excess rises with mass: upper minus lower is +86 ± 40 (2.2σ), against +154 predicted.
- The implied exponent is 0.54 ± 0.27, against eRASS1's 0.86 (0.73–1.00): consistent at about 1σ.
- The shortfall in the upper half is in the direction the depolarization selection predicts, which is stronger in massive clusters.
