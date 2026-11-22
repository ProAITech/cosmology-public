# Entry 176: how far from galaxy clusters the magnetised plasma shows

Radio sources seen through galaxy clusters scatter more in rotation measure (RM) than sources in the field. This entry measures that scatter against distance from the nearest eRASS1 cluster, in units of the cluster's R500, with SPICE-RACS DR2, and asks how far out the excess reaches.

The rule was frozen in the Workbench (section 176, "the rule, frozen first") in commit `68ba070` at 10:18 UTC on 21 November 2026 and published as Workbench version 190, before any cluster was matched to a rotation measure. The code and its synthetic test were committed (`17163fd`) before the run.

## Result: a clear excess inside R500, and none beyond it

9,050 eRASS1 clusters and groups (`PCONT` < 0.5, z ≥ 0.05; median z = 0.30, median M500 = 2.1 × 10¹⁴ M☉, median θ500 = 3.0′) and 92,169 SPICE-RACS sources at |b| ≥ 20°. V is the variance of ΔRM with the catalogue noise removed; the null is the same cluster list turned in Galactic longitude 500 times.

| range | sources | V (rad²/m⁴) | rotated clusters | t | excess variance | in quadrature |
|---|---|---|---|---|---|---|
| A: inside R500 | 589 | 179.3 | 73.7 ± 10.8 | **+9.77** | **+105.5** | 10.3 rad/m² |
| B: 1 to 2 R500 | 1,456 | 71.7 | 72.3 ± 7.1 | −0.09 | −0.6 (below +11.1 at 95%) | below 3.3 rad/m² |
| C: 2 to 4 R500 | 5,223 | 71.3 | 71.9 ± 4.7 | −0.13 | −0.6 (below +7.1 at 95%) | below 2.7 rad/m² |

- Frozen decision: **range A, the method sees cluster plasma (t ≥ +3); range B, no excess; range C, no excess; no claim beyond R500.**
- None of the 500 rotations reaches the real value in range A. The second null (ΔRM shuffled among sources within cells of latitude and RM error) gives t = +13.1, +0.1 and −0.3.

## Reported only (item 7 of the rule)

**The profile in finer ranges** (200 rotations):

| x = θ/θ500 | 0–0.5 | 0.5–1 | 1–1.5 | 1.5–2 | 2–3 | 3–4 | 4–6 | 6–8 |
|---|---|---|---|---|---|---|---|---|
| sources | 221 | 368 | 608 | 848 | 2,303 | 2,920 | 7,036 | 7,605 |
| the same, rotated clusters | 117 | 338 | 549 | 754 | 2,034 | 2,602 | 6,241 | 6,711 |
| excess variance (rad²/m⁴) | +174.6 | +79.7 | +5.6 | −3.2 | −2.3 | +0.1 | +0.0 | −1.7 |
| its error | 21.5 | 11.9 | 10.3 | 9.0 | 5.4 | 5.0 | 3.7 | 4.0 |
| t | +8.1 | +6.7 | +0.5 | −0.4 | −0.4 | +0.0 | +0.0 | −0.4 |

The excess falls from 175 to 80 between the two halves of R500 and is gone in the first half-R500 outside it. The inner range holds nearly twice as many sources as the rotated clusters do, so about half of them are the clusters' own radio galaxies.

| check | A: t (excess) | B: t | C: t |
|---|---|---|---|
| north Galactic hemisphere | +5.6 (+92) | +0.0 | +1.1 |
| south Galactic hemisphere | +6.9 (+107) | −0.1 | −1.1 |
| Faraday-simple sources only | +9.2 (+101) | −0.1 | +0.1 |
| more massive half of the clusters | +14.3 (+204) | −0.8 | +0.1 |
| less massive half | +3.5 (+60) | +0.5 | −0.2 |
| nearer half (z below 0.295) | +9.4 (+119) | −0.8 | −0.2 |
| farther half | +1.7 (+44; 96 sources) | +1.1 | −0.4 |
| 1° foreground | +12.6 (+85) | +1.7 | −0.9 |

**Noise and counts.** The mean e_RM² is 5.1, 4.6 and 4.8 rad²/m⁴ in ranges A, B and C, and 4.8 in the field: the noise does not differ. Sources far from any cluster (70,260) give V = 70.7.

## Exploratory, made after the result was seen (`explore.py`, `explore.json`, `figure1.svg`)

Not part of the rule.

| quarter of the clusters by M500 | median M500 (10¹³ M☉) | median z | excess inside R500 | t | excess at 1–2 R500 |
|---|---|---|---|---|---|
| 1 | 8.2 | 0.17 | +31 ± 21 | +1.5 | +5 ± 15 |
| 2 | 16.4 | 0.31 | +107 ± 20 | +5.3 | +7 ± 14 |
| 3 | 25.4 | 0.37 | +176 ± 22 | +8.0 | +4 ± 12 |
| 4 | 44.4 | 0.38 | +234 ± 18 | +12.7 | −9 ± 10 |

- **Mass.** The excess variance inside R500 rises with mass in every step, as M500 to the power 0.86 (0.73–1.00; χ² = 3.2 for 2). In the most massive quarter the inner half of R500 gives +717 ± 41 rad²/m⁴, 27 rad/m² in quadrature, from 41 sources.
- **No mass quarter shows anything at 1–2 R500.**
- **Tails.** Inside R500, 13.2% of sources have |ΔRM| above 30 rad/m², against 5.2 ± 1.2% for rotated clusters (t = +6.6). At 1–2 R500 it is 4.3%, and at 2–4 R500 4.4%, against 5.2% and 5.1%.

**Catalogue cuts** (`variants.py`, `variants.json`; 200 rotations each).

| clusters used | inside R500: t (excess) | 1–2 R500: t | 2–4 R500: t |
|---|---|---|---|
| `PCONT` < 0.1 | +11.1 (+130) | −0.0 | +0.0 |
| z from 0.03 | +10.8 (+124) | −0.2 | +0.3 |
| z from 0.1 | +7.4 (+100) | +0.5 | +0.2 |
| the 2,000 largest on the sky removed | +4.9 (+77; 219 sources) | +1.0 | −0.8 |

**A second sign of the same plasma: sources inside clusters are less polarized.** Plasma that scrambles rotation across the beam lowers a source's fractional polarization. The median fractional polarization is 3.2% within 0.5 R500 against 6.2 ± 0.6% for rotated clusters (t = −5.4). Sources there are also brighter, and brighter sources can be catalogued at lower fractional polarization, so each source was compared with field sources of its own flux decile: the ratio is 0.74 within 0.5 R500 (t = −4.1 against 1.00 ± 0.06), 0.90 at 0.5–1 R500 (t = −2.2), and 0.96, 1.00, 0.97, 0.99 in the ranges from 1 to 4 R500 (|t| ≤ 1.0). Depolarization, like the scatter, stops at R500.

## Cause and field strength (`model.py`, `model.json`), after the frozen result

Not a pre-registered test: the data were already open.

**The excess is not in the clusters' own radio galaxies.** A cluster's own galaxies add to the source count near it; plasma in front of unrelated sources does not. Against the rotated clusters, corrected for the 13% of them that leave the SPICE-RACS footprint:

| x = θ/θ500 | sources | expected without members | surplus | excess variance |
|---|---|---|---|---|
| 0–0.5 | 221 | 132 | +89 (+7.7σ) | +174.6 |
| 0.5–1 | 368 | 382 | −14 (−0.7σ) | +79.7 |
| 1–1.5 | 608 | 620 | −12 | +5.6 |
| 1.5–2 | 848 | 852 | −4 | −3.2 |

Members are confined inside 0.5 R500, where they are about 40% of the sources. Between 0.5 and 1 R500 there is no surplus of sources at all, yet the variance is doubled. That excess is therefore in sight lines to unrelated sources: it is made by the cluster's plasma in front of them.

**The field at 0.75 R500 is about 2 microgauss.** Each cluster is given a gas density n_e(r) = n_0/(1 + r²/r_c²) with r_c = 0.2 R500, normalised to its own gas mass from the eRASS1 catalogue (median n_0 = 2.5 × 10⁻³ cm⁻³; 1.6 × 10⁻⁴ at 0.75 R500), and a random field of reversal length Λ = 10 kpc and strength B ∝ n_e^η. Sources in front of the cluster get no rotation (their share follows an assumed redshift distribution of median 0.8). The model is pushed through the same statistic and fitted to the three ranges from 0.5 to 2 R500.

| η | B at 0.75 R500 (μG) | χ² (2 d.o.f.) | predicted excess at 0.5–1, 1–1.5, 1.5–2 | thermal over magnetic pressure |
|---|---|---|---|---|
| 0 (uniform field) | 1.60 (1.60–1.80); corrected 1.40 (1.40–1.60) | 4.7 | 71, 21, 9 | 17; corrected 23 |
| **0.5** | **2.02 (1.94–2.27); corrected 2.02 (1.70–2.35)** | **0.5** | 78, 10, 2 | **11** |
| 1 | 2.29 (2.13–2.62); corrected 2.38 (1.88–2.87) | 0.2 | 77, 5, 1 | 8 |

*Correction, 21 November 2026 (`model_errors.py`, `model_errors.json`).* The ranges above took their errors from the rotated clusters alone. Where the cluster term dominates, it adds its own sample variance, because path lengths and densities differ widely between sources. Entry 178's synthetic test showed that leaving this out makes intervals too narrow. Refitted with it (Gaussian likelihood, errors depending on the field, the same draws for every field), the 68% ranges are about twice as wide; the 95% range for η = 0.5 is 1.38–2.83 μG. The corrected values are marked "corrected".
| measured | | | 80 ± 12, 6 ± 10, −3 ± 9 | |

- With no sources in front, η = 0.5 gives 1.86 μG. B scales as (Λ/10 kpc)^−½, and the pressure ratio as Λ.
- A field that falls with the gas fits; a uniform field leaves too much beyond R500 (χ² = 4.7 for 2).
- **By mass** (range 0.5–1 R500, η = 0.5):

| quarter | median M500 (10¹³ M☉) | kT500 (keV) | B at 0.75 R500 (μG) | thermal over magnetic pressure |
|---|---|---|---|---|
| 1 | 8.2 | 1.8 | 1.36 (0.06–2.35); corrected 1.05 (0–2.29) | 7; corrected 12 |
| 2 | 16.4 | 2.9 | 1.66 (0.95–2.13); corrected 1.50 (0.87–2.21) | 13; corrected 16 |
| 3 | 25.4 | 4.0 | 2.27 (2.09–2.79); corrected 2.27 (1.66–3.14) | 11; corrected 11 |
| 4 | 44.4 | 5.8 | 2.59 (2.40–2.68); corrected 2.30 (1.82–3.07) | 15; corrected 19 |

  The field rises with mass while the ratio of thermal to magnetic pressure stays near 11–19 (corrected): the field's energy follows the gas pressure.
- **The inner range was left out of the fit, and the fit shows why.** A field that fits 0.5–2 R500 predicts an excess of 500 to 8,000 rad²/m⁴ inside 0.5 R500, where 175 is measured. The sources that would carry it are missing: behind the dense gas they depolarize and leave the catalogue, and those that remain are 0.74 times as polarized as field sources.

**How much the field depends on the model** (`model_variants.py`, `model_variants.json`; η = 0.5, Λ = 10 kpc, fitted to 0.5–2 R500).

| assumption changed | B at 0.75 R500 (μG) |
|---|---|
| standard (core 0.2 R500, β = 2/3, median source redshift 0.8, line of sight to 3 R500) | 2.10 (1.94–2.27) |
| core 0.1 R500 / 0.3 R500 | 2.30 / 1.94 |
| β = 0.6 / 0.75 | 1.83 / 2.26 |
| no sources in front / twice the share | 1.86 / 2.51 |
| median source redshift 0.5 / 1.2 | 2.51 / 1.94 |
| line of sight to 2 R500 / 5 R500 | 2.02 / 2.10 |

Every variant lies between 1.83 and 2.51 μG. In round numbers: **B = 2.0 ± 0.3 (statistical, corrected) ± 0.4 (model) μG × (Λ/10 kpc)^−½ at 0.75 R500.** The reversal length is the one thing the method cannot fix: between 5 and 30 kpc it moves the field from 3.0 to 1.2 μG.

**What the data say about the reversal length** (`coherence.py`, `coherence.json`). Between 0.5 and 1 R500 the rotation measure differs from source to source by 8.9 rad/m². If it varied as much across each source, the sources would keep only 17% of their polarization at these wavelengths (799.5–1087.5 MHz; corrected on 22 November 2026 from 744–1032 MHz, the band of the first data release, which gave 10%). They keep 91 ± 4% of the field value, and at most 14% of them are missing (2σ). So the screen is nearly uniform across a source: the variance inside a source is at most 14–19% of the variance between sources, which for a Gaussian-correlated screen and a disc-shaped source means **Λ is at least 0.75–0.9 times the projected diameter of the polarized emission**. At the clusters' median redshift that is at least 7–8 kpc for emission 2″ across, 17–20 kpc for 5″, and 34–41 kpc for 10″. The source sizes are not measured here, so this is a bound in units of the source size. It makes 10 kpc a floor if the polarized emission is 3″ across or more, and so makes 2.0 μG a ceiling in that case (2.4 μG for emission 2″ across): for Λ of 20 to 50 kpc the field is 1.5 to 0.9 μG and the thermal-to-magnetic pressure ratio 20 to 55.

**Source sizes, and pairs of components** (`coherence2.py`, `coherence2.json`). One more pass over the catalogue kept each component's identifier, size and rotation measure (246,508 components; the beam is 14″).

- *Sizes.* Sources at 0.5–1 R500 have deconvolved major axes of 4.5″, 7.1″ and 14.5″ at the quartiles (12% unresolved), the same as field sources: 8, 17 and 40 kpc at their clusters. If the polarized emission fills the component, the bound above becomes Λ ≥ 16–18 kpc for the median source. The data do not confirm that it does: the larger half of the sources is not more depolarized than the smaller half (0.98 ± 0.05 against 0.90 ± 0.04 of the field value), so the polarized emission may be more compact than the component.
- *Pairs.* 9,649 sources have two components at least 30″ apart, each with its own rotation measure. Their difference cancels the Galactic foreground: away from clusters its noise-removed variance is 12 rad²/m⁴, against 71 for single sight lines.

| range | pairs (expected without members) | separation | variance of the difference | rotated clusters | t |
|---|---|---|---|---|---|
| inside R500 | 80 (46) | 114 kpc | 189.8 | 11.6 ± 7.2 | +24.8 |
| 0.5–1 R500 | 40 (35) | 117 kpc | 120.2 | 11.9 ± 8.7 | +12.4 |
| 1–2 R500 | 158 (133) | 140 kpc | 15.5 | 11.4 ± 4.1 | +1.0 |

  The cluster signal is there again, with the foreground gone, and again stops at R500. None of 200 rotations reaches the real value in the first two ranges. The size of the excess is much less certain than its existence: resampling the pairs gives +164 ± 80 inside R500 and +79 ± 77 at 0.5–1 R500. Inside R500 about 40% of these double sources are cluster members; at 0.5–1 R500 there is no significant surplus.
- *An upper bound on the reversal length.* If the rotation at the two components were unrelated, the variance of the difference would rise by twice the single-sight-line excess of the same components; if the screen were uniform across the pair, it would not rise. The ratio is 0.92 at the median (above 0.38 at 95%), so sight lines 114 kpc apart are at most 62% correlated, and **Λ is below 117 kpc** for a Gaussian correlation.
- *Together.* The reversal length lies between the size of the polarized emission and 117 kpc. With the polarized emission at least 2″ across, that is 7 to 117 kpc, and **the field at 0.75 R500 is between 0.6 and 2.4 μG**, with thermal over magnetic pressure between about 8 and 130. If the polarized emission fills the median component (17 kpc), the field is at most 1.8 μG.

**Against single clusters.** For Coma, Bonafede et al. (2010, arXiv:1002.0594) fit a central field of 4.7 μG falling as n_e^0.5, with fluctuations on scales of 2 to 34 kpc. At the density where this entry measures (1.6 × 10⁻⁴ cm⁻³) that model gives about 1.0 μG if Coma's central density is 3.4 × 10⁻³ cm⁻³ (a value recalled, not checked today). This entry's 2.1 μG assumes Λ = 10 kpc; with Λ near the top of Coma's range, 30 kpc, it is 1.2 μG. The two agree within what Λ allows. Stacks of Planck clusters report mean fields of a few μG (Osinga et al. 2025, arXiv:2408.07178; read in search summary only).

## What it does and does not show

- Magnetised cluster plasma is detected at about ten standard deviations with one survey's rotation measures and one cluster catalogue, in both Galactic hemispheres, under two nulls, and it grows with cluster mass. Sharma, Ravi and Krause (2026, arXiv:2610.37646), stacking SPICE-RACS around DESI groups and ROSAT clusters, reported limits only; the gain here is the size of the eRASS1 catalogue and the scaling of each cluster by its own R500.
- The excess ends at R500 at this sensitivity. Beyond it, the added scatter is below 3.3 rad/m² (1–2 R500) and 2.7 rad/m² (2–4 R500) at 95%. Cluster outskirts and the gas between clusters do not show.
- That bears on Entries 171–174. Sight lines passing within 4 R500 of an eRASS1 cluster but outside R500 carry no extra scatter, so listed clusters' outskirts are not where a growth of scatter with redshift would come from.
- The size inside R500 is a lower bound on the true scatter. At 799.5–1087.5 MHz, sources behind the densest plasma depolarize and leave the catalogue, and the median-based scatter discounts the tail that remains.
- No redshifts are used. The field strength above rests on an assumed share of sources in front of the clusters, a β-model for the gas, and a reversal length of 10 kpc that this entry cannot measure.

## Files

- `cluster_rm.py`: the test (`synthetic`, `run`). The SPICE-RACS reader and merge are Entry 171's.
- `synthetic.json`, `results.json`.
- `explore.py`, `explore.json`, `variants.py`, `variants.json`, `figure.py`, `figure1.svg` (exploratory); `model.py`, `model.json` (members and field strength), `model_variants.py`, `model_variants.json` (its dependence on assumptions), `coherence.py`, `coherence.json`, `coherence2.py`, `coherence2.json` (the reversal length: depolarization, source sizes, component pairs); `draft.md`, a draft manuscript.

Data, git-ignored under `data/`:

- `data/erass1-clusters/erass1cl_primary_v3.2.fits`: eRASS1 clusters and groups (Bulbul et al. 2024, arXiv:2402.08452), https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/BulbulE_DR1/erass1cl_primary_v3.2.fits.tgz, 23,764,095 bytes, sha256 `f9e9e8f750a2831ac68dad1fb24a028b024c00d3b5879f7a5760af81fc9a40ca`, fetched 21 November 2026; unpacked 50,636,160 bytes.
- `data/spice-racs-dr2/spice-racs.dr2.fits.gz`: as in Entry 171.

Run from the repository root:

```
.venv/bin/python -I entries/176-cluster-faraday-profile/cluster_rm.py synthetic
.venv/bin/python -I entries/176-cluster-faraday-profile/cluster_rm.py run        # about 8 minutes
.venv/bin/python -I entries/176-cluster-faraday-profile/explore.py               # about 8 minutes
```

## The X-COP density profile (after the result)

`model_xcop.py`, `model_xcop.json`. The β-model is replaced by the piecewise power-law fit to the scaled density of the X-COP clusters (Ghirardini et al. 2019, Table 2). Everything else is as in `model_variants.py`'s standard case.
- **X-COP shape, normalised to each cluster's eRASS1 gas mass:** B = 2.25 μG at 0.75 R500, with median n_e 1.73 × 10⁻⁴ cm⁻³ there (β-model: 2.02 μG, 1.64 × 10⁻⁴).
- **X-COP absolute density:** 3.84 × 10⁻⁴ cm⁻³ at 0.75 R500, 2.3 times what the gas masses imply, as expected for clusters more massive than ours. It gives B = 0.99 μG.

These fits use the rotation-only errors, so only their central values are quoted.
