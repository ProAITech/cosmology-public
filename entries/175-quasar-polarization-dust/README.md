# Entry 175: quasar polarization alignments against Planck's dust polarization

Hutsemékers and colleagues reported that the optical polarization angles of quasars line up over regions about a gigaparsec across. Entry 132 reproduced the alignment in their 355 quasars, found the angles linked to the polarization of nearby stars, and pre-registered this test: "with Planck's 353 GHz dust polarization as the foreground template, the full sample staying above 3σ in 3D would reopen cosmic order". This entry runs it.

The rule was frozen in the Workbench (section 175, "the rule, frozen first") in commit `7f866dd` at 09:47 UTC on 21 November 2026 and published as Workbench version 189, before the dust map was opened. The code and its synthetic tests were committed (`939b11d`) before the run, and the convention check on stars (`stars_check.json`) before any quasar was corrected.

## Result: Galactic, by the frozen rule

The foreground is what Planck's dust map predicts for starlight, with no fitted number: (q, u) = −(Q, −U) × 287.45 / 5.4, turned to the equatorial frame. Its median is 0.16% (quartiles 0.10–0.30%), against a median quasar polarization of 1.37%.

| foreground removed | σ for n = 10, 20, 30, 50, 80 neighbours | global p (10,000 shuffles) | S |
|---|---|---|---|
| none (f = 0) | 3.20, 3.60, 4.24, 4.23, 3.33 | 0.0006 | 3.24σ |
| **Planck's prediction (f = 1)** | 1.77, 1.85, 2.38, 2.51, 1.74 | 0.032 | **1.85σ** |
| the level that ends the link (f* = 2.05) | 0.50, 0.73, 0.87, 0.60, 0.56 | 0.36 | 0.37σ |

- Frozen decision: **Galactic. The foreground Planck predicts removes the alignment (S(1) < 2σ).**
- The f = 0 row reproduces Entry 132 (3.2–4.2σ), from a Python port of its browser code.
- **The link.** Before any correction the quasar angles follow the direction Planck predicts for dust-polarized starlight: mean cos 2Δ = +0.230 ± 0.038, six standard deviations. It falls to +0.150 at f = 1 and crosses zero at f* = 2.05.

## The convention check on stars (item 3 of the rule), passed before the quasars

- The code's Galactic-to-equatorial turn reproduces Heiles's two angle columns: median difference 0.26°, within 2° for 99.9% of 8,065 stars.
- For 144 stars at |b| ≥ 30°, beyond 300 pc, detected at 3σ: mean cos 2Δ between the star's angle and Planck's predicted angle is **+0.72 ± 0.06**. The stars are polarized 1.2 times the prediction in the median, so the conversion holds at high latitude to about 20%.
- For 2,407 distant stars at low latitude: +0.72 ± 0.01.

## Reported only (item 7 of the rule)

- **Random directions.** The same vectors, each turned by a random angle before subtraction, 10 trials: at f = 1 the n = 30 alignment stays at 3.56–4.73σ (true direction: 2.38σ); at f* it stays at 2.45–4.11σ (true direction: 0.87σ). The drop comes from the foreground's direction, not from added noise.
- **Cut applied again.** After the f = 1 correction 29 quasars fall below p = 0.6%. Without them (326 quasars): 1.2–1.9σ, global p = 0.098 (1.30σ).
- **Sky only.** Neighbours on the sky, redshift ignored: 3.00σ with nothing removed, 1.47σ at f = 1, 0.78σ at f*.
- **Link by polarization.** +0.29 ± 0.07 for p = 0.6–1.0% (114 quasars), +0.23 ± 0.07 for 1.0–1.6% (95), +0.18 ± 0.09 for 1.6–3% (66), +0.19 ± 0.08 above 3% (80).

## Exploratory, made after the result was seen

Not part of the rule.

**The removal level scanned** (`explore.py`, `explore.json`, `figure1.svg`). The alignment falls smoothly as the predicted foreground is removed: n = 30 gives 4.26σ at f = 0, 3.37 at 0.5, 2.38 at 1, 1.12 at 1.5, 0.76 at 1.75, 0.85 at 2, then rises again (2.28 at f = 3) as over-subtraction lines the angles up the other way. Random-direction trials stay at 3.2–4.7σ up to f = 1.5. A second set of 10,000 shuffles puts the global significance at 3.06σ for f = 0 and 1.84σ for f = 1; the frozen run's values are 3.24σ and 1.85σ.

**How firmly the result sits below 2σ** (`variants.py`, `variants.json`). The averaging aperture does not matter: 0.25°, 0.5°, 1° and 2° give 1.83σ, 1.82σ, 1.83σ and 1.84σ at f = 1 (one common set of 10,000 shuffles). The conversion ratio does: it is known to about 10%, and f = 0.9, 1.0 and 1.1 give 2.04σ, 1.82σ and 1.61σ. So the rule's 2σ line is crossed at a foreground about 8% weaker than Planck's prediction. The stars point the other way: at high latitude they are polarized 1.2 times the prediction, and f = 1.2 gives 1.42σ.

**Can dust of the predicted size make all of this?** (`simulate.py`, `simulate.json`) At the 355 real positions, with the real predicted foreground and the real measurement errors, synthetic quasars get random angles and a log-normal intrinsic polarization, and are redrawn until they pass p ≥ 0.6%. There is no alignment in the simulation except what the foreground and the cut make. 120 samples per case:

| parent median p | foreground | link before | link after f = 1 | f* | alignment, n = 30 (σ) | after f = 1 (σ) | median observed p |
|---|---|---|---|---|---|---|---|
| **the data** | | **0.230** | **0.150** | **2.05** | **4.24** | **2.38** | **1.37%** |
| 0.5% | Planck's | 0.247 ± 0.039 | 0.129 | 1.72 ± 0.21 | 4.45 ± 1.73 | 1.48 ± 1.42 | 1.08% |
| 0.8% | Planck's | 0.184 ± 0.034 | 0.083 | 1.61 ± 0.25 | 2.74 ± 1.50 | 0.70 ± 1.11 | 1.32% |
| 0.3% | Planck's | 0.312 ± 0.034 | 0.184 | 1.84 ± 0.17 | 6.64 ± 1.68 | 2.85 ± 1.39 | 0.94% |
| 0.5% | twice Planck's | 0.409 ± 0.036 | 0.326 | 2.79 ± 0.20 | 9.56 ± 1.92 | 6.84 ± 1.80 | 1.14% |
| 0.5% | none | 0.001 ± 0.039 | −0.115 | 0.28 ± 0.21 | 0.01 ± 1.07 | 1.22 ± 1.23 | 1.08% |

- A foreground of Planck's size acting on a sample cut at 0.6% makes an alignment of 3–6σ and a link of 0.18–0.31 with no cosmic order at all. The data sit inside that range.
- The cut is why the link ends near twice the predicted level and why about 1.5–2.4σ remains after the correction: both appear in the simulation.
- A foreground twice Planck's is excluded: it would give a link of 0.41 and an alignment near 10σ.
- With the parent that matches the observed median polarization (0.8%), every simulated number is 1 to 1.8 standard deviations below the data. The data are a little more aligned than this simple parent gives; a parent with more weakly polarized quasars and a tail of strongly polarized ones would sit between the rows.

**The link among strongly polarized quasars** (`simulate_bands.py`, `simulate_bands.json`). The data keep a link of +0.19 ± 0.08 for the 80 quasars above 3%, which looked like more than an added 0.16% could give. The simulation says otherwise: in 400 samples each, the band above 3% shows a link of +0.14 ± 0.16 (parent median 0.5%), +0.12 ± 0.11 (0.8%) and +0.08 ± 0.09 (a two-part parent, 80% with median 0.4% and 20% with median 2.5%), and 37%, 25% and 12% of samples reach the data's value. The foreground is not 0.16% everywhere: where it is large it both lifts a quasar's observed polarization and sets its angle. The two-part parent reproduces all four bands (+0.29, +0.27, +0.17, +0.08 against +0.29, +0.23, +0.18, +0.19) and the overall link (0.228 against 0.230).

## What it does and does not show

- The gigaparsec alignment in this sample does not survive the removal of the Galactic foreground that Planck predicts. Entry 132 reached the same verdict from stars with a fitted level; here nothing is fitted, and all 355 quasars are used.
- A residual of 1.85σ remains at f = 1, and the link ends only at twice the predicted level. The simulation above shows the sample's cut does both. In words: The sample was cut at p ≥ 0.6%, and a weakly polarized quasar is more likely to pass the cut when the foreground adds to it than when it subtracts, so the survivors lean toward the foreground's angle even after the true foreground is removed. And a pencil beam can be more polarized than the 80′ average.
- The link of +0.19 ± 0.08 among quasars polarized above 3%, which Entry 132 also saw and which first looked unexplained, is within what the simulation gives (+0.08 to +0.14; 12–37% of samples reach the data). Nothing in this sample is left that the foreground and the cut do not account for.
- Pelgrims (2019, arXiv:1709.10271) used Planck's 353 GHz map on the same quasars and concluded that dust shows at 2σ in about 30% of them and cannot fully account for the alignments. The numbers here do not contradict the first part. On the second, what remains after the predicted foreground is removed is below 2σ, which is not evidence for an alignment. Only his abstract was read.
- The 2014 claim that quasar polarizations line up with the large quasar groups that host them is a different claim and is not tested.

## Files

- `qso_dust.py`: `synthetic`, `stars`, `run`.
- `synthetic.json`, `stars_check.json`, `results.json`.
- `explore.py`, `explore.json`, `figure1.svg`; `variants.py`, `variants.json`; `simulate.py`, `simulate.json`, `simulate_bands.py`, `simulate_bands.json` (exploratory).
- `draft.md`: a draft manuscript.

Data, all git-ignored under `data/`:

- `data/planck-pr3-gnilc-dust/COM_CompMap_IQU-thermaldust-gnilc-unires_2048_R3.00.fits`: Planck PR3 GNILC thermal dust at 353 GHz, 80′, https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/component-maps/foregrounds/COM_CompMap_IQU-thermaldust-gnilc-unires_2048_R3.00.fits, 1,811,949,120 bytes, sha256 `2d68599b73aa5034340c85eed89a5270b68aaa33dfd52a286b3a2c7320491f88`.
- `data/hutsemekers-2005/qsos.tsv`: VizieR J/A+A/441/915/qsos, https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=J/A%2BA/441/915/qsos&-out.max=unlimited&-out.all, 40,046 bytes, sha256 `19573d69c6e2dee160b1662d6c1d322c6601287a4eb723db9ca177e58a183692`.
- `data/heiles-2000/catalog.tsv`: VizieR II/226/catalog, same mirror, 1,981,311 bytes, sha256 `092fb3d15ae0485de090976cf63ba467609ce7069ae123c542662a94a0c22200`.

All three fetched 21 November 2026.

Run from the repository root:

```
.venv/bin/python -I entries/175-quasar-polarization-dust/qso_dust.py synthetic
.venv/bin/python -I entries/175-quasar-polarization-dust/qso_dust.py stars
.venv/bin/python -I entries/175-quasar-polarization-dust/qso_dust.py run        # about 2 minutes
```

## The original regions A1 and A3 (exploratory, after the result)

`regions.py`, `regions.json`. Hutsemékers et al. (2005, Sect. 3.1 and Table 1, read from arXiv:astro-ph/0507274) defined two regions and sectors:
- **A1:** B1950 RA 11h15m–14h29m, z 1.0–2.3, angles expected in 146°–46°;
- **A3:** 21h20m–24h00m, z 0.7–1.5, angles expected in 84°–164°.

They tested each with a binomial count. Their "new sample" leaves out the 1998 quasars from which the regions were chosen. The script reproduces all four of their counts:
- A1: 42/56 full, 27/40 new;
- A3: 24/29 full, 14/18 new.

The new sample is r_Pol ≥ 7, less the re-measured B1222+228 and B1246−057.

| region, sample | n | in sector: measured | f = 1 | f* | dust alone | dust-only share ≥ measured | V toward sector centre: measured, f = 1 |
|---|---|---|---|---|---|---|---|
| A1 new | 40 | 27 (P 0.0028) | 27 | 27 | 16.5–17.0 ± 3.1 | 0.1–0.2% | 1.99, 2.18 |
| A1 full | 56 | 42 (3.3e−6) | 41 | 41 | 23.0–23.6 ± 3.6 | 0 of 2,000 | 3.41, 3.62 |
| A3 new | 18 | 14 (0.0043) | 12 (0.049) | 11 (0.12) | 10.3–11.1 ± 2.0 | 6–11% | 2.62, 2.19 |
| A3 full | 29 | 24 (2.6e−5) | 22 | 21 | 16.4–17.4 ± 2.6 | 0.3–0.9% | 4.36, 3.84 |

- **A3 is Galactic.** The foreground (median 0.21%, mean angle 130°, 93% of its angles in the sector), the 17 Heiles stars in the patch (mean 132°) and the quasars (126–132°) share one direction. The correction brings the new sample down to what dust alone gives.
- **A1 is not Galactic.** The foreground (0.14%, mean 58°, 16% in the sector) and the 50 stars in the patch (mean 58°) point 50° away from the quasars (8–18°). Dust alone puts fewer quasars in the sector than uniform angles would. The correction leaves the new-sample count at 27/40 and raises the directional statistic from V = 1.99 to 2.18 (P = 0.015). A sector-free Rayleigh test gives P = 0.09.
- **Reading.** The frozen full-sample decision (Galactic) stands. Of the two regions defined in advance, one (A1, high redshift, where the 2005 paper saw the strongest alignment) is not explained by dust, at P ≈ 0.003 (binomial) to 0.015 (directional) on its independent quasars. That is short of 3σ.
- **A decisive test needs new data.** If the 27/40 rate holds, 60 new A1 quasars would reach 3σ with probability 0.71, and 100 with probability 0.95.
