# Entry 171: Faraday rotation against quasar redshift, repeated on SPICE-RACS DR2

Does the scatter of the rotation measure (RM) grow with quasar redshift? This is Entry 164's test, which used NVSS rotation measures, repeated on SPICE-RACS DR2. SPICE-RACS has RM errors of about 1–2 rad/m² against NVSS's several to tens, and about three times the density.
The rule was frozen in the Workbench (section 171, "the rule, frozen first") in commit `c53b0e3` at 05:14 UTC on 21 November 2026 and published as Workbench version 185, before any SPICE-RACS rotation measure was matched to a quasar.

## Result: the scatter grows with redshift, and the growth claim passes every deciding check

| | value |
|---|---|
| quasars kept (Milliquas Type Q, z > 0, SPICE-RACS match within 10″, \|b\| ≥ 20°) | 4,867 (5,264 matched, 397 at low latitude) |
| median redshift (range) | 1.40 (0.12–5.25) |
| slope b of the scatter against log₁₀(1+z) | **+5.08** rad/m² per unit (fit error ±1.33, χ² = 2.9 for 6) |
| shuffled-redshift null, 1,000 shuffles in 12 cells | +0.45 ± 1.22; largest of 1,000 below +5.08 (99.87% point +3.66) |
| t | **+3.81** (rule: t ≥ +3 means the scatter grows) |
| deciding checks, t (each must reach +1.5) | north +3.12, south +2.35, smaller-e_RM half +3.03, all classes +2.73, longitude sectors +3.89 |
| frozen decision | **grows with redshift; growth claimed: every check reaches t ≥ 1.5** |

Scatter σ (1.4826 × MAD of ΔRM, rad/m²) in the 8 equal-count bins, 608–609 quasars each:

| mean log₁₀(1+z) | 0.186 | 0.263 | 0.314 | 0.359 | 0.400 | 0.443 | 0.496 | 0.581 |
|---|---|---|---|---|---|---|---|---|
| median z | 0.57 | 0.83 | 1.06 | 1.29 | 1.51 | 1.77 | 2.14 | 2.71 |
| σ | 7.96 | 7.51 | 8.28 | 8.72 | 9.30 | 9.02 | 9.33 | 9.62 |
| bootstrap error | 0.46 | 0.49 | 0.50 | 0.47 | 0.52 | 0.43 | 0.58 | 0.40 |
| median e_RM | 0.72 | 0.98 | 1.10 | 1.08 | 1.27 | 1.30 | 1.32 | 1.42 |
| σ with the mean e_RM² removed | 7.83 | 7.34 | 8.09 | 8.55 | 9.12 | 8.84 | 9.13 | 9.43 |

**Size.** With the catalogue noise removed, the scatter rises from about 7.6 rad/m² at z ≈ 0.6–0.8 to about 9.3 at z ≈ 2.1–2.7. In quadrature that adds about 5.4 rad/m², 28 ± 10 rad²/m⁴ of variance, from first to last bin.

**Reported only (item 8 of the rule).** (a) 1° foreground radius, at least 5 others: 4,770 quasars, b = +5.51, null +0.73 ± 1.03, t = +4.62; the scatter falls to 6.1–7.9 rad/m² and the growth stays. (b) Faraday-simple sources only (no component with `complex_flag`): 4,628, b = +5.41, t = +3.84. (c) Chance matches with every quasar moved 5′ north: 92, 1.7% of the 5,264. (d) The 95% limit (7.08) is 1.05 times Entry 164's (6.75). That ratio is moot, since this test finds growth rather than a limit.

**Entry 164 agrees.** Its flat NVSS result allowed b up to 6.75 at 95%. +5.08 lies inside that, and NVSS's null was three times wider (±3.46 against ±1.22 here). The growth was below what NVSS could see.

## What it does and does not show

The rule says growth is what intergalactic fields or currents would make. It is also what intervening galaxies make: the gas around foreground galaxies, seen as Mg II absorbers, rotates polarization, and more distant quasars have more of them in front. In the same catalogue Malik et al. (2026, arXiv:2605.16924) find 612 DESI DR1 quasar sightlines with Mg II absorbers carrying √(σ²_MgII − σ²_control) = 4.13 ± 0.91 rad/m² more scatter than 1,871 without, about 17 rad²/m⁴ of variance per absorber sightline. Earlier, smaller samples found no growth in sightlines without absorbers (Malik et al. 2018, arXiv:1710.10396).

A rough budget, not part of the rule: 25% of DESI's sightlines carry a detected absorber, about 0.25 per unit redshift path. Between z ≈ 0.6 and 2.7 that adds about 0.5 absorbers per sightline, or about 8 rad²/m⁴ of variance, against the 28 ± 10 measured here. On that budget detected absorbers carry about a third of the growth. Weaker absorbers below DESI's limit could carry more. **The growth is established; its cause is not.** Entry 173 is meant to separate the two with the DESI DR1 Mg II absorbers already on disk.

## Exploratory checks, made after the result was seen (`explore.py`, `explore.json`)

None of these is part of the rule. Each repeats the statistic with 300 shuffles.

| check | n | b | t |
|---|---|---|---|
| single-component sources only | 2,033 | +6.01 | +2.84 |
| multi-component sources only | 2,834 | +3.91 | +2.19 |
| Stokes I above / below median | 2,434 / 2,433 | +6.65 / +3.81 | +2.13 / +2.61 |
| polarized fraction above / below median | 2,434 / 2,433 | +4.60 / +7.87 | +3.23 / +3.09 |
| \|b\| ≥ 40° | 3,202 | +6.95 | **+4.78** |
| shuffles within HEALPix nside-4 cells (about 15°) × e_RM tercile, 387 cells | 4,867 | +5.08 | +4.14 |
| shuffles within cells of Stokes I × polarized fraction × e_RM terciles × \|b\| band, 103 cells | 4,867 | +5.08 | +3.15 |

These exclude the obvious instrumental and selection causes:
- **Merged components.** Low-redshift sources more often merge several components (75% in the first bin, 44% in the last), but both kinds show the growth.
- **Noise.** RM errors grow with z, but too little to matter.
- **Flux and polarized fraction.** Both halves of each show the growth.
- **Sky position.** Finer sky cells keep it.
- **Galactic foreground.** The growth is strongest away from the Galactic plane and with the finer foreground, the opposite of what a residual Galactic foreground would do.

## Deviations from the frozen text

None in substance. Implementation choices the rule left open:
1. Random numbers come from numpy (seed 20261121), not Entry 164's JavaScript generator.
2. Equal-count bin edges use Entry 164's rounding.
3. e_RM terciles are taken over the kept sample, as in Entry 164.
4. Galactic coordinates use Entry 164's matrix and agree with the catalogue's own `l`, `b` to 7 × 10⁻⁶ degrees.
5. With 1,000 shuffles the one-sided p cannot go below 1/1001; t is the stated statistic.

The Milliquas file parsed fully: 1,021,800 records, none unparsed, 860,100 Type Q with z > 0, the ReadMe's count.

## Files and how to run

- `spice_faraday_z.py`: the test. `synthetic` runs the tests on synthetic skies and reads no data (`synthetic.json`: three null skies t = +0.09, +0.58, −0.39; planted growth t = +7.9, +10.4). `extract` streams the SPICE-RACS table once (30 s) into `work/` (git-ignored). `run` writes `results.json`, about 6 minutes. Run from the repository root with `.venv/bin/python -I`.
- `results.json` (14 kB): counts, the primary fit with bins and null histogram, the five deciding checks, the reported-only checks.
- `explore.py`, `explore.json` (4.5 kB): the exploratory checks above, with per-bin profiles of e_RM, flux, polarized fraction, component counts and latitude.

## Data

- SPICE-RACS DR2: `data/spice-racs-dr2/spice-racs.dr2.fits.gz`, CSIRO DAP csiro:64891 (https://data.csiro.au/collection/csiro:64891), 4,968,301,262 bytes, sha256 `a1a2f3c1d12a8a8b3aa8526d74b8e0115054f610e2e4186108f0f068005d0053`. Entry 166's quality cuts: 333,174 rows, merged into 134,400 sources, median source error 1.78 rad/m².
- Milliquas v8: `data/milliquas-v8/milliquas.zip`, https://quasars.org/milliquas.zip, 32,701,849 bytes, sha256 `9128d6aa91354d0ee1c83dde3ab8764a20d542aa82327e4a969a4357ba9ac984`, downloaded 21 November 2026 with the user's permission. It holds `milliquas.txt`, read from the zip; the ReadMe is https://quasars.org/Milliquas-ReadMe.txt.

## What the cloud session should write up

Entry 171 needs these pieces:
- **The result, stated plainly.** The scatter grows with redshift, t = +3.81, and every deciding check passes. Give the size: about 5 rad/m² added in quadrature from z ≈ 0.6 to 2.7. Say that Entry 164 agrees, its limit 6.75 sitting above b = 5.08.
- **The caution.** The growth's cause is not decided: intervening galaxies (Mg II absorbers, Malik et al. 2026) are the ordinary candidate, and the rough budget gives them about a third.
- **A figure.** Scatter against log₁₀(1+z), raw and with the noise removed, with the fit (`primary.bins`), beside the null histogram (`primary.null.hist`) with the real slope marked.
- **The exploratory table,** marked as after the fact.
- **The Gap Map.** Intergalactic Faraday rotation goes from "no growth, limit a factor of ten above expectation" to "growth detected at 3.8σ; source of the growth open (intervening galaxies vs intergalactic medium)".
- **A Ledger sentence,** and the section 171 rule marked as run.
