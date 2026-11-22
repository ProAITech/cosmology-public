# Entry 173: Faraday growth with and without Mg II absorbers

Entry 171 found the scatter of rotation measure (RM) growing with quasar redshift (t = +3.81, 4,867 Milliquas quasars behind SPICE-RACS DR2). Growth is what intergalactic fields or currents would make, and also what intervening galaxies make. This entry asks whether the growth is carried by sight lines with Mg II absorbers, using the DESI DR1 quasar catalogue and DESI's Mg II absorber catalogue.

The rule was frozen in the Workbench (section 173, "the rule, frozen first") in commit `a104830` at 08:02 UTC on 21 November 2026 and published as Workbench version 187. The code and its synthetic tests were committed (`01756bb`) before the run. The section's heading gives 08:10 UTC; the commit is eight minutes earlier.

## Result: unresolved by the frozen rule

| | sight lines | slope b | shuffled null | t |
|---|---|---|---|---|
| D, all DESI quasars behind a SPICE-RACS source | 3,842 | +3.41 | +0.84 ± 1.53 | **+1.68** |
| N, no Mg II absorber in DESI's catalogue | 2,869 | +2.90 | +0.36 ± 1.73 | **+1.47** |
| random subsets of D at N's redshifts (1,000) | 2,869 | +3.52 ± 1.06 | | |
| t_share = (b_N − subsets' mean) / sd | | | | **−0.58** |

- Drop in slope when absorber sight lines are removed: **+0.62 ± 1.06** (95% upper limit +2.36).
- u (test of "absorbers carry at least half of D's growth") = −0.51; the rule needs u ≤ −1.645 and t_D ≥ 2.5.
- Frozen decision: **unresolved**. Removing absorber sight lines does not lower the slope beyond chance (t_share > −2), but growth without them is not established (t_N < 3) and the half-share bound is not reached.

Why it is unresolved: the growth itself is weak in this sample. D gives t = +1.68, where Entry 171 gave +3.81. The two agree within errors (net slope 2.57 ± 1.53 here, 4.63 ± 1.22 there), but with so little growth in D there is little for the absorbers to carry or not carry.

Counts: 1,650,736 DESI quasars (2,182,309 rows, one kept per `TARGETID`), 204,019 sight lines with absorbers in DESI's catalogue. 3,905 quasars lie within 10″ of a qualifying SPICE-RACS component; 3,868 sources keep one quasar each; 26 fall below |b| = 20°; 3,842 remain, 973 with an absorber (25%) and 2,869 without.

Scatter σ (1.4826 × MAD of ΔRM, rad/m²) in the 8 equal-count redshift bins of D, 480–481 sight lines each:

| median z | 0.59 | 0.93 | 1.16 | 1.36 | 1.57 | 1.81 | 2.11 | 2.68 |
|---|---|---|---|---|---|---|---|---|
| share with an absorber | 7% | 15% | 18% | 20% | 28% | 31% | 41% | 44% |
| absorbers per sight line | 0.07 | 0.17 | 0.22 | 0.25 | 0.38 | 0.39 | 0.66 | 0.74 |
| σ, all | 8.56 ± 0.56 | 7.62 ± 0.51 | 8.65 ± 0.56 | 8.49 ± 0.54 | 9.05 ± 0.57 | 9.59 ± 0.57 | 9.20 ± 0.66 | 9.15 ± 0.54 |
| σ, no absorber | 8.15 ± 0.53 | 7.80 ± 0.55 | 8.72 ± 0.59 | 8.76 ± 0.57 | 8.96 ± 0.71 | 8.62 ± 0.59 | 9.45 ± 0.82 | 8.62 ± 0.71 |
| σ, with absorber | 13.6 ± 2.3 | 5.4 ± 1.4 | 8.3 ± 1.4 | 7.1 ± 1.2 | 9.1 ± 1.2 | 11.2 ± 1.1 | 9.2 ± 1.0 | 9.9 ± 0.7 |

## Reported only (item 6 of the rule)

- **(a) Absorber sight lines against the rest.** σ_A = 9.42 ± 0.50 against 8.84 ± 0.24 for N drawn to A's redshifts: an excess variance of +10.5 ± 10.3 rad²/m⁴ (3.2 rad/m² in quadrature), t = +1.45 against random relabelling. Malik et al. (2026, arXiv:2605.16924) report 4.13 ± 0.91 rad/m² from 612 absorber and 1,871 control sight lines, after cuts on Galactic H I and Hα that this rule does not make. The two agree; this sample alone does not detect the excess.
- **(b) Good spectra** (z ≥ 0.4 and `TSNR2_QSO` ≥ 29.8, the lower quartile of A's): 2,604 sight lines, 730 with absorbers. D: b = +7.35, null +1.01 ± 2.01, t = +3.15. N: b = +4.80, t = +1.82. t_share = −1.80, drop +2.18 ± 1.21.
- **(c) Strong absorbers only** (`EW_2796` ≥ 1 Å; 615 sight lines) as the class removed: the rest gives b = +2.95, t = +1.42; t_share = −1.13, drop +0.85 ± 0.76. Their scatter exceeds the rest's by +10.8 ± 14.2 rad²/m⁴.
- **(d) 1° foreground** (3,779 sight lines): D t = +1.89, N t = +1.09, t_share = −0.92, drop +0.81 ± 0.88.
- **(e) By number of absorbers:** none, σ = 8.73 ± 0.20 (2,869, median z 1.35); one, 9.26 ± 0.67 (686, median z 1.71); two or more, 9.34 ± 0.75 (287, median z 2.09). Slope in A alone: +6.41 against +2.30 ± 4.05, t = +1.01.
- **(f) Overlap:** 2,009 of the 3,842 (52%) are also in Entry 171's sample; their DESI and Milliquas redshifts agree within 0.05 for 99.0%.
- **(g)** The drop is 13% of Entry 171's net slope of 4.63, with a 95% upper limit of 51%.

Every variant leans the same way: the slope falls a little when absorber sight lines go (drop +0.6 to +2.2), never by more than 1.8 standard deviations.

## Exploratory, made after the result was seen (`explore.py`, `explore.json`)

Not part of the rule. 300 shuffles each, 1,000 for the union.

| sample | sight lines | b | null | t |
|---|---|---|---|---|
| DESI sight lines also in Entry 171 | 2,009 | +6.36 | +0.34 ± 2.02 | +2.99 |
| DESI-only sight lines | 1,833 | +2.08 | +1.56 ± 2.54 | +0.20 |
| D with `TSNR2_QSO` ≥ 29.8 | 2,643 | +6.22 | +1.04 ± 1.89 | +2.74 |
| D with `TSNR2_QSO` < 29.8 | 1,199 | +1.79 | +0.28 ± 2.95 | +0.52 |
| D with z ≥ 0.4 | 3,753 | +4.46 | +1.01 ± 1.71 | +2.02 |
| Entry 171's sight lines with a DESI quasar | 2,009 | +5.90 | +0.49 ± 2.04 | +2.65 |
| Entry 171's sight lines without one | 2,858 | +4.73 | +0.61 ± 1.59 | +2.59 |
| union of Entry 171 and DESI-only | 6,700 | +4.20 | +0.77 ± 1.14 | +3.02 |

- Entry 171's growth is in both halves of its own sample, with and without a DESI quasar.
- The 1,833 quasars that DESI adds show no growth of their own (net +0.5 ± 2.5), but their error is too wide to contradict Entry 171's 4.6: the difference is 1.6 standard deviations.
- On the union of both samples the growth stands at t = +3.02, down from +3.81. That number is post hoc.

## What it does and does not show

- It does not show that absorbers carry the growth: no variant reaches t_share ≤ −2.
- It does not show that they do not: growth without absorbers is t = +1.47, and the bound on their share is loose (below 51% of Entry 171's slope at 95%, if D shares that slope, which D does not establish).
- It weakens Entry 171 a little. An independent replication was the hope; the DESI-only sight lines neither confirm nor contradict it. The growth needs a test on sight lines that neither entry has used.

## Files

- `absorber_split.py`: the test (`synthetic`, `extract`, `run`). Entry 171's functions are copied in so the folder stands alone.
- `synthetic.json`: tests on synthetic skies, made before the rule was frozen.
- `results.json`: the run of 21 November 2026.
- `explore.py`, `explore.json`: the exploratory table above.

Data, all git-ignored under `data/`:

- `data/spice-racs-dr2/spice-racs.dr2.fits.gz`: SPICE-RACS DR2, https://data.csiro.au/collection/csiro:64891, 4,968,301,262 bytes, sha256 `a1a2f3c1…0053`.
- `data/desi-dr1-qso/QSO_cat_iron_cumulative_v0.fits`: DESI DR1 quasar catalogue, https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/QSO/iron/QSO_cat_iron_cumulative_v0.fits, 816,203,520 bytes, sha256 `a65604303cc29af8865db9193a23934b2540591b16c5afd2a4bea2d61464aef4`, fetched 21 November 2026.
- `data/desi-dr1-mgii/MgII-Absorbers-DR1.fits`: DESI DR1 Mg II absorbers v1.0, https://data.desi.lbl.gov/public/dr1/vac/dr1/mgii-absorber/v1.0/, 65,750,400 bytes.
- `data/milliquas-v8/milliquas.zip`: Milliquas v8, https://quasars.org/milliquas.zip, 32,701,849 bytes, sha256 `9128d6aa…c984` (overlap count and the exploratory table only).

Run from the repository root:

```
.venv/bin/python -I entries/173-spice-desi-absorbers/absorber_split.py synthetic
.venv/bin/python -I entries/173-spice-desi-absorbers/absorber_split.py run       # about 7 minutes
.venv/bin/python -I entries/173-spice-desi-absorbers/explore.py                  # about 5 minutes
```
