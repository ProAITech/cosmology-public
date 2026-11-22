# Entry 174: Faraday growth on sight lines not used before (Quaia)

Entry 171 found the scatter of rotation measure (RM) growing with quasar redshift, t = +3.81, on 4,867 Milliquas quasars behind SPICE-RACS DR2. Entry 173 found the growth weak in DESI's quasars (t = +1.68). This entry tests the growth on sight lines neither entry used: SPICE-RACS sources behind Quaia quasars (Gaia + unWISE, G < 20.5, all sky), after removing every source that carries a Milliquas object or a DESI DR1 quasar.

The rule was frozen in the Workbench (section 174, "the rule, frozen first") in commit `9a42e43` at 08:28 UTC on 21 November 2026 and published as Workbench version 188. The code and its synthetic tests were committed (`d00c57a`) before the run. The section's heading gives 08:50 UTC; the commit is earlier.

## Result: unresolved by the frozen rule, narrowly

| | value |
|---|---|
| Quaia quasars within 10″ of a qualifying SPICE-RACS component | 10,058, on 10,046 sources |
| removed: source carries a Milliquas object or a DESI quasar | 5,049 |
| removed: Galactic \|b\| < 20° | 1,248 |
| new sight lines kept | **3,749** (median z 1.20, median G 19.6) |
| slope b | **+3.02** (χ² = 9.8 for 6) |
| shuffled-redshift null, 1,000 shuffles | −0.09 ± 1.67 |
| t | **+1.87** (one-sided p = 0.031) |
| net slope, and its 95% upper limit | +3.11, +5.85 |
| Entry 171's net slope seen through Quaia redshifts (4.63 × β, β = 0.835) | 3.87 |
| frozen decision | **unresolved**: t < 2, and Entry 171's slope is not excluded |
| combination with Entry 171, (3.81 + 1.87)/√2 (reported, decides nothing) | 4.01 |

The rule needed t ≥ 2 to call the growth supported; the new sight lines give +1.87. Their slope has the sign of Entry 171's and four-fifths of its size (+3.11 against an expected 3.87), with an error of 1.67. So the independent sample neither confirms the growth by the rule nor contradicts it.

Scatter σ (1.4826 × MAD of ΔRM, rad/m²) in the 8 equal-count bins, 468–469 sight lines each:

| median z (Quaia) | 0.45 | 0.74 | 0.91 | 1.09 | 1.33 | 1.58 | 1.95 | 2.58 |
|---|---|---|---|---|---|---|---|---|
| σ | 8.44 | 10.99 | 8.69 | 9.98 | 9.23 | 9.65 | 9.88 | 10.27 |
| bootstrap error | 0.61 | 0.71 | 0.52 | 0.57 | 0.60 | 0.67 | 0.58 | 0.64 |

## Reported only (item 6 of the rule)

| check | sight lines | b | null | t |
|---|---|---|---|---|
| north Galactic hemisphere | 1,113 | +11.26 | +0.31 ± 3.54 | +3.09 |
| south Galactic hemisphere | 2,636 | +1.99 | +0.14 ± 1.90 | +0.97 |
| Dec < −20° | 2,306 | +2.42 | −0.44 ± 2.18 | +1.32 |
| Dec ≥ −20° | 1,443 | +2.08 | +0.40 ± 2.85 | +0.59 |
| smaller redshift-error half | 1,875 | +5.81 | −0.45 ± 2.18 | +2.88 |
| G < 20 | 2,632 | +1.26 | −0.47 ± 1.99 | +0.87 |
| 1° foreground, at least 5 others | 3,709 | +3.07 | −0.10 ± 1.32 | +2.41 |
| Entry 171's sight lines that Quaia has, Milliquas z | 3,659 | +6.13 | +0.33 ± 1.50 | +3.85 |
| the same sight lines, Quaia z | 3,659 | +4.61 | +0.06 ± 1.36 | +3.35 |
| **union: Entry 171 + DESI-only + new** | **10,449** | **+3.57** | **+0.16 ± 0.96** | **+3.54** |

- **Chance matches.** With every Quaia quasar moved 5′ north, 170 match a source, 153 of them on unused sources: about 3% of the 4,997 unused matches.
- **Quaia's redshifts.** On the 3,659 sight lines shared with Entry 171, 67.7% of Quaia's redshifts lie within 0.01 (1+z) of Milliquas's and 88.6% within 0.2 (1+z). Replacing Milliquas's redshifts by Quaia's on those sight lines lowers the net slope from 5.80 to 4.55, a factor 0.78, close to the β = 0.835 the rule uses.
- **The union.** Entry 171's 4,867 sight lines, the 1,833 DESI-only ones and the 3,749 new ones, each with its own redshift: 10,449 sight lines, slope +3.57 against +0.16 ± 0.96, t = +3.54; none of 1,000 shuffles reaches the real slope. The scatter rises bin by bin: 8.18, 8.89, 8.80, 8.98, 9.05, 9.26, 9.60, 9.68 rad/m² (errors about 0.33) from median z = 0.52 to 2.65. With the catalogue noise removed, the variance grows by about 19 rad²/m⁴ between the first two bins and the last two, 4.3 rad/m² in quadrature.

## What it does and does not show

- The growth is not confirmed by an independent sample at the level the rule asked for. Two independent additions now exist: DESI-only (net slope +0.5 ± 2.5, Entry 173, post hoc) and Quaia (+3.1 ± 1.7 here, or +3.7 ± 2.0 corrected for redshift error). Their weighted mean is +2.5 ± 1.6, against Entry 171's 4.63 ± 1.22.
- It is not contradicted either. All three samples agree within errors, and together they give t = +3.54 on 10,449 sight lines. The union's net slope, 3.41 ± 0.96, is smaller than Entry 171's, as expected when a first detection is followed up.
- The checks are uneven: the growth in the new sample sits mostly in the northern Galactic hemisphere (t = +3.09 on 1,113) and in the better redshifts (+2.88), and is weak in the south (+0.97) and for G < 20 (+0.87). None of these was a deciding check.
- The cause is still open (Entry 173).

## Files

- `quaia_growth.py`: the test (`synthetic`, `run`). Entry 171's functions are copied in so the folder stands alone.
- `synthetic.json`: tests on synthetic skies, made before the rule was frozen.
- `results.json`: the run of 21 November 2026.

Data, all git-ignored under `data/`:

- `data/quaia/quaia_G20.5.fits`: Quaia (Storey-Fisher et al. 2024, arXiv:2306.17749), Zenodo doi:10.5281/zenodo.10403370, https://zenodo.org/api/records/10403370/files/quaia_G20.5.fits/content, 171,020,160 bytes, md5 `98659ac4bd8a09da2c4ce653690d53df`, sha256 `9e357a760f2ff2db9417a6791fb71dc421ee30324301f22db898038c9e648584`, fetched 21 November 2026.
- `data/spice-racs-dr2/spice-racs.dr2.fits.gz`, `data/milliquas-v8/milliquas.zip`, `data/desi-dr1-qso/QSO_cat_iron_cumulative_v0.fits`: as in Entries 171 and 173.

Run from the repository root:

```
.venv/bin/python -I entries/174-quaia-replication/quaia_growth.py synthetic
.venv/bin/python -I entries/174-quaia-replication/quaia_growth.py run       # about 9 minutes
```
