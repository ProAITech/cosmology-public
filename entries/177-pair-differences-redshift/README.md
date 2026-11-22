# Entry 177: two components of one source, against redshift

Entries 171–174 found the scatter of rotation measure (RM) growing with quasar redshift at about 3.5σ, cause unknown. The difference between the RMs of two well-separated components of the same source cancels the Galactic foreground and anything else coherent across the source, and keeps twice the variance of any screen smaller than the separation. If small intervening screens make the growth, the variance of that difference must grow with redshift twice as fast as the single-sight-line variance does.

The rule was frozen in the Workbench (section 177, "the rule, frozen first") in commit `d154028` at 11:29 UTC on 21 November 2026 and published as Workbench version 191, before any difference was set against a redshift. The code and its synthetic test were committed (`3715167`) before the run.

## Result: unresolved by the frozen rule, narrowly

| | value |
|---|---|
| sources at \|b\| ≥ 20° with two components at least 30″ apart | 9,649 |
| of those with a redshift | **650** (Milliquas 430, DESI 74, Quaia 146) |
| slope of the difference's variance against log₁₀(1+z) | +17.2 |
| shuffled redshifts, 1,000 shuffles in 9 cells | +0.4 ± 17.5 |
| t | **+0.96** |
| slope unrelated small screens would give (2 β b_single) | 122.9 (from b_single = 63.8 ± 18.2) |
| net slope as a share of that | **0.14**; below 0.74 at 95% with the allowance for rotation at the source |
| u for "small screens carry all of the growth" | −1.55 (the rule needs −1.645 or less) |
| u for "at least half" | +0.04 |
| frozen decision | **unresolved**: no growth in the difference is established, and the whole of the single-sight-line growth is not excluded |

Variance of the difference in the six bins of redshift, 108–109 pairs each:

| median z | 0.13 | 0.35 | 0.55 | 0.80 | 1.12 | 1.82 |
|---|---|---|---|---|---|---|
| V (rad²/m⁴) | 19.5 | 11.8 | 9.7 | 27.1 | 19.5 | 21.8 |
| bootstrap error | 6.5 | 6.4 | 6.2 | 8.2 | 8.2 | 6.6 |

- The difference between two components never has a variance above 27 rad²/m⁴ in any bin, so screens smaller than the separation add at most about 14 rad²/m⁴ to a single sight line at any redshift reached here.
- The verdict falls short because of the allowance the rule makes for rotation at the source (45.5 in slope units, from the first bin's 19.5): without it u for "all" would be −2.7. The allowance is deliberately the largest possible.
- The pairs lie at lower redshift than the quasars of Entries 171–174 (median 0.66, quartiles 0.35–1.12), so they cover less of the range where the single-sight-line growth was seen.

## Reported only (item 7 of the rule)

| sample | pairs | slope | shuffles | t |
|---|---|---|---|---|
| Milliquas and DESI redshifts only | 504 | +2.7 | +1.6 ± 17.8 | +0.06 |
| separation below the median (39.5″) | 325 | +50.9 | +2.7 ± 34.0 | +1.42 |
| separation above the median | 325 | +34.5 | +1.2 ± 19.7 | +1.69 |
| smaller-error half | 325 | +49.3 | +5.3 ± 21.3 | +2.06 |

- With spectroscopic redshifts only, the share is 0.01 (below 0.64 at 95%) and u for "all" is −1.87; that variant was not the deciding one.
- The two sight lines are 128, 195 and 271 kpc apart (quartiles) at half the source's redshift; separations are 33″, 40″ and 54″.
- The whole sample has V = 19.3 rad²/m⁴ (raw scatter 5.0 rad/m², mean error² 5.8).

## What it does and does not show

- It does not show that small screens carry the growth: the difference does not grow significantly.
- It does not exclude it by the rule: "all of the growth" is disfavoured at u = −1.55, short of the 95% line.
- It agrees with Entry 173, which asked the same thing a different way. There, removing sight lines with detected Mg II absorbers lowered the single-sight-line slope by 18% of the union's net slope (below 69% at 95%). Here, screens smaller than a few hundred kpc account for 14% (below 74%). Two independent methods put the small-screen share near a sixth and neither excludes two thirds.
- If the single-sight-line growth is real and small screens carry a minority of it, the rest is coherent across a few hundred kpc: large-scale structure, or something common to both components. This test cannot tell those apart, nor from a single-sight-line growth that is a fluctuation.

## Files

- `pair_z.py`: `synthetic`, `extract`, `run`. Loaders and the FITS reader are copied from Entry 174 so the folder stands alone.
- `synthetic.json`, `results.json`.

Data, git-ignored under `data/`: as in Entries 171, 173 and 174 (SPICE-RACS DR2, Milliquas v8, the DESI DR1 quasar catalogue, Quaia). No new download.

Run from the repository root:

```
.venv/bin/python -I entries/177-pair-differences-redshift/pair_z.py synthetic
.venv/bin/python -I entries/177-pair-differences-redshift/pair_z.py extract    # about 1 minute
.venv/bin/python -I entries/177-pair-differences-redshift/pair_z.py run        # about 7 minutes
```
