# Entry 172: a current through nearby galaxies along their spin, in SPICE-RACS DR2

Entry 61 found a lean of about 2σ toward a current running through galaxies along their spin: +0.22 ± 0.13 × 10¹⁸ A from seven galaxies in POSSUM's first fields, and +0.20 ± 0.11 with Entry 60's northern stack. This entry repeats Entry 61's fit, galaxy by galaxy, on SPICE-RACS DR2, which covers the whole sky south of Dec +49°.
The rule was frozen in the Workbench (section 172, "the rule, frozen first") in commit `6875799` at 05:36 UTC on 21 November 2026 and published as Workbench version 186, before any rotation measure near a listed galaxy was fitted. The note's heading says 05:40 UTC; the commit came four minutes before that, and the run after both.

## Result: unresolved, with the stack closer to zero than Entry 61's

| | value |
|---|---|
| galaxies listed (LVHIS rotation fits, and PHANGS within 12 Mpc) | 57 |
| passing Entry 60's cuts (inclination ≥ 40°, \|b\| ≥ 20°, 1.5–80 Mpc) | 31 |
| also ≥ 30 sight lines at 30–200 kpc and every quadrant ≥ 10% | **20** |
| stack, weights 1/error² (10¹⁸ A, + along the spin) | **+0.086** |
| blank-sky null of the stack (200,000 draws) | mean +0.004, sd 0.136 |
| p_up (current along the spin needs ≤ 0.00135) | 0.25 |
| one-sided 95% upper bound U ("not confirmed" needs U < +0.20) | +0.30 |
| frozen decision | **unresolved** |

Stacks reported only:
- **Without Entry 61's galaxies:** 16 galaxies, +0.049 ± 0.152, p_up = 0.36, U = +0.29. Of Entry 61's seven, four pass every cut here; IC 1954, NGC 685 and NGC 1672 lie beyond 11 Mpc and have too few sight lines.
- **All well-covered galaxies, whatever the cuts:** 37 galaxies, −0.092 ± 0.105, U = +0.08. This includes low-latitude galaxies whose errors are understated (below).

Other settings, weighted as the primary (reported only; Entry 61's values for its seven galaxies in brackets):

| setting | stack | positive |
|---|---|---|
| foreground of order 1 | −0.008 [+0.155] | 10 of 20 |
| order 3 | +0.006 [+0.127] | 9 of 20 |
| from 50 kpc | +0.085 [+0.160] | 12 of 20 |
| to 300 kpc | −0.042 [+0.131] | 9 of 20 |
| 100–200 kpc only | **−0.270** [+0.654, all seven positive] | 8 of 20 |
| template turned by 90° | −0.047 | 10 of 20 |
| LVHIS second fits | +0.055 | 11 of 20 |
| step template, uniform gas (units of Entry 61's step) | +0.102 [+1.28] | 10 of 20 |

Entry 61's galaxies measured again:

| galaxy | here | Entry 61 |
|---|---|---|
| NGC 300 | +0.14 ± 0.43 | +0.23 ± 0.28 |
| ESO 154-G023 | +0.00 ± 0.32 | +0.13 ± 0.16 |
| NGC 1311 | +0.21 ± 0.34 | +0.16 ± 0.19 |
| IC 1959 | +0.72 ± 0.36 | +0.11 ± 0.18 |

The errors are the robust blank-sky spreads.

## What it means

The frozen outcome is *unresolved*: the stack is neither large enough for a current nor bounded below Entry 61's +0.20. Read beyond the rule, every sign points away from a current tied to spin:
- the stack fell from +0.22 to +0.09;
- the 16 new galaxies give +0.05;
- the positive counts sit at half;
- the 100–200 kpc pattern that was all positive in Entry 61 is negative here;
- the turned template gives as much as the real one.

The test did not gain precision. At 3.7 sources per square degree each galaxy's blank-sky error is 0.25–0.5, so twenty galaxies reach ±0.14, no better than Entry 61's seven at POSSUM's density. Entry 61's forecast of ±0.04 for 53 galaxies needs POSSUM's full survey, not SPICE-RACS.

Combining with Entry 61 is not clean: four galaxies are shared and the sky is the same. With only the 16 new galaxies added to Entry 61's north–south value, a weighted mean gives +0.15 ± 0.09.

## Notes

- Low-latitude galaxies give wild amplitudes, for example NGC 4945 (b = +13°) +4.0, ESO 215-G009 (+10°) −16.9 and Circinus (−4°) −3.9. Their blank-sky errors come from |b| ≥ 20° places and are far too small, as the rule's "where it can fail" warned. Entry 60's latitude cut keeps them out of the decision; they are why the "all well-covered" stack is not usable.
- Nine PHANGS galaxies at 7.9 Mpc or beyond, and AM 0605-341 and ESO 461-G036 have fewer than 30 sight lines. NGC 4826 enters from PHANGS (+84° latitude); Entry 60 noted its counter-rotating gas, and its value is −0.54 ± 0.34.
- Blank places: 400 per galaxy at 6 angles each, 2,400 amplitudes; every galaxy reached the cap.

## Deviations from the frozen text

None. Implementation choices:
1. Offsets are gnomonic (tangent-plane) coordinates times the distance. Entry 61's code took them from its own neighbour search, whose projection is not recorded.
2. A merged source takes its first row's position, as in Entries 166 and 171.
3. Blank-place candidates are a grid from Dec −90° to +49°.
4. The places are taken in a random order with seed 1234, and the stack draws use seed 1357.

## Files and how to run

- `galaxy_currents.py`: `synthetic` (injection test on a synthetic sky, reads no data; +0.5 injected returns +0.500 and +0.499), `extract` (30 s; `work/` is git-ignored), `run` (45 s), from the repository root with `.venv/bin/python -I`.
- `galaxies.txt`: the 57-galaxy table (LVHIS tables 9 and 2, Koribalski et al. 2018, VizieR J/MNRAS/478/1611; UNGC, Karachentsev et al. 2013, VizieR J/AJ/145/101; PHANGS, Lang et al. 2020, VizieR J/ApJ/897/122 table 1), read in browser memory on 21 November 2026 and typed in.
- `results.json` (31 kB): every galaxy's fit, coverage, blank-sky spread and percentile, the three stacks with their nulls, and the other settings.

## Data

SPICE-RACS DR2, `data/spice-racs-dr2/spice-racs.dr2.fits.gz`, CSIRO DAP csiro:64891 (https://data.csiro.au/collection/csiro:64891), 4,968,301,262 bytes, sha256 `a1a2f3c1d12a8a8b3aa8526d74b8e0115054f610e2e4186108f0f068005d0053`. Entry 166's quality cuts: 134,400 sources.

## What the cloud session should write up

Entry 172 needs these pieces:
- **The result.** Unresolved by the frozen rule: +0.09 ± 0.14 from 20 galaxies, U = +0.30. The 16 new galaxies give +0.05 ± 0.15, and every secondary sign points away from a current along the spin, including Entry 61's 100–200 kpc pattern, now −0.27.
- **A figure.** Galaxy amplitudes with blank-sky errors, Entry 61's seven beside them, and the stacks.
- **The Gap Map node.** "A current through other galaxies" changes from "a lean of about 2σ" to "the lean is not repeated in SPICE-RACS (+0.05 ± 0.15 from 16 new galaxies); settling it at ±0.04 needs POSSUM's full survey".
- **A Ledger sentence,** and section 172 marked as run.
