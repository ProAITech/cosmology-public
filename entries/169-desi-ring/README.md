# Entry 169: the Big Ring on DESI DR1 Mg II absorbers

Rule frozen in the Workbench, section 169 (commit f3374e3, published as Workbench version 184, 20 November 2026, 04:04 UTC), before any count in the ring's area.

## Outcome: not tested

DESI's first release does not reach the Big Ring. The ring's annulus (140–260 Mpc around RA 216.0°, Dec +50.0°) holds one absorber in the slice and one in the outer windows. The median across the lattice is 22 outer absorbers per 45° sector, so each sector needs 11; the ring's eight sectors hold 0, 0, 0, 0, 0, 0, 0 and 1. Under item 3 of the rule, the outcome is *not tested*. Checks (a) to (d) use the same area and carry no information. Their numbers are in `results.json` only for completeness.

Coverage, from all sample absorbers at any redshift with RA between 205° and 227°:

| Dec | absorbers |
|---|---|
| 34–36° | 1,298 |
| 36–38° | 883 |
| 38–40° | 529 |
| 40–42° | 219 |
| 42–44° | 164 |
| 44–46° | 164 |
| 46–48° | 99 |
| 48–50° | 38 (all at RA < 208°) |
| 50–62° | 0 |

Within 6° of the ring's nominal centre there are 75 absorbers at any redshift. Within 6° of the Arc's centre (Entry 167) there are 2,149. DR1's northern edge in this RA range lies at about Dec +48° to +50°, through the ring's southern half.

## What the run still gives

- **Precision a later release would have.** 40 lattice nodes passed coverage, the rule's minimum being 30. At a node's own centre the annulus excess E scatters by sd₁ = 0.137 (±14%, as the rule expected from Entry 167's ±12% band). The 81-centre search lifts the null to a mean E_max of +0.233, with sd 0.132. Testing the ring this way needs an absorber catalogue that covers Dec +45° to +55° at RA 205°–227°, such as DESI's second release or a later one.
- **The Giant Ring, check (e), southern half only.** The frame was inferred as centred 240 Mpc south of the ring's nominal centre, at (216.0°, +45.16°). The elliptical shell's northern part is empty in DR1: sectors 1–6 of 16, the angles from east through north to north-west, hold no outer absorber. In the part DR1 covers, which includes the Arc:
  - 40 Mpc shell (as published): 26 absorbers against 29.4 expected, E = −0.12 ± 0.19.
  - 120 Mpc shell: 83 against 80.9, E = +0.03 ± 0.13.

  Only 6 and 7 lattice nodes pass the shell's coverage rule, fewer than the 10 the rule needs, so there is no rank. This half-shell overlaps the Arc's band, which Entry 167 already found without excess, so it adds little and decides nothing.

## Deviations and implementation choices

- No deviation from the frozen rule.
- Implementation choices the rule left open:
  - Absorbers are projected within 13° of each centre (24° for the shell), more than the 10.5° the search and surround need.
  - For the Giant Ring null, nodes within 25° of the inferred frame centre are dropped. The rule gave a 15° exclusion only for the ring.
  - The shell's sixteen sectors are equal angles in the ellipse's own normalised coordinates.
  - The coverage median is taken over all lattice nodes, including the dropped ones, and the ring.
- The code was tested before the real run on a synthetic uniform table: no excess (E_max at the ring +0.29 against a null of +0.25 ± 0.12). A ring of 120 extra absorbers was then planted at (216.5°, +50.8°) and recovered at (216.8°, +51.0°) with z = 21. The test script is not kept.

## Lesson

The rule had a *not tested* branch, and it was needed. Before freezing a rule on a survey's catalogue, check the survey's published footprint map for the target's position. That needs no value from the table. Here it would have shown that DR1 stops near Dec +48° to +50° at this right ascension.

## Files

- `desi_ring.py`: the whole analysis. Run it from the repository root with `.venv/bin/python -I entries/169-desi-ring/desi_ring.py`. It takes about 10 seconds, is deterministic and uses no random numbers.
- `results.json`: sample, R, the ring's counts and sectors, the null summary and the 40 nodes as (RA, Dec, E_max, E_nom, disc E), the checks, and the Giant Ring shell.

## Data

`data/desi-dr1-mgii/MgII-Absorbers-DR1.fits`: DESI DR1 Mg II absorber value-added catalogue v1.0, https://data.desi.lbl.gov/doc/releases/dr1/vac/mgii-absorber/, 65,750,400 bytes, sha256 `38c2e6281ce983a2b34b5810fdf0e3021cebaad3848cfc1cca974f185901bddf`. Table MGII_ABSORBERS, 270,529 rows; 191,538 in the sample after Entry 167's cuts and de-duplication.

Claims tested:
- Lopez, Clowes & Williger (2024), *A Big Ring on the Sky*, JCAP 07, 055, [arXiv:2402.07591](https://arxiv.org/abs/2402.07591). The centre was read from its Figure 6.
- Lopez & Clowes (2026), *A Giant Ring on the Sky*, [arXiv:2604.17534](https://arxiv.org/abs/2604.17534). Ellipse: a = 588, b = 442 Mpc, θ = 0.278 rad, centre (5.12, 52.2) Mpc.
