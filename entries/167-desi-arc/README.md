# Entry 167: the Giant Arc on DESI DR1 Mg II absorbers

Does the band through the Giant Arc's place hold more Mg II absorbers than the sightlines predict, on DESI's much larger sample? Entry 66 found +34% (36 against 26.8, chance 5%) in a second SDSS-based catalogue, with ±19% counting noise.
Rule frozen in the Workbench (section 167, "the rule, frozen first") in commit `96074ac` on 14 November 2026, before any value in the table was read; the first read of values was the global-only summary the rule specifies, and it needed no amendment.

## Result: no excess at the Arc's place

By the frozen rule the outcome is **no excess**.

| | |
|---|---|
| absorbers in the band (1,000 × 200 Mpc, z = 0.742–0.862) | 78 |
| expected from absorbers at other redshifts on the same sightlines (R = 0.2696 × 307) | 82.8 |
| excess E | **−0.058 ± 0.119** (Poisson) |
| lattice null, 52 bands elsewhere in DESI | mean +0.013, sd 0.145; 97.5th percentile +0.25, maximum +0.28 |
| (E − mean)/sd | −0.49; 34 of the 52 null bands are higher |
| one-sided 95% limit on the excess | +0.18 |
| against Entry 66's +0.34 ± 0.19 | z = −1.66 |

The DESI band holds no more than expected, and an excess above +18% is excluded at 95%. The earlier +34% is 1.7 errors above this value, so the two are not strongly in conflict; the second catalogue's excess is not confirmed here.

## Checks (reported, not decisive)

- **(a) 300 Mpc wide.** 118 against 126.4 expected, E = −0.067 (Entry 66: +30%); rank p = 0.71.
- **(b) The neighbouring slices.** The slice just *below* (z = 0.630–0.742) holds **80 against 59 expected, E = +0.355**, null mean −0.018 ± 0.157, rank p = 0.037. The slice above (0.862–0.991): E = −0.143, rank p = 0.92. Entry 66 had seen the same pattern in the SDSS-based catalogue, 33 against 22.8 in the slice just below, and noted it might mean the feature lies deeper than one slice. Here it is one of two neighbouring slices looked at, p = 0.037 each; allowing for two, about 0.07. It is not part of the frozen decision and it is not a detection. The DESI and SDSS samples share an unknown fraction of sightlines, so the two slice-below excesses are not independent.
- **(c) Linking.** The test is uninformative here. With 610 sample absorbers in the field the link is 89 Mpc, 0.8 of the mean spacing, close to the percolation threshold (about 0.87), and the largest group holds 477 of the 610 (78%), spanning 1,537 Mpc, the width of the field, with its axis 0.3° from east–west and its centre 178 Mpc from the cluster. A group that fills the box is no match for the Arc by position. Its member count ranks at p = 0.28 among the 52 lattice fields (null median 214, 84th percentile 580, maximum 685). The DESI footprint near the field is not complete, which concentrates the absorbers and pushes the linking toward percolation. The paper's hull significance was not computed, as the rule said.

## What this test could and could not do

- The counting noise is ±12% per band, not the ±7% projected in the rule's introduction. The projection took the whole file (270,529 absorbers against the second catalogue's 35,629, 7.6 times) to apply to the Arc's field. In the field the sample holds 610 absorbers in the slice against the second catalogue's 336, 1.8 times as many, and the band's expected count is 83 against 26.8, 3.1 times. The null's scatter (0.145) is above the Poisson noise (about 0.12); the extra, about 0.08 in quadrature, is DESI's tile pattern and real clustering together, and I did not separate them.
- The test could have confirmed an excess of +0.34 only at about 2σ. It excludes much larger ones and puts the 95% limit at +0.18. It cannot show that no modest structure exists. In the standard model a band of this size differs from the mean by about 3%.
- Shared sightlines with the discovery were not removed; as the outcome is "no excess" that does not affect the result.

## Deviations and implementation choices

1. **Median over cells.** The rule says each of the band's five cells must hold at least half "the median number of outer absorbers per cell over the lattice". Taken over all lattice cells the median would be zero, as most of the lattice lies outside DESI. I took the median over cells with at least one outer absorber (43 per cell, threshold 21.5).
2. **Duplicates.** 30 of 191,568 rows were removed by the rule's duplicate step (65 rows shared a TARGETID and redshift before any cut; the others are the 1″ rule).
3. **Overlap rule.** A lattice node is dropped when its centre lies within |x| < 1,000 Mpc and |y| < 200 Mpc of the cluster in the Arc's tangent plane.
4. **Positions in the field.** Offsets use 2,837 Mpc for every absorber, as the rule says, and the third coordinate for linking is the comoving distance minus 2,837 Mpc.
5. **Lattice.** 218 nodes; 52 remain after the coverage and overlap rules. The nodes' fields (check c) overlap heavily, so those ranks are not independent.
6. The rule's introduction said the 7.6-fold increase "brings the band's noise to about ±7%". It did not (above).

## Files and how to run

- `desi_arc.py`: the whole analysis (4 s): `.venv/bin/python -I entries/167-desi-arc/desi_arc.py`. Deterministic.
- `results.json` (10 kB): sample counts, R, the Arc's counts and E, every null node (RA, Dec, N_band, N_oth, E), the checks and the linking result.
- Checked on synthetic catalogues before the real run: distances 2,665.05 / 2,837.05 / 3,003.27 Mpc against the rule's 2,665 / 2,837 / 3,003; null mean ≈ 0 with Poisson-sized scatter; an injected +35% in the band returned +31% ± 12%.

## Data

DESI DR1 Mg II absorber value-added catalogue v1.0 (https://data.desi.lbl.gov/doc/releases/dr1/vac/mgii-absorber/), `data/desi-dr1-mgii/MgII-Absorbers-DR1.fits`, 65,750,400 bytes, sha256 `38c2e6281ce983a2b34b5810fdf0e3021cebaad3848cfc1cca974f185901bddf` (verified against the published checksum file). Table MGII_ABSORBERS, 270,529 rows, 33 columns. Used: TARGETID, RA, DEC, ZWARN, Z_QSO, Z_MGII, LINE_SNR_MIN, TSNR2_QSO. Global summary read first: all rows SURVEY = main; ZWARN 0 for 267,484; LINE_SNR_MIN 25th percentile 3.15; 204,019 distinct TARGETIDs.

## What the cloud session should write up

Entry 167: **no excess** at the Arc's place on DESI (E = −0.06 ± 0.12; limit +0.18; 1.7σ below Entry 66's +0.34), the Arc's claim not confirmed at the sensitivity DESI gives in this field; the correction that the field's DESI sample is 1.8 times the second catalogue and the noise ±12%, not eight times and ±7%; the slice-below excess (+0.36, p = 0.037, about 0.07 for two slices) described as a pre-registered check and not a detection, with the note that Entry 66 saw the same in the SDSS-based catalogue; the linking check as uninformative; a figure of E for the Arc's band beside the 52 null bands (`primary.null.nodes`), with the slice-below and slice-above values marked; the Gap Map update (the Arc's band excess is not seen on DESI, with a 95% limit of +18%; Entry 66's linking result is not tested by DESI here, because the check was uninformative; independence from the discovery's sightlines was not needed for a null result); a Ledger sentence; and section 167 marked as run. What would settle the slice-below question: remove sightlines shared with SDSS (cross-match with the SDSS quasar catalogues in browser memory) and fix a rule for a redshift-resolved profile before looking.
