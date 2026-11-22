# Entry 165: the frozen deeper-slab void test, on Planck PR3 SMICA

The test, its catalogue, statistic, errors and thresholds were frozen in **Entry 148** (11 November 2026), before any CMB map was opened. This folder runs it.

**Choices the protocol left open are fixed below, 13 November 2026 21:42 UTC, before any void was evaluated on the map.** Before this file was written the map was opened only for its global properties: file format, 77.94% of the sky unmasked by the common mask, rms 107.4 μK over the unmasked pixels. Anything that changes after the stack is computed will be listed in the Workbench entry as a deviation.

## Data

| Item | Path | Source | Size (bytes) | sha256 |
|---|---|---|---|---|
| Planck PR3 SMICA, table COMP-MAP, columns I_STOKES and I_STOKES_INP | `data/planck-pr3-smica/COM_CMB_IQU-smica_2048_R3.00_full.fits` | https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/component-maps/cmb/COM_CMB_IQU-smica_2048_R3.00_full.fits | 2,013,312,960 | `60952c645eb33d151905ddf5837477e15ca02a3261feca4ceca3c3feece4f9ac` |
| Planck PR3 common intensity mask, column TMASK (1 = use) | `data/planck-pr3-masks/COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits` | https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/masks/COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits | 201,335,040 | `23f49a3479073408b83a7f1bb2cca490d9be7f934b2e0e891524b6a1d06a5fc4` |

The archive publishes no checksums; the byte counts equal the server's `content-length`. Both files are Nside 2048, NESTED, Galactic. `data/` is not in git: re-fetch from the URLs above.

Catalogue: `code/voids_slab1.txt` (36 voids) and `code/voids_slab2.txt` (62 voids, of which 43 fresh and 19 overlapping), as frozen in Entry 148. Spectrum for simulated skies: `code/cmbreion_tau054.json`, the notebook's TT spectrum scaled to a first peak of 5,733 μK², exactly as in `code/prereg.py`.

## Choices fixed here

1. **Resolution.** All decisions use Nside 256 (NESTED). A low-resolution pixel is used only if all 64 of its Nside-2048 sub-pixels are unmasked; its value is the mean of those 64 (`healpy.ud_grade(..., pess=True)`). The simulations run at this resolution, so the decision rules are applied to it. The same statistic on the native Nside-2048 map, using every unmasked pixel, is reported as a cross-check on the real map only.
2. **Map.** I_STOKES, not the inpainted column; K_CMB × 10⁶ gives μK. Masked pixels are dropped. Monopole and dipole are fitted by least squares on the used pixels (`healpy.fit_dipole`) and subtracted, identically for the real map and for every simulation.
3. **Regions.** Pixel centres at angular distance d ≤ θ/2 from the void centre form the inner disc, and θ < d ≤ 1.5θ the annulus, with θ the catalogue radius. Means are over used pixels, equal weight. A void enters a stack only if both regions hold at least one used pixel. Voids dropped are counted and listed.
4. **Samples**, as in Entry 148's table. "Largest quarter" means radius at or above the 75th percentile of its own sample, as in `code/prereg.py`: nearest slab 36 (largest quarter 9), fresh 43 (largest quarter 12), overlapping 19, all 62.
5. **Simulated skies.** 1,000 Gaussian skies at Nside 256, ℓ ≤ 700, from the notebook's TT spectrum (the basis of the forecast), with no beam or pixel window, the same mask, the same monopole and dipole removal, at the same void positions. Seed 20261113. A second 1,000, from SMICA's own spectrum (the inpainted map, degraded to Nside 256), is reported beside them. The notebook-spectrum set is the primary error.
6. **Rotations.** 1,000 random shifts Δl, uniform in [0°, 360°), applied to the whole catalogue at once (both slabs together), latitudes unchanged. Seeded.
7. **Decision.** The thresholds of Entry 148 apply to the raw stack in μK: the check is +10 μK; the decision is +20 μK or more (real effect), below +15 μK (chance, or the nearest structures only), and in between unresolved. The σ equivalents are reported but do not move the thresholds. The stack with the null mean subtracted is reported too.
8. **Not done here.** The repeat in NILC, SEVEM and Commander (protocol line "repeated in…") needs maps not yet downloaded, and is left until they are.

## The repeat in three more maps

**Fixed 21 November 2026, 09:55 UTC, before any of the three maps was opened** (only their sizes on the server were read). Entry 148's protocol says the stack is "repeated in NILC, SEVEM and Commander". `other_maps.py` does it.

1. **Maps.** Planck PR3 full-mission CMB maps from NILC, SEVEM and Commander, column `I_STOKES`, from the same IRSA directory as SMICA: `COM_CMB_IQU-nilc_2048_R3.00_full.fits` (1,610,660,160 bytes), `COM_CMB_IQU-sevem_2048_R3.00_full.fits` (2,013,278,400) and `COM_CMB_IQU-commander_2048_R3.00_full.fits` (1,610,660,160), under `data/planck-pr3-nilc/`, `-sevem/` and `-commander/`. Checksums are added below once the files are on disk.
2. **Everything else as for SMICA,** by calling the same functions: the common intensity mask, the Nside-256 map from fully unmasked pixels, monopole and dipole fitted and removed on the used pixels, the same regions, samples and void list.
3. **Errors.** The simulated-sky nulls of the SMICA run are reused, since they depend only on the mask and the void positions. The rotation null is recomputed on each map with the same 1,000 shifts.
4. **Decision.** Entry 148's thresholds, applied to each map separately: check, +10 μK for the largest quarter of the nearest slab; decision on the 43 fresh voids, +20 μK or more (real effect), below +15 μK (chance, or the nearest structures only), between (unresolved). The test is complete if all four maps fall on the same side; if they disagree, the entry reports the disagreement as a foreground-cleaning dependence and the outcome is unresolved.
5. Nothing here changes once a map is stacked.

### Result of the repeat (21 November 2026)

The section above was committed at 09:11 UTC (`4696044`), before the maps were opened; its heading gives a later time. Files on disk: NILC sha256 `8f9a8533bc40b46b51b6437b8624dc523559c14bff38aee5df85500dbd4f9436`, SEVEM `89cf9495c9ab289f0035d9dad0b187c93318f12d6201d014c4fbe27830515f87`, Commander `fe0332e61087ab445656ef272b31578ca4b1edcc7c3784a7662577684e2685a4`; sizes as listed in item 1.

| map | check: largest quarter of the nearest slab (μK) | decision: 43 fresh voids (μK) | all 62 deeper voids (μK) | +35 μK excluded at | step 2 |
|---|---|---|---|---|---|
| SMICA (13 November) | +12.05 | +8.30 ± 5.13 | +4.57 | 5.2σ | below +15 μK: chance, or the nearest structures only |
| NILC | +12.08 | +8.36 ± 5.13 | +4.65 | 5.2σ | below +15 μK: chance, or the nearest structures only |
| SEVEM | +11.78 | +8.26 ± 5.13 | +4.59 | 5.2σ | below +15 μK: chance, or the nearest structures only |
| Commander | +11.83 | +8.29 ± 5.13 | +4.53 | 5.2σ | below +15 μK: chance, or the nearest structures only |

- **The test is complete, and the four maps agree.** The fresh voids stack to +8.3 μK in every map, with the same error, below the +15 μK line. The nearest-slab check stays above +10 μK in every map.
- Void by void on the fresh sample, the three new maps differ from SMICA by 0.7 to 0.9 μK rms and by less than 0.07 μK in the mean. Foreground cleaning does not move the result.
- 42 of the 43 fresh voids are used in every map (one has its inner disc inside the mask, as with SMICA). Each map keeps 76.6% of the sky at Nside 256 and has an rms of 100 μK after the monopole and dipole are removed.
- Rotation nulls recomputed on each map give ±4.84 to ±4.85 μK, against ±5.13 from simulated skies.

Results: `results_other_maps.json`.

## Files

- `smica_voids.py`: the analysis. Stages `prep`, `null`, `real`, or `all`. Run with `.venv/bin/python entries/165-smica-voids/smica_voids.py all` from the repository root. Python 3.12 with numpy 2.5.4, scipy 1.18.1, healpy 1.20.1, astropy 8.0.2.
- `results.json`: per-void values, stacks, null widths and decisions (SMICA).
- `other_maps.py`, `results_other_maps.json`: the repeat in NILC, SEVEM and Commander.
- `work/`: intermediate arrays (git-ignored).
