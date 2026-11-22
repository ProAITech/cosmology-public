# Entry 187: is the CMB colder around spirals beyond redshift 0.02?

**Bottom line: no. By the frozen rule, the claimed cooling is excluded.** On the Planck SMICA map, 3,583 late-type spirals (Sb–Sd) at z = 0.02–0.05 imply a cooling inside 0.2 Mpc of +0.0 ± 3.0 μK, with t = +0.01 against 1,000 random rotations of the sample. The 95% bound is 4.9 μK, or 3.3 μK for the published profile. That excludes the −15 μK found near nearby spirals at 5.0σ and the −27 μK of spirals in dense filaments at 9.0σ. It is also below the 4–10 μK that the fine-scale CMB power allows any galaxy (Entry 47). NILC, SEVEM and Commander agree to within 0.4 μK. The larger half and the denser half also show nothing (+0.6 ± 3.8 and +0.2 ± 4.6 μK). Whatever makes the nearby signal is not a property of late-type spirals, and stops before z = 0.02.

The rule is in Workbench section 187. It was frozen at 19:49 UTC on 21 November 2026 (commit `7fcf22d`, on GitHub before any galaxy was read) and published as Workbench version 199.

## Question

The claim (Hansen et al. 2025, A&A; arXiv:2411.15307) is that the CMB is colder inside 0.2 Mpc of large late-type spirals at z = 0.004–0.02: −27 μK for 1,250 in dense filaments, 3–4σ after allowing for the choices made. Luparello et al. (2023) give about 15 μK for nearby spirals in general. Entry 47 showed that a universal effect is capped at 4–10 μK by the fine-scale CMB power, and named the test: the same stack on spirals farther away.

## Data

- **Planck PR3 CMB maps** (git-ignored, cited in Entry 165): SMICA, NILC, SEVEM and Commander, `data/planck-pr3-<map>/COM_CMB_IQU-<map>_2048_R3.00_full.fits` (columns I_STOKES and I_STOKES_INP, K_CMB), with the common mask `data/planck-pr3-masks/COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits` (TMASK). `prep_maps.py` writes them to `work/` (git-ignored) as int16 files.
- **Nair & Abraham (2010)**, VizieR `J/ApJS/186/427/table2`: 14,034 SDSS galaxies at 0.01 < z < 0.1, g < 16, with T-types assigned by eye. Read in browser memory (`encode.js`).
- **2MRS** (Huchra et al. 2012), VizieR `J/ApJS/199/26/table3`: the claim's own selection, used as a positive control. Read in browser memory.
- **Transfer.** The galaxies the rule names were encoded in the browser (HEALPix index at Nside 4096, redshift and two flags; 6 bytes each, with a checksum on every line). They were read out of the page and decoded locally (`decode.py`), in the manner of `code/raw_chunks.py`. Chrome could not reach a server on this machine's localhost, so the stack could not run in the browser.

## Files

- `power.py`: noise of the disc stack on random positions, run before the rule.
- `predict.py`, `predict.json`: the share of a 0.2 Mpc top hat, and of the published profile, that Planck's beam keeps.
- `healpix.js`: HEALPix ang2pix (RING), checked against healpy on 242,001 points.
- `prep_maps.py`: the maps as int16 files.
- `stack.py`: the stack, the rotation null, the injection test (`synthetic`) and the decision.
- `encode.js`, `decode.py`: the browser read-out and its decoding.

## Result

The equal-weight stack uses 1,000 uniform random rotations of the whole sample as its null. C is the disc mean minus the mean at 2–3 Mpc, and D the disc mean, the claim's own statistic. δ is the cooling implied for a top hat inside 0.2 Mpc (C excess divided by f̄ = 0.672).

| sample | map | used | C (μK) | null | t_C | δ (μK) | 95% bound | D (μK) | null | t_D |
|---|---|---|---|---|---|---|---|---|---|---|
| **P, late spirals z 0.02–0.05** | **SMICA** | **3,583** | **+0.22** | **+0.20 ± 2.01** | **+0.01** | **+0.04 ± 2.99** | **4.9** | +2.29 | −0.01 ± 5.28 | +0.43 |
| P | NILC | 3,583 | +0.05 | +0.20 ± 2.01 | −0.08 | −0.23 ± 2.99 | 5.1 | +2.02 | −0.38 ± 5.27 | +0.46 |
| P | SEVEM | 3,583 | +0.44 | +0.20 ± 2.01 | +0.12 | +0.36 ± 2.99 | 4.6 | +0.14 | −2.16 ± 5.27 | +0.44 |
| P | Commander | 3,583 | +0.28 | +0.19 ± 2.01 | +0.04 | +0.13 ± 2.99 | 4.8 | +2.83 | +0.88 ± 5.25 | +0.37 |
| P | SMICA, ℓ ≤ 1 removed | 3,583 | +0.22 | +0.20 ± 2.01 | +0.01 | +0.04 ± 2.99 | 4.9 | +3.02 | +0.49 ± 7.77 | +0.33 |
| larger half (Rp90 > 8.55 kpc) | SMICA | 2,124 | +0.59 | +0.23 ± 2.45 | +0.15 | +0.56 ± 3.77 | 5.6 | +0.68 | +0.04 ± 5.61 | +0.11 |
| denser half (rho5 > −0.238) | SMICA | 1,763 | +0.25 | +0.15 ± 3.11 | +0.03 | +0.15 ± 4.63 | 7.5 | +3.26 | −0.02 ± 6.38 | +0.51 |
| control: 2MRS, z 0.004–0.02 | SMICA | 630 | −2.15 | +0.02 ± 4.03 | −0.54 | −2.6 ± 4.8 | 10.5 | −5.83 | −0.08 ± 4.71 | −1.22 |
| control | SMICA, ℓ ≤ 1 removed | 630 | −2.14 | +0.02 ± 4.03 | −0.54 | −2.6 ± 4.8 | 10.5 | −9.76 | +0.62 ± 5.30 | −1.96 |

- **Decision (P, SMICA).** t = +0.01, not ≤ −3. δ = +0.04 is at least −15 + 3 × 2.99 = −6.0, so *the claimed cooling is excluded*. For the published profile, which fades to zero at 1.5 Mpc, δ = +0.03 ± 2.04 μK (95% bound 3.3 μK).
- **The claim's own statistic, D.** It carries 5.3 μK of scatter from the large-scale CMB the discs share, and stands at t = +0.43. On random positions it reached only t = −2.2 with a −15 μK signal injected (SYNTH). This is why the rule used C.
- **Control.** The claim's own selection in 2MRS at z < 0.02 leans cold on D: t = −1.2 with ℓ ≤ 5 removed, −2.0 with ℓ ≤ 1 removed (one-sided p = 0.10 and 0.025). On C it stands at −0.5. It counts 694 galaxies, where the authors have 2,550. So the frozen isophotal-radius cut is stricter than their size measure, and the control is weaker than planned.

## Deviations

1. **2MRS control.** It has 694 galaxies rather than about 2,550, as above. This is reported, not decided.
2. **Output bug.** The first run of `stack.py run` stopped while writing its JSON (numpy float32). The computation is seeded and unchanged, and was rerun after the fix.
3. **Reading order.** The larger-half variant was read before the primary result, because a loop passed the wrong arguments in zsh. That affects only the order of reading.
4. **Transfer route.** The stack ran locally, not in the browser, because Chrome could not reach a server on this machine. The rule's Data paragraph already named the checksummed read-out that was used instead. `samples/readout.txt` (37 lines, 4,340 galaxies, sha256 290b7eb2…c42829) decodes to the counts the browser reported: P 3,646, larger half 2,161, denser half 1,809, 2MRS 694.

## Reproduce

```
.venv/bin/python -I entries/187-cmb-spirals/prep_maps.py        # maps to work/ (about 500 MB)
.venv/bin/python -I entries/187-cmb-spirals/decode.py           # samples/readout.txt -> sample files
.venv/bin/python -I entries/187-cmb-spirals/stack.py run P smica_l5
.venv/bin/python -I entries/187-cmb-spirals/summarize.py; .venv/bin/python -I entries/187-cmb-spirals/figure.py
```

More files: `samples/readout.txt` and `samples/check.py` (the transfer and its check); `result_<sample>_<map>.json` and `results.json`; `figure.py` and `figure1.svg`.
