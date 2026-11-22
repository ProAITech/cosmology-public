# Entry 178: has the cluster field changed since redshift 0.4, at equal thermal pressure?

**Result: unresolved by the frozen rule (Workbench section 178, commit `043d1e6`, 17:03 UTC, 21 November 2026).** The farther half of the eRASS1 clusters (median z 0.46) shows no detectable rotation-measure excess between 0.5 and 2 R500. At equal thermal pressure, its field is ρ = 0.18 times that of the nearer half (68% range 0.10–0.53, 95% range 0.10–1.10). The 95% range reaches the bottom of the grid, which makes the rule return *unresolved*. Equal fields (ρ = 1) are disfavoured by Δχ² = 3.3, about 1.8σ. The data lean toward a field that was weaker, relative to the gas pressure, 4 billion years ago, but do not establish it.

## What was done

- Sources: Entry 176's 92,169 SPICE-RACS DR2 sources.
- Clusters: Entry 176's 9,050 eRASS1 clusters, split at the median redshift, 0.2953.
- Assignment: each source belongs to the half of its nearest cluster, in units of R500.
- Excess: the noise-removed variance of ΔRM in each half and range, minus its mean over 500 rotations of the whole cluster list.
- Model: each cluster's field is tied to its thermal pressure, B²/8π = 1.92 n_e kT/β_pl. It uses Entry 176's gas model (n_e from MGAS500, β model with r_c = 0.2 R500), isothermal kT from M500 and R500, a reversal length of 10 kpc, and the front and member shares of `model.py`. The model is fitted to 0.5–1, 1–1.5 and 1.5–2 R500.
- Ratio: ρ = B_far/B_near at equal pressure = √(β_near/β_far), from the profile likelihood.

| half | clusters | median z | median M500 (10¹³ M☉) | kT (keV) | excess at 0.5–1 / 1–1.5 / 1.5–2 R500 (rad²/m⁴) | sources | 1/β_pl (68%) | B at 0.75 R500 (μG) |
|---|---|---|---|---|---|---|---|---|
| near | 4,524 | 0.173 | 15.2 | 2.6 | +85.8 ± 13.5, −1.8 ± 11.2, −7.9 ± 9.0 | 309, 495, 685 | 0.093 (0.063–0.133) | 1.58 (1.29–1.89) |
| far | 4,526 | 0.458 | 25.1 | 4.2 | +28.6 ± 28.8, +1.6 ± 21.4, +5.5 ± 19.2 | 59, 113, 163 | 0.003 (0–0.025) | 0.47 (0–1.36) |

- Inside 0.5 R500 (not fitted): near +252 ± 23 (t = +10.8, 184 sources); far +15 ± 53 (37 sources).
- The near half alone gives β_pl = 11 (7.5–16), the value Entry 176 found for all clusters.
- Reported only: with the median redshift of background sources at 0.6 or 1.0 instead of 0.8, ρ = 0.22 (95% up to 1.29) and 0.19 (up to 1.07).

## Synthetic test (before the run)

`synthetic.json`. A first version took its errors from the rotated clusters alone. On a sky with equal fields it returned ρ = 0.66 (95% range 0.44–0.85), a false "weaker". The cluster term brings its own sample variance, and that variance now enters the errors. With the correction:

| synthetic sky | ρ returned (95% range) |
|---|---|
| equal fields, four skies | 0.85, 1.08, 0.62, 1.29 (every range covers 1) |
| true ρ = 0.5 | 0.45 (0.20–0.97) |
| true ρ = 1.5 | 1.49 (1.01–2.19) |
| no field | unresolved |

The same omission was in Entry 176's field fit and has been corrected there (`entries/176-cluster-faraday-profile/model_errors.py`, commit `475b47c`). The field is unchanged at 2.02 μG, but its 68% range is now 1.70–2.35 μG.

## Deviations

- The two reported-only variants first ran without the changed redshift distribution, a bug that gave them identical numbers. They were rerun with it. The primary result is deterministic and did not change.

## Reading

- The far half has a fifth as many sources between 0.5 and 1 R500 as the near half (59 against 309). Its clusters are smaller on the sky, more of the sources behind them are in front of them, and their rotation is diluted by (1+z)⁻⁴.
- At the near half's pressure ratio (1/β_pl = 0.093), the model predicts an excess of 103 rad²/m⁴ at 0.5–1 R500 in the far half. Across the near half's 68% range the prediction is 81–127. The far half measures 29 ± 29, with an error from the rotations alone; with a field, the error is larger.
- Four effects could make the far half look weaker without any change in the field:
  - a reversal length that was shorter in the past (ρ measures B√Λ);
  - centring errors, which matter more for clusters 2′ across;
  - photometric-redshift scatter, which blurs the split;
  - masses or gas masses biased with redshift.
- None of these was tested here.
- A decisive test needs about three times more sources behind z ≈ 0.4–0.6 clusters, which POSSUM's full survey should provide.

## Files

- `field_z.py`: `synthetic` and `run`.
- `cluster_rm.py`: Entry 176's, unchanged.
- `synthetic.json`, `results.json`.
- Data: as Entry 176 (SPICE-RACS DR2 through Entry 171's cache, `data/erass1-clusters/erass1cl_primary_v3.2.fits`).
