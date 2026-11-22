# Entry 183: can the Faraday growth be separated from radio luminosity? A power study, not run on the data

**Bottom line: not with this sample.** The growth of RM scatter with quasar redshift (Entries 171–174, t = +3.54 on 10,449 sight lines) might come from the sources rather than from space between: in a flux-limited sample the distant quasars are more luminous. The natural test shuffles redshifts at fixed radio luminosity and asks how much growth is left. On simulated skies it cannot tell the two causes apart. Luminosity and redshift correlate at 0.85, so a growth made by redshift and one made by luminosity leave t_L values that differ by about one standard deviation. A test that would almost always return *unresolved* was not frozen or run. No flux density was read.

## Design tested

`lum_growth.py`:
- **Luminosity.** L = 4π D_L² S (1+z)^(−1−α) at 888 MHz, α = −0.7, Planck 2018 distances.
- **Statistic.** Entry 174's union statistic with redshifts shuffled within 60 cells (latitude band × RM-error tercile × luminosity quintile) instead of 12.
- **Decision.** Simulated on the sample's own z and L: the observed t_L is compared with its distribution when the growth comes from redshift only (R_z) and from luminosity only (R_L), each matched to the observed rise of 19 rad² m⁻⁴.

## Synthetic skies (`synthetic.json`)

The skies have 10,000 sources with flux densities from N(>S) ∝ S^−1.5 above 8 mJy and log-normal redshifts (median 1.3), giving corr(log L, log(1+z)) = 0.85. There are three truths: growth 19 from luminosity, growth 19 from redshift, and none. Each was judged against 30 calibration skies per model, with 80 shuffles each.

| truth | t₁₂ | t_L | calibration, redshift growth: t_L, share kept | calibration, luminosity growth: t_L, share kept | decision |
|---|---|---|---|---|---|
| luminosity | 4.00 | 0.59 | 1.31 ± 1.03, 0.25 | 0.61 ± 1.02, 0.11 | unresolved |
| redshift | 2.94 | 2.46 | 1.61 ± 0.90, 0.30 | 0.38 ± 0.78, 0.08 | not a luminosity effect |
| none | 0.90 | 0.37 | 1.29 ± 0.81, 0.25 | 0.18 ± 1.30, 0.03 | (t₁₂ < 3) |

- **Absorption.** Shuffling at fixed luminosity takes up 70–75% of a genuine redshift growth, not just the luminosity part.
- **Overlap.** The two calibration distributions of t_L sit 0.7–1.2 apart, with widths of 0.8–1.3. The correct call for the redshift truth came from a high draw (2.46, about the 85th percentile of its own model).
- **Why.** This is not a defect of the method. Two predictors correlated at 0.85 inflate the error on each by 1/√(1 − 0.85²) ≈ 1.9. A 3.5σ growth then leaves at most about 1.8σ for redshift at fixed luminosity.

## What would separate them

- **Sources of very different luminosity at the same redshift.** A survey several times deeper at the same frequency (POSSUM's full survey, or MeerKAT) would widen the luminosity range at each redshift and lower the correlation.
- **Sources whose rotation cannot come from their surroundings:** pulsars and fast radio bursts at known distance. These test the intervening medium without any source term.
- **Rest-frame tests.** Rotation at the source is diluted by (1+z)⁻². The growth adds 4.3 rad m⁻² in quadrature at z = 2.65, which a source term could supply only with about 57 rad m⁻² in the rest frame. If the most luminous radio quasars carried that much, the few at z < 0.5 would show it directly, diluted by a factor of 2 at most.

## Files

- `lum_growth.py`: the design; the `run` stage was never executed.
- `quaia_growth.py`: Entry 174, copied unchanged.
- `synthetic.json`.
