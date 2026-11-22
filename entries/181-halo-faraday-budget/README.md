# Entry 181: how much of the Faraday growth can intervening haloes make?

**Bottom line: almost none.** Entries 171–174 found the rotation-measure scatter of quasar sight lines growing with redshift: 63.8 ± 18.2 rad²/m⁴ per unit log₁₀(1+z), about 19 rad²/m⁴ between z = 0.52 and 2.65. Entries 176–180 measured what one cluster or group does to a sight line that crosses it. Turned into a budget, those measurements show that haloes crossed on the way produce 0.3–3% of the growth for the evolution Entry 178 supports. Even with no dilution, and with the measured mass scaling extrapolated down to galaxy-mass haloes, they produce at most about a quarter. If the growth is real, it does not come from virialised haloes. It must come from structure outside them (filaments, the intergalactic medium), from the sources, or be a fluctuation.

This is a calculation from measured quantities. No new data were read and no rule was frozen: nothing is decided by it.

## Inputs

- **Response of the statistic.** ψ = dV/df is the rise of (1.4826 MAD)² − ⟨e²⟩ in a large sample when a small share f of its sight lines cross a halo inside R500. It is measured by mixing field sources with real sources inside R500 (Entry 180's samples):
  - all 669 sources around 10,151 clusters: ψ = 64.9 rad²/m⁴;
  - the 74 background quasars: ψ = 88.2 rad²/m⁴.
- **Weight of a crossing.** w = (M/M_eff)^0.86 × [(1 + z_ref)/(1 + z)]⁴, with M_eff = 2.3 × 10¹⁴ M☉ and z_ref = 0.12 for the clusters behind which ψ was measured.
  - The exponent 0.86 is Entry 176's mass relation.
  - The (1+z)⁻⁴ is the dilution of rotation made at z. Entry 178's far half matches it: 29 ± 29 rad²/m⁴ measured, about 35 expected from dilution alone.
  - A variant with no dilution serves as an upper bound.
- **Crossings.**
  - (a) Catalogued clusters (eRASS1, MCXC, PSZ2), counted on the positions of the 11,554 real quasar sight lines.
  - (b) All haloes above M_min, from the Tinker et al. (2008) mass function for M500c, using an Eisenstein & Hu (1998) no-wiggle spectrum with Planck 2015 parameters (σ8 = 0.816 reproduced).
  - Only haloes at z < z_quasar − 0.05 count.

## Results

**Catalogued clusters on the real sight lines.** 0.3–1.0% of quasar sight lines cross one inside R500, in all eight redshift bins. The added variance is 0.2–0.8 rad²/m⁴, with no trend: slope 0.18 rad²/m⁴ per unit log₁₀(1+z), 0.3% of the growth.

*Check of the mass function.* For M500 > 2 × 10¹⁴ M☉ at z < 0.3, the prediction is that 0.64% of sight lines cross a halo over the whole sky. The catalogues, covering about half the sky with incomplete MCXC coverage, give 0.25%.

**All haloes (mass function).** Expected crossings and added variance for sight lines to z = 0.52 and z = 2.65, with ψ = 64.9:

| M_min (M☉) | crossings to 0.52 / 2.65 | with dilution: V_add, slope, share | no dilution: V_add, slope, share |
|---|---|---|---|
| 10¹⁴ | 0.025 / 0.038 | 1.6 / 1.7; 0.4; 1% | 1.9 / 2.6; 2.1; 3% |
| 3 × 10¹³ | 0.083 / 0.152 | 2.6 / 2.9; 1.0; 2% | 3.1 / 5.0; 5.1; 8% |
| 10¹³ | 0.18 / 0.39 | 3.2 / 3.7; 1.4; 2% | 3.9 / 6.9; 8.0; 13% |
| 3 × 10¹² | 0.34 / 0.87 | 3.6 / 4.2; 1.7; 3% | 4.4 / 8.4; 10.6; 17% |
| 10¹² | 0.55 / 1.61 | 3.7 / 4.5; 1.9; 3% | 4.6 / 9.2; 12.1; 19% |

- With the background quasars' ψ (88.2), every share is 1.36 times larger. The extreme case (galaxy-mass haloes, no dilution) then reaches 26%.
- Haloes add 3–5 rad²/m⁴ to every extragalactic sight line, almost independent of redshift. That is a floor to subtract before reading the remaining scatter as intergalactic.

## Reading

- The cause of the growth was open after Entries 173 (Mg II absorbers, about 13–18%, below about half) and 177 (screens smaller than a few hundred kpc, about 14%, below 74%).
- This budget closes the halo route from the other side. The rotation clusters and groups are measured to produce is far too small, and too flat in redshift, to make the growth.
- Galaxy-mass haloes are an extrapolation of the cluster mass relation by two orders of magnitude, and they still give at most a quarter.
- What remains possible:
  - rotation in structure coherent over more than a few hundred kpc outside haloes (the cosmic web);
  - rotation at the sources that grows with redshift;
  - a statistical fluctuation (the growth is 3.5σ and not independently confirmed).

## Limits

- ψ is linear in f only for small shares. At the 1.6 crossings per sight line of the most extreme case the approximation is generous: it overstates the halo term.
- The mass function is good to about 30% (n(>10¹⁴ M☉) = 6 × 10⁻⁶ Mpc⁻³ at z = 0).
- The 0.86 mass exponent was measured over 8 × 10¹³ to 4 × 10¹⁴ M☉ only.

## Files

- `budget.py`
- `results.json`
- `cluster_rm.py`, `other_half.py` and `background.py`: Entries 176, 179 and 180, copied unchanged.
