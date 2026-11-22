# Entry 168: stretched waves with the collaboration's own selection

Entry 79 asked whether gravitational waves are stretched with distance, as light is: detector-frame masses read heavier by (1 + z)^α, α = 1 for expansion, α = 0 if the waves arrive unstretched (as any plasma or tired-light redshift would leave them). It found α = 0.9 from the open catalogue, but the answer moved with a detection threshold the open data fix only roughly, and the 104 binaries above 40 M☉ alone disfavoured expansion by 5.2 in log-likelihood. This entry keeps Entry 79's likelihood and replaces its selection with the LVK injection set: which of 1.7 million simulated signals, put into the real O1–O4a data, the searches found.
Rule frozen in the Workbench (section 168, "the rule, frozen first") in commit `9747176` on 15 November 2026, before any injection value was read; the global read it specifies needed no amendment.

## Result: stretch confirmed; the heavy-binary residual is gone

| | all 155 events (O1–O4a) | 65 heavy (M_det > 40 M☉) | 127 light (< 60 M☉) |
|---|---|---|---|
| best α (68%) | 1.4 (1.29–1.82) | 0.6 (0.23–1.29) | 1.0 (0.67–1.21) |
| unstretched (α = 0) worse by Δln L | 8.15 | 1.15 | 4.08 |
| expansion (α = 1) worse by Δln L | 1.58 | **0.08** (Entry 79: 5.2) | 0.00 |

- **Decision on the stretch: stretch confirmed.** T = ln L(α = 1) − ln L(α = 0) = **6.57**. 50 mocks with unstretched waves give −1.52 ± 1.70, maximum 1.67: the real value is above every one, 4.8 of their standard deviations above their mean. 50 mocks with stretched waves give 1.16 ± 3.46, minimum −2.70: the real value lies inside them (+1.6 sd).
- **Decision on the heavy binaries: the residual is gone.** h = ln L(best) − ln L(α = 1) for the heavy events is **0.08**. The stretched mocks give a median of 2.69 (84th percentile 5.18, maximum 9.66); 49 of 50 are at least as large, so p_h = 0.98.
- **Convergence: passed.** From the two halves of the found injections, T = 6.55 and 6.58, h = 0.088 and 0.070 (differences 0.03 and 0.02; the rule allows 1.0).

## Checks (reported, not decisive)

- **(a) Entry 79's selection on the same 155 events:** T = 3.6, best α 1.4 (0.85–1.52); heavy h = 2.83, best α −0.5 at the scan's edge. So of the heavy binaries' 5.2 in Entry 79, leaving out the 104 O4b events (39 of them heavy) accounts for about half (5.2 → 2.8), and the collaboration's selection removes the rest (2.8 → 0.08). With the new selection the overall preference for a stretch also strengthens (T 3.6 → 6.6).
- **(b) Smoothing of 0.16 instead of 0.08:** T = 7.65, h = 0.00; best α 1.7 (1.35–1.85) for all, 1.0 (0.27–1.40) for the heavy.
- **(c) By run group:** O1+O2 (10 events) α = 0.2, uninformative (68% −0.5 to 2.5); O3 (59) α = 1.0 (−0.07 to 1.63); O4a (86) α = 2.2 (1.54 to 2.5).

## What to watch

The best α for all events is 1.4, and α = 1 is 1.58 worse in log-likelihood, about 1.8σ formally; O4a alone gives 2.2. The stretched mocks show this is within what a true α = 1 produces with these distances and errors: 13 of 50 (26%) find α ≥ 1.4, and their T spreads widely (−2.7 to 21.8). It is not evidence for more stretch than light's, but O4a's high value is the first place to look if the selection for that run is imperfect (for example, events taken from the cumulative GWOSC list against a FAR cut over the searches in the injection file).

## Deviations and implementation choices

1. **Mocks draw each mass through the selection at the true distance**, as the frozen rule says. Entry 79's code drew them through the selection at the catalogue distance; its text was ambiguous. The fits of the mocks use the catalogue distance, as the real fit does.
2. **α grid.** Mocks are profiled on the same grid as the real data (−0.5 to 2.5 in 0.1), finer than Entry 79's mock grid (0.2).
3. **False-alarm-rate units.** The file has no unit attributes. The values (medians 0.16 to 55, many below 1/365, minima down to 10⁻⁹⁵) fit per-year units, which is how they are used; I could not confirm this from the file itself.
4. **Importance weights are very uneven overall** (an effective sample of 69 to 81 per run group over the whole plane), because the reference population differs from the draws in mass ratio, spin and distance. At the events themselves the kernel-weighted effective number is at least 82 (median 1,874; at least 370 for every heavy event), and the halves agree to 0.03, so the estimate is stable where it is used.
5. **P_det** is scaled to a maximum of 1 per run group and floored at 10⁻⁶ at an event's own point, as Entry 79 did. A constant per group cancels in the likelihood.
6. The likelihood is imported unchanged from `code/gw_redshift.py`, pointed at the 155-event subset (`work/gw_rows_o1_o4a.txt`, written by the script).

## Files and how to run

- `gw_injections.py`: the whole test (38 minutes, mostly the 100 mocks; seeds fixed): `.venv/bin/python -I entries/168-gw-injections/gw_injections.py`. Needs `h5py` (3.16.0, installed into `.venv/` for this entry with the user's permission).
- `results.json`: injection counts, the effective-sample diagnostics, the profiles for all, heavy and light (for a figure), the convergence halves, every mock's T, best α and (stretched) heavy h, the decisions and the checks.

## Data

- Events: `code/gw_rows.txt` (Entry 79; GWOSC cumulative GWTC event API, queried 11 November 2026), run codes 1–5 only: 155 binary black holes (O1+O2 10, O3 59, O4a 86).
- Injections: LVK GWTC-4.0 sensitivity estimates, Zenodo record 16740128 (https://zenodo.org/records/16740128), `data/lvk-gwtc4-injections/mixture-semi_o1_o2-real_o3_o4a-cartesian_spins_20250503134659UTC.hdf`, 624,783,664 bytes, sha256 `fdaa4743a169fef3a3ca4071b2fb939af3433df9de8fa9973305fad1186090f1`. Table `events`: 1,697,750 accepted of 1,130,738,485 generated; found under the rule: 369,995 (O1+O2), 221,324 (O3), 474,603 (O4a).

## What the cloud session should write up

Entry 168: **stretch confirmed** (T = 6.57, above every unstretched mock, 4.8 of their standard deviations) and **the heavy-binary residual gone** (h = 0.08 against a stretched-mock median of 2.69, p = 0.98); the decomposition from check (a), half of Entry 79's heavy residual from the O4b events and the rest from its selection model; the note that the best α is 1.4 (O4a 2.2), within the stretched mocks; a figure of the three likelihood profiles (`profiles`) beside the mock distributions of T (`mocks`); a correction note at the end of Entry 79 (its heavy-binary disagreement was the selection model and the O4b events, not the waves); the Gap Map update (gravitational waves are stretched as light is, now with the collaboration's own selection and mock-calibrated at 4.8σ against unstretched waves, which no plasma redshift can produce); a Ledger sentence; and section 168 marked as run. Limits to keep: single middle values for masses and distances, distances inferred assuming expansion, and O4b left out because no injections cover it.
