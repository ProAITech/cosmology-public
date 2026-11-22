# Entry 164: Faraday rotation against quasar redshift

Does the scatter of the rotation measure (RM) grow with quasar redshift, as intergalactic fields or currents would make it?
The rule was frozen in the Workbench (Entry 164, "the rule, frozen first") on 14 November 2026, in commit `c70e257`, before any analysis.

**Result: flat, the standard picture.** Slope of the residual scatter against log10(1+z): b = +1.06 rad/m² per unit log10(1+z); shuffled-redshift null +2.36 ± 3.46; t = −0.38 (rule: |t| < 2 is flat). One-sided 95% limit b + 1.645 sd = 6.75, so the scatter added between z = 0 and z = 1 is at most 0.301 × 6.75 = 2.0 rad/m², against a scatter of about 17.8 rad/m² in every redshift bin. All five checks agree (t between −0.65 and +0.81).

## Files

- `faraday_redshift.js`: the whole analysis. Run in the browser on the VizieR CfA page (`https://vizier.cfa.harvard.edu/viz-bin/VizieR`); data are held in browser memory and nothing is written to disk. Paste the file, call `F164.start()`, poll `F164.state` (about 9 minutes), then `F164.compact(F164.state.result)` gives `results.json`.
- `results.json`: the results (4.5 kB). Binned scatter and the fit for the primary sample and each check; the shuffled-slope null as a histogram (bin width 1 rad/m² per unit log10(1+z), from `lo`), its mean, sd and quantiles (0.13, 2.5, 50, 97.5, 99.87%); t, one-sided p-values, the 95% limit; sample counts.

## Data (VizieR, loaded through `/viz-bin/votable`; no local copy)

- NVSS rotation measures, Taylor, Stil & Sunstrum (2009): `J/ApJ/702/1230/catalog`, constraint `RM=-10000..10000`, columns RAJ2000, DEJ2000, RM, e_RM: 37,543 sources, all usable.
- Million Quasars catalogue, Milliquas: `VII/294`, columns RAJ2000, DEJ2000, Name, Type, z; constraints `z>0`, `DEJ2000>-41`, `Type=Q*` (then A*, B*, K*, N*), in 36 slices of 10° of RA, 180 requests. Candidates kept (first letter of Type upper case): Q 854,823; A 45,546; B 2,009; K 5,909; N 38,121.
- Both catalogues are as served by VizieR on 14 November 2026. Positions come back in decimal degrees (`-oc.form=d`).

## Sample

Type Q with z > 0, nearest NVSS RM source within 10″ (one quasar per NVSS source, the closest): 2,946 matched. Galactic |b| ≥ 20°: 2,688 kept (258 dropped); every one had at least 10 other NVSS sources within 3°. All five classes together: 3,744 matched, 3,410 kept. The median redshift is 1.34 (range 0.12–5.25).

## Deviations from the frozen text

1. **Sample size.** The rule said "about 5,000 are expected". The rule as written gives 2,688. The expectation came from a lost partial run (2,356 matched after half the sky); I did not reproduce it and do not know what selection it used. Nothing was changed to reach 5,000.
2. **Milliquas fetched one Type letter at a time.** A single constraint `Type=Q*|A*|B*|K*|N*` returned no rows at VizieR, so each class has its own request. The result is the same set.
3. **One quasar per NVSS source.** The rule did not say what to do when two Milliquas objects fall within 10″ of one NVSS source; the closest is kept.
4. **Check calibrations use 300 shuffles each.** The rule fixed 1,000 shuffles for the primary calibration and did not give a count for the checks.
5. **Before the freeze**, only three header rows of each catalogue were read, to learn column names and type codes; and one trial query that returned no rows (the combined Type constraint of item 2).

## Notes on reading it

- The null mean (+2.4) is not zero. I did not investigate why. The rule compares the real slope with the null mean, so t is measured from it.
- The scatter in every bin is 16–19 rad/m². I did not separate how much is catalogue noise (e_RM) and how much is Galactic structure that survives the 3° median. A field of 1 nG at the mean baryon density gives about 0.2 rad/m² per gigaparsec, so the test is blind to it by a factor of about ten. It bounds only fields or currents about ten times stronger.
- Chance matches, estimated from densities and not measured: about 26 Milliquas Q objects per square degree and a 10″ radius give about 23 chance matches among the 2,946, under 1%.

## What the cloud session should write up

Entry 164: the result above (flat; t = −0.38; limit 2.0 rad/m² at z = 1), the deviations, a figure of scatter against log10(1+z) with the fit (`primary.bins`) beside the null histogram (`primary.null.hist`) with the real slope marked, a Gap Map update (intergalactic Faraday rotation: no growth, but the limit is a factor ~10 above the standard expectation), a Ledger sentence, and the pre-registered rule marked as run.
