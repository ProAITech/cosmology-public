# Entry 180: the cluster Faraday signal on quasars behind the clusters

**Result: seen, by the frozen rule (Workbench section 180, commit `73d4030`, 17:39 UTC, 21 November 2026).** On sight lines to quasars at least 0.1 in redshift behind the cluster, the rotation-measure scatter is raised inside R500 at 4.8σ, and not beyond. These sources cannot be cluster members, so the excess of Entries 176 and 179 comes from cluster gas in front of the sources.

- Sources: Entry 176's 92,169 SPICE-RACS DR2 sources, of which 11,554 have a quasar redshift (Milliquas 6,063, DESI DR1 1,742, Quaia 3,749).
- Clusters: 10,151 (Entry 176's 9,050 eRASS1 clusters and Entry 179's 1,101 MCXC and PSZ2 clusters).
- A quasar is used only if its redshift exceeds that of its nearest cluster (in units of R500) by more than 0.1.

| range | background quasars (rotated) | V (rad²/m⁴) | null | t | excess |
|---|---|---|---|---|---|
| A: x < 1 | 74 (64) | 216.8 | 79.9 ± 28.6 | **+4.78** | +136.9 |
| B: 1–2 | 171 (184) | 82.4 | 80.1 ± 16.3 | +0.14 | +2.3 |
| C: 2–4 | 688 (671) | 77.3 | 78.6 ± 8.5 | −0.16 | −1.3 |

Decision: A *seen*; B and C *no excess*.

Reported only:

| range | all sources (rotated) | t | excess |
|---|---|---|---|
| x < 1 | 669 (545) | +12.25 | +117.1 |
| 1–2 | 1,700 (1,560) | +0.26 | +1.7 |
| 2–4 | 6,124 (5,650) | −0.41 | −1.5 |

- The table above is every source around all 10,151 clusters on both halves of the sky.
- Per source, the background-quasar excess inside R500 is 1.17 times the all-source excess. A share of the all-source sample sits in front of the clusters or among them and dilutes the excess. The background quasars have neither problem.

**Reading.** Together with Entry 176 (no surplus of sources at 0.5–1 R500) and Entry 179 (the same signal around other clusters on the other half of the sky), this settles where the rotation arises: in magnetised gas inside R500 that lies in front of background sources.

## Files

- `background.py`: `synthetic` and `run`.
- `cluster_rm.py` (Entry 176) and `other_half.py` (Entry 179), copied unchanged.
- `synthetic.json`: gas in front of all sources gives t = +3.0; rotation in members only gives −0.9; no effect gives −0.4.
- `results.json`.
- Data: SPICE-RACS DR2, Milliquas v8, DESI DR1 quasars, Quaia, eRASS1, MCXC and PSZ2, all cited in the script header.
