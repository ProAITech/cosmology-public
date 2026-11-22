# A cosmology inquiry: analysis code and derived results

This repository holds the analysis scripts, the derived results and the small data tables of an inquiry that tests claims made for plasma cosmology, and related open questions, against public data. A run of tests on public catalogues (Entries 164–187) is the basis of six manuscripts in preparation. It is a snapshot taken on 22 November 2026 from the working repository, which is private; its history is not included.

## Layout

| Path | What it holds |
|---|---|
| `entries/NNN-<slug>/` | One folder per larger test (Entries 164–187): its scripts, its results files (`.json`), figures, and a README naming its data, what was fixed in advance, and what it found. |
| `code/` | Scripts and small data tables of the earlier entries (1–163) and of the smaller later ones (for example Entry 170, `code/neutron_li*.py`). `code/INDEX.md` lists them by entry. |
| `DATA.md` | Every raw catalogue and map the tests read: source, size and checksum, so the inputs can be fetched again. |
| `requirements.txt` | The Python packages used. |

## Entries 164–187

| Entry | Question | Outcome | Code | Manuscript |
|---|---|---|---|---|
| 164 | Does the RM scatter of NVSS quasars grow with redshift? | flat | `entries/164-faraday-redshift/` |  |
| 165 | Do 43 voids at z = 0.03–0.06 that nobody had looked at warm the CMB (the test frozen in Entry 148)? | no hot-void effect | `entries/165-smica-voids/` |  |
| 166 | Do the RM steps across Orion A and Perseus follow the clouds? | not decided | `entries/166-spice-clouds/` |  |
| 167 | Is there an Mg II excess at the Giant Arc in DESI DR1? | no excess | `entries/167-desi-arc/` | Giant Arc test |
| 168 | Are gravitational waves stretched, with the LVK injection selection? | stretch confirmed | `entries/168-gw-injections/` |  |
| 169 | Is the Big Ring in DESI DR1 absorbers? | not tested | `entries/169-desi-ring/` |  |
| 170 | Can nucleons from evaporating black holes cut lithium (Eur. Phys. J. C 85, 1362; arXiv:2606.03316)? | no (theory check) | `code/neutron_li.py` | Lithium comment |
| 171 | Does the RM scatter of SPICE-RACS quasars grow with redshift? | growth claimed | `entries/171-spice-faraday-redshift/` | Faraday rotation and redshift |
| 172 | Does a current thread other galaxies (SPICE-RACS)? | unresolved | `entries/172-spice-galaxy-currents/` |  |
| 173 | Do Mg II absorber sight lines carry the growth (DESI DR1)? | unresolved | `entries/173-spice-desi-absorbers/` | Faraday rotation and redshift |
| 174 | Does the growth appear on sight lines not used before (Quaia)? | unresolved, narrowly | `entries/174-quaia-replication/` | Faraday rotation and redshift |
| 175 | Do quasar polarizations align once Planck’s dust is removed? | Galactic | `entries/175-quasar-polarization-dust/` | Quasar polarization alignments |
| 176 | How far from clusters does Faraday rotation show their plasma? | inside R500 only | `entries/176-cluster-faraday-profile/` | Cluster magnetic fields |
| 177 | Does the RM difference within double sources grow with redshift? | unresolved, narrowly | `entries/177-pair-differences-redshift/` | Faraday rotation and redshift |
| 178 | Has the cluster field changed since z ≈ 0.46? | unresolved | `entries/178-cluster-field-evolution/` | Cluster magnetic fields |
| 179 | Does the cluster signal repeat on the other half of the sky? | replicated | `entries/179-cluster-faraday-other-half/` | Cluster magnetic fields |
| 180 | Is the cluster rotation in front of the sources? | seen | `entries/180-cluster-faraday-background/` | Cluster magnetic fields |
| 181 | Can intervening haloes make the Faraday growth? | no (budget) | `entries/181-halo-faraday-budget/` | Faraday rotation and redshift |
| 182 | Does region A1 repeat on quasars measured after 2005? | not repeated | `entries/182-a1-new-quasars/` | Quasar polarization alignments |
| 183 | Can source luminosity be separated from redshift in the growth? | not with this sample (power study) | `entries/183-faraday-luminosity/` | Faraday rotation and redshift |
| 184 | Is Entry 113’s helium neutron clock right? | confirmed (check) | `code/neutron_clock_prymordial.py` |  |
| 186 | Is the cluster Faraday signal seen at 1.4 GHz (NVSS)? | no decision | `entries/186-nvss-clusters/` |  |
| 187 | Is the CMB colder around spirals at z = 0.02–0.05? | excluded | `entries/187-cmb-spirals/` || CMB around spirals |

Entry 170 also uses `code/neutron_li_primat.py`, `code/neutron_li_shapes.py`, `code/neutron_li_prymordial.py` and `code/neutron_li_protons.py`. Entry 184 is `code/neutron_clock_prymordial.py`. Entry 185 is a summary table and has no code.

## Running

- Python 3.12 with the packages in `requirements.txt` (numpy, scipy, healpy, astropy, h5py). The scripts in `code/` need only numpy and scipy.
- Fetch the raw data listed in `DATA.md` into `data/`.
- Scripts in `entries/` run from the repository root, for example `python3 entries/176-cluster-faraday-profile/cluster_rm.py`. Many keep intermediate files in `entries/*/work/`, which git ignores; each README gives the order in which its scripts run. Scripts in `code/` run from inside `code/`.
- A script prints the numbers its entry quotes and writes its results file beside itself. Running it again overwrites that file with the same numbers.
- The `.js` files were run in a browser tab against live catalogue services and are kept as a record of the query.
- A few earlier scripts read the inquiry's pages (`forecasts2.py` and `scoreboard92.py` read the Gap Map), which are not part of this repository.

## Pre-registration

Most tests from Entry 164 on were decided by a rule (the statistic, the null test and the thresholds) fixed before the data were combined, in a dated commit to the private working repository. Each entry's README states what was fixed and what had already been read. The dated records are not part of this repository; they are available from the corresponding author on request.

## Contact

Tony Tong, Pro AI Technology. Write through the issues of this repository.
