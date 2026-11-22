# Scripts and data tables of the earlier entries

This index was written for the earlier entries (1–163), whose write-ups are in the inquiry's Workbench, a page not included here; "page section" names the Workbench section of each entry. The entries from 164 on are in `entries/` and are listed in the top-level README.

### How to run

- Python 3 with `numpy` and `scipy`. No other packages.
- Run a script from inside the `code/` folder, for example `cd code` then `python3 s8_board.py`. Scripts read their data tables from the folder they sit in, and some import a neighbour (`cosmo.py`, `dark_energy_fit.py`, `bbn_net.py`, `static_hubble.py`, `pantheon_rows.py`, `deuterium_yield.py`, `gravity_beta.py`).
- A script prints the numbers quoted in its entry. Most also write a `.json` file of results next to themselves.
- The `.js` files were run in a browser tab against live public catalogues (VizieR, SDSS SkyServer, Gaia, GWOSC, LoTSS, POSSUM and others). They are kept as a record of what was computed, and need those catalogue services to run again. Each one's header names the service and the table.
- The data tables were typed in or fetched from published tables. Most carry a first line saying what the columns are and where they come from.

### Run order and arguments

Eight results files are read by another script, so run the first before the second:

| First | Then |
|---|---|
| `deuterium_yield.py` and `lithium_thick_target.py` | `deuterium_budget.py` |
| `distance_test_repeat.py` | `dimming_dilemma.py` |
| `loop_flux.py` | `flux_budget.py` |
| `clumping_fit.py` | `forecasts.py` |
| `dynamical_mass_test.py` | `trend_checks.py` |
| `static_distance_law.py` | `gas_fraction_test.py` |
| `cmbpol.py` (with the arguments below) | `pol_compare.py` |
| `cmbreion.py`, three runs (below) | `reion_compare.py` |

The results files these second scripts need are included, so each can also be run on its own.

Scripts that take arguments:

- `python3 cmbpol.py 0.0222 0.2975 0.6851 0.965 noncmb`: baryon density, matter share, h and spectral index, all from measurements that do not use the CMB, then a label. It writes `cmbpol_noncmb.json`, which `pol_compare.py` reads.
- `TAU_RE=0.054 python3 cmbreion.py 0.0222 0.2975 0.6851 0.965 tau054`: the same parameters, with the optical depth to reionization set by `TAU_RE`. Entry 94 uses three runs, labelled `tau0`, `tau054` and `tau09` (`TAU_RE` = 0, 0.054 and 0.09); each writes `cmbreion_<label>.json` (and a 5 MB `srcR_<label>.npz` of intermediate sources, which nothing else reads), and `reion_compare.py` reads all three.
- `python3 gw_redshift.py gw_rows.txt [n_mocks] > gw_redshift.json`: the event table, then the number of mock catalogues (30 by default).
- `python3 spt_xmm_distance.py data` and `python3 selection_sim_xmm.py data`: the folder that holds `bulbul_spt.txt`.
- `distance_test_repeat.py`, `ratio_selection_sim.py` and `sz_xray_distance.py` take an optional data folder (default: this one).
- `python3 deuterium_yield.py quick` runs a shorter version; a number instead sets the samples per energy (2,500 by default).
- `python3 raw_chunks.py <folder>` and `python3 raw_chunks3.py <folder>`: these rebuild six data tables from the rows as first typed in, and write them to the folder you name. The output is identical to the tables here.

Pages one level up:

- `forecasts2.py` and `scoreboard92.py` read `gapmap.html` from the folder above `code/`, so they count today's Gap Map.
- `corrections_check.py` and `corrections_de.py` read nine small tables from `code/corr/`.

### Running time

Most scripts finish in seconds. These take minutes:

- Two and a half to five minutes each: `dark_energy_fit.py`, `deuterium_yield.py`, `followup.py`, `hull_significance.py`, `pantheon_plus_test.py`, `plasma_checks.py`, `sn_shift_scan.py`, `tail_scan.py`; about 3 minutes, `ccc_tl_test.py`; about 2, `ccctl_bbn.py`.
- About 6 minutes: `cmbpol.py`. About 8 minutes for each `cmbreion.py` run.
- `gw_redshift.py` took 4 minutes with 2 mock catalogues; the published run used the default 30, so allow roughly 40 minutes to an hour, or pass a smaller number.
- `tail_robust.py` imports `tail_scan.py`, which repeats that scan, so allow ten minutes or more.

Timings are on one core of a small cloud machine.

### Scripts by entry

Each entry's files, in the order the entry links them. The description is the first lines of the file's own header.


**Entry 15 · Internal review and corrections** (page section `#review`)

- `code/hybrid_redshift.py`: Hybrid redshift: 1+z = (1+z_c)(1+z_t). z_c from expansion (FLRW, matter + Lambda), z_t from photon energy loss. Family A: constant log-fraction, ln(1+z_t) = f ln(1+z).
- `code/wide_binaries_model.js`: Forward model for the Gaia wide-binary test (Plasma Cosmology Workbench, Entry 15).
- `code/filament_foreground_check.js`: Filament-current test without the Galactic foreground map (Plasma Cosmology Workbench, Entry 15).

**Entry 16 · Precession, Sirius and the force budget** (page section `#precession`)

- `code/precession.py`: Precession: the standard torque account versus motion around a current filament or a Sirius companion (Workbench Entry 16).

**Entry 17 · Arp's ejected quasars: a motion test** (page section `#arp`)

- `code/arp_proper_motion.js`: Test of the claim that quasars are ejected from nearby galaxies (Arp; Burbidge; Bell).
- `code/arp_model.py`: What ejection of quasars from nearby galaxies predicts for Gaia, and what the measured outward drift allows.

**Entry 18 · The Milky Way's halo circuit: what drives 10¹⁷ A** (page section `#circuit`)

- `code/halo_circuit.py`: The Milky Way's halo circuit. What drives ~1e17 A, and is it fed from outside?

**Entry 19 · Is there a universal current direction? A correction and a first count** (page section `#sign`)

- `code/halo_sign_stack.js`: Does the rotation measure of background radio sources seen through the halos of nearby inclined galaxies carry the four-quadrant sign pattern that outward axial currents would give (Myserlis & Contopoulos 2021)?
- `code/sign_test.py`: Is there a universal direction for galactic currents?

**Entry 20 · Do young stars leave the gas behind? A test of magnetic rotation curves** (page section `#gasstars`)

- `code/gas_vs_stars.js`: If magnetic forces, which act on gas but not on stars, supplied part of the centripetal force that keeps the Galaxy's gas on its orbit, stars would leave the gas behind as soon as they formed: they would stream outward, lag the gas, and scatter.
- `code/gas_vs_stars_model.py`: Can a force that acts on gas but not on stars (magnetic support) hold up the Galaxy's rotation curve?

**Entry 21 · A variable speed of light: timing passed, geometry failed** (page section `#vsl`)

- `code/vsl_shape.py`: A variable speed of light as the cause of redshift (Unzicker's Dirac-Sciama-Dicke cosmology, arXiv:0708.3518).

**Entry 22 · Is intergalactic space a radio fog?** (page section `#fog`)

- `code/radio_fog.js`: Is intergalactic space a radio fog? Lerner's reading of the CMB needs the space between galaxies to absorb or scatter radio waves.
- `code/sz_shadows.js`: Workbench Entry 22, Part 7: do the shadows that galaxy clusters cast on the CMB fade with distance? A cluster's hot gas leaves a cool spot in the CMB (the Sunyaev-Zeldovich effect) at 98-150 GHz.
- `code/fog_model.py`: Is intergalactic space a radio fog? Derived numbers. Lerner's reading of the CMB needs the intergalactic medium to absorb or scatter radio waves between roughly 150 MHz and 200 GHz.
- data: `code/fog_results.json`

**Entry 23 · Could stars have made the helium?** (page section `#helium`)

- `code/helium_budget.py`: Could stars have made the helium? An energy account.

**Entry 24 · Were galaxy clusters denser in the past? Temperature against shadow depth** (page section `#denser`)

- `code/cluster_density_test_v2.py`: Were galaxy clusters denser in the past? A test that uses no distances.
- `code/cluster_density_test.py`: Workbench Entry 24, second version: were galaxy clusters denser in the past?
- `code/selection_sim.py`: Entry 24, second version: how much of the trend is the selection of the sample?
- `code/sz_shadows.js`: Workbench Entry 22, Part 7: do the shadows that galaxy clusters cast on the CMB fade with distance? A cluster's hot gas leaves a cool spot in the CMB (the Sunyaev-Zeldovich effect) at 98-150 GHz.
- data: `code/clusters_y_T_v2.txt`, `code/clusters_y_T.txt`, `code/temps_matched.txt`, `code/act_signals.txt`, `code/first_sample_act.txt`, `code/act_q_curve.txt`

**Entry 25 · Deuterium and lithium from cosmic rays: an objection that does not hold** (page section `#lithium`)

- `code/lithium_thick_target.py`: Lithium and deuterium from cosmic rays, thin target against thick target. The classic objection to making deuterium with cosmic rays: the same cosmic rays make far too much lithium (alpha + alpha -> 6Li, 7Li, 7Be).
- `code/gamma_albedo.py`: Workbench Entry 25, Part 3: do stars hide the gamma rays? Making deuterium with cosmic rays also makes neutral pions, which decay to gamma rays.

**Entry 26 · El Gordo: the conflict sits in a fifth of the collision speed** (page section `#elgordo`)

- `code/elgordo_infall.py`: El Gordo.

**Entry 27 · The shadow and the glow: a distance test with 66 clusters** (page section `#shadowglow`)

- `code/sz_xray_distance.py`: The shadow-plus-X-ray distance test. A cluster's hot gas is seen twice. Its shadow on the CMB (ACT DR5, Hilton et al. 2021) measures the pressure summed along the line of sight and does not care how far away the cluster is.
- `code/raw_chunks.py`: Raw strings copied from the browser session (VizieR tab), in the order they were read out. Joined and written to data files.
- data: `code/sanders_act_profiles.txt`, `code/act_signals.txt`, `code/temps_matched.txt`

**Entry 28 · The distance test repeated with two other pairs of telescopes** (page section `#repeat`)

- `code/distance_test_repeat.py`: The shadow-against-glow distance test, repeated with two other pairs of instruments. Entry 27 compared ACT shadows with Chandra profiles (Sanders et al. 2018) for 66 clusters.
- `code/ratio_selection_sim.py`: Does the way the clusters were chosen bend the shadow-against-glow test?
- `code/raw_chunks3.py`: Raw strings read out of the browser session (VizieR), round 21. Joined, checked and written to data files.
- data: `code/sanders_spt.txt`, `code/planck_lovisari.txt`

**Entry 29 · Deuterium by neutron capture: the model's yield holds** (page section `#dyield`)

- `code/deuterium_yield.py`: How many deuterons does a cosmic ray make when it stops in dense hydrogen and helium?
- `code/deuterium_budget.py`: Workbench Entry 29, second part: variants for a pure-hydrogen target, other spectra, and what the yield means for the model's budget.

**Entry 30 · The nitrogen that comes with stellar helium** (page section `#nitrogen`)

- `code/nitrogen_check.py`: The nitrogen that comes with helium made in intermediate-mass stars. The helium-from-stars model (Lerner's GOLE paper, section 4) makes the helium in stars of 4-12 solar masses, which return it to the gas in winds.

**Entry 31 · The static universe's Hubble law against the supernovae** (page section `#snlaw`)

- `code/static_hubble.py`: The static universe's Hubble law against the supernovae.
- `code/dimming_dilemma.py`: Workbench Entry 31, part 3: is there any dimming law that saves a static distance?
- `code/pantheon_check.py`: Workbench Entry 31, check: the static Hubble laws against the Pantheon supernova sample (Scolnic et al. 2018; 1,048 supernovae, z = 0.010 to 2.26).
- `code/pantheon_rows.py`: The Pantheon supernova sample (Scolnic et al. 2018), read from VizieR J/ApJ/859/101/fullz in eleven redshift slices, sorted by redshift.

**Entry 32 · The distance test to redshift 1.5, with a second X-ray telescope** (page section `#third`)

- `code/spt_xmm_distance.py`: The shadow-against-glow distance test with a third X-ray data set, reaching redshift 1.5.
- `code/selection_sim_xmm.py`: Does selection bend the SPT-with-XMM test? The simulation of Entry 28, rerun for this sample. The 57 clusters reach down to SPT significance 4.7, lower than the Chandra sample's 6.8, so noise-boosted detections matter more.
- data: `code/data/bulbul_spt.txt`

**Entry 33 · The static universe's distance law, measured instead of assumed** (page section `#measured`)

- `code/static_distance_law.py`: The static universe's distance law, measured instead of assumed. Earlier entries tested particular static distance laws (d = cz/H0, d = (c/H0) ln(1+z)).
- `code/static_hubble.py`: The static universe's Hubble law against the supernovae.
- `code/pantheon_rows.py`: The Pantheon supernova sample (Scolnic et al. 2018), read from VizieR J/ApJ/859/101/fullz in eleven redshift slices, sorted by redshift.
- data: `code/entry28_points.json`, `code/entry32_points.json`

**Entry 34 · Gas fractions and radio cores: two more tests of the last static model** (page section `#fgas`)

- `code/gas_fraction_test.py`: A second test of the one static model left by Entry 33, with cluster gas fractions.
- `code/radio_core_ruler.py`: Workbench Entry 34, Part 3: compact radio cores as a ruler (120 quasars, Cao et al. 2017). theta x D_M(z) should go as (1+z)^x: x = 1 in expanding space if the cores have one physical size, x = 0 in static space with the measured distance law.
- data: `code/mantz2022_fgas.txt`, `code/mantz2014_fgas.txt`, `code/cao2017_quasars.txt`

**Entry 35 · Alfvén and Klein's antimatter universe, in numbers** (page section `#antimatter`)

- `code/alfven_klein.py`: Alfven and Klein's antimatter universe, in numbers.
- `code/turnaround_density.py`: Workbench Entry 35, Part 4: the Alfven-Klein cloud at the turnaround density its 1962 paper states (about 1e-2 atoms per cm3).

**Entry 36 · Is dark energy changing? DESI's signal rebuilt and taken apart** (page section `#darkenergy`)

- `code/dark_energy_fit.py`: Is dark energy changing? DESI's signal rebuilt from published numbers, and taken apart. Models: flat LCDM (w = -1), and w(a) = w0 + wa (1 - a).
- `code/cosmo.py`: Distances and likelihoods for Workbench Entry 36: flat models with dark energy w(a) = w0 + wa (1 - a), photons, and neutrinos (N_eff = 3.044, one of 0.06 eV).
- `code/static_hubble.py`: The static universe's Hubble law against the supernovae.
- `code/pantheon_rows.py`: The Pantheon supernova sample (Scolnic et al. 2018), read from VizieR J/ApJ/859/101/fullz in eleven redshift slices, sorted by redshift.
- `code/clumping_fit.py`: Workbench Entry 36, Part 4: the free factor on the sound horizon replaced by the clumping model of Entry 13.
- data: `code/clump_track.json`

**Entry 37 · The dark-energy signal: one galaxy sample, and the Hubble constant as referee** (page section `#referee`)

- `code/followup.py`: Two follow-ups on the dark-energy signal of Entry 36. 1. The galaxies at redshift 0.71 carry half of it.
- `code/cosmo.py`: Distances and likelihoods for Workbench Entry 36: flat models with dark energy w(a) = w0 + wa (1 - a), photons, and neutrinos (N_eff = 3.044, one of 0.06 eV).
- `code/dark_energy_fit.py`: Is dark energy changing? DESI's signal rebuilt from published numbers, and taken apart. Models: flat LCDM (w = -1), and w(a) = w0 + wa (1 - a).

**Entry 38 · What to measure next: the nearby supernovae, and five forecasts** (page section `#next`)

- `code/sn_shift_scan.py`: Entry 38, block 5 in finer steps: dark-energy preference against a brightness shift of the supernovae nearer than z = 0.1.
- `code/sn_shift_check.py`: Two checks on the supernova scan. 1. The matter share a supernova-only fit returns for each shift, to compare with the published supernova-only values. 2.
- `code/sn_seam.py`: Does the scan depend on where "nearby" ends? The same shift applied below z = 0.05 and below z = 0.2, and the offset the nearby supernovae show against the curve the distant ones prefer.
- `code/forecasts.py`: What to measure next. Forecasts from the notebook's own fits. Each block takes a result of Entries 31-37, adds hypothetical data, and asks how much is needed to decide the question it left open.
- `code/cosmo.py`: Distances and likelihoods for Workbench Entry 36: flat models with dark energy w(a) = w0 + wa (1 - a), photons, and neutrinos (N_eff = 3.044, one of 0.06 eV).
- `code/dark_energy_fit.py`: Is dark energy changing? DESI's signal rebuilt from published numbers, and taken apart. Models: flat LCDM (w = -1), and w(a) = w0 + wa (1 - a).

**Entry 39 · The supernova seam, measured in Pantheon+: dark energy and the static law** (page section `#seam`)

- `code/pantheon_plus_test.py`: The supernova test of Entry 38 on the Pantheon+ table itself. Data: 1,590 light curves with zHD > 0.01 from Brout et al. (2022), table 7 (VizieR J/ApJ/938/110): zHD, corrected magnitude, and the diagonal error.
- `code/pantheon_plus_scan.py`: The nearby-supernova scan of Entry 38 run on Pantheon+ in finer steps, and both samples binned against one curve.
- `code/pantheon_plus_ps1.py`: Entry 39, part 3: the seam inside one telescope. Foundation (nearby) and the Pan-STARRS Medium Deep supernovae (distant) were both observed with the Pan-STARRS1 camera.
- `code/static_law_plus.py`: Entry 39, part 3: Entry 31's test of the static Hubble laws, rerun on Pantheon+ (diagonal errors, repeats averaged). Every law gets a free additive constant.
- data: `code/pantheon_plus_rows.txt`, `code/pantheon_plus_foundation.txt`, `code/pantheon_plus_ps1md.txt`

**Entry 40 · Pantheon+ with no survey trusted against another** (page section `#ties`)

- `code/pantheon_plus_surveys.py`: Pantheon+ with no survey's brightness level trusted against another's. Each survey (or each camera) gets its own free level, so only the shape of the Hubble diagram inside a survey is used.
- `code/pantheon_plus_expected.py`: What each curve predicts for the power q measured inside one survey. A curve that bends (the standard one) gives a different q in each redshift range; a pure power law (Lerner's static law, q = 1/2) gives the same q everywhere.
- data: `code/survey_labels.txt`

**Entry 41 · The same-camera tie with its corrections undone** (page section `#undo`)

- `code/corrections_check.py`: Could the corrections applied to the supernova magnitudes create the same-camera step that sinks the static law?
- `code/corrections_de.py`: Entry 41, second part: the same-camera pair with and without the model corrections, in the two fits of Entry 40.
- data: `code/lightcurve_foundation.txt`, `code/lightcurve_panstarrs.txt`

**Entry 42 · A fifth pair of instruments: clusters not chosen by their shadows** (page section `#fifth`)

- `code/boxsz_test.py`: The shadow-against-glow distance test with a fifth pair of instruments, and clusters not chosen by their shadows. Data: Czakon et al.
- data: `code/boxsz_table3.txt`

**Entry 43 · What moved the anchor: the same supernovae, five hundredths of a magnitude** (page section `#moved`)

- `code/anchor_shift.py`: What moved the nearby anchor between the 2018 supernova compilation and Pantheon+? Supernovae in both tables are matched by name.
- data: `code/near2018.txt`, `code/nearplus.txt`, `code/far2018.txt`, `code/farplus.txt`

**Entry 45 · The line-width test on 71 galaxies, with a nearby anchor** (page section `#widths`)

- `code/dynamical_mass_test.py`: Stellar against dynamical mass in 71 quiescent galaxies at z = 0.8 to 2.2, under expanding and static distances. A galaxy's velocity dispersion comes straight from its spectrum and needs no distance.
- `code/trend_checks.py`: Entry 45, checks on the trend inside the distant sample: what drives it, and does galaxy size account for it? Run dynamical_mass_test.py first; it writes the json file this script reads.
- data: `code/vds2013_table4.txt`

**Entry 46 · Which frame is at rest? Supernovae, the CMB dipole and the count dipoles** (page section `#frame`)

- `code/sn_frame.py`: Which frame is at rest? The Sun's motion measured from nearby supernovae alone, against the CMB dipole.
- data: `code/pantheon_plus_nearby_coords.txt`

**Entry 47 · A fog at one part in 100,000? The claimed cooling of the CMB near nearby galaxies** (page section `#cooling`)

- `code/cooling_test.py`: The claimed cooling of the CMB near nearby galaxies and filaments, and energy loss in plasma.

**Entry 48 · Counting the numbers: a shorter sound horizon against changing dark energy** (page section `#knobs`)

- `code/knob_count.py`: Counting the numbers. Published and our own fits that ease the mismatch between DESI's ruler and the CMB, either by changing recombination (a shorter sound horizon) or by changing dark energy.

**Entry 49 · The line-width test on a second sample: LEGA-C at redshift 0.6 to 1** (page section `#legac`)

- `code/legac_test.py`: The line-width test on a second, independent sample at redshift 0.6 to 1.
- data: `code/legac_box.txt`, `code/cosmos_box.txt`

**Entry 50 · Fifty entries: the scoreboard, and the measurements that would decide what is left** (page section `#fifty`)

- `code/forecasts2.py`: Two new forecasts and the scoreboard. 1. How many galaxies a stack needs to find, or rule out, a 27 uK cooling of the CMB around spirals at redshift 0.05 to 0.3.

**Entry 51 · The last opening for static sizes: can older stars hide the mass?** (page section `#agecap`)

- `code/age_split.py`: Can older stars hide the mass that static distances need? The published stellar masses cap each galaxy's age at the age of the expanding universe at its redshift.

**Entry 52 · How often do random points make a Giant Arc?** (page section `#arcs`)

- `code/random_arcs.py`: How often do random points make a 'Giant Arc'?
- `code/hull_significance.py`: Workbench Entry 52, second part: the published 'convex hull' significance, applied to random points. The method is that of Clowes et al.
- `code/slab_check.py`: Workbench Entry 52, check: can plain random points reproduce the simulation paper's own numbers? Sawala et al.
- `code/reply_check.py`: Workbench Entry 52, check: the reply's 'nowhere to be seen', on random points.
- `code/breakup.py`: Workbench Entry 52, check: what happens to a random chain of about 44 when the linking length is changed?

**Entry 53 · Lithium: can the primordial plasma fix it?** (page section `#tails`)

- `code/bbn_net.py`: A small Big Bang nucleosynthesis network (Workbench Entry 53). Eight species (n, p, d, t, 3He, 4He, 7Li, 7Be), the neutron-proton weak rates, and the eleven reactions that carry the light elements.
- `code/tail_scan.py`: Can a non-thermal tail in the primordial plasma fix lithium? Charged nuclei must tunnel, so their reaction rates are set by the few nuclei several kT up the tail.
- `code/tail_robust.py`: Workbench Entry 53, robustness of the tail scan: flat S-factors instead of the adopted shapes, and the lithium cut that deuterium allows.
- `code/plasma_checks.py`: Workbench Entry 53, checks. (1) How thermal is the plasma during nucleosynthesis? Coulomb collision rates against nuclear reaction rates and the expansion rate.

**Entry 54 · The Galaxy's vertical field: where does the flux come from?** (page section `#vflux`)

- `code/flux_budget.py`: The Milky Way's vertical magnetic flux, and what could supply it. Entry 18 showed that the halo's field rings are the Galaxy's own vertical field wound up by rotation, and left the source of that vertical field open.
- `code/loop_flux.py`: Workbench Entry 54, check: how much net magnetic flux does a tangled field put through a large loop? A field made by causal processes in the early universe is tangled on a small scale lam.

**Entry 55 · The lopsided counts: an erratum, a new measurement, and the slice test corrected** (page section `#slices`)

- `code/slices.py`: The count dipole after one erratum and one new measurement, and what redshift slices would show. (A) The Quaia quasar dipole after its authors' erratum.

**Entry 56 · The Hubble constant as referee, with the 2026 lens result** (page section `#referee2`)

- `code/referee2.py`: The Hubble constant as referee, updated. Entry 37 compared what each reading of the DESI-CMB mismatch predicts for the Hubble constant with four direct measurements.

**Entry 57 · A current through each galaxy: the Faraday test, run on what exists** (page section `#through`)

- `code/through_current_stack.js`: Does the rotation measure of background radio sources behind the halos of nearby inclined galaxies carry the pattern that a current running THROUGH each galaxy along its spin axis would give (Peratt's filament picture)?
- `code/forecast_through.py`: How big is the rotation-measure signal of a current running through a galaxy along its spin axis, and what would it take to see it? 1.

**Entry 58 · Age-corrected supernovae and the logarithmic law: a claim from August 2026** (page section `#loglaw`)

- `code/age_bias_referee.py`: A 2026 paper (Nguyen, arXiv:2608.00806) reports that supernovae corrected for a progenitor-age bias (Son et al. 2025) prefer the logarithmic law d_L = (c/H0)(1+z) ln(1+z) over flat LCDM.
- `code/des_dovekie_fits.js`: The DES five-year supernovae with the 2025 "Dovekie" calibration (Popovic et al., arXiv:2511.07517), read in a browser tab from the collaboration's public table github.com/des-science/DES-SN5YR  4_DISTANCES_COVMAT/DES-Dovekie_Metadata.csv   (1,820 rows,…
- `code/ruler_by_survey.py`: Workbench Entry 58, check added 11 November 2026 Does the clustering ruler exclude the logarithmic law (steady expansion) in each survey on its own?

**Entry 59 · Do JWST galaxies keep one size? Two static-space claims and a test with line widths** (page section `#jwstsize`)

- `code/line_width_jwst.py`: Two 2025-26 papers (Raikov, Tsymbal & Lovyagin, arXiv:2507.19651; Tsymbal, Raikov & Lovyagin, arXiv:2604.27867) find that JWST galaxies keep one linear size at all redshifts if distance follows static tired light, and ask for a test of sizes and masses…
- `code/astrodeep_sizes.js`: The size-redshift result of Raikov, Tsymbal & Lovyagin (arXiv:2507.19651) re-derived from the catalogue they used: ASTRODEEP-JWST (Merlin et al. 2024), VizieR J/A+A/691/A240, seven fields.

**Entry 60 · A current tied to spin: the signed test on 111 spirals** (page section `#spinsign`)

- `code/signed_stack.js`: The through-current test of Entry 57 with each galaxy's spin sense added, so that the pattern can be stacked with its sign.
- `code/through_current_stack.js`: Does the rotation measure of background radio sources behind the halos of nearby inclined galaxies carry the pattern that a current running THROUGH each galaxy along its spin axis would give (Peratt's filament picture)?

**Entry 61 · The southern survey's first fields: the spin test, galaxy by galaxy** (page section `#possum`)

- `code/possum_through.js`: The through-current test of Entries 57 and 60, galaxy by galaxy, on the first public POSSUM rotation measures (Anderson et al. 2024, VizieR J/MNRAS/533/4068: 22,817 measurements in 35 ASKAP fields, about 30 per square degree, typical error 1 rad/m2).
- data: `code/possum_results.json`

**Entry 62 · Line widths, audited: where they rule out static space and where they do not** (page section `#accel`)

- `code/local_relations_z5.py`: In a universe that is static and does not evolve, galaxies at any redshift are the same population as nearby ones, so they must lie on the relations that nearby galaxies obey.
- data: `code/local_relations_z5.json`

**Entry 63 · The quiescent test at redshift 3 to 4.7: static laws stay off, and the standard reading drifts** (page section `#quiescent3`)

- `code/quiescent_z3.py`: Entry 62 audited the line-width test against the acceleration relation of nearby galaxies and said where the audit should bite hardest: quiescent galaxies beyond redshift 3, which are denser than those at redshift 1 to 2, so that a static distance law…
- data: `code/quiescent_z3.json`

**Entry 64 · Pole currents: is the Milky Way's halo circuit there in other galaxies?** (page section `#poles`)

- `code/pole_currents.js`: Entry 61 fitted, galaxy by galaxy, the pattern of a current running THROUGH a galaxy along its spin axis (the same ring sense above and below the disc).
- `code/possum_through.js`: The through-current test of Entries 57 and 60, galaxy by galaxy, on the first public POSSUM rotation measures (Anderson et al. 2024, VizieR J/MNRAS/533/4068: 22,817 measurements in 35 ASKAP fields, about 30 per square degree, typical error 1 rad/m2).
- data: `code/pole_currents.json`

**Entry 65 · CCC+TL after its press release: the shape test it has not been given** (page section `#ccctl`)

- `code/ccc_tl_test.py`: The CCC+TL model (Gupta 2023, MNRAS 524, 3385): part of the redshift comes from expansion with "covarying coupling constants", part from tired light:  1 + z = (1 + z_x)(1 + z_t).
- data: `code/ccc_tl_test.json`

**Entry 66 · The Giant Arc's recipe on real absorbers: a second catalogue and 53 other volumes** (page section `#arcreal`)

- `code/arc_fields.js`: The Giant Arc's own recipe, run on real absorbers from an independent catalogue, in the Arc's field and in 53 other volumes of the same sky.
- `code/band_variance.py`: How much should the density of a volume the size of the Giant Arc's band differ from the mean, in the standard model, and how does that compare with the counting noise of the absorber catalogue?
- data: `code/arc_fields.json`

**Entry 67 · CCC+TL and the light elements: a change of clock, and two branches** (page section `#cccbbn`)

- `code/ccctl_bbn.py`: CCC+TL's nucleosynthesis claim, and its second branch on the shape test.
- `code/bbn_net.py`: A small Big Bang nucleosynthesis network (Workbench Entry 53). Eight species (n, p, d, t, 3He, 4He, 7Li, 7Be), the neutron-proton weak rates, and the eleven reactions that carry the light elements.
- `code/ccc_tl_test.py`: The CCC+TL model (Gupta 2023, MNRAS 524, 3385): part of the redshift comes from expansion with "covarying coupling constants", part from tired light:  1 + z = (1 + z_x)(1 + z_t).
- data: `code/ccctl_bbn.json`

**Entry 68 · Filament currents on a second, denser grid** (page section `#filspice`)

- `code/filament_spice.js`: The filament sign-step test of Entries 13 and 15, rerun on the SPICE-RACS first release.
- `code/combine_filaments.py`: The node-oriented filament current from two independent rotation-measure grids, in units of 1e17 A (n_e = 1e-5 cm^-3 within 1 Mpc of the spine).
- data: `code/filament_spice.json`

**Entry 69 · Does the gas density decide the filament test? A correction** (page section `#fildens`)

- `code/density_free.py`: Does the filament gas density decide the filament-current test of Entry 68?
- data: `code/density_free.json`

**Entry 71 · Is the Sun powered from outside? The electric Sun's power budget** (page section `#esun`)

- `code/electric_sun.py`: Could the Sun be powered from outside by an electric current (Juergens 1972, 1982; the 'electric Sun')?
- data: `code/electric_sun.json`

**Entry 72 · Electric comets against Rosetta** (page section `#ecomet`)

- `code/electric_comet.py`: The 'electric comet' (a rocky body discharging in the Sun's field, its water made from solar-wind protons and oxygen stripped from surface minerals) against Rosetta's measurements of comet 67P at perihelion (13 August 2015, 1.24 au).
- data: `code/electric_comet.json`

**Entry 73 · Currents along star-forming filaments? What the field reversals allow** (page section `#clouds`)

- `code/cloud_currents.py`: If the line-of-sight field reversals across Orion A, California and Perseus were the toroidal field of an axial current, how large would the current be, is it physically allowed, and what measurement would tell a current from a bent ambient field?
- data: `code/cloud_currents.json`

**Entry 75 · The radio background, and a synchrotron microwave background from currents** (page section `#rsb`)

- `code/radio_background.py`: Peratt's synchrotron microwave background, and whether current-carrying cosmic filaments could supply the measured excess radio synchrotron background (RSB).  cgs units.
- data: `code/radio_background.json`

**Entry 76 · Quantized quasar redshifts: Karlsson's periodicity in 749,742 quasars** (page section `#karlsson`)

- `code/karlsson_dr16q.js`: Is there a Karlsson periodicity, Δlog10(1+z) = 0.089 with peaks at z = 0.06, 0.30, 0.60, 0.96, 1.41, 1.96, in the redshifts of 749,742 SDSS DR16 quasars (Lyke et al. 2020, VizieR VII/289/dr16q, column z, QSO = 1)?
- data: `code/karlsson_dr16q.json`

**Entry 77 · Quantized galaxy redshifts: 37.6 and 72 km/s in ALFALFA** (page section `#tifft`)

- `code/tifft_alfalfa.js`: Galaxy redshift periodicity (Tifft's 72 km/s; Guthrie & Napier's 37.6 km/s in galactocentric velocities of nearby spirals) tested on ALFALFA HI velocities (Haynes et al. 2018, VizieR J/ApJ/861/49/table2, Vhel < 3000 km/s, HI code 1).
- data: `code/tifft_alfalfa.json`

**Entry 78 · Are companion galaxies redshifted? 199,443 members of SDSS groups** (page section `#companions`)

- `code/companion_redshifts.js`: Are companion galaxies redshifted relative to their dominant galaxy (Arp 1994: all 22 major companions in the Local Group and M81 group)?
- data: `code/companion_redshifts.json`

**Entry 79 · Are gravitational waves stretched too? 259 black-hole mergers** (page section `#gwz`)

- `code/gw_fetch.js`: The event table behind gw_redshift.py. Run in a browser tab on https://gwosc.org/eventapi/html/GWTC/ (same-origin fetch); queried 11 November 2026.
- `code/gw_redshift.py`: Are gravitational waves redshifted? A binary's chirp is a clock. If the waves are stretched on the way, as expansion predicts, the masses read off the waveform (detector frame) are (1+z) times the true ones.
- data: `code/gw_rows.txt`, `code/gw_redshift.json`

**Entry 80 · The quasar Hubble diagram: does it favour a static universe? 2,421 quasars** (page section `#quasarhd`)

- `code/quasar_hd.js`: Does the quasar Hubble diagram favour a static universe? Data: 2,421 quasars with rest-frame 2500 A and 2 keV fluxes (Lusso et al. 2020, VizieR J/A+A/642/A150).
- data: `code/quasar_hd.json`

**Entry 81 · Are distant galaxies younger? The 4000 Å break at fixed velocity dispersion** (page section `#dbreak`)

- `code/d4000_static.js`: Are distant galaxies younger? The 4000 A break at fixed velocity dispersion.
- data: `code/d4000_static.json`

**Entry 82 · Does the Lyman-α forest thicken with distance?** (page section `#lyaforest`)

- `code/lya_static.py`: Does the Lyman-alpha forest thicken with distance?
- data: `code/lya_static.json`

**Entry 83 · Radio source counts: too many faint sources for a universe that does not change** (page section `#rcounts`)

- `code/radio_counts.js`: Radio source counts and a universe that does not change. Data: NVSS 1.4 GHz catalogue (Condon et al. 1998), VizieR VIII/65/nvss on the CfA mirror, queried 11 November 2026.
- data: `code/radio_counts.json`

**Entry 84 · Is the galaxy distribution a fractal of dimension 2? 144,568 SDSS galaxies** (page section `#fractal`)

- `code/fractal_sdss.js`: Is the galaxy distribution a fractal of dimension 2? Run in a browser tab on https://skyserver.sdss.org/dr18/SearchTools/sql (same-origin SqlSearch), 11 November 2026.
- data: `code/fractal_sdss.json`

**Entry 85 · Black holes or plasmoids? The sizes of two shadows against masses weighed by stars** (page section `#bhshadow`)

- `code/shadow_scale.py`: Black holes or plasmoids? Ring sizes of Sgr A* and M87* against masses weighed by stars.
- data: `code/shadow_scale.json`

**Entry 86 · Where does the CMB come from? Its lensing by galaxies out to redshift 4, and radio from redshift 6** (page section `#cmblens`)

- `code/cmb_lens.py`: Is the CMB lensed from behind redshift 4?
- data: `code/cmb_lens.json`

**Entry 87 · Is cluster gas held by gravity or by currents? Galaxies and gas in one well** (page section `#clheat`)

- `code/sdss_sigma.js`: Galaxy velocity dispersions of X-ray clusters from SDSS DR18 spectra. Run in the browser on skyserver.sdss.org (same-origin SQL search).
- `code/gravity_beta.py`: Gravity-only expectation for beta_spec = mu m_p sigma_los^2 / kT. One NFW well holds both: galaxies as isotropic tracers (Jeans equation), gas in hydrostatic equilibrium.
- `code/cluster_heat.py`: Is cluster gas held and heated by gravity, or by fields and currents? beta_spec = mu m_p sigma^2 / kT compares the galaxies' line-of-sight velocity dispersion (galaxies feel gravity only) with the gas temperature.
- data: `code/cluster_heat.json`

**Entry 88 · Are double radio galaxies young? Peratt's sequence against the stars around them** (page section `#radiohost`)

- `code/radio_hosts.py`: Are double radio galaxies young? Peratt's sequence against the stars of their hosts.
- data: `code/radio_hosts.json`

**Entry 89 · Does gas get poorer in metals with distance? Damped Lyman-α absorbers to redshift 5** (page section `#dlametal`)

- `code/dla_metals.py`: Does the gas get poorer in metals with distance? Damped Lyman-alpha absorbers (DLAs) are gas clouds with N(HI) >= 2e20 cm^-2 seen against quasars.
- data: `code/dla_metals.json`

**Entry 90 · Does the corona redden sunlight? Plasma redshift tested at the Sun** (page section `#plasmaz`)

- `code/solar_plasma_z.py`: Does the corona redden sunlight?
- data: `code/solar_plasma_z.json`

**Entry 91 · Is matter less clumped than the CMB predicts? The S₈ measurements of 2026** (page section `#s8`)

- `code/s8_board.py`: Is matter less clumped than the CMB predicts? The S8 measurements of 2026. S8 = sigma8 (Omega_m/0.3)^0.5, all in flat LCDM, as published.
- data: `code/s8_board.json`

**Entry 92 · Is the CMB polarized like sound in a plasma? E-modes predicted without the CMB** (page section `#cmbpol`)

- `code/cmbpol.py`: Semi-analytic CMB temperature AND E-mode polarization spectra from primordial-plasma physics (Entry 92: extends the Entry 12 code with the polarization source of Zaldarriaga & Seljak 1997).
- `code/pol_compare.py`: Is the CMB polarized the way sound in the primordial plasma requires? Measured: positions and heights of the peaks and troughs of the TT, TE and EE spectra (Planck 2018 results I, Table 5).
- data: `code/pol_compare.json`

**Entry 93 · Ninety-two entries: the scoreboard again, and the measurements still to come** (page section `#score92`)

- `code/scoreboard92.py`: The scoreboard after 92 entries, counted from the Gap Map's nodes.
- data: `code/scoreboard92.json`

**Entry 94 · The second plasma era: does the CMB see the universe re-ionize?** (page section `#reion`)

- `code/cmbreion.py`: Semi-analytic CMB temperature AND E-mode polarization spectra from primordial-plasma physics (Entry 92: extends the Entry 12 code with the polarization source of Zaldarriaga & Seljak 1997).
- `code/reion_compare.py`: The second plasma era.
- data: `code/reion_compare.json`, `code/cmbreion_tau0.json`, `code/cmbreion_tau054.json`, `code/cmbreion_tau09.json`

### Data tables and results files

Tables typed in or fetched from published sources, and the results files kept. "Made by" names the script that writes a results file; "used by" the scripts that read the table or file. A results file no other script reads is quoted by its entry.

| File | Made by | Used by |
|---|---|---|
| `code/act_q_curve.txt` | typed in or fetched | cluster_density_test.py, raw_chunks.py, selection_sim.py |
| `code/act_signals.txt` | typed in or fetched | cluster_density_test.py, raw_chunks.py, sz_xray_distance.py |
| `code/age_bias_referee.json` | age_bias_referee.py | – |
| `code/arc_fields.json` | arc_fields.js | – |
| `code/band_variance.json` | band_variance.py | – |
| `code/boxsz_table3.txt` | typed in or fetched | boxsz_test.py |
| `code/cao2017_quasars.txt` | typed in or fetched | radio_core_ruler.py |
| `code/ccc_tl_test.json` | ccc_tl_test.py | – |
| `code/ccctl_bbn.json` | ccctl_bbn.py | – |
| `code/cloud_currents.json` | cloud_currents.py | – |
| `code/clump_track.json` | kept from its entry's run | clumping_fit.py |
| `code/cluster_heat.json` | cluster_heat.py | – |
| `code/clusters_y_T.txt` | typed in or fetched | cluster_density_test.py, cluster_density_test_v2.py, selection_sim.py, sz_shadows.js |
| `code/clusters_y_T_v2.txt` | typed in or fetched | cluster_density_test_v2.py |
| `code/cmb_lens.json` | cmb_lens.py | – |
| `code/cmbpol_noncmb.json` | cmbpol.py | pol_compare.py |
| `code/cmbreion_tau0.json` | kept from its entry's run | – |
| `code/cmbreion_tau054.json` | kept from its entry's run | – |
| `code/cmbreion_tau09.json` | kept from its entry's run | – |
| `code/companion_redshifts.json` | companion_redshifts.js | – |
| `code/corr/found_a.txt` | typed in or fetched | corrections_check.py, corrections_de.py |
| `code/corr/found_b.txt` | typed in or fetched | corrections_check.py, corrections_de.py |
| `code/corr/found_c.txt` | typed in or fetched | corrections_check.py, corrections_de.py |
| `code/corr/found_d.txt` | typed in or fetched | corrections_check.py, corrections_de.py |
| `code/corr/ps1_a.txt` | typed in or fetched | corrections_check.py, corrections_de.py |
| `code/corr/ps1_b.txt` | typed in or fetched | corrections_check.py, corrections_de.py |
| `code/corr/ps1_c.txt` | typed in or fetched | corrections_check.py, corrections_de.py |
| `code/corr/ps1_d.txt` | typed in or fetched | corrections_check.py, corrections_de.py |
| `code/corr/ps1_e.txt` | typed in or fetched | corrections_check.py, corrections_de.py |
| `code/cosmos_box.txt` | typed in or fetched | age_split.py, legac_test.py |
| `code/d4000_static.json` | d4000_static.js | – |
| `code/data/bulbul_spt.txt` | typed in or fetched | selection_sim_xmm.py, spt_xmm_distance.py |
| `code/data/sanders_spt.txt` | typed in or fetched | distance_test_repeat.py, ratio_selection_sim.py, raw_chunks3.py, spt_xmm_distance.py |
| `code/decia18.txt` | typed in or fetched | dla_metals.py |
| `code/density_free.json` | density_free.py | – |
| `code/des_dovekie_results.json` | des_dovekie_fits.js | – |
| `code/dla_metals.json` | dla_metals.py | – |
| `code/electric_comet.json` | electric_comet.py | – |
| `code/electric_sun.json` | electric_sun.py | – |
| `code/entry28_points.json` | kept from its entry's run | boxsz_test.py, forecasts.py, spt_xmm_distance.py, static_distance_law.py |
| `code/entry32_points.json` | kept from its entry's run | boxsz_test.py, static_distance_law.py |
| `code/far2018.txt` | typed in or fetched | anchor_shift.py |
| `code/farplus.txt` | typed in or fetched | anchor_shift.py |
| `code/filament_spice.json` | filament_spice.js | combine_filaments.py |
| `code/first_sample_act.txt` | typed in or fetched | cluster_density_test.py |
| `code/fog_results.json` | radio_fog.js | fog_model.py |
| `code/fractal_sdss.json` | fractal_sdss.js | – |
| `code/fricat.txt` | typed in or fetched | radio_hosts.py |
| `code/friicat.txt` | typed in or fetched | radio_hosts.py |
| `code/gw_redshift.json` | gw_redshift.py | – |
| `code/gw_rows.txt` | typed in or fetched | gw_fetch.js, gw_redshift.py |
| `code/karlsson_dr16q.json` | karlsson_dr16q.js | – |
| `code/legac_box.txt` | typed in or fetched | age_split.py, legac_test.py |
| `code/lightcurve_foundation.txt` | typed in or fetched | – |
| `code/lightcurve_panstarrs.txt` | typed in or fetched | – |
| `code/line_width_jwst.json` | line_width_jwst.py | – |
| `code/local_relations_z5.json` | local_relations_z5.py | – |
| `code/lya_static.json` | lya_static.py | – |
| `code/mantz2014_fgas.txt` | typed in or fetched | gas_fraction_test.py |
| `code/mantz2022_fgas.txt` | typed in or fetched | forecasts.py, gas_fraction_test.py |
| `code/maughan2008_matched.txt` | typed in or fetched | cluster_density_test_v2.py |
| `code/migkas_cand.txt` | typed in or fetched | cluster_heat.py, sdss_sigma.js |
| `code/near2018.txt` | typed in or fetched | anchor_shift.py |
| `code/nearplus.txt` | typed in or fetched | anchor_shift.py |
| `code/pantheon_plus_foundation.txt` | typed in or fetched | pantheon_plus_ps1.py |
| `code/pantheon_plus_nearby_coords.txt` | typed in or fetched | sn_frame.py |
| `code/pantheon_plus_ps1md.txt` | typed in or fetched | pantheon_plus_ps1.py |
| `code/pantheon_plus_rows.txt` | typed in or fetched | age_bias_referee.py, anchor_shift.py, ccc_tl_test.py, corrections_check.py, corrections_de.py, pantheon_plus_expected.py, pantheon_plus_ps1.py, pantheon_plus_scan.py, pantheon_plus_surveys.py, pantheon_plus_test.py, static_law_plus.py |
| `code/planck_lovisari.txt` | typed in or fetched | distance_test_repeat.py, raw_chunks3.py |
| `code/pol_compare.json` | pol_compare.py | – |
| `code/pole_currents.json` | pole_currents.js | – |
| `code/possum_results.json` | possum_through.js | – |
| `code/quasar_hd.json` | quasar_hd.js | – |
| `code/quiescent_z3.json` | quiescent_z3.py | – |
| `code/radio_background.json` | radio_background.py | – |
| `code/radio_counts.json` | radio_counts.js | – |
| `code/radio_hosts.json` | radio_hosts.py | – |
| `code/rafelski12.txt` | typed in or fetched | dla_metals.py |
| `code/reion_compare.json` | reion_compare.py | – |
| `code/ruler_by_survey.json` | ruler_by_survey.py | – |
| `code/s8_board.json` | s8_board.py | – |
| `code/sanders_act_profiles.txt` | typed in or fetched | raw_chunks.py, sz_xray_distance.py |
| `code/sanders_spt.txt` | typed in or fetched | distance_test_repeat.py, ratio_selection_sim.py, raw_chunks3.py, spt_xmm_distance.py |
| `code/scoreboard92.json` | scoreboard92.py | – |
| `code/sdss_hist_conc.txt` | typed in or fetched | radio_hosts.py |
| `code/sdss_hist_zoo.txt` | typed in or fetched | radio_hosts.py |
| `code/sdss_sigma.txt` | typed in or fetched | cluster_heat.py |
| `code/shadow_scale.json` | shadow_scale.py | – |
| `code/signed_stack.json` | signed_stack.js | – |
| `code/solar_plasma_z.json` | solar_plasma_z.py | – |
| `code/static_distance_law.json` | static_distance_law.py | gas_fraction_test.py |
| `code/survey_labels.txt` | typed in or fetched | corrections_de.py, pantheon_plus_expected.py, pantheon_plus_surveys.py |
| `code/temps_matched.txt` | typed in or fetched | cluster_density_test.py, raw_chunks.py, sz_xray_distance.py |
| `code/through_results.json` | forecast_through.py | – |
| `code/tifft_alfalfa.json` | tifft_alfalfa.js | – |
| `code/vds2013_table4.txt` | typed in or fetched | ccc_tl_test.py, dynamical_mass_test.py, local_relations_z5.py, trend_checks.py |
| `code/wu98_raw.txt` | typed in or fetched | cluster_heat.py |


### Added since Round 89 (Entries 95–163)

168 files were added after the first version of this bundle; the earlier sections above list only what existed at Entry 94, so for newer entries see the Workbench section for that entry (each names its code files in its Limits paragraph).

Test of the new Python scripts, run from `code/` with a 120 s limit each: 60 of 65 ran cleanly. Three exceed 120 s (`biref.py`, `frb_shape.py`, `wb36.py`) and were not run to the end. Two need results files made earlier in the same entry's chain and stop at once when run alone (`polvoid.py`, `prereg.py`). Browser scripts (`.js`) were run against live catalogues and need those services to run again.

New files: age_edge.json, age_edge.py, apogee_check.js, astero.js, astero.json, astero.py, astero_raw.json, biref.json, biref.py, bubbles.json, bubbles.py, civ_bands.json, civ_bands.py, cloud_pin.js, cloud_pin.json, cloud_pin.py, cloud_step.js, cloud_step.json, cloud_step.py, clump_limits.json, clump_limits.py, cmbreion_flash09.json, cosmic_heat.json, cosmic_heat.py, cr_circuits.json, cr_circuits.py, decider.json, decider.py, drift.json, drift.py, dsph_plasma.json, dsph_plasma.py, eb_feh.txt, ebz.json, ebz.py, efe.json, efe.py, eker_pairs.txt, estars.py, estars_results.json, flare_circuits.js, flare_circuits_results.json, flash_compare.json, flash_compare.py, flows.json, flows.py, frb93.txt, frb_shape.json, frb_shape.py, gaia_vesc.js, gas_fraction_test.json, gc_filaments.json, gc_filaments.py, grb_circuits.js, grb_circuits_results.json, grz.json, grz.py, grz_stack.js, grz_tempel.json, grz_tempel.py, grz_tempel_raw.json, grz_why.js, grz_why.json, grz_why.py, grz_why_raw.json, halo_loop.json, halo_loop.py, helium.json, helium.py, hot_voids.json, hot_voids.py, hubble_routes.json, hubble_routes.py, jet_currents.json, jet_currents.py, lens_static.json, lens_static.py, li3gc.json, li3gc.py, li6397.json, li6397.py, li_ism.json, li_ism.py, lrd_5000K.json, lrd_5000K.py, lrd_lags.json, lrd_lags.py, lrd_scatter.json, lrd_scatter.py, m30.csv, massflux.json, massflux.py, mw_mass.json, mw_mass.py, n6397.csv, n6752.csv, neutron_bbn.json, neutron_bbn.py, nvss_dipole.js, nvss_dipole.json, nvss_dipole.py, oldest.json, oldest.py, one_field.json, one_field.py, ou24.txt, outflow_spin.json, outflow_spin.py, plens.json, plens.py, pmf_21cm.json, pmf_21cm.py, pmf_halos.json, pmf_halos.py, pn_pinch.json, pn_pinch.py, polvoid.json, polvoid.py, poster_models.json, poster_models.py, prereg.json, prereg.py, progenitors.json, progenitors.py, pta_plasma.json, pta_plasma.py, pulsars.json, pulsars.py, qso_align.js, qso_align_results.json, quasar_dilation.js, quasar_dilation_results.json, s4tm.txt, scaling.json, scaling.py, scoreboard126.json, scoreboard126.py, scoreboard142.json, scoreboard142.py, scoreboard158.json, scoreboard158.py, seis.csv, seisz.json, seisz.py, slacs.txt, spec100.json, static_hubble.json, subgiant_hist.json, tau21.json, tau21.py, tau_from_xhi.json, tau_from_xhi.py, timescape.json, timescape.py, tout.py, vesc_apogee.json, vesc_apogee.py, vesc_results.json, void_fields.json, void_fields.py, void_finder.js, voids_slab1.txt, voids_slab2.txt, wb36.json, wb36.py, wb_2026.json, xrism_balance.json, xrism_balance.py
