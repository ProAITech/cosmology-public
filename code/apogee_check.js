// Entry 160, run in the browser on the ESA Gaia archive page (VizieR CfA allows cross-origin reads).
// 1. APOGEE DR17 (VizieR III/286/catalog), 36 slices of 10 degrees in RA, plx > 0.333 mas, columns
//    RAJ2000, DEJ2000, GaiaEDR3, Tel, HRV, s_HRV, e_HRV, Nvis, SNR, plx, e_plx, pmRA, e_pmRA, pmDE, e_pmDE, rgeo, RV, e_RV, [M/H], Teff, logg.
//    496,454 rows; duplicates removed by Gaia source number, keeping the highest SNR: 433,064 stars.
// 2. Entry 151's asynchronous ADQL query of gaiadr3.gaia_source (parallax_over_error > 20, ruwe < 1.4, radial velocity,
//    parallax > 0.333, Galactocentric speed > 300 km/s, Sun at (12.9, 245.6, 7.78) km/s), adding grvs_mag and rv_method_used: 36,308 stars.
// 3. Galactocentric speed with the ICRS-to-Galactic matrix of Entry 151 and distance 1/parallax, once with Gaia's and once with APOGEE's velocity.
// APOGEE velocity accepted if e_HRV < 0.5 km/s, SNR >= 20, Nvis >= 2 and s_HRV < 1 km/s (sample B); sample A also allows single visits.
const ok = a => isFinite(a[4]) && isFinite(a[6]) && a[6] < 0.5 && a[8] >= 20 && a[7] >= 2 && a[5] < 1;
// for each Gaia star with an APOGEE match passing ok(): bright = grvs_mag < 12; count |RV_gaia - HRV| < 10 km/s; count speeds > 350 with each velocity.
// Results: bright 124 checked, 123 agree, 26 above 350 with either velocity; faint 100 checked, 88 agree, 8 above 350 with either.
