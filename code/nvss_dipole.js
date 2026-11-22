// Entry 163, run in the browser on the VizieR CfA page.
// NVSS (VIII/65/nvss), all sources with S1.4 >= 15 mJy, fetched in nine flux slices: 422,541 sources (positions and fluxes only).
// Sky cut into 100 bands of equal width in sin(Dec) x 100 in RA: 10,000 equal-area cells of 4.1 square degrees.
// Cells kept if their centre lies more than 1.5 deg beyond the masks (|b| > bcut, Dec > dmin).
// Least squares: count_i = N0 (1 + d . n_i) (optionally + 5 quadrupole terms); errors from var(count) = N0, then scaled by
// sqrt(1.297), the measured variance-to-mean ratio of the cells after removing the dipole (S >= 20 mJy, |b| > 15, Dec > -30).
// Count slope x = ln[N(>S)/N(>1.5S)] / ln 1.5 over |b| > 15 deg. CMB dipole direction RA 167.94, Dec -6.94.
