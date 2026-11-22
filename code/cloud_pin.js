// Entry 162, run in the browser on the VizieR CfA page, with the data, spines and source selection of Entry 161 (cloud_step.js).
// Scatter sig fixed from the fit without a step (RM = c0 + c1 d + c2 s), then for each step position d0 = -6..6 deg (0.25 steps):
//   RM = c0 + c1 d + c2 s + A sgn(d - d0);  profile = chi2(no step) - chi2(step at d0).
// Null: the same scan for the spines shifted by dl = +-10..60 deg and db = -4, 0, +4 deg (N > 60): best d0 and its delta-chi2.
// Extent: steps confined to |d| < w for w = 1, 2, 3, 5, 10 deg; two-zone fit with separate amplitudes inside 2 deg and 2-10 deg.
// Orientation: spine rotated by -80..90 deg (10 deg steps) about its midpoint; best delta-chi2 for |d0| <= 2 deg.
