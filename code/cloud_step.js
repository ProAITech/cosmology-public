// Entry 161, run in the browser on the VizieR CfA page.
// Data: J/ApJ/702/1230/catalog (Taylor, Stil & Sunstrum 2009), all 37,543 NVSS rotation measures (query needs a constraint, RM=-10000..10000);
// sexagesimal positions converted to Galactic with the ICRS-to-Galactic matrix of Entry 151.
// Spines (l, b): Orion A (208.6, -19.1) to (214.5, -19.6); Perseus (158.35, -20.5) [NGC 1333] to (160.57, -16.8) [Barnard 5].
// Local flat frame at the spine midpoint: x = dl cos b0, y = db; s = position along the spine, d = signed distance (toward higher b positive).
// Sources kept: |s| <= half-length + 1 deg, 0.05 < |d| <= 10 deg, |RM - median| < 5 MAD.
// Weighted least squares, weights 1/(e_RM^2 + sig^2), sig iterated from the residual variance (6 passes):
//   RM = c0 + c1 d + c2 s + A_in sgn(d) [|d| < 1.5] + A_out sgn(d) [|d| >= 1.5]
// Near/far variant: window 5 deg (A_in = 1.5-5 deg is then 0-5 deg; A_out = 5-10 deg).
// Null: the same spine shifted by dl = +-10..60 deg (step 5) and db = -4, 0, +4 deg, refitted where N > 60; z = A_out / its formal error.
