// Entry 151: the local escape speed from Gaia DR3, run in the browser on the ESA Gaia archive page (TAP, ADQL, synchronous).
// Selection: parallax_over_error > 20, ruwe < 1.4, a radial velocity, parallax > 0.333 mas (within 3 kpc),
// and Galactocentric speed above 300 km/s, computed inside the query with the Sun's motion (12.9, 245.6, 7.78) km/s
// (Drimmel & Poggio 2018) and the ICRS-to-Galactic rotation matrix.
const T=[[-0.0548755604162154,-0.8734370902348850,-0.4838350155487132],[0.4941094278755837,-0.4448296299600112,0.7469822444972189],[-0.8676661490190047,-0.1980763734312015,0.4559837761750669]];
const SUN=[12.9,245.6,7.78];
function adqlVgc2(){const a="RADIANS(ra)",d="RADIANS(dec)",vr="radial_velocity",va="(4.74047*pmra/parallax)",vd="(4.74047*pmdec/parallax)";
 const ex=`(${vr}*COS(${d})*COS(${a}) - ${va}*SIN(${a}) - ${vd}*SIN(${d})*COS(${a}))`,ey=`(${vr}*COS(${d})*SIN(${a}) + ${va}*COS(${a}) - ${vd}*SIN(${d})*SIN(${a}))`,ez=`(${vr}*SIN(${d}) + ${vd}*COS(${d}))`;
 const g=[0,1,2].map(i=>`(${T[i][0]}*${ex} + ${T[i][1]}*${ey} + ${T[i][2]}*${ez} + ${SUN[i]})`); return `(POWER(${g[0]},2)+POWER(${g[1]},2)+POWER(${g[2]},2))`;}
const q=`SELECT source_id, ra, dec, l, b, parallax, parallax_error, pmra, pmra_error, pmdec, pmdec_error, radial_velocity, radial_velocity_error, phot_g_mean_mag FROM gaiadr3.gaia_source WHERE parallax_over_error > 20 AND ruwe < 1.4 AND radial_velocity IS NOT NULL AND parallax > 0.333 AND ${adqlVgc2()} > 90000`;
const csv=await (await fetch('/tap-server/tap/sync',{method:'POST',body:new URLSearchParams({REQUEST:'doQuery',LANG:'ADQL',FORMAT:'csv',QUERY:q})})).text();
// per star: Galactocentric speed, its scatter from 60 Monte Carlo draws of parallax, proper motions and radial velocity (correlations ignored)
// Leonard & Tremaine (1990) tail, P(v) = (k+1)(v_esc - v)^k / (v_esc - v_min)^(k+1) for v_min < v < v_esc,
// plus a 1% outlier floor uniform to 900 km/s; likelihood on a grid of v_esc (400-800, step 2) and k (0.5-7, step 0.1);
// stars kept if their speed error is below a cut (15 or 25 km/s). Marginal v_esc for flat k priors on chosen ranges.
