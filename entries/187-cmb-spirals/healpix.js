/* Entry 187: HEALPix ang2pix in RING order (Gorski et al. 2005), ported from the HEALPix C++ loc2pix. Validated against healpy. */
function ang2pixRing(nside, z, phi) {          // z = cos(theta), phi in radians
  const za = Math.abs(z), twopi = 2 * Math.PI;
  let tt = phi % twopi; if (tt < 0) tt += twopi; tt /= (0.5 * Math.PI);
  const ncap = 2 * nside * (nside - 1), npix = 12 * nside * nside;
  if (za <= 2 / 3) {
    const t1 = nside * (0.5 + tt), t2 = nside * z * 0.75;
    const jp = Math.floor(t1 - t2), jm = Math.floor(t1 + t2);
    const ir = nside + 1 + jp - jm, kshift = 1 - (ir & 1);
    let ip = Math.floor((jp + jm - nside + kshift + 1) / 2); ip = ((ip % (4 * nside)) + 4 * nside) % (4 * nside);
    return ncap + (ir - 1) * 4 * nside + ip;
  }
  const tp = tt - Math.floor(tt), tmp = nside * Math.sqrt(3 * (1 - za));
  const jp = Math.floor(tp * tmp), jm = Math.floor((1 - tp) * tmp), ir = jp + jm + 1;
  let ip = Math.floor(tt * ir); ip = ((ip % (4 * ir)) + 4 * ir) % (4 * ir);
  return z > 0 ? 2 * ir * (ir - 1) + ip : npix - 2 * ir * (ir + 1) + ip;
}
if (typeof module !== 'undefined') module.exports = { ang2pixRing };
