"""Gravity-only expectation for beta_spec = mu m_p sigma_los^2 / kT.
One NFW well holds both: galaxies as isotropic tracers (Jeans equation), gas in
hydrostatic equilibrium. T is averaged over 0.2-0.5 r500 (as in Migkas et al. 2020);
sigma_los over galaxies projected inside an aperture R_ap (in units of r200).
No normalisation is free: both scale with G M200 / r200.
"""
import numpy as np
def m(y): return np.log(1+y) - y/(1+y)
def model(c=4.0, cg=4.0, slope=2.2, Rap=1.0, rmax=3.0, n=1500):
    r = np.logspace(-3, np.log10(rmax), n)            # in r200
    M = m(c*r)/m(c)                                    # M(<r)/M200
    g = M/r**2                                         # G=M200=r200=1 units, V200=1
    nu = 1/((cg*r)*(1+cg*r)**2)                        # tracer density (NFW shape, conc. cg)
    # sigma_r^2 nu = int_r^rmax nu g dr'
    integrand = nu*g
    cum = np.concatenate([[0], np.cumsum(0.5*(integrand[1:]+integrand[:-1])*np.diff(r))])
    snu = cum[-1]-cum                                  # int from r to rmax
    sr2 = snu/nu
    # projection: Sigma(R) = 2 int nu dz ; Sigma sigma_los^2 = 2 int nu sr2 dz (isotropic)
    Rg = np.linspace(0.005, Rap, 200)
    num = 0; den = 0
    for R in Rg:
        rr = r[r > R]; z = np.sqrt(rr**2-R**2)
        nn = np.interp(rr, r, nu); ss = np.interp(rr, r, sr2)
        dz = np.gradient(z)
        S = np.sum(nn*dz); SS = np.sum(nn*ss*dz)
        w = 2*np.pi*R
        num += w*SS; den += w*S
    sig2 = num/den
    # r500 in units of r200: mean density 500 vs 200 times critical
    rs = np.linspace(0.3, 1.0, 20000); f = m(c*rs)/m(c)/rs**3
    r500 = rs[np.argmin(abs(f-2.5))]
    rr = np.linspace(0.2*r500, 0.5*r500, 50)
    kT = np.mean((m(c*rr)/m(c)/rr)/slope)              # kT/(mu m_p) = G M(r)/r / slope
    return sig2/kT, r500
if __name__ == '__main__':
    for c in (3, 4, 6):
        for cg in (c/2, c):
            for s in (2.0, 2.2, 2.5):
                print(c, cg, s, [round(model(c, cg, s, Rap)[0], 3) for Rap in (0.4, 0.6, 1.0)], round(model(c, cg, s)[1], 3))
