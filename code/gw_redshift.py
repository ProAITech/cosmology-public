"""Workbench Entry 79: are gravitational waves redshifted?

A binary's chirp is a clock. If the waves are stretched on the way, as expansion predicts, the masses read off the
waveform (detector frame) are (1+z) times the true ones. If the cosmological redshift were energy lost by light to
plasma or anything else light meets, the waves, which meet nothing, would arrive unstretched and the detector-frame
masses would be the true ones.  Parameterise  M_true = M_det / (1+z(D))**alpha:  alpha = 1 expansion, alpha = 0 static.

Data: 259 binary black holes with FAR < 1/yr from the cumulative GWTC on GWOSC (gw_rows.txt).
Method: likelihood of each detector-frame chirp mass GIVEN its gravitational-wave distance D, so the number of events
at each distance (rates, volumes, evolution of the rate) never enters:
    p(M_i | D_i) = f(M_i/(1+z_i)^a) (1+z_i)^-a  Pdet(M_i, D_i)  /  integral of the same over M
f is one mass distribution shared by all distances (log-density on 21 knots, fitted); Pdet is the chance that a binary
of detector-frame chirp mass M at catalogue distance D is detected.  alpha is profiled over f.
Selection, as finally used (run_emp): the catalogue distance of an event is set by its SNR and loudness,
D = kappa x loudness(M) / SNR, where kappa carries the orientation the distance estimate assumed. kappa is measured for
every event; its spread is narrow (0.23 in ln). So in (mass, catalogue distance) the chance of detection is
    Pdet(M, D) = mean over the events' kappa of  Phi((kappa x loudness(M)/D - rho50)/width)
with the soft threshold rho50 = 9.1, width 0.7 fitted to the observed SNRs (a rho^-4 parent; tail slope measured -4.2).
Nothing about orientation or detector reach is assumed. A first version (run_case, kept for comparison) used the
orientation-averaged chance of detection at the catalogue distance, from a phenomenological IMR amplitude, an analytic
Advanced LIGO noise curve and isotropic orientations, with each run's reach fitted to its events' SNRs; it is softer
than the selection the catalogue distances actually carry, and its answer moved with the assumed reach.
Loudness against mass: the design-curve waveform corrected to the events ('IMR+data', a quadratic in ln M fitted so
that the median kappa is the same at every mass), the uncorrected waveform ('IMR'), or the inspiral law M^(5/6).
Run:  python3 -I gw_redshift.py gw_rows.txt [n_mocks]  > gw_redshift.json
"""
import sys, json, math
import numpy as np
from scipy.optimize import minimize
from scipy.integrate import quad

MSUN_S = 4.925491e-6           # G Msun / c^3 in seconds
rows = [list(map(float, r.split(','))) for l in open(sys.argv[1]) if not l.startswith('#') for r in l.strip().split(';') if r]
A = np.array(rows)
RUN, MDET, D, DLO, DHI, SNR, Q, ZLVK = A.T
RUN = RUN.astype(int)
N = len(A)

# ---------- cosmology used only to turn D into z for the expanding case (Planck 2015, as LVK use) ----------
H0, OM = 67.74, 0.3075
c = 299792.458
zg = np.linspace(0, 3, 30001)
Ez = np.sqrt(OM*(1+zg)**3 + 1-OM)
dc = np.concatenate([[0], np.cumsum(0.5*(1/Ez[1:]+1/Ez[:-1])*np.diff(zg))])*c/H0
dl = (1+zg)*dc
z_of_D = lambda d: np.interp(d, dl, zg)

# ---------- orientation factor Theta in [0,1] for one interferometer, isotropic ----------
rng = np.random.default_rng(1)
n = 400000
ct, ph, ps, ci = rng.uniform(-1, 1, n), rng.uniform(0, 2*np.pi, n), rng.uniform(0, np.pi, n), rng.uniform(-1, 1, n)
Fp = 0.5*(1+ct**2)*np.cos(2*ph)*np.cos(2*ps) - ct*np.sin(2*ph)*np.sin(2*ps)
Fx = 0.5*(1+ct**2)*np.cos(2*ph)*np.sin(2*ps) + ct*np.sin(2*ph)*np.cos(2*ps)
TH = np.sort(np.sqrt(Fp**2*(1+ci**2)**2/4 + Fx**2*ci**2))
def surv(x):                      # P(Theta > x)
    return 1 - np.searchsorted(TH, x)/n
thg = np.linspace(0, 1, 201)
thpdf = np.histogram(TH, bins=thg, density=True)[0]
def pdf(x):
    i = np.clip((x*200).astype(int), 0, 199)
    return np.where((x >= 0) & (x <= 1), thpdf[i], 0.0)

# ---------- optimal SNR shape against detector-frame chirp mass ----------
def Sh(f, f0=215.0):              # Advanced LIGO analytic fit (Sathyaprakash & Schutz 2009)
    x = f/f0
    return 1e-49*(x**-4.14 - 5*x**-2 + 111*(1 - x**2 + 0.5*x**4)/(1 + 0.5*x**2))
COEF = {'merg': (2.9740e-1, 4.4810e-2, 9.5560e-2), 'ring': (5.9411e-1, 8.9794e-2, 1.9111e-1),
        'sig': (5.0801e-1, 7.7515e-2, 2.2369e-2), 'cut': (8.4845e-1, 1.2848e-1, 2.7299e-1)}
def rho1(mc, q=0.7, flow=20.0, imr=True):
    if not imr:
        return mc**(5/6)
    eta = q/(1+q)**2
    mt = mc*eta**-0.6*MSUN_S
    fr = {k: (a*eta**2+b*eta+cc)/(math.pi*mt) for k, (a, b, cc) in COEF.items()}
    fm, fR, sg, fc = fr['merg'], fr['ring'], fr['sig'], fr['cut']
    w = math.pi*sg/2*(fR/fm)**(-2/3)
    def amp2(f):
        if f < fm: a = f**(-7/6)
        elif f < fR: a = fm**(-7/6)*(f/fm)**(-2/3)
        else: a = fm**(-7/6)*w*(sg/(2*math.pi))/((f-fR)**2 + sg**2/4)
        return a*a/Sh(f)
    if flow >= fc: return 0.0
    edges = [flow] + [e for e in (fm, fR) if flow < e < fc] + [fc]
    s = sum(quad(amp2, edges[k], edges[k+1], limit=200)[0] for k in range(len(edges)-1))
    return mc**(5/6)*math.sqrt(s)

def make_rho(**kw):
    g = np.exp(np.linspace(np.log(1.5), np.log(600), 400))
    v = np.array([rho1(m, **kw) for m in g])
    v = v/np.interp(30.0, g, v)                    # normalise to 1 at a 30 Msun chirp mass
    return lambda m: np.interp(np.log(m), np.log(g), v)

# ---------- the fit ----------
LMIN, LMAX = np.log(2.0), np.log(400.0)
KNOTS = np.linspace(np.log(3.0), np.log(200.0), 21)   # spacing 0.21 in ln(mass); 8 and 11 knots were too coarse to resolve the 8-10 Msun peak
MG = np.exp(np.linspace(LMIN, LMAX, 700))
DLM = np.diff(np.log(MG)).mean()

def logf(lm, th):                  # log-density in ln(source chirp mass), linear between knots, steep fall outside
    v = np.interp(lm, KNOTS, th)
    v = v - 8.0*np.clip(KNOTS[0]-lm, 0, None)/0.1 - 8.0*np.clip(lm-KNOTS[-1], 0, None)/0.1
    return v

def setup(rho, rth, Dv, run, K):
    """Pdet on the mass grid for every event, and at the event's own mass."""
    rg = rho(MG)
    PG = surv(rth*Dv[:, None]/(K[run][:, None]*rg[None, :]))
    PE = surv(rth*Dv/(K[run]*rho(MDET)))
    return PG, PE

def _hat(x):
    """indices and weights of the two knots bracketing x (clamped at the ends), for value and gradient"""
    k = np.clip(np.searchsorted(KNOTS, x) - 1, 0, len(KNOTS) - 2)
    t = np.clip((x - KNOTS[k])/(KNOTS[k+1] - KNOTS[k]), 0.0, 1.0)
    return k, t

def _out(x):                       # the fixed fall-off outside the knots
    return 8.0*np.clip(KNOTS[0]-x, 0, None)/0.1 + 8.0*np.clip(x-KNOTS[-1], 0, None)/0.1

D2 = np.diff(np.eye(len(KNOTS)), 2, axis=0)
PEN = D2.T @ D2/0.25**2           # smoothness on the log-density (second differences)

def nll(th, a, zz, PG, PE, sel):
    """negative log-likelihood and its gradient with respect to the knot values"""
    nk = len(KNOTS)
    sh = a*np.log1p(zz[sel])
    u = np.log(MDET[sel]) - sh
    ku, tu = _hat(u)
    lnum = (1-tu)*th[ku] + tu*th[ku+1] - _out(u) + np.log(np.maximum(PE[sel], 1e-300))
    v = np.log(MG)[None, :] - sh[:, None]
    kv, tv = _hat(v)
    lf = (1-tv)*th[kv] + tv*th[kv+1] - _out(v)
    dens = np.exp(lf)*PG[sel]
    S = dens.sum(1)
    lden = np.log(S*DLM + 1e-300)
    q = dens/np.maximum(S, 1e-300)[:, None]
    g = -(np.bincount(ku, 1-tu, nk) + np.bincount(ku+1, tu, nk))
    g += np.bincount(kv.ravel(), (q*(1-tv)).ravel(), nk) + np.bincount((kv+1).ravel(), (q*tv).ravel(), nk)
    f = -np.sum(lnum - lden) + 0.5*th @ PEN @ th
    return f, g + PEN @ th

def profile(alphas, zz, PG, PE, sel):
    out, th = [], np.zeros(len(KNOTS))
    for a in alphas:
        r = minimize(nll, th, args=(a, zz, PG, PE, sel), method='L-BFGS-B', jac=True)
        th = r.x
        out.append(r.fun)
    out = np.array(out)
    return out - out.min()

def calibrate(rho, rth, Dv):
    """Sensitivity K of each run: maximise the likelihood of the observed SNRs given each event's mass and distance,
       SNR = K rho(M) Theta / D, with Theta truncated by the threshold."""
    K = {}
    for r in sorted(set(RUN)):
        s = RUN == r
        kap = SNR[s]*Dv[s]/rho(MDET[s])
        def nl(lk):
            k = math.exp(lk)
            p = pdf(kap/k)/k
            P = surv(rth*Dv[s]/(k*rho(MDET[s])))
            if np.any(p <= 0) or np.any(P <= 0): return 1e30
            return -np.sum(np.log(p) - np.log(P))
        grid = np.linspace(np.log(kap.max()*1.0001), np.log(kap.max()*8), 400)
        vals = [nl(g) for g in grid]
        K[r] = math.exp(grid[int(np.argmin(vals))])
    return K

def emp_rho(rho, deg=2):
    """loudness corrected to the events: a polynomial in ln(M/30) fitted to ln(SNR x D / loudness) after removing
       each run's median; the orientation factor implied by the distance estimates has the same median at every mass"""
    y = np.log(SNR*D/rho(MDET))
    for r in set(RUN): y[RUN == r] -= np.median(y[RUN == r])
    c = np.polyfit(np.log(MDET/30.0), y, deg); c[-1] = 0.0
    return (lambda m: rho(m)*np.exp(np.polyval(c, np.log(m/30.0)))), [round(float(x), 4) for x in c]

def loudness(name, **kw):
    if name == 'IMR': return make_rho(**kw), None
    if name == 'IMR+data': return emp_rho(make_rho(**kw))
    if name == 'inspiral': return make_rho(imr=False), None

AL = np.round(np.arange(-0.5, 2.51, 0.1), 2)
def summarise(prof, alphas=AL):
    i = int(np.argmin(prof)); fine = np.linspace(alphas[0], alphas[-1], 3001); pf = np.interp(fine, alphas, prof)
    inside = fine[pf <= 0.5]
    d0 = float(np.interp(0.0, alphas, prof)); d1 = float(np.interp(1.0, alphas, prof))
    return {'alpha_best': float(alphas[i]), 'alpha_68': [round(float(inside.min()), 2), round(float(inside.max()), 2)],
            'dlnL_alpha0': round(d0, 2), 'dlnL_alpha1': round(d1, 2), 'sigma_alpha0': round(math.sqrt(2*max(d0, 0)), 2),
            'sigma_alpha1': round(math.sqrt(2*max(d1, 0)), 2)}

def run_case(name, loud='IMR+data', scale=1.0, rth=9.0, rho_kw={}, Dscale=1.0, sel=None):
    rho, coef = loudness(loud, **rho_kw)
    Dv = D*Dscale
    K = calibrate(rho, rth, Dv)
    Karr = np.zeros(7)
    for k, v in K.items(): Karr[k] = v*scale
    PG, PE = setup(rho, rth, Dv, RUN, Karr); PE = np.maximum(PE, 1e-6)
    zz = z_of_D(Dv)
    sel = np.ones(N, bool) if sel is None else sel
    prof = profile(AL, zz, PG, PE, sel)
    res = {'case': name, 'loudness': loud, 'scale': scale, 'n': int(sel.sum())}
    res.update(summarise(prof))
    res['K_by_run_Mpc'] = {str(k): round(v*scale, 0) for k, v in K.items()}
    if coef is not None: res['loudness_correction_coef'] = coef
    return res, (prof, rho, Karr, zz, PG, PE)

from scipy.special import ndtr
RHO50, WIDTH = 9.1, 0.7

def emp_selection(rho, r50=RHO50, w=WIDTH):
    kap = SNR*D/rho(MDET)
    med = {r: float(np.median(kap[RUN == r])) for r in set(RUN)}
    m12 = float(np.median(kap[RUN <= 2])); med[1] = m12; med[2] = m12      # O1 and O2 hold 3 and 7 events
    u = kap/np.array([med[r] for r in RUN])                                 # orientation factor relative to the run
    def pdet(Dv, M, runs):
        # Dv, M, runs broadcast together; returns the chance of SNR above the soft threshold
        need = (np.vectorize(med.get)(runs)*rho(M)/Dv)
        sh = need.shape; x = need.reshape(-1)
        out = np.empty(x.size)
        for a in range(0, x.size, 20000):
            out[a:a+20000] = ndtr((x[a:a+20000, None]*u[None, :] - r50)/w).mean(1)
        return out.reshape(sh)
    return pdet, med, u

def run_emp(name, loud='IMR+data', r50=RHO50, w=WIDTH, Dscale=1.0, sel=None):
    rho, coef = loudness(loud)
    pdet, med, u = emp_selection(rho, r50, w)
    Dv = D*Dscale
    # a common scale on all distances cancels in kappa x loudness / D, so selection uses the catalogue values; only z moves
    PG = pdet(D[:, None], MG[None, :], RUN[:, None])
    PE = np.maximum(pdet(D, MDET, RUN), 1e-6)
    zz = z_of_D(Dv)
    sel = np.ones(N, bool) if sel is None else sel
    prof = profile(AL, zz, PG, PE, sel)
    res = {'case': name, 'loudness': loud, 'rho50': r50, 'width': w, 'n': int(sel.sum())}
    res.update(summarise(prof))
    if coef is not None: res['loudness_correction_coef'] = coef
    return res, (prof, rho, pdet, zz, PG, PE)

if __name__ == '__main__':
    out = {'n_events': N, 'runs': {str(r): int((RUN == r).sum()) for r in sorted(set(RUN))}}
    # the soft threshold, fitted to the observed SNRs with a rho^-4 parent
    g = np.linspace(5, 200, 40000)
    def thr_nll(p, s, beta=4.0):
        r50, w = p
        if w < 0.05: return 1e9
        norm = np.trapezoid(g**-beta*ndtr((g-r50)/w), g)
        return -np.sum(-beta*np.log(s) + np.log(ndtr((s-r50)/w) + 1e-300) - np.log(norm))
    thr = {}
    for beta in (4.0, 4.5):
        r = min((minimize(thr_nll, [a, b], args=(SNR, beta), method='Nelder-Mead') for a in (8, 9, 10) for b in (0.5, 1, 2)), key=lambda r: r.fun)
        thr['beta%.1f' % beta] = [round(float(x), 2) for x in r.x]
    st = SNR[SNR > 12]
    thr['tail_slope_above_12'] = round(float(-(1 + len(st)/np.sum(np.log(st/12)))), 2)
    out['threshold_fit'] = thr
    base, (prof, rho, pdet, zz, PG, PE) = run_emp('base')
    out['base'] = base
    out['base_profile'] = [[float(a), round(float(p), 3)] for a, p in zip(AL, prof)]
    kap = SNR*D/rho(MDET)
    out['kappa_spread_ln'] = round(float(np.std(np.log(kap/np.array([np.median(kap[RUN == r]) for r in RUN])))), 3)
    chk = {}
    for nm in ('IMR', 'IMR+data', 'inspiral'):
        r_, _ = loudness(nm); y = np.log(SNR*D/r_(MDET))
        for r in set(RUN): y[RUN == r] -= np.median(y[RUN == r])
        edges = [0, 11, 16, 25, 35, 45, 60, 80, 200]
        chk[nm] = [[edges[k], edges[k+1], int(((MDET >= edges[k]) & (MDET < edges[k+1])).sum()),
                    round(float(np.median(y[(MDET >= edges[k]) & (MDET < edges[k+1])])), 3)] for k in range(len(edges)-1)]
    out['loudness_check'] = chk
    cases = []
    for nm, kw in [('threshold 8.6', dict(r50=8.6)), ('threshold 9.6', dict(r50=9.6)), ('sharper threshold, width 0.3', dict(w=0.3)),
                   ('loudness: design waveform, uncorrected', dict(loud='IMR')), ('loudness: inspiral law', dict(loud='inspiral')),
                   ('distances 20% larger', dict(Dscale=1.2)), ('distances 20% smaller', dict(Dscale=0.8)),
                   ('2023-2025 runs only', dict(sel=RUN >= 5)), ('2015-2020 runs only', dict(sel=RUN <= 4)),
                   ('detector chirp mass below 60', dict(sel=MDET < 60)), ('detector chirp mass above 15', dict(sel=MDET > 15)),
                   ('detector chirp mass above 40', dict(sel=MDET > 40))]:
        r_, _ = run_emp(nm, **kw); cases.append(r_)
    out['variations'] = cases
    # the first version, orientation-averaged selection with a fitted reach, for comparison
    cmp_ = []
    for sc in (0.9, 1.0, 1.2):
        r_, _ = run_case('orientation-averaged selection, reach %.1f' % sc, scale=sc); cmp_.append(r_)
    out['orientation_model'] = cmp_
    # figure bands, base case
    fig = {'events': [[int(D[i]), round(float(MDET[i]), 2), int(RUN[i])] for i in range(N)]}
    knots = {}
    for a in (0.0, 1.0):
        th = minimize(nll, np.zeros(len(KNOTS)), args=(a, zz, PG, PE, np.ones(N, bool)), method='L-BFGS-B', jac=True).x
        knots[a] = th; band = []
        for d in np.linspace(200, 7000, 35):
            zd = float(z_of_D(d)); w = np.zeros_like(MG)
            for r in sorted(set(RUN)):
                P = pdet(np.full_like(MG, d), MG, np.full(MG.shape, r)); ww = np.exp(logf(np.log(MG) - a*np.log1p(zd), th))*P
                if ww.sum() > 0: w += ww/ww.sum()*(RUN == r).sum()
            cdf = np.cumsum(w)/w.sum()
            band.append([int(d)] + [round(float(np.interp(p, cdf, MG)), 2) for p in (0.16, 0.5, 0.84)])
        fig['alpha%d' % int(a)] = band
    out['figure'] = fig
    # mock catalogues: each event keeps its catalogue distance and run; its true distance is drawn from the quoted
    # 90% interval (log-normal, sigma = half-width/1.645); its mass is drawn at that true distance under each hypothesis,
    # through the selection at the catalogue distance, then scattered by 5%. Fitted exactly as the real data.
    nmock = int(sys.argv[2]) if len(sys.argv) > 2 else 30
    sig = 0.5*(np.log(D+DHI) - np.log(np.maximum(D-DLO, 1)))/1.645
    ALM = np.round(np.arange(-0.6, 2.21, 0.2), 2)
    PGrid = PG
    mocks = {}
    for a in (0.0, 1.0):
        th = knots[a]; Ts, Ab = [], []
        for k in range(nmock):
            gg = np.random.default_rng(1000*int(a)+k)
            Dt = D*np.exp(gg.normal(0, 1, N)*sig)
            W = np.exp(logf(np.log(MG)[None, :] - a*np.log1p(z_of_D(Dt))[:, None], th))*PGrid
            cdf = np.cumsum(W, 1); cdf /= cdf[:, -1:]
            Mm = np.array([np.interp(gg.random(), cdf[i], MG) for i in range(N)])*np.exp(gg.normal(0, 0.05, N))
            keep = MDET.copy(); MDET[:] = Mm
            PEm = np.maximum(pdet(D, MDET, RUN), 1e-6)
            pm = profile(ALM, zz, PGrid, PEm, np.ones(N, bool))
            MDET[:] = keep
            Ts.append(round(float(np.interp(0.0, ALM, pm) - np.interp(1.0, ALM, pm)), 2)); Ab.append(float(ALM[int(np.argmin(pm))]))
        mocks['alpha%d' % int(a)] = {'T': Ts, 'alpha_best': Ab}
    mocks['observed_T'] = round(base['dlnL_alpha0'] - base['dlnL_alpha1'], 2)
    out['mocks'] = mocks
    print(json.dumps(out))
