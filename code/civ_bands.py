"""Entry 130: Alfven-Arrhenius critical-velocity bands against the modern solar system and exoplanets.

Gamma = GM/r (km^2/s^2), the gravitational potential per unit mass at a body's orbit.
Hypothesis: infalling gas is ionized where its free-fall speed reaches v_c = sqrt(2 e V_ion / m),
i.e. at Gamma_ion = v_c^2 / 2; grains then settle at 2/3 of that distance, Gamma_orb = 0.75 v_c^2.
Bodies should therefore sit in bands of Gamma fixed by atomic constants, the same for every central body.
"""
import json, math
import numpy as np
from scipy.stats import ks_2samp

e, u = 1.602176634e-19, 1.66053907e-27
# first ionization energy (eV, NIST) and mean atomic mass (u)
EL = {'H': (13.598, 1.008), 'He': (24.587, 4.0026), 'C': (11.260, 12.011), 'N': (14.534, 14.007),
      'O': (13.618, 15.999), 'Ne': (21.565, 20.180), 'Na': (5.139, 22.990), 'Mg': (7.646, 24.305),
      'Al': (5.986, 26.982), 'Si': (8.152, 28.085), 'S': (10.360, 32.06), 'Ca': (6.113, 40.078),
      'Fe': (7.902, 55.845), 'Ni': (7.640, 58.693), 'Cs': (3.894, 132.905), 'U': (6.194, 238.03)}
vc = {k: math.sqrt(2 * V * e / (A * u)) / 1e3 for k, (V, A) in EL.items()}
GROUPS = {'H': ['H'], 'He': ['He'], 'CNO': ['C', 'N', 'O', 'Ne'],
          'metals': ['Na', 'Mg', 'Al', 'Si', 'S', 'Ca', 'Fe', 'Ni']}

def zones(k):
    return {g: (math.log10(k * min(vc[x] for x in L) ** 2), math.log10(k * max(vc[x] for x in L) ** 2))
            for g, L in GROUPS.items()}

VAR = {'A': 0.75, 'B': 0.5}   # A: with the 2/3 settling law; B: at the ionization radius

GM = {'Sun': 1.32712440018e11, 'Jupiter': 1.266865319e8, 'Saturn': 3.79312062e7, 'Uranus': 5.7939513e6,
      'Neptune': 6.8350995e6, 'Earth': 398600.4, 'Mars': 42828.37, 'Pluto': 975.5}
RAD = {'Sun': 695700, 'Jupiter': 71492, 'Saturn': 60268, 'Uranus': 25559, 'Neptune': 24764,
       'Earth': 6378.1, 'Mars': 3396.2, 'Pluto': 1188.3}
MASS_KG = {k: v / 6.6743e-20 for k, v in GM.items()}
AU = 1.495978707e8

# (primary, name, a_km, discovery year); regular = prograde, e < 0.2, i < 20 deg (JPL SSD mean elements)
B = [('Sun', 'Mercury', 0.38710 * AU, 0), ('Sun', 'Venus', 0.72333 * AU, 0), ('Sun', 'Earth', 1.0 * AU, 0),
     ('Sun', 'Mars', 1.52368 * AU, 0), ('Sun', 'Jupiter', 5.2026 * AU, 0), ('Sun', 'Saturn', 9.5549 * AU, 0),
     ('Sun', 'Uranus', 19.2184 * AU, 1781), ('Sun', 'Neptune', 30.1104 * AU, 1846),
     ('Jupiter', 'Io', 421800, 1610), ('Jupiter', 'Europa', 671100, 1610), ('Jupiter', 'Ganymede', 1070400, 1610),
     ('Jupiter', 'Callisto', 1882700, 1610), ('Jupiter', 'Amalthea', 181400, 1892), ('Jupiter', 'Thebe', 221900, 1979),
     ('Jupiter', 'Adrastea', 129000, 1979), ('Jupiter', 'Metis', 128000, 1979),
     ('Saturn', 'Mimas', 186000, 1789), ('Saturn', 'Enceladus', 238400, 1789), ('Saturn', 'Tethys', 295000, 1684),
     ('Saturn', 'Dione', 377700, 1684), ('Saturn', 'Rhea', 527200, 1672), ('Saturn', 'Titan', 1221900, 1655),
     ('Saturn', 'Hyperion', 1481500, 1848), ('Saturn', 'Iapetus', 3561700, 1671), ('Saturn', 'Helene', 377600, 1980),
     ('Saturn', 'Telesto', 295000, 1980), ('Saturn', 'Calypso', 295000, 1980), ('Saturn', 'Methone', 194700, 2004),
     ('Saturn', 'Polydeuces', 377600, 2004), ('Saturn', 'Janus', 151500, 1966), ('Saturn', 'Epimetheus', 151400, 1966),
     ('Saturn', 'Atlas', 137700, 1980), ('Saturn', 'Prometheus', 139400, 1980), ('Saturn', 'Pandora', 141700, 1980),
     ('Saturn', 'Pan', 133600, 1990), ('Saturn', 'Pallene', 212300, 2004), ('Saturn', 'Daphnis', 136500, 2005),
     ('Saturn', 'Anthe', 198100, 2007), ('Saturn', 'Aegaeon', 167500, 2008), ('Saturn', 'S/2009 S 1', 117061, 2009),
     ('Uranus', 'Ariel', 190929, 1851), ('Uranus', 'Umbriel', 265986, 1851), ('Uranus', 'Titania', 436298, 1787),
     ('Uranus', 'Oberon', 583511, 1787), ('Uranus', 'Miranda', 129846, 1948), ('Uranus', 'Puck', 86004, 1985),
     ('Uranus', 'Cordelia', 49755, 1986), ('Uranus', 'Ophelia', 53765, 1986), ('Uranus', 'Bianca', 59170, 1986),
     ('Uranus', 'Cressida', 61770, 1986), ('Uranus', 'Desdemona', 62663, 1986), ('Uranus', 'Juliet', 64362, 1986),
     ('Uranus', 'Portia', 66101, 1986), ('Uranus', 'Rosalind', 69930, 1986), ('Uranus', 'Belinda', 75258, 1986),
     ('Uranus', 'Perdita', 76418, 1999), ('Uranus', 'Mab', 97737, 2003), ('Uranus', 'Cupid', 74396, 2003),
     ('Uranus', 'S/2025 U 1', 57844, 2025),
     ('Neptune', 'Naiad', 48300, 1989), ('Neptune', 'Thalassa', 50100, 1989), ('Neptune', 'Despina', 52600, 1989),
     ('Neptune', 'Galatea', 62000, 1989), ('Neptune', 'Larissa', 73600, 1981), ('Neptune', 'Proteus', 117700, 1989),
     ('Neptune', 'Hippocamp', 105300, 2013)]
SMALL = [('Earth', 'Moon', 384400, 0), ('Mars', 'Phobos', 9375, 1877), ('Mars', 'Deimos', 23457, 1877),
         ('Pluto', 'Charon', 19600, 1978), ('Pluto', 'Styx', 43200, 2012), ('Pluto', 'Nix', 49300, 2005),
         ('Pluto', 'Kerberos', 58300, 2011), ('Pluto', 'Hydra', 65200, 2005)]
RINGS = {'Saturn A–C rings': ('Saturn', 74658, 136775), 'Main asteroid belt': ('Sun', 2.1 * AU, 3.3 * AU),
         'Classical Kuiper belt': ('Sun', 42 * AU, 48 * AU)}

def lg(p, a): return math.log10(GM[p] / a)
bodies = [dict(p=p, n=n, a=a, y=y, lg=lg(p, a), rR=a / RAD[p]) for p, n, a, y in B]
old = [b for b in bodies if b['y'] < 1976]
new = [b for b in bodies if b['y'] >= 1976]

def dist(x, Z):
    return min(max(0.0, lo - x, x - hi) for lo, hi in Z.values())

def frac_in(xs, Z, w, d=0.0):
    Zs = {k: (lo + d, hi + d) for k, (lo, hi) in Z.items()}
    return float(np.mean([dist(x, Zs) <= w for x in xs]))

DELTA = np.round(np.arange(-0.6, 0.6001, 0.002), 3)
out = {'vc': {k: round(v, 2) for k, v in vc.items()}, 'zones': {}, 'enrich': {}, 'wmin': {}, 'gap': {}}
for V, k in VAR.items():
    Z = zones(k); out['zones'][V] = {g: [round(10 ** lo, 1), round(10 ** hi, 1)] for g, (lo, hi) in Z.items()}
    for label, S in (('old', old), ('new', new), ('all', bodies)):
        xs = [b['lg'] for b in S]
        for w in (0.05, 0.1, 0.15, 0.2):
            fo = frac_in(xs, Z, w); fn = np.array([frac_in(xs, Z, w, d) for d in DELTA])
            out['enrich'][f'{V}_{label}_{w}'] = dict(n=len(xs), f=round(fo, 3), null=round(float(fn.mean()), 3),
                                                     ratio=round(fo / fn.mean(), 2), p=round(float(np.mean(fn >= fo - 1e-12)), 3))
    # smallest common half-width that puts every pre-1976 body in a band, and what that leaves uncovered
    dd = sorted(((dist(b['lg'], Z), b['n']) for b in old), reverse=True)
    w0 = dd[0][0]
    lo_r, hi_r = min(b['lg'] for b in old), max(b['lg'] for b in old)
    grid = np.linspace(lo_r, hi_r, 20001)
    cov = float(np.mean([dist(x, Z) <= w0 + 1e-9 for x in grid]))
    gaps = []; inside = False
    for x in grid:
        g_ = dist(x, Z) > w0 + 1e-9
        if g_ and not inside: start = x; inside = True
        if not g_ and inside: gaps.append((10 ** start, 10 ** x)); inside = False
    if inside: gaps.append((10 ** start, 10 ** hi_r))
    out['wmin'][V] = dict(w=round(w0, 3), factor=round(10 ** w0, 2), worst=[(n, round(d, 3)) for d, n in dd[:5]],
                          coverage=round(cov, 3), gaps=[[round(a, 1), round(b, 1)] for a, b in gaps])
    # what sits in the leftover gaps
    occ = [b['n'] + f" ({b['p']}, {b['y'] if b['y'] else 'ancient'})" for b in bodies
           if any(a <= 10 ** b['lg'] <= c for a, c in gaps)]
    ring_occ = [r for r, (p, a1, a2) in RINGS.items() if any(not (GM[p] / a1 < a or GM[p] / a2 > c) for a, c in gaps)]
    out['gap'][V] = dict(bodies=occ, rings=ring_occ)

# scaling: do giant-planet satellite systems line up in Gamma (atomic bands) or in planetary radii (tides)?
sc = {}
for p in ('Jupiter', 'Saturn', 'Uranus', 'Neptune'):
    S = [b for b in bodies if b['p'] == p]
    inner = min(S, key=lambda b: b['a'])
    sc[p] = dict(n=len(S), med_G=10 ** float(np.median([b['lg'] for b in S])),
                 med_rR=float(np.median([b['rR'] for b in S])), inner=inner['n'],
                 inner_G=10 ** inner['lg'], inner_rR=inner['rR'], GM_R=GM[p] / RAD[p])
def sd(v): return float(np.std(np.log10(v), ddof=1))
out['scaling'] = dict(systems={k: {kk: (round(vv, 3) if isinstance(vv, float) else vv) for kk, vv in v.items()} for k, v in sc.items()},
                      sd_medG=round(sd([v['med_G'] for v in sc.values()]), 3),
                      sd_medrR=round(sd([v['med_rR'] for v in sc.values()]), 3),
                      sd_innerG=round(sd([v['inner_G'] for v in sc.values()]), 3),
                      sd_innerrR=round(sd([v['inner_rR'] for v in sc.values()]), 3))
ks = {'G': [], 'rR': []}
P4 = ['Jupiter', 'Saturn', 'Uranus', 'Neptune']
for i in range(4):
    for j in range(i + 1, 4):
        A = [b for b in bodies if b['p'] == P4[i]]; Bq = [b for b in bodies if b['p'] == P4[j]]
        ks['G'].append(ks_2samp([b['lg'] for b in A], [b['lg'] for b in Bq]).statistic)
        ks['rR'].append(ks_2samp([math.log10(b['rR']) for b in A], [math.log10(b['rR']) for b in Bq]).statistic)
out['scaling']['ks_mean'] = {k: round(float(np.mean(v)), 3) for k, v in ks.items()}
# slope of innermost Gamma against GM/R (tidal scaling predicts 1, a fixed band predicts 0)
x = np.log10([v['GM_R'] for v in sc.values()]); y = np.log10([v['inner_G'] for v in sc.values()])
out['scaling']['inner_slope'] = round(float(np.polyfit(x, y, 1)[0]), 3)

# forbidden zones
floorA = {k: round(0.75 * vc[k] ** 2, 2) for k in ('Ni', 'Fe', 'Cs', 'U')}
floorB = {k: round(0.5 * vc[k] ** 2, 2) for k in ('Ni', 'Fe', 'Cs', 'U')}
out['floor'] = dict(A=floorA, B=floorB, H_line=dict(A=round(0.75 * vc['H'] ** 2), B=round(0.5 * vc['H'] ** 2)))
out['small'] = [dict(n=n, p=p, G=round(GM[p] / a, 4), y=y, rR=round(a / RAD[p], 2)) for p, n, a, y in SMALL]
didy = 6.6743e-20 * 5.6e11 / 1.19
out['didymos_G'] = didy
out['below_floor_regular'] = [(b['n'], round(10 ** b['lg'], 2)) for b in bodies if 10 ** b['lg'] < floorA['Fe']]

# exoplanets: NASA Exoplanet Archive pscomppars, queried 11 Nov 2026 (6,455 with host mass and orbit);
# histograms of log10 Gamma from -2.0 in 0.1-dex bins (last bin open-ended)
z = lambda n: [0] * n
EX = {'Transit': z(42) + [1, 0, 0, 3, 2, 5, 10, 12, 35, 52, 87, 102, 147, 239, 268, 399, 445, 525, 494, 511, 533, 420, 202, 135, 89, 60, 15, 7],
      'RV': z(33) + [1, 1, 1, 13, 11, 15, 15, 21, 33, 26, 39, 42, 52, 77, 55, 79, 64, 85, 59, 56, 48, 53, 50, 50, 33, 51, 47, 30, 33, 30, 21, 8, 4, 1, 1, 0, 0],
      'Imaging': [0, 0, 0, 1, 0, 0, 1, 3, 2, 4, 1, 1, 1, 1, 0, 1, 3, 2, 2, 2, 3, 3, 2, 6, 3, 3, 0, 5, 3, 4, 1, 2, 4, 2, 2, 2, 5, 4, 3, 1, 4, 2, 0, 0, 1] + z(25),
      'Microlensing': z(32) + [1, 4, 2, 2, 5, 7, 10, 36, 35, 50, 72, 43, 10, 3, 3, 2, 3] + z(21),
      'Other': z(32) + [1, 0, 0, 0, 0, 1, 4, 1, 4, 2, 1, 5, 2, 2, 5, 2, 3, 1, 4, 1, 0, 2, 5, 3, 5, 5, 2, 0, 0, 0, 0, 0, 3, 3, 0, 2, 2, 3]}
for k, v in EX.items(): assert len(v) == 70, (k, len(v))
for k, nk in (('Transit', 4798), ('RV', 1205), ('Imaging', 90), ('Microlensing', 288), ('Other', 74)): assert sum(EX[k]) == nk, (k, sum(EX[k]))
tot = np.sum([v for v in EX.values()], axis=0)
out['exo'] = dict(n=int(tot.sum()), above_A=5102, above_B=5328, giants=1893, giants_above_A=919, hot_jupiters=717,
                  rv_giants=781, rv_giants_above_A=132, below_U_floor=45, max_G=282269, min_G=0.0209,
                  usp=[('TOI-2431 b', 2.08, 93078), ('TOI-2260 b', 2.22, 90542), ('TOI-2109 b', 2.27, 71971)],
                  hist=tot.tolist())
out['fig'] = dict(bodies=[[b['p'], b['n'], round(b['lg'], 3), b['y']] for b in bodies],
                  small=[[p, n, round(math.log10(GM[p] / a), 3), y] for p, n, a, y in SMALL],
                  logM={k: round(math.log10(v), 3) for k, v in MASS_KG.items()},
                  rings={k: [p, round(math.log10(GM[p] / a2), 3), round(math.log10(GM[p] / a1), 3)] for k, (p, a1, a2) in RINGS.items()},
                  zonesA={g: [round(lo, 3), round(hi, 3)] for g, (lo, hi) in zones(0.75).items()},
                  zonesB={g: [round(lo, 3), round(hi, 3)] for g, (lo, hi) in zones(0.5).items()},
                  rR={p: [[b['n'], round(math.log10(b['rR']), 3), round(b['lg'], 3), b['y']] for b in bodies if b['p'] == p] for p in P4})
json.dump(out, open('civ_bands.json', 'w'), ensure_ascii=False, indent=1)
print(json.dumps({k: v for k, v in out.items() if k != 'fig' and k != 'exo'}, ensure_ascii=False, indent=1))
print({k: v for k, v in out['exo'].items() if k != 'hist'})
