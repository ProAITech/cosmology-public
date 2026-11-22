"""Workbench Entry 170, check: the neutron-injection scan redone with an independent set of nuclear rates.

neutron_li.py used the network of Entry 53 (bbn_net.py), whose thermonuclear rates are Kawano's 1993 fits.  Here the same scan is run
with the tabulated rates of the PRIMAT compilation (Pitrou, Coc, Uzan & Vangioni 2018, Phys. Rep. 754, 1), and again with the NACRE II
set, both read from the plain-text tables distributed with the public code PRyMordial (Burns, Tait & Valli 2024, arXiv:2307.07061):
  data/prymordial/PRyMordial-main.zip   https://github.com/vallima/PRyMordial (archive of the main branch), 1,653,114 bytes,
      sha256 e7d41823af285241d8bb28c19f40f772132176781b723d99e583bebce81d8074, fetched 21 November 2026, unpacked beside the zip.
      Tables used: PRyMrates/nuclear/key_primat_rates/*.txt and key_nacreii_rates/*.txt (twelve reactions each) and
      other_nucl_rates/Be7naa.txt; three columns: T9, N_A<sigma v> in cm3 s-1 mol-1, uncertainty factor; T9 = 0.001 to 10.
Only the tables are read.  PRyMordial itself is not run: this is a check of the nuclear inputs, not a rerun in an independent code.
The thermodynamics, the weak rates and the solver are still Entry 53's.

Three variants:
  kawano     Entry 170 as published (the baseline, rerun here at the lowest-deuterium epochs only)
  primat     the eleven reactions of Entry 53 with PRIMAT's rates, plus two reactions Entry 53 lacks:
             t(p,gamma)4He and 7Be(n,alpha)4He (a second neutron channel that destroys 7Be)
  nacreii    the same with the NACRE II set for the twelve key rates
Reverse rates keep Entry 53's detailed-balance factors.  Above T9 = 10, where the tables end and no nucleus yet survives, Kawano's
fits are used.

Run from code/:  ../.venv/bin/python -I neutron_li_primat.py      (about 25 minutes)
"""
import os, sys, json, inspect, textwrap, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import bbn_net as bn

TABLES = os.path.join(HERE, '..', 'data', 'prymordial', 'PRyMordial-main', 'PRyMrates', 'nuclear')
MAP = dict(npdg='npdg', dpg='dpHe3g', ddn='ddHe3n', ddp='ddtp', tdn='tdan', tag='taLi7g', hent='He3ntp', hdp='He3dap', hag='He3aBe7g',
           benp='Be7nLi7p', lipa='Li7paa', tpag='tpag')
ETA, T9F = 6.124, 0.006
D_OBS, D_ERR = 2.527e-5, 0.030e-5

# ---- the network with a neutron source (as in neutron_li.py) and two more reactions
src = textwrap.dedent(inspect.getsource(bn.run))
def swap(a, b):
    global src
    assert src.count(a) == 1, (src.count(a), a)
    src = src.replace(a, b)
swap("def run(eta10=6.10,", "def run_x(A_inj=0.0, T9c=0.3, wid=0.15, extra=True, eta10=6.10,")
swap("        w=lnp*Y[N]-lpn*Y[P]; d=np.zeros(9); d[0]=dT\n",
     "        r12=R*f['benaa']*Y[B]*Y[N] if extra else 0.0            # 7Be + n -> 2 4He\n"
     "        r13=R*f['tpag']*Y[T]*Y[P] if extra else 0.0             # t + p -> 4He + gamma\n"
     "        w=lnp*Y[N]-lpn*Y[P]; d=np.zeros(9); d[0]=dT\n")
swap("d[1+N]=-w-r1+r3-r5+r6-r10",
     "d[1+N]=-w-r1+r3-r5+r6-r10-r12+A_inj*np.exp(-0.5*(np.log(T9/T9c)/wid)**2)/(np.sqrt(2*np.pi)*wid)*(-dT/T9)")
swap("d[1+P]= w-r1-r2+r4+r5+r7+r10-r11", "d[1+P]= w-r1-r2+r4+r5+r7+r10-r11-r13")
swap("d[1+T]= r4+r5-r6-r9", "d[1+T]= r4+r5-r6-r9-r13")
swap("d[1+A]= r6+r7-r8-r9+2*r11", "d[1+A]= r6+r7-r8-r9+2*r11+2*r12+r13")
swap("d[1+B]= r8-r10", "d[1+B]= r8-r10-r12")
exec(src, bn.__dict__)
run_x = bn.run_x
KAWANO = bn.rates


def load(folder, name):
    a = np.loadtxt(os.path.join(TABLES, folder, name + '.txt'))
    return np.log(a[:, 0]), np.log(a[:, 1])


def table_rates(key_folder):
    tabs = {k: load(key_folder, v) for k, v in MAP.items()}
    tabs['benaa'] = load('other_nucl_rates', 'Be7naa')
    lo, hi = np.exp(tabs['npdg'][0][0]), np.exp(tabs['npdg'][0][-1])
    def rates(t9):
        f = KAWANO(t9)
        f.setdefault('benaa', 0.0); f.setdefault('tpag', 0.0)
        if lo <= t9 <= hi:
            lt = np.log(t9)
            for k, (x, y) in tabs.items():
                f[k] = float(np.exp(np.interp(lt, x, y)))
        return f
    return rates


def kawano_rates(t9):
    f = KAWANO(t9); f['benaa'] = 0.0; f['tpag'] = 0.0
    return f


def checked(r):
    assert r['t_end'] > 1e4, 'integration stopped early'
    return r


def solve_A(T9c, base, extra, target=1 / 3.):
    """Neutrons per baryon that bring total A = 7 to target x standard (neutron_li.py's search)."""
    lo, hi, A, hist = None, None, 1e-5, []
    for it in range(14):
        r = checked(run_x(A_inj=A, T9c=T9c, extra=extra, eta10=ETA, T9f=T9F))
        ratio = r['Li7'] / base['Li7']
        hist.append((A, ratio, r))
        if abs(np.log(ratio / target)) < 0.01:
            return A, r, hist
        if ratio > target:
            lo = (A, ratio)
        else:
            hi = (A, ratio)
        if lo and hi:
            la, lb = np.log(lo[0]), np.log(hi[0]); ra, rb = np.log(lo[1]), np.log(hi[1])
            A = float(np.exp(la + (np.log(target) - ra) * (lb - la) / (rb - ra)))
        elif lo:
            A *= 10.0 if ratio > 0.9 else 3.0
        else:
            A /= 3.0
        if A > 0.01:
            return None, r, hist
    return A, r, hist


def scan(label, rates, extra, epochs, t0):
    bn.rates = rates
    base = checked(run_x(A_inj=0.0, extra=extra, eta10=ETA, T9f=T9F))
    out = dict(standard=dict(Yp=base['Yp'], D=base['D'], He3=base['He3'], Li7=base['Li7'], Be_frac=base['Be_frac']), rows=[])
    print(label, 'standard', {k: (round(v, 5) if k in ('Yp', 'Be_frac') else float('%.4g' % v)) for k, v in out['standard'].items()}, flush=True)
    for T9c in epochs:
        A, r, hist = solve_A(T9c, base, extra)
        if A is None:
            r = min(hist, key=lambda h: h[1])[2]
        row = dict(T9c=T9c, keV=round(T9c * 86.17, 2), solved=A is not None, A_per_baryon=A, Li_ratio=r['Li7'] / base['Li7'], D=r['D'],
                   D_rise=r['D'] / base['D'] - 1.0, D_pull_vs_obs=(r['D'] - base['D']) / D_ERR, He3_rise=r['He3'] / base['He3'] - 1.0,
                   dYp=r['Yp'] - base['Yp'], runs=len(hist))
        out['rows'].append(row)
        print('%s T9c %.2f (%5.1f keV): A %.3g per baryon, Li x%.3f, D %+.1f%% (%.0f sigma), 3He %+.1f%%, Yp %+.4f  [%d runs, %d s]' % (
            label, T9c, T9c * 86.17, A or float('nan'), row['Li_ratio'], 100 * row['D_rise'], row['D_pull_vs_obs'], 100 * row['He3_rise'],
            row['dYp'], len(hist), time.time() - t0), flush=True)
    solved = [x for x in out['rows'] if x['solved']]
    if solved:
        best = min(solved, key=lambda x: x['D_rise'])
        out['minimum'] = dict(T9c=best['T9c'], keV=best['keV'], D_rise=best['D_rise'], D=best['D'], A_per_baryon=best['A_per_baryon'],
                              D_pull_vs_obs=best['D_pull_vs_obs'])
    return out


def main():
    t0 = time.time()
    out = dict(eta10=ETA, T9_end=T9F, tables='PRyMordial main branch, PRyMrates/nuclear (see the header)',
               note='A check of the nuclear inputs only: thermodynamics, weak rates and solver are bbn_net.py (Entry 53).')
    out['kawano'] = scan('kawano', kawano_rates, False, (0.7, 0.5, 0.4), t0)
    out['primat'] = scan('primat', table_rates('key_primat_rates'), True, (1.0, 0.7, 0.5, 0.4, 0.35, 0.3, 0.25, 0.15, 0.05), t0)
    out['nacreii'] = scan('nacreii', table_rates('key_nacreii_rates'), True, (0.7, 0.5, 0.4, 0.3), t0)
    out['seconds'] = round(time.time() - t0)
    json.dump(out, open(os.path.join(HERE, 'neutron_li_primat.json'), 'w'), indent=1)
    print('minima', {k: out[k].get('minimum') for k in ('kawano', 'primat', 'nacreii')})
    print('written neutron_li_primat.json', out['seconds'], 's')


if __name__ == '__main__':
    main()
