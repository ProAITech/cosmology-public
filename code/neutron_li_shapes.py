"""Workbench Entry 170, check: does the shape of the neutron source change the least deuterium excess?

neutron_li.py injected neutrons as a Gaussian in ln T of width 0.15.  Here the same question (how many neutrons cut A = 7 to a third,
and what does that do to D/H?) is asked for narrower and wider bursts and for the source a decaying particle would make,
dYn/dt = (A / tau_X) exp(-t / tau_X).  The network is neutron_li_primat.py's: Entry 53's with PRIMAT's tabulated rates and the two
added reactions, t(p,gamma)4He and 7Be(n,alpha)4He.  Data: the rate tables cited in neutron_li_primat.py.

Run from code/:  ../.venv/bin/python -I neutron_li_shapes.py        (about 40 minutes)
                 ../.venv/bin/python -I neutron_li_shapes.py wide   (about 15 minutes: adds hotter centres for the two wide bursts,
                                                                     whose least-deuterium centre lay above the first grid, and
                                                                     merges them into neutron_li_shapes.json)
"""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import bbn_net as bn
import neutron_li_primat as P

src = P.src
a = "def run_x(A_inj=0.0,"
b = "+A_inj*np.exp(-0.5*(np.log(T9/T9c)/wid)**2)/(np.sqrt(2*np.pi)*wid)*(-dT/T9)"
assert src.count(a) == 1 and src.count(b) == 1
src = src.replace(a, "def run_s(tauX=0.0, A_inj=0.0,")
src = src.replace(b, "+(A_inj*np.exp(-t/tauX)/tauX if tauX>0 else A_inj*np.exp(-0.5*(np.log(T9/T9c)/wid)**2)/(np.sqrt(2*np.pi)*wid)*(-dT/T9))")
exec(src, bn.__dict__)
run_s = bn.run_s
bn.rates = P.table_rates('key_primat_rates')
ETA, T9F, D_ERR = P.ETA, P.T9F, P.D_ERR


def solve(base, target=1 / 3., **kw):
    lo, hi, A, hist = None, None, 1e-5, []
    for it in range(14):
        r = P.checked(run_s(A_inj=A, eta10=ETA, T9f=T9F, **kw))
        ratio = r['Li7'] / base['Li7']
        hist.append((A, ratio, r))
        if abs(np.log(ratio / target)) < 0.01:
            return A, r, len(hist)
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
        if A > 0.02:
            return None, min(hist, key=lambda h: h[1])[2], len(hist)
    return A, r, len(hist)


def main():
    t0 = time.time()
    wide = len(sys.argv) > 1 and sys.argv[1] == 'wide'
    base = P.checked(run_s(A_inj=0.0, eta10=ETA, T9f=T9F))
    out = dict(standard=dict(Yp=base['Yp'], D=base['D'], He3=base['He3'], Li7=base['Li7']), gaussian=[], decay=[])
    if wide:
        out = json.load(open(os.path.join(HERE, 'neutron_li_shapes.json')))
    def row(kind, A, r, n, **kw):
        d = dict(kw, solved=A is not None, A_per_baryon=A, Li_ratio=r['Li7'] / base['Li7'], D=r['D'], D_rise=r['D'] / base['D'] - 1.0,
                 D_pull_vs_obs=(r['D'] - base['D']) / D_ERR, He3_rise=r['He3'] / base['He3'] - 1.0, dYp=r['Yp'] - base['Yp'], runs=n)
        out[kind].append(d)
        print(kind, {k: (float('%.4g' % v) if isinstance(v, float) else v) for k, v in d.items()}, '%d s' % (time.time() - t0), flush=True)
    grid = [(wid, T9c) for wid in (0.05, 0.3, 0.6) for T9c in (0.6, 0.5, 0.45, 0.4)]
    if wide:
        grid = [(0.3, 0.8), (0.3, 0.7), (0.3, 0.65), (0.6, 1.0), (0.6, 0.8), (0.6, 0.7)]
    for wid, T9c in grid:
        A, r, n = solve(base, T9c=T9c, wid=wid)
        row('gaussian', A, r, n, wid=wid, T9c=T9c, keV=round(T9c * 86.17, 1))
    for tau in (() if wide else (100., 200., 300., 450., 600., 800., 1100., 1500., 2500., 5000., 20000.)):
        A, r, n = solve(base, tauX=tau)
        row('decay', A, r, n, tau_s=tau)
    for kind in ('gaussian', 'decay'):
        ok = [x for x in out[kind] if x['solved']]
        if ok:
            out['minimum_' + kind] = min(ok, key=lambda x: x['D_rise'])
    out['minimum_by_width'] = {str(wd): min((x for x in out['gaussian'] if x['solved'] and x['wid'] == wd), key=lambda x: x['D_rise'])
                               for wd in sorted({x['wid'] for x in out['gaussian'] if x['solved']})}
    out['seconds'] = out.get('seconds', 0) + round(time.time() - t0) if wide else round(time.time() - t0)
    json.dump(out, open(os.path.join(HERE, 'neutron_li_shapes.json'), 'w'), indent=1)
    print('minima', out.get('minimum_gaussian'), out.get('minimum_decay'))


if __name__ == '__main__':
    main()
