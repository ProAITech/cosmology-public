"""Entry 187: figure1.svg - the cooling inside 0.2 Mpc implied by each stack (delta, top-hat normalisation, +-1 sigma from 1,000
rotations), against the published levels (-15 uK, Luparello et al. 2023; -27 uK, Hansen et al. 2025, dense filaments) and the
4-10 uK ceiling of Entry 47.  Reads results.json.  Run from the repository root."""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
A = json.load(open(os.path.join(HERE, 'results.json')))['all']
rows = [('P|smica_l5', 'SMICA (decides)', 'P'), ('P|nilc_l5', 'NILC', 'P'), ('P|sevem_l5', 'SEVEM', 'P'), ('P|commander_l5', 'Commander', 'P'),
        ('P_large|smica_l5', 'larger half', 'S'), ('P_dense|smica_l5', 'denser half', 'S'), ('A_2MRS|smica_l5', '2MRS, z &lt; 0.02 (control)', 'A')]
inj = (-14.09, 3.26)
W, H, L, T, PW = 640, 420, 170, 50, 430
rowh = 36; PH = rowh * (len(rows) + 1)
xmin, xmax = -35, 12
X = lambda v: L + PW * (v - xmin) / (xmax - xmin)
o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {T + PH + 66}" font-family="Helvetica, Arial, sans-serif" font-size="12">', '<rect width="100%" height="100%" fill="#fff"/>']
o.append(f'<rect x="{X(-10):.1f}" y="{T}" width="{X(-4) - X(-10):.1f}" height="{PH}" fill="#f1e6c8"/>')
o.append(f'<text x="{X(-7):.1f}" y="{T - 8}" text-anchor="middle" fill="#8a6d1f" font-size="11">ceiling, Entry 47</text>')
for v in range(-35, 11, 5):
    o.append(f'<line x1="{X(v):.1f}" x2="{X(v):.1f}" y1="{T}" y2="{T + PH}" stroke="#e6e6e6"/>')
    lab = ('0' if v == 0 else f'{v:+d}').replace('-', '\u2212')
    o.append(f'<text x="{X(v):.1f}" y="{T + PH + 16}" text-anchor="middle" fill="#333">{lab}</text>')
o.append(f'<line x1="{X(0):.1f}" x2="{X(0):.1f}" y1="{T}" y2="{T + PH}" stroke="#888"/>')
for v, lab in ((-15, '−15 μK, nearby spirals'), (-27, '−27 μK, dense filaments')):
    o.append(f'<line x1="{X(v):.1f}" x2="{X(v):.1f}" y1="{T}" y2="{T + PH}" stroke="#b5462f" stroke-width="1.5" stroke-dasharray="5 4"/>')
    o.append(f'<text x="{X(v):.1f}" y="{T - 8}" text-anchor="middle" fill="#b5462f" font-size="11">{lab.split(",")[0]}</text>')
for i, (k, lab, kind) in enumerate(rows + [('inj', 'injected −15 μK (test)', 'I')]):
    y = T + rowh * (i + 0.75)
    if kind == 'I': d, s = inj
    else: d, s = A[k]['delta_top'], A[k]['sd_delta_top']
    col = {'P': '#2f5bb7', 'S': '#2f5bb7', 'A': '#7a7a7a', 'I': '#9a9a9a'}[kind]
    o.append(f'<line x1="{X(d - s):.1f}" x2="{X(d + s):.1f}" y1="{y:.1f}" y2="{y:.1f}" stroke="{col}" stroke-width="2"/>')
    o.append(f'<circle cx="{X(d):.1f}" cy="{y:.1f}" r="{6 if i == 0 else 4.5}" fill="{col if kind != "I" else "#fff"}" stroke="{col}" stroke-width="1.5"/>')
    o.append(f'<text x="{L - 10}" y="{y + 4:.1f}" text-anchor="end" fill="#222"{" font-weight=\"bold\"" if i == 0 else ""}>{lab}</text>')
o.append(f'<text x="{L + PW / 2:.1f}" y="{T + PH + 36}" text-anchor="middle" fill="#333">cooling inside 0.2 Mpc implied by the stack (μK, ±1σ)</text>')
o.append(f'<text x="{L + PW / 2:.1f}" y="{T + PH + 54}" text-anchor="middle" fill="#b5462f" font-size="11">dashed: published coolings, −15 μK for nearby spirals and −27 μK for spirals in dense filaments</text>')
o.append(f'<text x="{L - 160}" y="22" fill="#222" font-size="13" font-weight="bold">3,583 late-type spirals at z = 0.02–0.05, Planck PR3</text>')
o.append('</svg>')
open(os.path.join(HERE, 'figure1.svg'), 'w', encoding='utf-8').write('\n'.join(o)); print('figure1.svg written')
