"""Entry 186: figure1.svg - excess RM variance against distance from the nearest eRASS1 cluster, NVSS (1.4 GHz, this entry) beside
SPICE-RACS (0.9 GHz, Entry 176), with the field model's lowest prediction inside 0.5 R500.  Reads results.json.
Run from the repository root:  .venv/bin/python -I entries/186-nvss-clusters/figure.py"""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
R = {r['range']: r for r in json.load(open(os.path.join(HERE, 'results.json')))['ranges']}
spice = {'0-0.5': (174.6, 21.5), '0.5-1': (79.7, 11.9), '1-2': (-0.6, 7.1), '2-4': (-0.6, 4.7)}
cats = ['0-0.5', '0.5-1', '1-2', '2-4']; labels = ['&lt; 0.5', '0.5–1', '1–2', '2–4']
W, H, L, T, PW, PH = 620, 380, 70, 40, 520, 280
ymin, ymax = -300, 600
Y = lambda v: T + PH * (ymax - v) / (ymax - ymin)
X = lambda i, off: L + PW * (i + 0.5) / 4 + off
o = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" font-family="Helvetica, Arial, sans-serif" font-size="12">' % (W, H), '<rect width="100%" height="100%" fill="#fff"/>']
for v in range(ymin, ymax + 1, 150):
    o.append('<line x1="%d" x2="%d" y1="%.1f" y2="%.1f" stroke="#e4e4e4"/>' % (L, L + PW, Y(v), Y(v)))
    o.append('<text x="%d" y="%.1f" text-anchor="end" fill="#333">%d</text>' % (L - 6, Y(v) + 4, v))
o.append('<line x1="%d" x2="%d" y1="%.1f" y2="%.1f" stroke="#888"/>' % (L, L + PW, Y(0), Y(0)))
o.append('<line x1="%.1f" x2="%.1f" y1="%.1f" y2="%.1f" stroke="#b5462f" stroke-width="2" stroke-dasharray="6 4"/>' % (X(0, -50), X(0, 50), Y(500), Y(500)))
o.append('<text x="%.1f" y="%.1f" fill="#b5462f">field model, lowest prediction</text>' % (X(0, -48), Y(500) - 6))
for i, c in enumerate(cats):
    for (val, err), off, col, mk in ((spice[c], -12, '#2f5bb7', 'c'), ((R[c]['excess'], R[c]['null_sd']), 12, '#d08a1e', 's')):
        o.append('<line x1="%.1f" x2="%.1f" y1="%.1f" y2="%.1f" stroke="%s" stroke-width="2"/>' % (X(i, off), X(i, off), Y(max(val - err, ymin)), Y(min(val + err, ymax)), col))
        if mk == 'c':
            o.append('<circle cx="%.1f" cy="%.1f" r="5" fill="%s"/>' % (X(i, off), Y(val), col))
        else:
            o.append('<rect x="%.1f" y="%.1f" width="10" height="10" fill="%s"/>' % (X(i, off) - 5, Y(val) - 5, col))
    o.append('<text x="%.1f" y="%d" text-anchor="middle" fill="#333">%s</text>' % (X(i, 0), T + PH + 18, labels[i]))
o.append('<text x="%.1f" y="%d" text-anchor="middle" fill="#333">distance from the nearest eRASS1 cluster (R500)</text>' % (L + PW / 2, T + PH + 38))
o.append('<text x="18" y="%.1f" transform="rotate(-90 18 %.1f)" text-anchor="middle" fill="#333">excess RM variance (rad² m⁻⁴)</text>' % (T + PH / 2, T + PH / 2))
o.append('<circle cx="%d" cy="22" r="5" fill="#2f5bb7"/><text x="%d" y="26" fill="#333">SPICE-RACS, 0.9 GHz (Entry 176)</text>' % (L + 300, L + 310))
o.append('<rect x="%d" y="9" width="10" height="10" fill="#d08a1e"/><text x="%d" y="18" fill="#333">NVSS, 1.4 GHz (this entry)</text>' % (L + 5, L + 20))
o.append('</svg>')
open(os.path.join(HERE, 'figure1.svg'), 'w').write('\n'.join(o))
print('written figure1.svg')
