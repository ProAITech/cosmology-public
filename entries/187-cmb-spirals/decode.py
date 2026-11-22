"""Entry 187: decode the browser read-out of encode.js into the sample files stack.py reads.
Input: samples/readout.txt, lines "<line>:<count>:<checksum>:<base64>" exactly as the page showed them (encode.js header).
Checks every line's count and checksum and stops on any mismatch.  Writes samples/P.txt, P_large.txt, P_dense.txt, A_2MRS.txt:
Galactic unit vector (HEALPix RING Nside 4096 pixel centre) and redshift, one galaxy per line.
Sources: Nair & Abraham (2010), VizieR J/ApJS/186/427/table2; 2MRS, Huchra et al. (2012), VizieR J/ApJS/199/26/table3.
Run from the repository root: .venv/bin/python -I entries/187-cmb-spirals/decode.py"""
import base64, sys
import numpy as np, healpy as hp
from pathlib import Path
S = Path(__file__).resolve().parent / 'samples'
recs = []
for ln in open(S / 'readout.txt').read().split():
    i, n, cs, b = ln.split(':'); raw = base64.b64decode(b); assert len(raw) == 6 * int(n), ('length', i)
    c = 0
    for k in range(int(n)):
        V = int.from_bytes(raw[6 * k:6 * k + 6], 'big'); pix, z5, fl = V >> 20, (V >> 4) & 0xFFFF, V & 15
        c = (c + pix % 999983 + 7 * z5 + 13 * fl) % 1000000007; recs.append((pix, z5, fl))
    if c != int(cs): sys.exit(f'checksum mismatch on line {i}')
R = np.array(recs); v = np.array(hp.pix2vec(4096, R[:, 0])).T; z = R[:, 1] / 1e5; fl = R[:, 2]
A = (fl & 4) > 0
fmt = '%.7f %.7f %.7f %.5f'
for name, m in (('P', ~A), ('P_large', ~A & ((fl & 1) > 0)), ('P_dense', ~A & ((fl & 2) > 0)), ('A_2MRS', A)):
    np.savetxt(S / f'{name}.txt', np.c_[v[m], z[m]], fmt=fmt, header=f'Entry 187 sample {name}: Galactic unit vector x y z, redshift (decode.py)')
    print(name, int(m.sum()))
