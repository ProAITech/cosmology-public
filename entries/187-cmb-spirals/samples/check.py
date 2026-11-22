"""Entry 187: check the count and checksum of each read-out line in the files given (see encode.js)."""
import base64, sys
for f in sys.argv[1:]:
    for ln in open(f).read().split():
        i, n, cs, b = ln.split(':')
        try: raw = base64.b64decode(b, validate=True)
        except Exception as e: print(i, 'base64 error', e); continue
        c = 0
        for k in range(int(n)):
            V = int.from_bytes(raw[6*k:6*k+6], 'big'); pix, z5, fl = V >> 20, (V >> 4) & 0xFFFF, V & 15
            c = (c + pix % 999983 + 7*z5 + 13*fl) % 1000000007
        print(i, 'ok' if (len(raw) == 6*int(n) and c == int(cs)) else 'MISMATCH')
