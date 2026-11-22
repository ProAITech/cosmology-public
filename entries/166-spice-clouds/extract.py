"""Entry 166, step 1: stream the SPICE-RACS DR2 table once and keep only what the frozen rule needs.

Source: data/spice-racs-dr2/spice-racs.dr2.fits.gz (CSIRO DAP csiro:64891, 4,968,301,262 bytes), one BINTABLE of 9,294,225 rows x 1030 bytes,
125 columns, read in chunks straight from the gzip stream (no decompressed copy on disk). Columns used: l, b, rm, rm_err, source_id, the four quality
flags of Entry 68 (snr_flag, leakage_flag, channel_flag, stokesI_fit_flag) and the survey identifiers catalog_name, telescope, dataref, type.
Nothing here looks at a spine, a window or a fit.  It writes (git-ignored) work/compact.npz and the small, tracked meta.json (value counts of the
categorical columns and of the four flags), which are what the survey selection and the flag convention are read from.
Run from the repository root:  .venv/bin/python -I entries/166-spice-clouds/extract.py
"""
import gzip, json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
PATH = os.path.join(ROOT, 'data', 'spice-racs-dr2', 'spice-racs.dr2.fits.gz')
WORK = os.path.join(HERE, 'work')
CHUNK_ROWS = 100_000


def read_header(f):
    cards = []
    while True:
        blk = f.read(2880)
        if len(blk) < 2880:
            raise SystemExit('short read in header')
        for i in range(0, 2880, 80):
            c = blk[i:i + 80].decode('ascii', errors='replace')
            cards.append(c)
            if c.startswith('END     '):
                return cards


def card(cards, key):
    for c in cards:
        if c[:8].strip() == key:
            return c[10:].split('/')[0].strip().strip("'").strip()
    return None


def table_dtype(h):
    nf = int(card(h, 'TFIELDS'))
    names, fmts = [], []
    for i in range(1, nf + 1):
        names.append(card(h, 'TTYPE%d' % i))
        fmts.append(card(h, 'TFORM%d' % i))
    fields, off = {}, 0
    for n, t in zip(names, fmts):
        if t.endswith('A'):
            code, size = 'S%d' % int(t[:-1] or 1), int(t[:-1] or 1)
        else:
            code, size = {'D': ('>f8', 8), 'E': ('>f4', 4), 'K': ('>i8', 8), 'J': ('>i4', 4), 'L': ('S1', 1)}[t]
        fields[n] = (code, off)
        off += size
    rowbytes = int(card(h, 'NAXIS1'))
    assert off == rowbytes, (off, rowbytes)
    return np.dtype({'names': list(fields), 'formats': [v[0] for v in fields.values()], 'offsets': [v[1] for v in fields.values()], 'itemsize': rowbytes}), rowbytes


def main():
    t0 = time.time()
    keep = {k: [] for k in ('l', 'b', 'rm', 'rm_err', 'source_id', 'catalog_name', 'flags')}
    cats = {k: {} for k in ('catalog_name', 'telescope', 'dataref', 'type', 'notes', 'flux_type', 'rm_method')}
    flagcounts = {k: {'T': 0, 'F': 0, 'other': 0} for k in ('snr_flag', 'leakage_flag', 'channel_flag', 'stokesI_fit_flag')}
    flagnames = list(flagcounts)
    with gzip.open(PATH, 'rb') as f:
        read_header(f)                       # primary HDU (no data)
        h = read_header(f)                   # the table header
        dt, rowbytes = table_dtype(h)
        nrows = int(card(h, 'NAXIS2'))
        done = 0
        while done < nrows:
            n = min(CHUNK_ROWS, nrows - done)
            buf = f.read(n * rowbytes)
            if len(buf) != n * rowbytes:
                raise SystemExit('short read in data at row %d' % done)
            a = np.frombuffer(buf, dtype=dt)
            for k in cats:
                u, c = np.unique(a[k], return_counts=True)
                for v, cc in zip(u, c):
                    s = v.decode('ascii', errors='replace').strip()
                    cats[k][s] = cats[k].get(s, 0) + int(cc)
            fl = np.zeros((n, 4), dtype=np.uint8)
            for j, k in enumerate(flagnames):
                x = a[k]
                fl[:, j] = np.where(x == b'T', 1, np.where(x == b'F', 0, 2))
                for sym, code in (('T', 1), ('F', 0), ('other', 2)):
                    flagcounts[k][sym] += int((fl[:, j] == code).sum())
            keep['l'].append(a['l'].astype('<f8'))
            keep['b'].append(a['b'].astype('<f8'))
            keep['rm'].append(a['rm'].astype('<f8'))
            keep['rm_err'].append(a['rm_err'].astype('<f8'))
            keep['source_id'].append(a['source_id'].copy())
            keep['catalog_name'].append(a['catalog_name'].copy())
            keep['flags'].append(fl)
            done += n
            if (done // CHUNK_ROWS) % 10 == 0:
                print('%d / %d rows, %.0f s' % (done, nrows, time.time() - t0), flush=True)
    os.makedirs(WORK, exist_ok=True)
    out = {k: np.concatenate(v) for k, v in keep.items()}
    np.savez(os.path.join(WORK, 'compact.npz'), **out)
    meta = {'rows': nrows, 'row_bytes': rowbytes, 'columns': len(dt.names), 'categorical_counts': cats, 'flag_counts': flagcounts,
            'seconds': round(time.time() - t0)}
    json.dump(meta, open(os.path.join(HERE, 'meta.json'), 'w'), indent=1)
    print('done', meta['seconds'], 's')


if __name__ == '__main__':
    main()
