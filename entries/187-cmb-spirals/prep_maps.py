"""Entry 187: the Planck PR3 CMB maps as compact files for the browser stack (stack.html).
For each map: I_STOKES in micro-K_CMB, RING order, with the multipoles l <= LCUT removed (a_lm of I_STOKES_INP, the
inpainted map, up to LCUT); pixels outside the Planck common temperature mask set to -32768; others stored as int16 in
units of 0.1 micro-K (clipped at +-3276.7).  Files: work/<map>_l<LCUT>.i16 (little-endian, 12*2048^2 values).
Data (git-ignored, cited in README.md):
  data/planck-pr3-{smica,nilc,sevem,commander}/COM_CMB_IQU-<map>_2048_R3.00_full.fits   columns I_STOKES, I_STOKES_INP (K_CMB)
  data/planck-pr3-masks/COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits                     column TMASK (1 = use)
Run from the repository root: .venv/bin/python -I entries/187-cmb-spirals/prep_maps.py"""
import numpy as np, healpy as hp, json
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]; W = HERE / 'work'; W.mkdir(exist_ok=True)
mask = hp.read_map(ROOT / 'data/planck-pr3-masks/COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits', field=0) > 0.5
info = {'used_fraction': float(mask.mean())}
for name, cuts in (('smica', (5, 1)), ('nilc', (5,)), ('sevem', (5,)), ('commander', (5,))):
    f = ROOT / f'data/planck-pr3-{name}/COM_CMB_IQU-{name}_2048_R3.00_full.fits'
    T = hp.read_map(f, field=0) * 1e6; Ti = hp.read_map(f, field=5) * 1e6
    for L in cuts:
        low = hp.alm2map(hp.map2alm(Ti, lmax=L, iter=3), 2048, lmax=L)
        x = T - low
        q = np.clip(np.round(x * 10), -32767, 32767).astype('<i2'); q[~mask] = -32768
        q.tofile(W / f'{name}_l{L}.i16')
        info[f'{name}_l{L}'] = {'rms_used_uK': float(x[mask].std()), 'clipped': int((np.abs(x[mask]) > 3276.7).sum())}
        print(name, L, info[f'{name}_l{L}'], flush=True)
json.dump(info, open(W / 'maps.json', 'w'), indent=1)
