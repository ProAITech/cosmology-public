"""Entry 187: the claim's signal as the Planck map would show it.  For a galaxy at redshift z the disc has radius
theta_d = 0.2 Mpc / D_A(z) (H0 = 70, Om = 0.3; 0.25 deg at z = 0.01, as Hansen et al. 2025 state).  Two profiles of
depth 1: 'top' = 1 inside 0.2 Mpc, 0 outside (the least the claim implies); 'fade' = 1 inside 0.2 Mpc, falling linearly
to 0 at 1.5 Mpc (the published profile 'remains low out to 1-2 Mpc').  Each is smoothed by the SMICA beam (5 arcmin FWHM)
and the Nside-2048 pixel (sigma 1.72/sqrt(12) arcmin), then averaged inside theta_d (f_disc) and at 2-3 Mpc (f_ann).
Writes predict.json: z grid, f_disc and f_ann for both profiles.  Run from the repository root."""
import numpy as np, json
from pathlib import Path
HERE = Path(__file__).resolve().parent
def DA(z, H0=70.0, Om=0.3):
    zz = np.linspace(0, z, 400); E = np.sqrt(Om * (1 + zz) ** 3 + 1 - Om)
    return 299792.458 / H0 * np.trapezoid(1 / E, zz) / (1 + z)
sig = np.hypot(5.0 / 2.35482, 1.72 / np.sqrt(12))      # arcmin
out = {'z': [], 'theta_d_arcmin': [], 'top': {'f_disc': [], 'f_ann': []}, 'fade': {'f_disc': [], 'f_ann': []}, 'sigma_arcmin': sig}
for z in np.round(np.arange(0.004, 0.1001, 0.002), 3):
    mpc = 1 / DA(z) * 180 / np.pi * 60                   # arcmin per Mpc
    td = 0.2 * mpc; r3 = 3.0 * mpc
    h = min(0.1, td / 25); L = r3 + 4 * sig; n = int(np.ceil(2 * L / h)); n += n % 2
    x = (np.arange(n) - n / 2 + 0.5) * h; X, Y = np.meshgrid(x, x); R = np.hypot(X, Y)
    k = np.fft.fftfreq(n, d=h); KX, KY = np.meshgrid(k, k); G = np.exp(-2 * (np.pi * sig) ** 2 * (KX ** 2 + KY ** 2))
    disc, ann = R < td, (R >= 2 * mpc) & (R < 3 * mpc)
    for name, prof in (('top', (R < td).astype(float)), ('fade', np.where(R < td, 1.0, np.clip((1.5 * mpc - R) / (1.3 * mpc), 0, 1)))):
        s = np.real(np.fft.ifft2(np.fft.fft2(prof) * G))
        out[name]['f_disc'].append(round(float(s[disc].mean()), 4)); out[name]['f_ann'].append(round(float(s[ann].mean()), 5))
    out['z'].append(float(z)); out['theta_d_arcmin'].append(round(td, 3))
json.dump(out, open(HERE / 'predict.json', 'w'))
for i in range(0, len(out['z']), 6):
    print(out['z'][i], out['theta_d_arcmin'][i], out['top']['f_disc'][i], out['fade']['f_disc'][i], out['fade']['f_ann'][i])
