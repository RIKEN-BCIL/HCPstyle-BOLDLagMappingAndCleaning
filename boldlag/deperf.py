"""drDeperf_hcp_seed / drDeperf_longTR -- removal of the perfusion lag structure.

For every lag value present in the (resliced) lag map, the sLFO that tracked that
lag (``Seeds.mat`` of the lag-mapping folder, time shifted accordingly) is
regressed out of the high-passed original data within that lag region; the
temporal mean is restored.  Output ``<basename>_dep.nii.gz`` and ``sLFO.mat``.
"""
import os
import numpy as np
import nibabel as nib
from scipy.io import savemat, loadmat
from scipy.signal import resample_poly
from .filters import bptf, hp_sigma, regfilt
from . import progress


def _round(x):
    """MATLAB round (half away from zero)."""
    return np.sign(x) * np.floor(np.abs(x) + 0.5)


def load_seeds(lagdir):
    p = os.path.join(lagdir, 'Seeds.npy')
    if os.path.exists(p):
        return np.load(p)
    try:
        return loadmat(os.path.join(lagdir, 'Seeds.mat'))['Seeds']
    except NotImplementedError:                # MATLAB v7.3
        import h5py
        with h5py.File(os.path.join(lagdir, 'Seeds.mat'), 'r') as h:
            return h['Seeds'][()].T


def shifted_seeds(Seeds, MaxLag):
    """``Motodata`` of drDeperf: column p (lag L = -MaxLag..MaxLag) is the seed used at
    that lag, shifted in time so that it is aligned with the voxels of that lag."""
    n = Seeds.shape[1]
    cols = []
    for p in range(1, MaxLag + 1):
        t = Seeds[:, n - p]
        cols.append(np.r_[np.zeros(MaxLag - p + 1), t[:len(t) - MaxLag - 1 + p]])
    cols.append(Seeds[:, MaxLag])
    for p in range(1, MaxLag + 1):
        t = Seeds[:, MaxLag - p]
        cols.append(np.r_[t[p:], np.zeros(p)])
    return np.stack(cols, 1)


def deperf(orig_vols, lag_nii, TR, section, Nruns, lagdir=None, reso=None, outdir='.', hp_hz=0.008, span=(0.0, 1.0)):
    """Deperfusion one run.

    orig_vols  original 4D run (full resolution)
    lag_nii    lag map resliced onto the grid of ``orig_vols`` (``rLagMap.nii``; NaN or 10000 outside)
    TR         repetition time (s)
    section    1-based run number within the concatenation used for lag mapping
    Nruns      number of concatenated runs
    lagdir     folder containing ``Seeds.mat`` (default: folder of ``lag_nii``)
    reso       tracking step used for lag mapping (None = TR)
    Returns the path of ``<outdir>/<basename>_dep.nii.gz``.
    """
    lagdir = lagdir or os.path.dirname(os.path.abspath(lag_nii))
    step = float(reso) if reso else float(TR)
    limg = nib.load(lag_nii)
    Lag = limg.get_fdata()
    Mask = 10000.0 * (Lag < 100)
    nib.save(nib.Nifti1Image(Mask.astype(np.float32), limg.affine), os.path.join(lagdir, 'Mask.nii'))
    Lag[Lag > 100] = np.nanmin(np.where(Lag > 100, np.nan, Lag))    # NaN stays NaN (no region)
    Lag = _round(Lag / step)

    Seeds = load_seeds(lagdir)
    MaxLag = (Seeds.shape[1] - 1) // 2
    Moto = shifted_seeds(Seeds, MaxLag)
    if reso and abs(step - TR) > 1e-9:
        Moto = resample_poly(Moto, int(round(reso * 100)), int(round(TR * 100)), axis=0, window=('kaiser', 5.0))
    img = nib.load(orig_vols)
    Nvols = img.shape[3]
    if Moto.shape[0] < Nruns * Nvols:
        raise ValueError(f'Seeds length {Moto.shape[0]} < {Nruns} x {Nvols} volumes')
    Moto = Moto[(section - 1) * Nvols:section * Nvols]
    savemat(os.path.join(outdir, 'sLFO.mat'), {'Motodata': Moto})

    sp = progress.Span('deperf: filtering', *span)
    sp(0.05)
    print('Filtering...', flush=True)
    Y = np.asanyarray(img.dataobj).astype(np.float32)
    Tmean = Y.mean(3)
    low = bptf(Y, hp_sigma(hp_hz, TR), -1)
    del Y
    out = np.repeat(Tmean[..., None], Nvols, axis=3).astype(np.float32)
    print('Cleaning images...', flush=True)
    LL = np.arange(-MaxLag, MaxLag + 1)
    sp.stage = 'deperf: regressing lag regions'
    for p, L in enumerate(LL):
        sp(0.4 + 0.55 * p / len(LL))
        m = (Lag == L) & (Tmean != 0)
        if not m.any():
            continue
        out[m] += regfilt(low[m], Moto[:, p], use_automask=False)
    base = os.path.basename(orig_vols)
    base = base[:-7] if base.endswith('.nii.gz') else os.path.splitext(base)[0]
    dst = os.path.join(outdir, base + '_dep.nii.gz')
    hdr = img.header.copy(); hdr.set_data_dtype(np.float32)
    nib.save(nib.Nifti1Image(out, img.affine, hdr), dst)
    sp(1.0)
    return dst
