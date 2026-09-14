"""drMerge4D -- concatenate fMRI runs for lag mapping.

Each run: NaN -> 0, temporal mean stored, high-pass (``hp_hz``, FSL bptf semantics:
mean removed), 5 % linear taper at both ends; runs are concatenated and the average
of the run means is added back.
"""
import os
import numpy as np
import nibabel as nib
from .filters import bptf, hp_sigma, taper_window


def merge4d(name, TR, vols, outdir='.', hp_hz=0.008):
    """Write ``<outdir>/<name>.nii.gz`` from the list of 4D files ``vols``; returns its path."""
    segs, Tmean = [], 0.0
    for f in vols:
        img = nib.load(f)
        Y = np.nan_to_num(np.asanyarray(img.dataobj).astype(np.float32), nan=0.0)
        N = Y.shape[3]
        Tmean = Tmean + Y.mean(3, dtype=np.float64)
        print(f'...filtering {os.path.basename(f)} ({N} vols)', flush=True)
        Y = bptf(Y, hp_sigma(hp_hz, TR), -1) * taper_window(N).astype(np.float32)
        segs.append(Y)
    Tmean = (Tmean / len(vols)).astype(np.float32)
    All = np.concatenate(segs, axis=3) + Tmean[..., None]
    out = os.path.join(outdir, name + '.nii.gz')
    hdr = img.header.copy()
    hdr.set_data_dtype(np.float32)
    nib.save(nib.Nifti1Image(All, img.affine, hdr), out)
    return out
