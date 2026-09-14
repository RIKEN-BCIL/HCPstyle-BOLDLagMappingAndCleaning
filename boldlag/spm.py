"""SPM12 operations used by the MATLAB scripts, re-implemented with scipy.

* :func:`smooth`  -- ``spm_smooth`` (Gaussian FWHM in mm, implicit NaN mask, ``im=1``)
* :func:`reslice` -- ``spm.spatial.coreg.write`` (interp 0 = nearest, 1 = trilinear; wrap 0, mask 0)
"""
import numpy as np
from scipy.special import erf
from scipy.ndimage import convolve1d, map_coordinates


def spm_kernel(fwhm_vox):
    """1-D kernel of spm_smooth (integrated Gaussian, see spm_smooth.m > smooth1)."""
    s1 = fwhm_vox / np.sqrt(8 * np.log(2))
    n = int(round(6 * s1))
    x = np.arange(-n, n + 1, dtype=float)
    s = s1 ** 2 + np.finfo(float).eps
    w1, w2, w3 = 0.5 * np.sqrt(2 / s), -0.5 / s, np.sqrt(s / 2 / np.pi)
    k = 0.5 * (erf(w1 * (x + 1)) * (x + 1) + erf(w1 * (x - 1)) * (x - 1) - 2 * erf(w1 * x) * x) \
        + w3 * (np.exp(w2 * (x + 1) ** 2) + np.exp(w2 * (x - 1) ** 2) - 2 * np.exp(w2 * x ** 2))
    k[k < 0] = 0
    return k / k.sum()


def smooth(vol, fwhm_mm, vox, chunk=200):
    """Smooth a 3D or 4D (x,y,z[,t]) array like SPM's batch ``spatial.smooth`` with
    ``im = 1``: NaNs are treated as 0 during the separable convolution and restored
    afterwards; at the volume boundary the truncated kernel is renormalised (as
    spm_conv_vol does).  ``vox`` = voxel size (mm) per axis.  Verified against SPM12
    (max difference 1e-5)."""
    fwhm = np.broadcast_to(np.asarray(fwhm_mm, float), (3,))
    kern = [spm_kernel(fwhm[a] / vox[a]) for a in range(3)]
    v = np.asarray(vol)
    is3d = v.ndim == 3
    if is3d:
        v = v[..., None]
    # edge normaliser: separable convolution of an all-ones volume (zero padded)
    norm = 1.0
    for a in range(3):
        n1 = convolve1d(np.ones(v.shape[a]), kern[a], mode='constant', cval=0.0)
        shape = [1, 1, 1, 1]; shape[a] = v.shape[a]
        norm = norm * n1.reshape(shape)
    out = np.empty(v.shape, np.float32)
    for t in range(0, v.shape[3], chunk):
        y = v[..., t:t + chunk].astype(np.float64)
        nanm = ~np.isfinite(y)
        y[nanm] = 0.0
        for a in range(3):
            y = convolve1d(y, kern[a], axis=a, mode='constant', cval=0.0)
        y /= norm
        y[nanm] = np.nan
        out[..., t:t + chunk] = y
    return out[..., 0] if is3d else out


def reslice(src, src_affine, ref_affine, ref_shape, order=0):
    """Reslice ``src`` (3D array with ``src_affine``) onto the grid defined by
    ``ref_affine``/``ref_shape``.  ``order`` 0 = nearest neighbour, 1 = trilinear.
    Voxels falling outside the source FOV become NaN (SPM ``mask=0``, ``wrap=0``)."""
    ref_shape = tuple(int(s) for s in ref_shape[:3])
    ijk = np.indices(ref_shape, dtype=float).reshape(3, -1)
    M = np.linalg.inv(src_affine) @ ref_affine
    xyz = M[:3, :3] @ ijk + M[:3, 3:4]
    n = np.asarray(src.shape[:3], float)[:, None]
    outside = ((xyz < 0) | (xyz > n - 1)).any(0)         # SPM: NaN if the (unrounded) coordinate leaves the FOV
    if order == 0:
        xyz = np.floor(xyz + 0.5)
    out = map_coordinates(np.asarray(src, np.float64), xyz, order=order, mode='nearest')
    out[outside] = np.nan
    return out.reshape(ref_shape)
