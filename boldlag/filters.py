"""Pure-numpy re-implementations of the FSL operations used by the MATLAB scripts.

* :func:`bptf`         -- ``fslmaths -bptf hp lp``  (FSL >= 5.0.7 / niimath: mean removed)
* :func:`regfilt`      -- ``fsl_regfilt``           (OLS removal of design columns, mean restored)
* :func:`subsamp2offc` -- ``fslmaths -subsamp2offc`` (2x2x2 non-centred box down-sampling)
* :func:`robust_range` -- ``fslstats -r``            (histogram based 2 % / 98 % limits)
* :func:`thrp`         -- ``fslmaths -thrp N``
* :func:`taper_window` -- the 5 % linear end taper used by drMerge4D

All were checked against the FSL 5.0.9 / 6.0.5 / niimath binaries on real data
(differences at float32 rounding level only).  Time is always the *last* axis.
"""
import numpy as np
from scipy.ndimage import correlate1d


def hp_sigma(hz, TR):
    """FSL sigma (volumes) for a cut-off in Hz, as used in the MATLAB scripts: ``1/(hz*2.35*TR)``."""
    return 1.0 / (hz * 2.35 * TR)


def bptf(y, hp_sig, lp_sig=-1, chunk=20000):
    """``fslmaths -nan -bptf hp_sig lp_sig`` on ``y`` (..., T); returns float32.

    High-pass: for every time point a Gaussian-weighted (sigma ``hp_sig``,
    half-width ``int(3*hp_sig)``) local linear fit is evaluated and subtracted
    (FSL's non-linear high-pass); the temporal mean is then removed.
    Low-pass: Gaussian smoothing (sigma ``lp_sig``, half-width ``int(20*lp_sig)+2``)
    renormalised at the edges.  A sigma <= 0 skips that filter.  NaN -> 0.
    """
    y = np.asarray(y)
    shp = y.shape
    T = shp[-1]
    Y = y.reshape(-1, T)
    out = np.empty(Y.shape, np.float32)
    ones = np.ones(T)
    if hp_sig > 0:
        m = int(hp_sig * 3)
        dt = np.arange(-m, m + 1, dtype=float)
        w = np.exp(-0.5 * dt * dt / (hp_sig * hp_sig))
        # correlate1d(x, k)[t] = sum_j k[j] x[t + j - m]   (dt = j - m; zero padded = window clipped)
        N = correlate1d(ones, w, mode='constant')
        A = correlate1d(ones, w * dt, mode='constant')
        C = correlate1d(ones, w * dt * dt, mode='constant')
        den = C * N - A * A
        ok = den != 0
        den_safe = np.where(ok, den, 1.0)
    if lp_sig > 0:
        ml = int(lp_sig * 20) + 2
        dl = np.arange(-ml, ml + 1, dtype=float)
        wl = np.exp(-0.5 * dl * dl / (lp_sig * lp_sig))
        S = correlate1d(ones, wl, mode='constant')
    for i in range(0, Y.shape[0], chunk):
        x = np.nan_to_num(Y[i:i + chunk].astype(np.float64), nan=0.0)
        if hp_sig > 0:
            B = correlate1d(x, w, axis=1, mode='constant')
            D = correlate1d(x, w * dt, axis=1, mode='constant')
            c = np.where(ok, (B * C - A * D) / den_safe, 0.0)
            x = x - c
            x -= x.mean(1, keepdims=True)
        if lp_sig > 0:
            x = correlate1d(x, wl, axis=1, mode='constant') / S
        out[i:i + chunk] = x
    return out.reshape(shp)


def automask(mean_img, frac=0.01):
    """fsl_regfilt / MELODIC automatic mask: ``mean >= min + frac*(max-min)``."""
    m = np.asarray(mean_img, np.float64)
    return m >= m.min() + frac * (m.max() - m.min())


def regfilt(y, design, filt=None, mask=None, use_automask=True, chunk=20000):
    """``fsl_regfilt -i y -d design -f filt``: remove the ``filt`` columns (0-based,
    default all) of ``design`` from ``y`` (..., T) by OLS on the full de-meaned
    design, restoring the voxel mean.  Voxels outside ``mask`` (default: the
    fsl_regfilt automask of the temporal mean) are set to 0 like the binary does."""
    y = np.asarray(y)
    shp = y.shape
    T = shp[-1]
    X = np.atleast_2d(np.asarray(design, float))
    if X.shape[0] != T:
        X = X.T
    X = X - X.mean(0)
    filt = np.arange(X.shape[1]) if filt is None else np.asarray(filt)
    P = np.linalg.pinv(X)                          # (K, T)
    Xf, Pf = X[:, filt], P[filt]
    Y = y.reshape(-1, T)
    out = np.zeros(Y.shape, np.float32)
    if mask is None:
        mask = automask(Y.mean(1)) if use_automask else np.ones(Y.shape[0], bool)
    mask = np.asarray(mask, bool).reshape(-1)
    idx = np.flatnonzero(mask)
    for i in range(0, idx.size, chunk):
        ii = idx[i:i + chunk]
        x = Y[ii].astype(np.float64)
        mu = x.mean(1, keepdims=True)
        xd = x - mu
        out[ii] = xd - (xd @ Pf.T) @ Xf.T + mu
    return out.reshape(shp)


def robust_range(x, nbins=1000):
    """``fslstats -r``: 2 % / 98 % limits from a 1000-bin histogram over [min, max]
    (lower edge of the 2 % bin, upper edge of the 98 % bin)."""
    x = np.asarray(x, np.float32).ravel()
    mn, mx = np.float32(x.min()), np.float32(x.max())
    if mx == mn:
        return float(mn), float(mx)
    h, e = np.histogram(x, bins=nbins, range=(mn, mx))
    c = np.cumsum(h)
    lo = np.searchsorted(c, 0.02 * x.size)
    hi = np.searchsorted(c, 0.98 * x.size)
    w = (mx - mn) / np.float32(nbins)
    return float(mn + np.float32(lo) * w), float(mn + np.float32(hi + 1) * w)


def thrp(x, pct):
    """``fslmaths -thrp pct``: zero everything below ``pct`` % of the robust range."""
    rmin, rmax = robust_range(x)
    thr = rmin + pct / 100.0 * (rmax - rmin)
    x = np.asarray(x, np.float32)
    return np.where(x < thr, 0, x).astype(np.float32)


def subsamp2offc(vol):
    """``fslmaths -subsamp2offc``: average 2x2x2 blocks (voxels 2i, 2i+1) along the
    first three axes; a trailing odd voxel is kept as is.  3D or 4D input; float32 out."""
    v = np.asarray(vol, np.float64)
    for ax in range(3):
        n = v.shape[ax]
        sl = [slice(None)] * v.ndim
        sl[ax] = slice(0, 2 * (n // 2), 2)
        a = v[tuple(sl)]
        sl[ax] = slice(1, 2 * (n // 2), 2)
        b = v[tuple(sl)]
        s = (a + b) / 2
        if n % 2:
            sl[ax] = slice(n - 1, n)
            s = np.concatenate([s, v[tuple(sl)]], axis=ax)
        v = s
    return v.astype(np.float32)


def subsamp2offc_affine(aff):
    """Affine of the ``subsamp2offc`` output (voxel size x2, centre at the midpoint of old voxels 0/1)."""
    aff = np.asarray(aff, float).copy()
    R = aff[:3, :3]
    aff[:3, 3] = aff[:3, 3] + R @ np.array([0.5, 0.5, 0.5])
    aff[:3, :3] = R * 2
    return aff


def taper_window(N):
    """drMerge4D end taper ``[0:1/((N-1)*.05):1 ones(1,ceil(N*.9)) 1:-1/((N-1)*.05):0 0](1:N)``."""
    st = 1 / ((N - 1) * .05)
    nup = int(np.floor(1 / st + 1e-10))
    up = np.arange(nup + 1) * st
    down = 1 - np.arange(nup + 1) * st
    W = np.concatenate([up, np.ones(int(np.ceil(N * .9))), down, [0.0]])
    return W[:N]
