"""Quick-look montage of a lag map (matplotlib, optional dependency)."""
import numpy as np
import nibabel as nib


def lagmap_montage(lagmap_nii, underlay_nii=None, out_png=None, lim=4.0, nslices=12, title=None):
    """Axial montage of ``lagmap_nii`` (jet, +-lim s) over ``underlay_nii`` (e.g. Tmean.nii).
    Returns the PNG path (default: next to the lag map)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    L = nib.load(lagmap_nii).get_fdata()
    U = nib.load(underlay_nii).get_fdata() if underlay_nii else None
    zs = np.flatnonzero(np.isfinite(L).any((0, 1)))
    zs = np.linspace(zs[0], zs[-1], nslices + 2)[1:-1].round().astype(int)
    ncol = 6
    nrow = int(np.ceil(len(zs) / ncol))
    fig, ax = plt.subplots(nrow, ncol, figsize=(2.2 * ncol, 2.5 * nrow), squeeze=False)
    for k, a in enumerate(ax.ravel()):
        a.axis('off')
        if k >= len(zs):
            continue
        z = zs[k]
        if U is not None:
            a.imshow(U[:, :, z].T, cmap='gray', origin='lower')
        im = a.imshow(np.ma.masked_invalid(L[:, :, z].T), cmap='jet', vmin=-lim, vmax=lim, origin='lower', alpha=0.9)
        a.set_title(f'z={z}', fontsize=8)
    fig.colorbar(im, ax=ax, shrink=0.6, label='lag (s)')
    fig.suptitle(title or lagmap_nii, fontsize=9)
    out_png = out_png or lagmap_nii.replace('.nii', '') + '_montage.png'
    fig.savefig(out_png, dpi=100)
    plt.close(fig)
    return out_png


def _lag_step(lagdir):
    """Tracking step (s) of a lag-map folder (params.json, else 1)."""
    import json, os
    p = os.path.join(lagdir, 'params.json')
    if os.path.exists(p):
        with open(p) as f:
            return float(json.load(f).get('step', 1.0))
    return 1.0


def lag_structure_plot(lagdir, out_png=None, t0=0, n=300, lim=4.0, shifted=True, title=None):
    """Rainbow plot of the sLFO time courses of a lag-map folder (``Seeds``), one line per
    lag coloured like the lag map (jet, +-lim s).  ``shifted=True`` plots the seeds
    shifted in time to their lag (what deperfusioning regresses out, cf. drDeperf's
    ``Motodata``); ``False`` plots ``Seeds.mat`` as stored (drPlotRainbow).
    ``t0``/``n`` select a window of samples.  Returns the PNG path."""
    import os
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import cm, colors
    from .deperf import load_seeds, shifted_seeds
    S = load_seeds(lagdir)
    MaxLag = (S.shape[1] - 1) // 2
    step = _lag_step(lagdir)
    if shifted:
        M = shifted_seeds(S, MaxLag)
        lags = np.arange(-MaxLag, MaxLag + 1) * step             # column p <-> lag -MaxLag+p
    else:
        M = S
        lags = (MaxLag - np.arange(S.shape[1])) * step           # Seeds column q <-> lag MaxLag-q
    t = np.arange(M.shape[0]) * step
    sl = slice(int(t0), int(t0) + int(n))
    norm = colors.Normalize(-lim, lim)
    fig, ax = plt.subplots(figsize=(14, 3.6))
    ax.set_facecolor((0.5, 0.5, 0.5))
    for k in np.argsort(np.abs(lags))[::-1]:                    # draw large |lag| first, lag 0 on top
        ax.plot(t[sl], M[sl, k], color=cm.jet(norm(lags[k])), lw=1.2)
    ax.set_xlim(t[sl][0], t[sl][-1])
    ax.set_xlabel('time (s)')
    ax.set_ylabel('sLFO (a.u.)')
    fig.colorbar(cm.ScalarMappable(norm=norm, cmap='jet'), ax=ax, pad=0.01, label='lag (s)')
    ax.set_title(title or f"{'shifted sLFO (regressors)' if shifted else 'Seeds.mat'} - {os.path.basename(lagdir)}", fontsize=9)
    fig.tight_layout()
    out_png = out_png or os.path.join(lagdir, 'sLFO_shifted.png' if shifted else 'Seeds.png')
    fig.savefig(out_png, dpi=100)
    plt.close(fig)
    return out_png


def slice_image(lagmap_nii, underlay_nii=None, axis=2, index=None, lim=4.0, out_png=None):
    """One slice of the lag map over the underlay (axis 0/1/2 = sagittal/coronal/axial).
    Returns (png path, number of slices along ``axis``)."""
    import os
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    L = nib.load(lagmap_nii).get_fdata()
    U = nib.load(underlay_nii).get_fdata() if underlay_nii else None
    nsl = L.shape[axis]
    index = nsl // 2 if index is None else int(np.clip(index, 0, nsl - 1))
    take = lambda a: np.take(a, index, axis=axis).T
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.axis('off')
    if U is not None:
        ax.imshow(take(U), cmap='gray', origin='lower')
    im = ax.imshow(np.ma.masked_invalid(take(L)), cmap='jet', vmin=-lim, vmax=lim, origin='lower', alpha=0.9)
    fig.colorbar(im, ax=ax, shrink=0.7, label='lag (s)')
    ax.set_title(f"{['sagittal', 'coronal', 'axial'][axis]} {index}/{nsl - 1}", fontsize=9)
    out_png = out_png or os.path.join(os.path.dirname(lagmap_nii), '_slice.png')
    fig.savefig(out_png, dpi=100, bbox_inches='tight')
    plt.close(fig)
    return out_png, nsl
