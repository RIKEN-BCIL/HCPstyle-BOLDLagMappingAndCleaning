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
