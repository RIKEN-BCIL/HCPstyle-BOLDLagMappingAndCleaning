"""Einsteining -- the whole pipeline on one subject (HCP-style directory layout):

1. per run: (optional) 2x down-sampling (``subsamp2offc``), DVARS spike detection,
   regression of 6 motion parameters + their forward/backward differences + FD +
   spike regressors (``fsl_regfilt`` equivalent)            -> ``mreg_z<run>.nii.gz``
2. concatenation with high-pass and end taper (:func:`merge4d`) -> ``REST<n>run.nii.gz``
3. lag mapping (:func:`lag4d`)                                  -> ``Lag_fix_..._cat<n>/``
4. reslicing of ``LagMap.nii`` onto the original grid (nearest) -> ``rLagMap.nii``
5. per run: deperfusioning (:func:`deperf`)                     -> ``<run>_dep.nii.gz``
6. ``<Results>/<run>_dep/`` folders with symbolic links (HCP layout)

Everything is written to ``<Results>/Lag_concat_scrub`` (an existing one is renamed
``old_Lag_concat_scrub``).  Human (cerebral seed mask) and monkey (whole brain,
0.5 s tracking step) data are handled by the same code through ``seed_mask``,
``reso`` and ``mask_pct``.
"""
import os, re, glob, json, shutil
import numpy as np
import nibabel as nib
from .filters import subsamp2offc, subsamp2offc_affine, regfilt
from .merge4d import merge4d
from .lag4d import lag4d
from .deperf import deperf
from .spm import reslice


def run_basename(f):
    b = os.path.basename(f)
    return b[:-7] if b.endswith('.nii.gz') else os.path.splitext(b)[0]


def framewise_displacement(movement_txt):
    """drFD: HCP ``Movement_Regressors.txt`` -> (FD (N,), first 12 columns).
    Rotations are scaled by 50 (radius 50 mm), as in the MATLAB code."""
    MP = np.loadtxt(movement_txt)[:, :12]
    tr = MP[:, 3:6] * 180 / np.pi * (2 * 50 * np.pi / 360)
    FD = np.r_[0, np.abs(np.diff(np.c_[MP[:, :3], tr], axis=0)).sum(1)]
    return FD, MP


def dvars_spikes(Y, thr=1.5):
    """DVARS of a 4D array and the spike volumes (0-based).

    Fixed version of the Einsteining criterion (GitHub issue #2): DVARS is the RMS of
    the frame-to-frame difference of the *raw* intensity within the brain (voxels
    that are never zero), expressed in percent of the mean brain signal; a volume
    (and the one before it) is a spike when DVARS > ``thr`` x median(DVARS)."""
    brain = (Y != 0).all(3) & (Y.mean(3) > 0)
    Ym = Y[brain].astype(np.float64)
    dvars = np.r_[0.0, np.sqrt(np.mean(np.diff(Ym, axis=1) ** 2, axis=0))] / Ym.mean() * 100
    spike = np.flatnonzero(dvars > thr * np.median(dvars))
    spike = np.unique(np.r_[spike, spike - 1])
    return dvars, spike[spike >= 0]


def scrub_run(run, outdir, downsample=True, movement_txt=None, spike_thr=1.5):
    """Step 1 of the pipeline for one run; returns the path of ``mreg_z<run>.nii.gz``.
    All regressors (6 motion + forward/backward differences + FD + spikes) are removed
    (the MATLAB versions before the issue #2 fix silently skipped the spike columns)."""
    fn = run_basename(run)
    prefix = 'z' if downsample else 'H'
    img = nib.load(run)
    Y = np.asanyarray(img.dataobj).astype(np.float32)
    aff = img.affine
    if downsample:
        print(f'Reducing resolution {fn}', flush=True)
        Y = subsamp2offc(Y)
        aff = subsamp2offc_affine(img.affine)
        hdr = img.header.copy(); hdr.set_data_dtype(np.float32)
        hdr.set_zooms(tuple(np.array(hdr.get_zooms()[:3]) * 2) + (hdr.get_zooms()[3],))
        nib.save(nib.Nifti1Image(Y, aff, hdr), os.path.join(outdir, f'z{fn}.nii'))
    print(f'Motion scrubbing {fn}', flush=True)
    dvars, spike = dvars_spikes(Y, spike_thr)
    np.savetxt(os.path.join(outdir, f'dvars_{fn}.txt'), dvars, fmt='%.4f')
    print(f'  DVARS median {np.median(dvars):.2f} %, max {dvars.max():.2f} %, {len(spike)} spike volumes', flush=True)
    N = Y.shape[3]
    Sreg = np.zeros((N, len(spike)))
    Sreg[spike, np.arange(len(spike))] = 1
    movement_txt = movement_txt or os.path.join(os.path.dirname(run), 'Movement_Regressors.txt')
    FD, MP = framewise_displacement(movement_txt)
    DP = np.diff(MP[:, :6], axis=0)
    X = np.c_[MP[:, :6], np.r_[np.zeros((1, 6)), DP], np.r_[DP, np.zeros((1, 6))], FD, Sreg]
    np.savetxt(os.path.join(outdir, f'regout_{fn}.txt'), X, fmt='%.6f', delimiter='\t')
    Y = regfilt(Y, X)
    out = os.path.join(outdir, f'mreg_{prefix}{fn}.nii.gz')
    hdr = nib.Nifti1Header(); hdr.set_data_dtype(np.float32)
    im = nib.Nifti1Image(Y, aff, hdr)
    im.header.set_zooms(tuple(np.abs(np.diag(aff)[:3])) + (img.header.get_zooms()[3],))
    nib.save(im, out)
    return out


def find_runs(subject_dir, pattern=r'_REST[12]', exclude=r'dep|@', nvols=None):
    """Runs ``<subject_dir>/MNINonLinear/Results/<d>/<d>.nii(.gz)`` whose folder name
    matches ``pattern`` and not ``exclude`` (and has ``nvols`` volumes if given)."""
    res = os.path.join(subject_dir, 'MNINonLinear', 'Results')
    runs = []
    for d in sorted(os.listdir(res)):
        if not os.path.isdir(os.path.join(res, d)) or not re.search(pattern, d) or re.search(exclude, d):
            continue
        for ext in ('.nii.gz', '.nii'):
            f = os.path.join(res, d, d + ext)
            if os.path.exists(f):
                if nvols is None or nib.load(f).shape[3] == nvols:
                    runs.append(f)
                break
    if not runs:
        raise FileNotFoundError('no runs found in ' + res)
    return runs


def reslice_lagmap(lagmap_nii, ref_nii, out_nii):
    """drReslice_Lag: NaN -> 10000, nearest-neighbour reslice onto ``ref_nii``'s grid."""
    li = nib.load(lagmap_nii)
    L = li.get_fdata()
    L[np.isnan(L)] = 10000
    ref = nib.load(ref_nii)
    r = reslice(L, li.affine, ref.affine, ref.shape[:3], order=0)      # NaN outside the source FOV, like SPM
    nib.save(nib.Nifti1Image(r.astype(np.float32), ref.affine), out_nii)
    return out_nii


def einsteining(runs, TR, PosiMax, THR=0.2, FIXED=1, Sm=8, only_lag=False, downsample=True,
                reso=None, seed_mask='hcp', mask_pct=10, lp_hz=None, workname='Lag_concat_scrub',
                results_dir=None, spike_thr=1.5):
    """Run the whole pipeline.  ``runs`` = list of 4D run files (HCP layout:
    ``.../MNINonLinear/Results/<run>/<run>.nii.gz`` with ``Movement_Regressors.txt``
    and ``*SBRef.nii.gz`` next to them).  Other arguments as in :func:`lag4d`.
    Returns the lag-map folder."""
    runs = [os.path.abspath(r) for r in runs]
    results_dir = results_dir or os.path.dirname(os.path.dirname(runs[0]))
    ref = glob.glob(os.path.join(os.path.dirname(runs[0]), '*SBRef.nii*'))
    if not ref:
        raise FileNotFoundError('no SBRef next to ' + runs[0])
    ref = ref[0]
    wd = os.path.join(results_dir, workname)
    if os.path.exists(wd):
        old = os.path.join(results_dir, 'old_' + workname)
        if os.path.exists(old):
            shutil.rmtree(old)
        os.rename(wd, old)
    os.makedirs(wd)
    with open(os.path.join(wd, 'Runs.json'), 'w') as f:
        json.dump(dict(Runs=runs, Ref=ref, TR=TR, PosiMax=PosiMax, THR=THR, FIXED=FIXED, Sm=Sm, reso=reso,
                       seed_mask=seed_mask, mask_pct=mask_pct, downsample=downsample, spike_thr=spike_thr), f, indent=1)
    n = len(runs)
    z = [scrub_run(r, wd, downsample, spike_thr=spike_thr) for r in runs]
    merged = merge4d(('' if downsample else 'hres') + f'REST{n}run', TR, z, wd)
    lagdir = lag4d(f'cat{n}', TR, merged, PosiMax, THR, FIXED, Sm, reso=reso, seed_mask=seed_mask,
                   mask_pct=mask_pct, lp_hz=lp_hz, cwd=wd)
    if only_lag:
        return lagdir
    rlag = reslice_lagmap(os.path.join(lagdir, 'LagMap.nii'), ref, os.path.join(wd, 'rLagMap.nii'))
    shutil.copy(rlag, lagdir)
    for r, run in enumerate(runs, 1):
        print(f'Deperfusioning {run_basename(run)}', flush=True)
        deperf(run, rlag, TR, r, n, lagdir=lagdir, reso=reso, outdir=wd)
    for run in runs:
        b = run_basename(run)
        d = os.path.join(results_dir, b + '_dep')
        os.makedirs(d, exist_ok=True)
        for link, target in [('Movement_Regressors.txt', f'../{b}/Movement_Regressors.txt'),
                             (f'{b}_dep_SBRef.nii.gz', f'../{b}/{b}_SBRef.nii.gz'),
                             (f'{b}_dep.nii.gz', f'../{workname}/{b}_dep.nii.gz')]:
            p = os.path.join(d, link)
            if os.path.lexists(p):
                os.unlink(p)
            os.symlink(target, p)
    return lagdir
