"""drLag4Drev7 -- lag mapping of 4D BOLD data (T. Aso, RIKEN-BDR), Python port.

One function covers the three MATLAB variants:

* ``drLag4Drev7``        (HCP, tracking step = TR, cerebral seed mask): ``reso=None, seed_mask=<mask>``
* ``drLag4Drev7_longTR`` (resampled to 1 s):                            ``reso=1.0``
* ``drLag4Drev7_monkey`` (resampled to 0.5 s, whole-brain seed):       ``reso=0.5``

Outputs (in ``<cwd>/Lag_{fix|rec}_{MaxLag}{TR|sec}_thr{10*THR}_sm{Sm}_{name}/``):
``LagOrig.nii`` (raw lag, s), ``LagMap.nii`` (eroded/filled, masked), ``e1LagOrig.nii``,
``MaxR.nii`` (peak cross-correlation), ``Seeds.mat`` (sLFO time courses used at each
lag, MATLAB compatible), ``RawSeed.mat``, ``InitSeed.mat``, ``Ysd.mat``, ``params.json``.
Intermediate files ``SD.nii``, ``Tmean.nii``, ``Mask.nii``, ``nanmask.nii`` and
``sm{Sm}_{MaxLag}{unit}.nii`` are written to ``cwd`` and reused if present.
"""
import os, json
import numpy as np
import nibabel as nib
from scipy.io import savemat
from scipy.signal import resample_poly

from .filters import bptf, hp_sigma, thrp
from .spm import smooth, reslice
from . import progress

HCP_SEED_MASK = os.path.join(os.path.dirname(__file__), 'data', 'BrainMask_lag_subsamp2offc.nii')


def _save(path, data, affine):
    nib.save(nib.Nifti1Image(np.asarray(data, np.float32), affine), path)


def _load4d(path, rng=None):
    img = nib.load(path)
    if rng is None:
        Y = np.asanyarray(img.dataobj).astype(np.float32)
    else:
        rng = np.asarray(rng)
        st = np.unique(np.diff(rng))
        if rng.size == 1 or (st.size == 1 and st[0] > 0):        # regular -> slice (partial read)
            step = int(st[0]) if st.size else 1
            Y = np.asanyarray(img.dataobj[..., int(rng[0]):int(rng[-1]) + 1:step]).astype(np.float32)
        else:
            Y = np.asanyarray(img.dataobj).astype(np.float32)[..., rng]
    return Y, img


def parse_range(rng):
    """MATLAB-style ``'a:b'`` / ``'a:s:b'`` (1-based, inclusive) or a python slice/array -> index array or None."""
    if rng is None or rng == '' or rng == []:
        return None
    if isinstance(rng, str):
        p = [int(float(v)) for v in rng.split(':')]
        if len(p) == 2:
            return np.arange(p[0] - 1, p[1])
        if len(p) == 3:
            return np.arange(p[0] - 1, p[2], p[1])
        raise ValueError('range must be a:b or a:s:b (MATLAB 1-based)')
    return rng


def prepare(vols, TR, MaxLag_sec, Sm, mask_pct=10, lp_hz=None, hp_hz=0.008, out_sm=None, cwd='.'):
    """Data preparation stage of drLag4D: SD, Tmean, Mask (``-thrp mask_pct``),
    percent signal change within the mask (NaN outside), SPM smoothing (``Sm`` mm
    FWHM, implicit mask), band-pass ``hp_hz`` .. ``lp_hz`` (default 0.9/MaxLag_sec).
    Returns the path of the filtered 4D file."""
    ULfreq = lp_hz if lp_hz else 1.0 / MaxLag_sec * .9
    Y, img = _load4d(vols)
    aff = img.affine
    Mean = Y.mean(3)
    _save(os.path.join(cwd, 'SD.nii'), Y.std(3, ddof=0), aff)
    _save(os.path.join(cwd, 'Tmean.nii'), Mean, aff)
    Mask = thrp(Mean, mask_pct)
    _save(os.path.join(cwd, 'Mask.nii'), Mask, aff)
    M = np.where(Mask == 0, np.nan, 1.0).astype(np.float32)
    _save(os.path.join(cwd, 'nanmask.nii'), M, aff)
    Y = 100 * Y / (Mean * M)[..., None]
    if Sm > 0:
        print('Smoothing...', flush=True)
        Y = smooth(Y, Sm, img.header.get_zooms()[:3])
    print('Filtering...', flush=True)
    Y = bptf(Y, hp_sigma(hp_hz, TR), hp_sigma(ULfreq, TR))
    _save(out_sm, Y, aff)
    return out_sm


def _corr_peak(seed, Y, Lim):
    """Cross-correlogram peak of ``seed`` (T,) against every column of ``Y`` (T,V)
    for shifts Lim..-Lim.  Returns (R, I) with I 1-based like MATLAB ``max``."""
    X = seed[Lim:len(seed) - Lim]
    xx = np.sqrt((X * X).sum())
    CC = []
    for Sft in range(Lim, -Lim - 1, -1):
        YY = Y[Lim + Sft:len(Y) - Lim + Sft]
        with np.errstate(invalid='ignore', divide='ignore'):
            CC.append((X @ YY) / (xx * np.sqrt((YY * YY).sum(0))))
    CC = np.stack(CC, -1)
    allnan = np.isnan(CC).all(1)
    CCf = np.where(np.isnan(CC), -np.inf, CC)
    I = CCf.argmax(1) + 1
    R = CCf.max(1)
    R[allnan] = np.nan
    I[allnan] = 1
    return R, I


def track(Y, Bmask, roi, THR, limit, FIXED=1, on_step=None, span=None):
    """Recursive / fixed-seed lag tracking (core of drLag4D).

    Y      (T, V) float64 filtered data (already amplitude-gated, resampled to the tracking step)
    Bmask  (V,) bool   voxels eligible for a lag value
    roi    (V,) float  weights (NaN outside) of the region giving the initial global signal
    Returns dict(Lag (V,) in steps, maxR, Seeds (T, 2*limit+1), RawSeed, InitSeed, Ysd).
    """
    T, V = Y.shape
    Y0 = Y                                   # un-normalised data (percent signal), kept for the region means
    with np.errstate(invalid='ignore', divide='ignore'):
        sY = Y * roi[None]
        RawSeed = np.nanmean(sY, 1)
        sYsd = np.nanstd(sY, 0, ddof=1)
        sYsd[sYsd == 0] = np.nan
        Seed = np.nanmean(sY / sYsd, 1)
        InitSeed = Seed.copy()
        Ysd = np.nanstd(Y, 0, ddof=1)
        Ysd[Ysd == 0] = np.nan
        Y = Y / Ysd
    if np.nanmax(Seed) == 0:
        raise RuntimeError('empty seed')
    print('Extracting sLFO...', flush=True)
    Lag = np.full(V, np.nan)
    Lag[~Bmask] = 100
    R, I = _corr_peak(Seed, Y, 2)
    I[R < THR] = 0
    I[~Bmask] = 0
    Lag[I == 3] = 0
    Seed0 = np.nanmean(Y[:, I == 3], 1)
    Seeds = [Seed0]
    maxR = R.copy()
    SeedU = Seed0.copy(); SeedD = Seed0.copy()
    Down = Y; Up = Y
    print('Tracking cross-correlogram peak...', flush=True)
    for p in range(1, limit + 1):
        Down = np.vstack([Down[1:], np.nanmean(Down, 0)[None]])
        R, I = _corr_peak(SeedD, Down, 1)
        maxR = np.fmax(maxR, R)
        I[R < THR] = 0
        if FIXED == 0:
            SeedD = Down[:, I == 2].mean(1)
        I[~np.isnan(Lag)] = 0
        Lag[I == 2] = -p
        Seeds.append(SeedD.copy())
        Up = np.vstack([np.nanmean(Up, 0)[None], Up[:-1]])
        R, I = _corr_peak(SeedU, Up, 1)
        maxR = np.fmax(maxR, R)
        I[R < THR] = 0
        if FIXED == 0:
            SeedU = Up[:, I == 2].mean(1)
        I[~np.isnan(Lag)] = 0
        Lag[I == 2] = p
        Seeds.insert(0, SeedU.copy())
        if on_step:
            on_step(p, Lag)
        if span:
            span(p / limit)
    Lag[Lag == 100] = np.nan
    # mean percent-signal time course of the voxels of each lag (-limit..limit): the sLFO with its real amplitude
    RegionMean = np.full((T, 2 * limit + 1), np.nan)
    with np.errstate(all='ignore'):
        for k, L in enumerate(range(-limit, limit + 1)):
            m = Lag == L
            if m.any():
                RegionMean[:, k] = np.nanmean(Y0[:, m], 1)
    return dict(Lag=Lag, maxR=maxR, Seeds=np.stack(Seeds, 1), RawSeed=RawSeed, InitSeed=InitSeed, Ysd=Ysd, RegionMean=RegionMean)


def _fill_nan_once(Y, strict=False):
    nb = np.stack([np.roll(Y, s, a) for a in range(3) for s in (-1, 1)], -1)
    with np.errstate(all='ignore'):
        return nb.mean(-1) if strict else np.nanmean(nb, -1)


def erode_lag(Lag, Brain, MaxLag=None):
    """drErode_Lag: drop |lag| >= MaxLag, fill every NaN iteratively with the mean of
    the 6 neighbours (circular), then mask with ``Brain``."""
    Y = np.array(Lag, np.float64)
    if MaxLag is not None:
        Y[np.abs(Y) >= MaxLag] = np.nan
    while np.isnan(Y).any():
        nanmask = np.isnan(Y)
        newY = _fill_nan_once(Y)
        Y[nanmask] = newY[nanmask]
        if np.isnan(Y).all():
            break
    Y[~(np.asarray(Brain) != 0)] = np.nan
    return Y


def erode1(Lag):
    """drErode1: one nanmean fill pass, then fill lone holes with the plain neighbour mean."""
    Y = np.array(Lag, np.float64)
    nanmask = np.isnan(Y)
    newY = _fill_nan_once(Y)
    Y[nanmask] = newY[nanmask]
    nanmask = np.isnan(Y)
    newY = _fill_nan_once(Y, strict=True)
    Y[nanmask] = newY[nanmask]
    return Y


def lag4d(name, TR, vols, PosiMax, THR=0.3, FIXED=1, Sm=8, rng=None, reso=None,
          seed_mask=None, mask_pct=10, lp_hz=None, amp_gate=4.0, cwd=None, overwrite=False, span=(0.0, 1.0)):
    """Lag mapping of a 4D BOLD file (drLag4Drev7 / _longTR / _monkey).

    name      string appended to the result folder name
    TR        repetition time (s)
    vols      4D NIfTI (motion corrected, normalised; e.g. output of :func:`merge4d`)
    PosiMax   tracking range +-PosiMax, in **TR** if ``reso`` is None, else in **seconds**
    THR       minimum cross-correlogram peak (0.3 human default, 0.2 in the pipelines)
    FIXED     1 = fixed-seed tracking, 0 = recursive tracking
    Sm        spatial smoothing FWHM (mm); 0 = none
    rng       time points to use (MATLAB 'a:b' / 'a:s:b' string, slice, or index array)
    reso      tracking step in seconds (data are resampled TR -> reso); None = one TR
    seed_mask NIfTI whose voxels (> 0.1 after trilinear reslice) give the initial global
              signal ('hcp' = bundled cerebral mask for HCP-style 2 mm -> subsamp2offc data);
              None = whole brain mask
    mask_pct  brain mask threshold, % of the robust range of the mean image (10 human, 15 monkey)
    lp_hz     low-pass cut-off; default 0.9 / (2*PosiMax [s])
    amp_gate  voxels whose |percent signal| exceeds this are discarded
    cwd       working directory (intermediates + result folder); default current dir
    Returns the result directory.
    """
    cwd = os.path.abspath(cwd or os.getcwd())
    TR = float(TR); PosiMax = float(PosiMax); THR = float(THR); FIXED = int(FIXED); Sm = float(Sm)
    if os.path.exists(vols):
        from .einsteining import check_tr
        check_tr(vols, TR)
    MaxLag = PosiMax * 2
    if reso:
        step, unit, limit = float(reso), 'sec', int(np.ceil(PosiMax / reso))
        MaxLag_sec = MaxLag
    else:
        step, unit, limit = TR, 'TR', int(np.ceil(PosiMax))
        MaxLag_sec = MaxLag * TR
    Smooth = f'{Sm:g}'
    tag = f'{MaxLag:g}{unit}_thr{int(round(10 * THR))}_sm{int(round(Sm))}_{name}'
    outdir = os.path.join(cwd, ('Lag_rec_' if FIXED == 0 else 'Lag_fix_') + tag)
    if os.path.exists(os.path.join(outdir, 'MaxR.nii')) and not overwrite:
        print('..use existing LagMap', outdir)
        return outdir

    p0, p1 = span
    sub = lambda a, b: (p0 + (p1 - p0) * a, p0 + (p1 - p0) * b)
    sm_file = os.path.join(cwd, f'sm{Smooth}_{MaxLag:g}{unit}.nii')
    if not os.path.exists(sm_file):
        print('Preparing the data...', flush=True)
        progress.report('lag4d: preparing (smoothing / filtering)', sub(0, 0.4)[0])
        prepare(vols, TR, MaxLag_sec, Sm, mask_pct, lp_hz, out_sm=sm_file, cwd=cwd)

    mask_img = nib.load(os.path.join(cwd, 'Mask.nii'))
    Mask = mask_img.get_fdata()
    if seed_mask:
        sm_path = HCP_SEED_MASK if seed_mask == 'hcp' else seed_mask
        sm_img = nib.load(sm_path)
        ROI = reslice(sm_img.get_fdata(), sm_img.affine, mask_img.affine, Mask.shape, order=1)
        ROI = np.nan_to_num(ROI)
        _save(os.path.join(cwd, 'rBrainMask.nii'), ROI, mask_img.affine)
        Bmask = (Mask + ROI) != 0
        ROI[ROI < 0.1] = np.nan
    else:
        Bmask = Mask != 0
        ROI = np.where(Bmask, 1.0, np.nan)

    progress.report('lag4d: reading volumes', sub(0.4, 0.5)[0])
    print('Reading volumes...', flush=True)
    rng = parse_range(rng)
    Y, img = _load4d(sm_file, rng)
    shp = Y.shape[:3]
    MAX = np.nanmax(np.abs(Y), 3)
    Y = Y * (MAX <= amp_gate)[..., None]
    Y = Y.reshape(-1, Y.shape[3]).T.astype(np.float64)          # (T, V)
    if reso and abs(TR - reso) > 1e-9:
        Y = resample_poly(Y, int(round(TR * 100)), int(round(reso * 100)), axis=0, window=('kaiser', 5.0))

    os.makedirs(outdir, exist_ok=True)
    aff = img.affine
    def on_step(p, Lag):
        L = Lag.copy(); L[L > 99] = np.nan
        _save(os.path.join(outdir, 'LagOrig_temp.nii'), L.reshape(shp) * step, aff)
        print(f'  +-{p * step:.2f} s', flush=True)
    res = track(Y, Bmask.reshape(-1), ROI.reshape(-1), THR, limit, FIXED, on_step,
                span=progress.Span('lag4d: tracking', *sub(0.5, 0.95)))
    del Y

    savemat(os.path.join(outdir, 'RawSeed.mat'), {'Seed': res['RawSeed'][:, None]})
    savemat(os.path.join(outdir, 'InitSeed.mat'), {'Seed': res['InitSeed'][:, None]})
    savemat(os.path.join(outdir, 'Ysd.mat'), {'Ysd': np.nan_to_num(res['Ysd']).reshape(shp).ravel(order='F')[None]})   # MATLAB voxel order
    savemat(os.path.join(outdir, 'Seeds.mat'), {'Seeds': res['Seeds']})
    np.save(os.path.join(outdir, 'Seeds.npy'), res['Seeds'])
    savemat(os.path.join(outdir, 'RegionMean.mat'), {'RegionMean': res['RegionMean']})
    np.save(os.path.join(outdir, 'RegionMean.npy'), res['RegionMean'])

    Lag1 = res['Lag'].reshape(shp) * step
    _save(os.path.join(outdir, 'LagOrig.nii'), Lag1, aff)
    LagMap = erode_lag(Lag1, Bmask, np.nanmax(Lag1))
    _save(os.path.join(outdir, 'LagMap.nii'), LagMap, aff)
    _save(os.path.join(outdir, 'e1LagOrig.nii'), erode1(Lag1), aff)
    _save(os.path.join(outdir, 'MaxR.nii'), res['maxR'].reshape(shp), aff)
    with open(os.path.join(outdir, 'params.json'), 'w') as f:
        json.dump(dict(name=name, TR=TR, vols=os.path.abspath(vols), PosiMax=PosiMax, THR=THR, FIXED=FIXED, Sm=Sm,
                       range=None if rng is None else np.asarray(rng).tolist(), reso=reso, step=step,
                       seed_mask=seed_mask, mask_pct=mask_pct, lp_hz=lp_hz, amp_gate=amp_gate), f, indent=1)
    progress.report('lag4d: finished', sub(1, 1)[0])
    print('Finished', outdir, flush=True)
    return outdir
