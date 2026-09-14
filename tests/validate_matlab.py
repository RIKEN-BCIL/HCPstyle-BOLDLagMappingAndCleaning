"""Stage-wise comparison of the Python port against MATLAB outputs.

Results on the author's data (2026-09): scrub, reslice, track (human) and monkey stages are
identical (lag values bit-exact, float32 rounding elsewhere); merge differs by <= 2.5e-5
relative (FSL float32 arithmetic); prep cannot be compared with the stored reference (its
REST4run.nii.gz was regenerated after sm8_22TR.nii) -> use the 'hcpfull' stage with a
reference produced by the released drLag4Drev7.m; deperf differs by ~1e-4 relative because
the reference was made with niimath, whose -bptf deviates from FSL's by up to ~1 % of the
filtered signal, plus one voxel zeroed by fsl_regfilt's automask.
usage: python tests/validate_matlab.py <stage> <valdir>   stage in scrub|merge|prep|track|hcpfull|reslice|deperf|monkey
Paths of the reference data are those on the author's system (edit HCP / NHP below)."""
import sys, os, time, shutil, glob
import numpy as np, nibabel as nib
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from boldlag import filters, spm
from boldlag import lag4d, merge4d, deperf, einsteining

HCP = '/mnt/hdd/data/HCP-data/HCPdata-retest/103818/MNINonLinear/Results'
LCS = HCP + '/Lag_concat_scrub'
LAG = LCS + '/Lag_fix_22TR_thr0_sm8_cat4'
NHP = '/mnt/pub/PROJ/NHP_NNP/MacaqueRhesus/rfmri-awake/sub-A03462/MNINonLinear/Results/Lag_concat_scrub'
TR = 0.72

def ld(f, dtype=np.float32):
    return np.asanyarray(nib.load(f).dataobj).astype(dtype)

def cmp(name, a, b):
    a = np.asarray(a, np.float64); b = np.asarray(b, np.float64)
    fa, fb = np.isfinite(a), np.isfinite(b)
    both = fa & fb
    d = np.abs(a[both] - b[both])
    scale = np.abs(b[both]).max() if both.any() else 1
    print(f'  {name:26s} shape={a.shape==b.shape} finite(py/ml)={fa.sum()}/{fb.sum()} nanpattern_agree={np.mean(fa==fb):.6f} '
          f'maxabs={d.max() if d.size else 0:.4g} (rel {d.max()/scale if d.size and scale else 0:.2g}) exact={np.mean(d==0) if d.size else 0:.4f} '
          f'r={np.corrcoef(a[both], b[both])[0,1] if both.sum()>2 and a[both].std()>0 else float("nan"):.6f}', flush=True)

def main(stage, V):
    os.makedirs(V, exist_ok=True); t0 = time.time()
    if stage == 'scrub':
        run = HCP + '/rfMRI_REST1_LR/rfMRI_REST1_LR.nii.gz'
        out = einsteining.scrub_run(run, V, downsample=True)
        cmp('z (subsamp2offc)', ld(V + '/zrfMRI_REST1_LR.nii'), ld(LCS + '/zrfMRI_REST1_LR.nii'))
        zi, zm = nib.load(V + '/zrfMRI_REST1_LR.nii'), nib.load(LCS + '/zrfMRI_REST1_LR.nii')
        print('  z affine equal:', np.allclose(zi.affine, zm.affine), zi.header.get_zooms(), zm.header.get_zooms())
        X, Xm = np.loadtxt(V + '/regout_rfMRI_REST1_LR.txt'), np.loadtxt(LCS + '/regout1.txt')
        print('  regout shape', X.shape, Xm.shape, ' maxabs', np.abs(X - Xm).max() if X.shape == Xm.shape else 'SHAPE MISMATCH')
        cmp('mreg_z (regfilt)', ld(out), ld(LCS + '/mreg_zrfMRI_REST1_LR.nii.gz'))
    elif stage == 'merge':
        z = [LCS + f'/mreg_zrfMRI_REST{a}_{b}.nii.gz' for a, b in [(1, 'LR'), (1, 'RL'), (2, 'LR'), (2, 'RL')]]
        out = merge4d.merge4d('REST4run', TR, z, V)
        cmp('REST4run', ld(out), ld(LCS + '/REST4run.nii.gz'))
    elif stage == 'prep':
        lag4d.prepare(LCS + '/REST4run.nii.gz', TR, 22 * TR, 8, mask_pct=15, out_sm=V + '/sm8_22TR.nii', cwd=V)
        for f in ['SD.nii', 'Tmean.nii', 'Mask.nii']:
            cmp(f, ld(V + '/' + f), ld(LCS + '/' + f))
        cmp('sm8_22TR (subset t<600)', ld(V + '/sm8_22TR.nii')[..., :600], np.asanyarray(nib.load(LCS + '/sm8_22TR.nii').dataobj[..., :600]))
    elif stage == 'track':
        for f in ['Mask.nii', 'Tmean.nii', 'sm8_22TR.nii']:
            if not os.path.exists(V + '/' + f): os.symlink(LCS + '/' + f, V + '/' + f)
        d = lag4d.lag4d('cat4', TR, LCS + '/REST4run.nii.gz', 11, THR=0, FIXED=1, Sm=8, seed_mask='hcp', mask_pct=15, cwd=V, overwrite=True)
        cmp('rBrainMask', ld(V + '/rBrainMask.nii'), ld(LCS + '/rBrainMask_lag_subsamp2offc.nii'))
        for f in ['LagOrig.nii', 'LagMap.nii', 'e1LagOrig.nii', 'MaxR.nii']:
            cmp(f, ld(d + '/' + f), ld(LAG + '/' + f))
        import h5py
        def h5(f, k):
            with h5py.File(f, 'r') as h: return h[k][()].T
        S, Sm_ = np.load(d + '/Seeds.npy'), h5(LAG + '/Seeds.mat', 'Seeds')
        cmp('Seeds', S, Sm_)
        from scipy.io import loadmat
        for f, k in [('InitSeed.mat', 'Seed'), ('RawSeed.mat', 'Seed'), ('Ysd.mat', 'Ysd')]:
            cmp(f, loadmat(d + '/' + f)[k].ravel(), h5(LAG + '/' + f, k).ravel())
    elif stage == 'hcpfull':
        # consistent reference: released drLag4Drev7('cat4','0.72',REST4run,'11','0','1','8') run with MATLAB+SPM in REF
        REF = os.environ.get('REF', '/home/aso/lagref_hcp')
        d = lag4d.lag4d('cat4', TR, REF + '/REST4run.nii.gz', 11, THR=0, FIXED=1, Sm=8, seed_mask='hcp', mask_pct=10, cwd=V, overwrite=True)
        for f in ['SD.nii', 'Tmean.nii', 'Mask.nii']:
            cmp(f, ld(V + '/' + f), ld(REF + '/' + f))
        print('  Mask binarised differing voxels:', int((ld(V + '/Mask.nii') > 0).sum() - (ld(REF + '/Mask.nii') > 0).sum()), (( ld(V + '/Mask.nii') > 0) != (ld(REF + '/Mask.nii') > 0)).sum())
        a, b = nib.load(V + '/sm8_22TR.nii'), nib.load(REF + '/sm8_22TR.nii')
        for t0_ in (0, 2400):
            cmp(f'sm8_22TR t={t0_}..+400', np.asanyarray(a.dataobj[..., t0_:t0_ + 400]), np.asanyarray(b.dataobj[..., t0_:t0_ + 400]))
        R = REF + '/Lag_fix_22TR_thr0_sm8_cat4'
        for f in ['LagOrig.nii', 'LagMap.nii', 'e1LagOrig.nii', 'MaxR.nii']:
            cmp(f, ld(d + '/' + f), ld(R + '/' + f))
        from scipy.io import loadmat
        def anymat(f, k):
            try:
                return loadmat(f)[k]
            except NotImplementedError:
                import h5py
                with h5py.File(f, 'r') as h: return h[k][()].T
        cmp('Seeds', np.load(d + '/Seeds.npy'), anymat(R + '/Seeds.mat', 'Seeds'))
        for f, k in [('InitSeed.mat', 'Seed'), ('RawSeed.mat', 'Seed'), ('Ysd.mat', 'Ysd')]:
            cmp(f, loadmat(d + '/' + f)[k].ravel(), anymat(R + '/' + f, k).ravel())
    elif stage == 'reslice':
        out = einsteining.reslice_lagmap(LAG + '/LagMap.nii', HCP + '/rfMRI_REST1_LR/rfMRI_REST1_LR_SBRef.nii.gz', V + '/rLagMap.nii')
        cmp('rLagMap', ld(out), ld(LCS + '/rLagMap.nii'))
    elif stage == 'deperf':
        ld_ = V + '/lagdir'; os.makedirs(ld_, exist_ok=True)
        for f in ['rLagMap.nii', 'Seeds.mat']:
            if not os.path.exists(ld_ + '/' + f): shutil.copy(LAG + '/' + f, ld_)
        out = deperf.deperf(HCP + '/rfMRI_REST1_LR/rfMRI_REST1_LR.nii.gz', ld_ + '/rLagMap.nii', TR, 1, 4, lagdir=ld_, outdir=V)
        import h5py
        with h5py.File(LCS + '/sLFO.mat', 'r') as h: M = h['Motodata'][()].T     # written by the last run (4) in MATLAB
        cmp('sLFO Motodata (section 4)', deperf.shifted_seeds(deperf.load_seeds(ld_), 11)[3 * 1200:4 * 1200], M)
        cmp('Mask.nii (deperf)', ld(ld_ + '/Mask.nii'), ld(LAG + '/Mask.nii'))
        a = nib.load(out); b = nib.load(LCS + '/rfMRI_REST1_LR_dep.nii.gz')
        for t0_ in (0, 600):
            cmp(f'dep t={t0_}..+200', np.asanyarray(a.dataobj[..., t0_:t0_ + 200]), np.asanyarray(b.dataobj[..., t0_:t0_ + 200]))
    elif stage == 'monkey':
        for f, g in [('Mask.nii', 'Mask.nii'), ('cat4_sm4_6sec_9.nii', 'sm4_6sec.nii')]:
            if not os.path.exists(V + '/' + g): os.symlink(NHP + '/' + f, V + '/' + g)
        for name, rng in [('day9_2half', '1559:3116'), ('cat4_9', None)]:
            d = lag4d.lag4d(name, 0.755, 'unused', 3, THR=0.2, FIXED=1, Sm=4, rng=rng, reso=0.5, seed_mask=None, mask_pct=15, cwd=V, overwrite=True)
            ref = NHP + f'/Lag_fix_6sec_thr2_sm4_{name}'
            for f in ['LagOrig.nii', 'LagMap.nii', 'MaxR.nii']:
                cmp(name + ' ' + f, ld(d + '/' + f), ld(ref + '/' + f))
    print(f'[{stage}] done in {time.time() - t0:.0f} s', flush=True)

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
