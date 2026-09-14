# BOLDLagMapping

## Extraction and removal of the sLFO with its time-lag structure in 4D blood oxygenation level dependent (BOLD) signal MRI data

### Introduction to lag mapping
```
Lim = 2;
YY = [];
for Sft = Lim:-1:-Lim
	YY = cat( 3, YY, Y( Lim+Sft+1:end-Lim+Sft, :));
end

disp('Tracking cross-correlogram peak...          ')
XX = repmat( Seed( Lim+1:end-Lim), [ 1 size( YY,2) size( YY,3)]);
CC = sum( XX.*YY, 1)./( sum( XX.*XX, 1).^.5 .* sum( YY.*YY, 1).^.5);
[ R, I] = max( CC, [], 3);
```

This is the core logic of lag mapping picked up from the part that extracts the first sLFO: Y is the original data made two-dimensional, time x voxel.
From this, we create the 3-dimensional array YY with third dimension of lag. They are time-shifted versions of the original data.
"Seed" is the whole-brain signal for this stage; it is repeated by "repmat" to a matrix of the same size to compute correlation with YY, and CC becomes the correlation.
CC is a 3-dimensional matrix 1 x number of voxels x lag, so we find the maximum value (R) in the third dimension direction. The location of maximum gives the lag value. Here Lim=2, so YY and XX have 5 layers in the third dimension, which means that if I=3, the correlation is maximum at the very phase of the whole-brain signal (if the phase is shifted, the correlation drops). Those voxels are determined to have lag = zero, and their average is the sLFO time course. We use the sLFO as the "Seed" to trace up- and downstream recursively.

There are dozens of regressors for removing sLFOs, like “regressor of a group of voxels whose lag is -6TR”. At first I tried to use the average time course of these voxels. But the smaller the voxel group of the lag, the further away from whole-brain variability. Sometimes the task response comes in, and the SN simply drops with small number of voxels.
For this reason, we take the safe approach of “removing the sLFO used to track that lag”. We regress out a time-shifted version of the first sLFO for each region. The next safest thing to do is to use the sLFO obtained during the "recursive lag tracking" where sLFO is updated in each step, but it has a problem similar to the above.
Since the sLFOs are distributed throughout the brain with various phases, their weighted average accounts for (a significant portion of) the low-frequency component of the whole-brain signal. However, if we regress that whole brain signal uniformly from all voxels, the correlation structure due to phase differences remains. This is why we need "deperfusioning".


## Python version (no MATLAB / SPM / FSL needed)

The `boldlag` package is a line-by-line port of the MATLAB scripts (drLag4Drev7,
drMerge4D, drDeperf, Einsteining).  Dependencies: `numpy`, `scipy`, `nibabel`
(`h5py` only to read MATLAB v7.3 `Seeds.mat`).  The FSL operations used by the
MATLAB code (`fslmaths -bptf / -subsamp2offc / -thrp`, `fsl_regfilt`) and the SPM
ones (`spm_smooth` with implicit mask, coregistration-reslice) are re-implemented in
numpy and were verified against the FSL 5.0.9 / 6.0.5 / niimath binaries and against
MATLAB outputs of a full HCP subject and of macaque data (see `tests/validate_matlab.py`;
lag maps are identical when the MATLAB run used fslmaths; with niimath, whose `-bptf`
deviates from FSL's by ~1e-3 of the signal, 0.06 % of voxels change by one lag step).

```
pip install .            # gives the `boldlag` command (or use `python -m boldlag`)
```

Human and monkey data are handled by the same code.  The only conceptual
difference is the region providing the initial global signal (`--seed-mask`):

| MATLAB script | Python |
|---|---|
| `drLag4Drev7( name, TR, vols, PosiMax, THR, FIXED, Smooth, range)` (HCP, step = TR, cerebral seed) | `boldlag lag4d name TR vols PosiMax --thr THR --fixed FIXED --sm Smooth --seed-mask hcp [--range 1:500]` |
| `drLag4Drev7_longTR` (step 1 s, PosiMax in s) | `... --reso 1` |
| `drLag4Drev7_monkey` (step 0.5 s, whole brain seed, mask 15 %) | `... --reso 0.5 --mask-pct 15` |
| `drMerge4D` | `boldlag merge4d REST4run TR run1.nii.gz run2.nii.gz ...` |
| `drDeperf_hcp_seed( vols, rLagMap, TR, section, Nruns)` | `boldlag deperf vols rLagMap.nii TR section Nruns` |
| `Einsteining_v07( Runs, TR, Nvols, MaxLag, MinR, Fixed, Sm, onlyLag)` | `boldlag einsteining TR MaxLag --runs run1.nii.gz run2.nii.gz --thr MinR --fixed Fixed --sm Sm [--only-lag]` |
| `Einsteining_v06( Sdir, ...)` (scan `MNINonLinear/Results` for `_REST1/_REST2`) | `boldlag einsteining TR MaxLag --subject-dir Sdir --pattern '_REST[12]' --nvols 1200 ...` |
| `Einsteining_v06_monkey` | `boldlag einsteining TR MaxLag --subject-dir Sdir --pattern 'BOLD_' --reso 0.5 --mask-pct 15 --seed-mask '' ...` |
| `Einsteining_v06_hireso_monkey` (no down-sampling) | `... --no-downsample` |

Example (HCP subject, as in the MATLAB release notes):

```
boldlag einsteining 0.72 9 --thr 0.2 --sm 8 \
    --runs /data/subject1/MNINonLinear/Results/rfMRI_REST1_PA/rfMRI_REST1_PA.nii.gz \
           /data/subject1/MNINonLinear/Results/rfMRI_REST1_AP/rfMRI_REST1_AP.nii.gz
```

Outputs go to `MNINonLinear/Results/Lag_concat_scrub/` with the same names as the
MATLAB version (`z*.nii`, `mreg_z*.nii.gz`, `REST<n>run.nii.gz`, `sm8_18TR.nii`,
`Lag_fix_18TR_thr2_sm8_cat<n>/{LagOrig,LagMap,MaxR}.nii`, `Seeds.mat`, `rLagMap.nii`,
`<run>_dep.nii.gz`, `sLFO.mat`) plus `<run>_dep/` folders with symbolic links.

Python API: `boldlag.lag4d.lag4d(...)`, `boldlag.merge4d.merge4d(...)`,
`boldlag.deperf.deperf(...)`, `boldlag.einsteining.einsteining(...)`; the building
blocks (`boldlag.filters.bptf/regfilt/subsamp2offc`, `boldlag.spm.smooth/reslice`,
`boldlag.lag4d.track`) can be used on numpy arrays directly.

Notes on differences from the MATLAB scripts
* The low-pass cut-off of the band-pass filter is `0.9 / (2*PosiMax)` Hz with
  PosiMax in seconds (`drLag4Drev7`, `_longTR`).  `drLag4Drev7_monkey` divided by TR
  once more; use `--lp-hz` to reproduce that if needed.
* `fsl_regfilt`'s automatic intensity mask is reproduced in the motion-scrubbing
  step, but not in the deperfusioning step, where it acted on high-passed
  (zero-mean) data and could zero out a few voxels arbitrarily.
* `Seeds.npy`/`params.json` are written next to the MATLAB-compatible `.mat` files.
* Motion scrubbing follows the fix of [issue #2](https://github.com/aso-toshihiko/BOLDLagMapping_Deperfusioning/issues/2):
  DVARS is computed on the raw intensity within the brain, in percent of the mean
  brain signal, and a volume (plus the previous one) is a spike when DVARS exceeds
  `--spike-thr` (default 1.5) x median; all regressors including the spike columns
  are removed (`--no-despike` to omit them).  The same fix is applied to the MATLAB scripts in `matlab/`.

（最後に日本語あり）

### Scripts: `matlab/` (MATLAB, issue #2 fixed) and `boldlag/` (Python port, see below). Older versions in Releases.
contact: Toshihiko ASO aso.toshihiko@gmail.com / https://www.researchgate.net/profile/Toshihiko_Aso

![lagmaps](https://github.com/RIKEN-BCIL/BOLDLagMapping/blob/master/LagMaps.jpg)
![lagmap_anim](https://github.com/RIKEN-BCIL/BOLDLagMapping/blob/master/lagmap_anim.gif)
![sLFO_anim](https://github.com/RIKEN-BCIL/BOLDLagMapping/blob/master/Lag_model_anim100.gif)

### Dependencies
For Linux/Mac. MATLAB scripts call [FSL][] commands and [SPM12] functions.
FSL6 + niimath or FSL5's fslmaths needed for resampling in "Einsteining" scripts.
Install FSL & MATLAB then evoke MATLAB from the shell.

[FSL]: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki "FSL"
[SPM12]: https://www.fil.ion.ucl.ac.uk/spm/software/spm12/

### See Releases for usage

**drLag4D** for tracking and **drDeperf** for deperfusioning.
**Einsteining** is the pipeline script.

![smoothnoisestructure](https://upload.wikimedia.org/wikipedia/commons/thumb/9/9c/Hybrid_image_decomposition.jpg/256px-Hybrid_image_decomposition.jpg)

BOLD deperfusioning is extracting Einstein (neurovascular coupling) by removing Marilyn Monroe (perfusion structure) from this image. For this purpose, first we smooth the original image to enhance the Marilyn.


### References

Recursive tracking

[Aso, T., Urayama, S., Hidenao, F., & Murai, T. (2019). Axial variation of deoxyhemoglobin density as a source of the low-frequency time lag structure in blood oxygenation level-dependent signals. PLoS ONE.](https://doi.org/10.1371/journal.pone.0222787) [(Correction here)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0225489)

[Nishida, S., Aso, T., Takaya, S., Takahashi, Y., Kikuchi, T., Funaki, T., … Miyamoto, S. (2018). Resting-state Functional Magnetic Resonance Imaging Identifies Cerebrovascular Reactivity Impairment in Patients With Arterial Occlusive Diseases: A Pilot Study. Neurosurgery, 85(5), 680-688.](https://doi.org/10.1093/neuros/nyy434)

[Aso, T., Jiang, G., Urayama, S. I., & Fukuyama, H. (2017). A resilient, non-neuronal source of the spatiotemporal lag structure detected by bold signal-based blood flow tracking. Frontiers in Neuroscience, 11(MAY), 1-13.](https://doi.org/10.3389/fnins.2017.00256)

Fixed-seed tracking

[Aso, T., Sugihara, G., Murai, T., Ubukata, S., Urayama, S., Ueno, T., Fujimoto, G., Thuy, D., Fukuyama, H., & Ueda, K. (2020). A venous mechanism of ventriculomegaly shared between traumatic brain injury and normal ageing. Brain, 143(JUN)](https://doi.org/10.1093/brain/awaa125)

[Satow, T., Aso, T., Nishida, S., Komuro, T., Ueno, T., Oishi, N., … Fukuyama, H. (2017). Alteration of venous drainage route in idiopathic normal pressure hydrocephalus and normal aging. Frontiers in Aging Neuroscience, 9(NOV), 1–10.](https://doi.org/10.3389/fnagi.2017.00387)


下記は最初のsLFOを作る部分です。Yは縦が時間、横が全ボクセルの2次元にした元データです。
ここからYYという3次元データを作り、この三次元目がラグになります。要するに縦（時間）が一個ずつズレていくだけです。
Seedは、この段階では全脳信号です。YYと相関を計算するために同じサイズの行列にrepmatで増やし、CCが相関になります。
CCは縦が１，横がボクセル数、三次元目がラグの3次元行列なので3次元目方向に最大値（下ではR）を求め、どこで最大になるか（下では変数I）がラグ値になります。
ここではLim=2なので、YY、XXは3次元目が５の大きさを持ち、つまりI=3であれば全脳信号の元の位相のときに相関が最大ということになります（逆に言うと位相をずらしたら相関が下がる）。そのボクセルたちをラグ＝ゼロと決定し、平均をsLFOとします。
あとはsLFOをSeedにして上流、下流に同じ方法でたどっていきます。

```
Lim = 2;
YY = [];
for Sft = Lim:-1:-Lim
	YY = cat( 3, YY, Y( Lim+Sft+1:end-Lim+Sft, :));
end

disp('Tracking cross-correlogram peak...          ')
XX = repmat( Seed( Lim+1:end-Lim), [ 1 size( YY,2) size( YY,3)]);
CC = sum( XX.*YY, 1)./( sum( XX.*XX, 1).^.5 .* sum( YY.*YY, 1).^.5);
[ R, I] = max( CC, [], 3);
```

sLFOを除去する際のregressorは、「ラグが-6TRであるボクセル群のregressor」、という風に何十個とあります。最初は本当にこれらのボクセルの平均タイムコースを使ってみました。しかし、そうするとボクセルが少ないラグほど全脳変動から離れていき、タスク変動も入ってきたり、単純にボクセル数に応じてSNも下がります。
このため、「そのラグをトラックしたときに使ったsLFOを除去する」という安全な方法をとっています。要するに最初に作ったsLFOを単に時間的にシフトしたものです。その次に安全なのが、再帰的なsLFOの更新を用いたトラッキングでのsLFOを使うことですが、上記と近い問題があります。
sLFOがいろんな位相で脳内に分布してるので、それらの加重平均が全脳信号の低周波成分（のかなりの部分）を占めています。しかし、その全脳信号を全体でregressionしてしまうと、位相差による相関構造が残る。このためラグ構造の全体を除去したほうがいいわけです。
