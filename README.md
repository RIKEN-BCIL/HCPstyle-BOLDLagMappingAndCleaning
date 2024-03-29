# HCPstyle-BOLDLagMappingAndCleaning

## Try the new scripts from "Release"

Feedback appreciated: Toshihiko ASO aso.toshihiko@gmail.com / https://www.researchgate.net/profile/Toshihiko_Aso

## **Pipeline for extraction and removal of the perfusion lag structure within 4D BOLD-fMRI data in HCP-style protocols**

![lagmaps](https://github.com/RIKEN-BCIL/BOLDLagMapping/blob/master/LagMaps.jpg)
![lagmap_anim](https://github.com/RIKEN-BCIL/BOLDLagMapping/blob/master/lagmap_anim.gif)
![smoothnoisestructure](https://upload.wikimedia.org/wikipedia/commons/thumb/9/9c/Hybrid_image_decomposition.jpg/256px-Hybrid_image_decomposition.jpg)

## BOLD deperfusioning is extracting Einstein (local neurovascular coupling) by removing smooth Marilyn Monroe (perfusion structure), so that the fMRI result becomes sharp and precise.

## Dependencies

For Linux/Mac. MATLAB scripts call [FSL][] commands and [SPM12] functions. 
- Install FSL & MATLAB then evoke MATLAB from the shell.
- **niimath** from https://github.com/rordenlab/niimath is preferred because **fslmaths** in FSL6 may not funciton with the option "-subsamp2offc". FSL5 works fine, but no longer actively maintained.


[FSL]: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki "FSL"
[SPM12]: https://www.fil.ion.ucl.ac.uk/spm/software/spm12/

## Usage

## Main script **Einsteining_vXX.m** calls **drLag4DrevXX.m** for lag mapping and **drDeperf_hcp_seed_niimath.m** for deperfusioning.

- **Einsteining** is a hard working process to pick up Einstein face (neural signal) by (1) extracting perfusion lag structure, a major component of physiological noise for fMRI, and (2) removing it. 
- Application **before** ICA-FIX is recommended (see Aso & Hayashi, ISMRM2022)
- Modify "Fdir" to point to your **FSL5** installation.
- For drLag4D scripts, note that tracking range should now be specified in TR, not seconds. Arguments are given in text at present.

## Examples

 For input data _BOLD_REST1_AP, BOLD_REST1_PA, BOLD_REST2_AP, BOLD_REST2_PA_ 
 
 output directories will be _dep_BOLD_REST1_AP, dep_BOLD_REST1_PA, dep_BOLD_REST2_AP, dep_BOLD_REST2_PA_

### Example usage:
 
 _Einsteining_v02( Sdir, 0.8, 478, 8, 0.2, 1, 8, 1)_

- Sdir: Subject directory
- TR = 0.8s
- Number of volumes per run = 375
- MaxLag = 8 -> lag mapping up to -8TR / +8TR
- MinR = 0.2 -> Cross-correlogram peak < 0.2 will not be used
- Fixed = 0  -> Recursive (flexible) LFO instead of Fixed-sLFO algorithm
- Spatial smoothing at 8mm FWHM
- onlyLag = 0 -> Lag mapping + deperfusioning
 

![sLFO_anim](https://github.com/RIKEN-BCIL/BOLDLagMapping/blob/master/Lag_model_anim100.gif)

# Older scripts and descriptions for more general use:
https://github.com/RIKEN-BCIL/BOLDLagMapping
