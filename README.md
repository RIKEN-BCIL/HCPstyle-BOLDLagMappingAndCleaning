# HCPstyle-BOLDLagMappingAndCleaning

contact: Toshihiko ASO aso.toshihiko@gmail.com / https://www.researchgate.net/profile/Toshihiko_Aso

## **Pipeline for extraction and removal of the perfusion lag structure within 4D BOLD-fMRI data in HCP-style protocols**

![lagmap_anim](https://github.com/RIKEN-BCIL/BOLDLagMapping/blob/master/lagmap_anim.gif)
![smoothnoisestructure](https://upload.wikimedia.org/wikipedia/commons/thumb/9/9c/Hybrid_image_decomposition.jpg/256px-Hybrid_image_decomposition.jpg)
![sLFO_anim](https://github.com/RIKEN-BCIL/BOLDLagMapping/blob/master/Lag_model_anim100.gif)

#### - BOLD deperfusioning is extracting Einstein (local neurovascular coupling) by removing smooth Marilyn Monroe (perfusion structure) from this image, so that the fMRI result becomes sharp and precise. 

### Dependencies
For Linux/Mac. MATLAB scripts call [FSL][] commands and [SPM12] functions. 
Install FSL & MATLAB then evoke MATLAB from the shell.

(fslmaths in FSL6 may cause errors with the option "-subsamp2offc".)

[FSL]: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki "FSL"
[SPM12]: https://www.fil.ion.ucl.ac.uk/spm/software/spm12/

### Usage

**NEW!** HCP-style pipeline **Einsteining_v01.m** (Einsteining version 0.1) released
with
**drLag4Drev4_hcp.m** for tracking and **drDeperf_hcp_seed.m** for deperfusioning.

- **Einsteining** is a hard working process to remove the perfusion lag structure, a major component of physiological noise for fMRI. 
- Application **before** ICA-FIX is recommended (see Aso & Hayashi, ISMRM2022)
- Modify "Fdir" to point to your **FSL5** installation.

-~-~-~-~-~-~-~-~-~-~

Older scripts for more general use:
https://github.com/RIKEN-BCIL/BOLDLagMapping
