"""boldlag -- Python port of the BOLDLagMapping / Deperfusioning MATLAB scripts (T. Aso).

Modules
-------
lag4d       drLag4Drev7 (human / longTR / monkey variants unified)
merge4d     drMerge4D
deperf      drDeperf_hcp_seed / drDeperf_longTR
einsteining Einsteining pipeline (scrubbing -> merge -> lag mapping -> deperfusioning)
filters     numpy re-implementations of fslmaths -bptf/-subsamp2offc/-thrp and fsl_regfilt
spm         SPM12 smoothing (implicit mask) and coregistration-reslice
"""
__version__ = '0.1.0'
from . import filters, spm, lag4d, merge4d, deperf, einsteining
__all__ = ['filters', 'spm', 'lag4d', 'merge4d', 'deperf', 'einsteining']
