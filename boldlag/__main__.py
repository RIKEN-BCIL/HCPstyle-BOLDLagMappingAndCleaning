"""Command line interface: ``boldlag <command> ...`` (also ``python -m boldlag``)."""
import argparse, sys, os


def _common_lag_args(p):
    p.add_argument('--thr', type=float, default=0.2, help='min cross-correlogram peak (default 0.2)')
    p.add_argument('--fixed', type=int, default=1, help='1 = fixed-seed (default), 0 = recursive tracking')
    p.add_argument('--sm', type=float, default=8, help='smoothing FWHM mm (default 8)')
    p.add_argument('--reso', type=float, default=None, help='tracking step in s (data resampled); default = 1 TR')
    p.add_argument('--seed-mask', default=None, help="'hcp' (bundled cerebral mask), a NIfTI file, or omit for whole brain")
    p.add_argument('--mask-pct', type=float, default=10, help='brain mask: %% of robust range of the mean image (10; 15 for monkey)')
    p.add_argument('--lp-hz', type=float, default=None, help='low-pass cut-off Hz (default 0.9/(2*PosiMax s))')


def main(argv=None):
    ap = argparse.ArgumentParser(prog='boldlag', description=__doc__)
    sub = ap.add_subparsers(dest='cmd', required=True)

    p = sub.add_parser('lag4d', help='lag mapping of one 4D file (drLag4Drev7)')
    p.add_argument('name'); p.add_argument('TR', type=float); p.add_argument('vols')
    p.add_argument('PosiMax', type=float, help='tracking range; in TR (default) or in s with --reso')
    _common_lag_args(p)
    p.add_argument('--range', default=None, help="time points, MATLAB style e.g. 1:500 or 1:2:500")
    p.add_argument('--cwd', default=None)

    p = sub.add_parser('merge4d', help='concatenate runs (drMerge4D)')
    p.add_argument('name'); p.add_argument('TR', type=float); p.add_argument('vols', nargs='+')
    p.add_argument('--outdir', default='.')

    p = sub.add_parser('deperf', help='deperfusion one run (drDeperf)')
    p.add_argument('vols'); p.add_argument('lag', help='rLagMap.nii on the grid of vols'); p.add_argument('TR', type=float)
    p.add_argument('section', type=int); p.add_argument('Nruns', type=int)
    p.add_argument('--lagdir', default=None, help='folder with Seeds.mat (default: folder of lag)')
    p.add_argument('--reso', type=float, default=None); p.add_argument('--outdir', default='.')

    p = sub.add_parser('einsteining', help='whole pipeline on one subject')
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--runs', nargs='+', help='4D run files (HCP layout)')
    g.add_argument('--subject-dir', help='subject dir; runs found under MNINonLinear/Results by --pattern')
    p.add_argument('--pattern', default=r'_REST[12]', help='regex for run folders (default _REST[12]; e.g. BOLD_)')
    p.add_argument('--nvols', type=int, default=None, help='keep only runs with this many volumes')
    p.add_argument('TR', type=float); p.add_argument('PosiMax', type=float)
    _common_lag_args(p)
    p.set_defaults(seed_mask='hcp')
    p.add_argument('--only-lag', action='store_true'); p.add_argument('--no-downsample', action='store_true')
    p.add_argument('--workname', default='Lag_concat_scrub')
    p.add_argument('--results-dir', default=None, help='where to write <workname>/ and <run>_dep/ (default: the Results folder of the runs)')
    p.add_argument('--spike-thr', type=float, default=1.5, help='spike if DVARS > thr x median (default 1.5)')
    p.add_argument('--no-despike', action='store_true', help='no spike regressors (motion + FD only)')

    a = ap.parse_args(argv)
    from boldlag import lag4d, merge4d, deperf, einsteining
    if a.cmd == 'lag4d':
        print(lag4d.lag4d(a.name, a.TR, a.vols, a.PosiMax, a.thr, a.fixed, a.sm, a.range, a.reso, a.seed_mask,
                          a.mask_pct, a.lp_hz, cwd=a.cwd))
    elif a.cmd == 'merge4d':
        print(merge4d.merge4d(a.name, a.TR, a.vols, a.outdir))
    elif a.cmd == 'deperf':
        print(deperf.deperf(a.vols, a.lag, a.TR, a.section, a.Nruns, a.lagdir, a.reso, a.outdir))
    elif a.cmd == 'einsteining':
        runs = a.runs or einsteining.find_runs(a.subject_dir, a.pattern, nvols=a.nvols)
        print(einsteining.einsteining(runs, a.TR, a.PosiMax, a.thr, a.fixed, a.sm, a.only_lag, not a.no_downsample,
                                      a.reso, a.seed_mask, a.mask_pct, a.lp_hz, a.workname, a.results_dir, spike_thr=None if a.no_despike else a.spike_thr))


if __name__ == '__main__':
    main()
