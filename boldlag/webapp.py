"""Browser front end (Streamlit): ``boldlag-web`` or ``python -m boldlag.webapp``
(= ``streamlit run boldlag/webapp.py``).  Requires ``pip install streamlit``.

Runs on the machine holding the data; paths are typed (or found from a subject
directory).  Jobs run inside the Streamlit process with live log, progress bar
and a lag-map montage at the end.  Settings can be downloaded / uploaded as JSON.
"""
import os, sys, io, json, contextlib

if __package__ in (None, ''):            # run as a plain script by `streamlit run`
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    import boldlag
    __package__ = 'boldlag'


def _run_streamlit():
    from streamlit.web import cli
    sys.argv = ['streamlit', 'run', os.path.abspath(__file__), '--server.headless', 'true',
                '--client.toolbarMode', 'viewer'] + sys.argv[1:]      # viewer mode hides the Deploy/settings menu
    sys.exit(cli.main())


def main():
    _run_streamlit()


class _StreamLog(io.TextIOBase):
    """stdout replacement that appends to a Streamlit placeholder."""
    def __init__(self, box, limit=200):
        self.box, self.lines, self.limit, self.buf = box, [], limit, ''
    def write(self, s):
        self.buf += s
        while '\n' in self.buf:
            line, self.buf = self.buf.split('\n', 1)
            if line.strip():
                self.lines.append(line)
        self.box.code('\n'.join(self.lines[-self.limit:]) or ' ')
        return len(s)
    def flush(self):
        pass


def app():
    import streamlit as st
    from . import progress
    st.set_page_config(page_title='boldlag', layout='wide')
    st.title('BOLD lag mapping / deperfusioning')
    S = st.session_state
    S.setdefault('runs', [])

    # ---------- sidebar: mode + settings file
    mode = st.sidebar.radio('Function', ['Pipeline (Einsteining)', 'Lag map only (drLag4D)', 'Deperfusioning (drDeperf)'])
    up = st.sidebar.file_uploader('Load settings (JSON)', type='json')
    if up is not None and S.get('_loaded') != up.name + str(up.size):
        d = json.load(up)
        for k, v in d.items():
            S[k] = v
        if 'runs' in d:
            S['runs_txt'] = '\n'.join(d['runs'])
        S['_loaded'] = up.name + str(up.size)
        st.rerun()

    def lag_params(prefix):
        c1, c2, c3 = st.columns(3)
        p = {}
        p['TR'] = c1.number_input('TR (s)', value=float(S.get(prefix + 'TR', 0.72)), format='%.4f', key=prefix + 'TR')
        p['PosiMax'] = c2.number_input('PosiMax (tracking range, TR or s with a tracking step)', value=float(S.get(prefix + 'PosiMax', 9.0)), key=prefix + 'PosiMax')
        p['THR'] = c3.number_input('Min cross-correlogram peak (THR; 0 = accept all)', value=float(S.get(prefix + 'THR', 0.0)), key=prefix + 'THR')
        p['FIXED'] = c1.selectbox('Tracking', ['fixed', 'recursive'], index=0 if S.get(prefix + 'FIXED', 'fixed') == 'fixed' else 1, key=prefix + 'FIXED')
        p['Sm'] = c2.number_input('Smoothing FWHM (mm): up to 8 mm; 8 for human, 4 for monkey, 0 = none', value=float(S.get(prefix + 'Sm', 8.0)), key=prefix + 'Sm')
        p['reso'] = c3.text_input('Tracking step (s; empty = one TR, 0.5 for monkey)', value=S.get(prefix + 'reso', ''), key=prefix + 'reso')
        p['seed'] = c1.text_input("Seed mask ('hcp' = bundled cerebral mask, a NIfTI path, or empty = whole brain)", value=S.get(prefix + 'seed', 'hcp'), key=prefix + 'seed')
        p['mask_pct'] = c2.number_input('Brain mask % of robust range (10 human, 15 monkey)', value=float(S.get(prefix + 'mask_pct', 10.0)), key=prefix + 'mask_pct')
        p['lp_hz'] = c3.text_input('Low-pass Hz (empty = 0.9/(2 PosiMax s))', value=S.get(prefix + 'lp_hz', ''), key=prefix + 'lp_hz')
        return p

    def kw(p):
        f = lambda s: float(s) if str(s).strip() else None
        return dict(TR=p['TR'], PosiMax=p['PosiMax'], THR=p['THR'], FIXED=1 if p['FIXED'] == 'fixed' else 0, Sm=p['Sm'],
                    reso=f(p['reso']), seed_mask=p['seed'].strip() or None, mask_pct=p['mask_pct'], lp_hz=f(p['lp_hz']))

    job = None
    if mode.startswith('Pipeline'):
        st.subheader('Runs (HCP layout: <run>/<run>.nii.gz with Movement_Regressors.txt and *SBRef)')
        c1, c2, c3 = st.columns([3, 1, 1])
        subj = c1.text_input('Subject directory (…/<subject>, containing MNINonLinear/Results)', value=S.get('subj', ''), key='subj')
        pattern = c2.text_input('run folder regex', value=S.get('pattern', '_REST[12]'), key='pattern')
        nvols = c3.text_input('volumes per run (optional)', value=S.get('nvols', ''), key='nvols')
        if st.button('Find runs'):
            from .einsteining import find_runs
            try:
                S['runs'] = find_runs(subj, pattern, nvols=int(nvols) if nvols.strip() else None)
                S['runs_txt'] = '\n'.join(S['runs'])          # keyed widgets ignore value= once they have state
            except Exception as e:
                st.error(str(e))
        S.setdefault('runs_txt', '\n'.join(S.get('runs', [])))
        runs_txt = st.text_area('Run files (one per line)', height=120, key='runs_txt')
        runs = [r.strip() for r in runs_txt.splitlines() if r.strip()]
        st.subheader('Lag mapping')
        p = lag_params('p.')
        st.subheader('Pipeline options')
        c1, c2 = st.columns(2)
        results_dir = c1.text_input('Output folder (empty = Results folder of the runs; an existing work folder is renamed old_*)', value=S.get('results_dir', ''), key='results_dir')
        workname = c2.text_input('Work folder name', value=S.get('workname', 'Lag_concat_scrub'), key='workname')
        c1, c2, c3, c4 = st.columns(4)
        downsample = c1.checkbox('2x down-sampling (subsamp2offc)', value=S.get('downsample', True), key='downsample')
        despike = c2.checkbox('Spike regressors', value=S.get('despike', True), key='despike')
        spike_thr = c3.number_input('spike if DVARS > thr x median', value=float(S.get('spike_thr', 1.5)), key='spike_thr')
        only_lag = c4.checkbox('Lag mapping only', value=S.get('only_lag', False), key='only_lag')
        settings = dict(mode=mode, runs=runs, **{'p.' + k: v for k, v in p.items()}, subj=subj, pattern=pattern, nvols=nvols,
                        results_dir=results_dir, workname=workname, downsample=downsample, despike=despike, spike_thr=spike_thr, only_lag=only_lag)
        if runs:
            k = kw(p)
            def job():
                from .einsteining import einsteining
                d = einsteining(runs, k['TR'], k['PosiMax'], k['THR'], k['FIXED'], k['Sm'], only_lag, downsample, k['reso'],
                                k['seed_mask'], k['mask_pct'], k['lp_hz'], workname, results_dir.strip() or None,
                                spike_thr=spike_thr if despike else None)
                return d, os.path.join(os.path.dirname(d), 'Tmean.nii')
    elif mode.startswith('Lag'):
        c1, c2 = st.columns([3, 1])
        vols = c1.text_input('4D file (motion corrected, normalised)', value=S.get('lag_vols', ''), key='lag_vols')
        name = c2.text_input('Name tag', value=S.get('lag_name', 'run1'), key='lag_name')
        rng = c2.text_input('Time range (MATLAB style, e.g. 1:500)', value=S.get('lag_range', ''), key='lag_range')
        cwd = c1.text_input('Work folder (empty = folder of the 4D file)', value=S.get('lag_cwd', ''), key='lag_cwd')
        p = lag_params('l.')
        settings = dict(mode=mode, lag_vols=vols, lag_name=name, lag_range=rng, lag_cwd=cwd, **{'l.' + k: v for k, v in p.items()})
        if vols.strip():
            k = kw(p)
            def job():
                from .lag4d import lag4d
                wd = cwd.strip() or os.path.dirname(os.path.abspath(vols))
                d = lag4d(name, k['TR'], vols, k['PosiMax'], k['THR'], k['FIXED'], k['Sm'], rng.strip() or None,
                          k['reso'], k['seed_mask'], k['mask_pct'], k['lp_hz'], cwd=wd)
                return d, os.path.join(wd, 'Tmean.nii')
    else:
        c1, c2 = st.columns([3, 1])
        vols = c1.text_input('Original 4D run', value=S.get('dep_vols', ''), key='dep_vols')
        lag = c1.text_input('rLagMap.nii (on the grid of the run)', value=S.get('dep_lag', ''), key='dep_lag')
        lagdir = c1.text_input('Lag folder with Seeds.mat (empty = folder of rLagMap)', value=S.get('dep_lagdir', ''), key='dep_lagdir')
        out = c1.text_input('Output folder', value=S.get('dep_out', '.'), key='dep_out')
        TR = c2.number_input('TR (s)', value=float(S.get('dep_TR', 0.72)), format='%.4f', key='dep_TR')
        sec = c2.number_input('Run number', value=int(S.get('dep_sec', 1)), min_value=1, key='dep_sec')
        n = c2.number_input('Number of runs', value=int(S.get('dep_n', 1)), min_value=1, key='dep_n')
        reso = c2.text_input('Tracking step (s; empty = TR)', value=S.get('dep_reso', ''), key='dep_reso')
        settings = dict(mode=mode, dep_vols=vols, dep_lag=lag, dep_lagdir=lagdir, dep_out=out, dep_TR=TR, dep_sec=sec, dep_n=n, dep_reso=reso)
        if vols.strip() and lag.strip():
            def job():
                from .deperf import deperf
                f = deperf(vols, lag, TR, int(sec), int(n), lagdir.strip() or None, float(reso) if reso.strip() else None, out.strip() or '.')
                print('written', f)
                return None, None

    st.sidebar.download_button('Save settings (JSON)', json.dumps(settings, indent=1, default=str), file_name='boldlag_settings.json')

    if st.button('Run', type='primary', disabled=job is None):
        bar = st.progress(0.0, text='starting')
        logbox = st.empty()
        progress.set_callback(lambda stage, frac: bar.progress(frac, text=f'{stage} ({100 * frac:.0f} %)'))
        log = _StreamLog(logbox)
        try:
            with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
                lagdir, under = job()
            bar.progress(1.0, text='finished')
            S['last'] = (lagdir, under)
        except Exception as e:
            st.exception(e)
        finally:
            progress.set_callback(None)
    if S.get('last') and S['last'][0]:
        lagdir, under = S['last']
        try:
            from .viewer import lagmap_montage
            png = lagmap_montage(os.path.join(lagdir, 'LagMap.nii'), under if under and os.path.exists(under) else None)
            st.image(png, caption=lagdir)
        except Exception as e:
            st.warning(f'montage: {e}')


def _in_streamlit():
    try:
        from streamlit.runtime.scriptrunner import get_script_run_ctx
        return get_script_run_ctx() is not None
    except Exception:
        return False


if _in_streamlit():          # executed by `streamlit run` / AppTest
    app()
elif __name__ == '__main__':  # `python boldlag/webapp.py` / `python -m boldlag.webapp`
    _run_streamlit()
