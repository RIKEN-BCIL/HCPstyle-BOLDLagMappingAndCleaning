"""Tkinter front end for boldlag (``boldlag-gui`` / ``python -m boldlag.gui``).

Three tabs -- the whole Einsteining pipeline, lag mapping of a single 4D file, and
deperfusioning of one run -- run in a background thread with their console output
shown in the log pane; the resulting lag map is displayed as a montage when done.
No dependencies beyond the package itself (matplotlib for the montage).
"""
import os, sys, glob, json, threading, queue, traceback
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from . import progress


class _QueueWriter:
    def __init__(self, q):
        self.q = q
    def write(self, s):
        if s:
            self.q.put(s)
    def flush(self):
        pass


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title('BOLD lag mapping / deperfusioning (boldlag)')
        self.geometry('1100x780')
        self.q = queue.Queue()
        self.worker = None
        self.result = None
        menu = tk.Menu(self); self.config(menu=menu)
        fm = tk.Menu(menu, tearoff=0); menu.add_cascade(label='File', menu=fm)
        fm.add_command(label='Load settings (JSON)...', command=self.load_settings)
        fm.add_command(label='Save settings (JSON)...', command=self.save_settings)
        fm.add_separator(); fm.add_command(label='Quit', command=self.destroy)
        nb = ttk.Notebook(self)
        nb.pack(fill='x', padx=8, pady=6)
        self.tab_pipe = ttk.Frame(nb); self.tab_lag = ttk.Frame(nb); self.tab_dep = ttk.Frame(nb)
        nb.add(self.tab_pipe, text=' Pipeline (Einsteining) ')
        nb.add(self.tab_lag, text=' Lag map only (drLag4D) ')
        nb.add(self.tab_dep, text=' Deperfusioning (drDeperf) ')
        self.nb = nb
        self._build_pipeline()
        self._build_lag()
        self._build_dep()
        bar = ttk.Frame(self); bar.pack(fill='x', padx=8)
        self.run_btn = ttk.Button(bar, text='Run', command=self.run)
        self.run_btn.pack(side='left')
        ttk.Button(bar, text='Show lag map', command=self.show_result).pack(side='left', padx=6)
        ttk.Button(bar, text='Clear log', command=lambda: self.log.delete('1.0', 'end')).pack(side='left')
        self.status = ttk.Label(bar, text='idle'); self.status.pack(side='right')
        pf = ttk.Frame(self); pf.pack(fill='x', padx=8, pady=(4, 0))
        self.pvar = tk.DoubleVar(value=0.0)
        self.pbar = ttk.Progressbar(pf, variable=self.pvar, maximum=1.0); self.pbar.pack(side='left', fill='x', expand=True)
        self.pstage = ttk.Label(pf, text='', width=45); self.pstage.pack(side='left', padx=6)
        progress.set_callback(lambda stage, frac: self.q.put(('progress', stage, frac)))
        self.log = tk.Text(self, height=18, wrap='word', font=('TkFixedFont', 9))
        self.log.pack(fill='both', expand=True, padx=8, pady=6)
        self.after(200, self._poll)

    # ---------- widgets helpers
    def _entry(self, parent, label, var, row, col=0, width=12, browse=None, tip=None):
        ttk.Label(parent, text=label).grid(row=row, column=col, sticky='e', padx=4, pady=2)
        e = ttk.Entry(parent, textvariable=var, width=width)
        e.grid(row=row, column=col + 1, sticky='we', padx=2, pady=2)
        if browse:
            ttk.Button(parent, text='...', width=3, command=lambda: self._browse(var, browse)).grid(row=row, column=col + 2)
        if tip:
            ttk.Label(parent, text=tip, foreground='gray').grid(row=row, column=col + 3, sticky='w', padx=4)
        return e

    def _browse(self, var, kind):
        if kind == 'file':
            p = filedialog.askopenfilename(filetypes=[('NIfTI', '*.nii *.nii.gz'), ('all', '*')])
        elif kind == 'dir':
            p = filedialog.askdirectory()
        else:
            p = ''
        if p:
            var.set(p)

    def _lag_params(self, parent, row0, monkey_note=True):
        """Common lag-mapping parameters; returns dict of tk variables."""
        v = dict(TR=tk.StringVar(value='0.72'), PosiMax=tk.StringVar(value='9'), THR=tk.StringVar(value='0'),
                 FIXED=tk.StringVar(value='fixed'), Sm=tk.StringVar(value='8'), reso=tk.StringVar(value=''),
                 seed=tk.StringVar(value='hcp'), mask_pct=tk.StringVar(value='10'), lp_hz=tk.StringVar(value=''))
        self._entry(parent, 'TR (s)', v['TR'], row0, tip='repetition time')
        self._entry(parent, 'Tracking range ±PosiMax', v['PosiMax'], row0 + 1, tip='in TR; in seconds when a tracking step is set')
        self._entry(parent, 'Min peak r (THR)', v['THR'], row0 + 2, tip='cross-correlogram peaks below this are ignored (0 = accept all)')
        ttk.Label(parent, text='Tracking').grid(row=row0 + 3, column=0, sticky='e', padx=4)
        ttk.Combobox(parent, textvariable=v['FIXED'], values=['fixed', 'recursive'], width=10, state='readonly').grid(row=row0 + 3, column=1, sticky='w')
        self._entry(parent, 'Smoothing FWHM (mm)', v['Sm'], row0 + 4, tip='up to 8 mm: 8 for human, 4 for monkey; 0 = none')
        self._entry(parent, 'Tracking step (s)', v['reso'], row0 + 5, tip='empty = one TR (human); 0.5 = monkey (data resampled)')
        self._entry(parent, 'Seed mask', v['seed'], row0 + 6, browse='file', tip="'hcp' = bundled cerebral mask, a NIfTI file, or empty = whole brain")
        self._entry(parent, 'Brain mask %', v['mask_pct'], row0 + 7, tip='% of robust range of the mean image (10 human, 15 monkey)')
        self._entry(parent, 'Low-pass (Hz)', v['lp_hz'], row0 + 8, tip='empty = 0.9/(2*PosiMax s)')
        info = ttk.Label(parent, text='', foreground='#0044aa')
        info.grid(row=row0 + 9, column=0, columnspan=4, sticky='w', padx=4, pady=4)
        def upd(*_):
            try:
                P, TR = float(v['PosiMax'].get()), float(v['TR'].get())
                if v['reso'].get().strip():
                    r = float(v['reso'].get())
                    info.config(text=f'= ±{P:g} s in steps of {r:g} s (resampled from TR {TR:g} s); band-pass 0.008–{0.9 / (2 * P):.3f} Hz')
                else:
                    info.config(text=f'= ±{P:g} TR = ±{P * TR:.2f} s in steps of one TR ({TR:g} s); band-pass 0.008–{0.9 / (2 * P * TR):.3f} Hz')
            except (ValueError, ZeroDivisionError):
                info.config(text='')
        for k in ('PosiMax', 'TR', 'reso'):
            v[k].trace_add('write', upd)
        upd()
        return v

    def _lag_kwargs(self, v):
        f = lambda s: float(s) if s.strip() else None
        seed = v['seed'].get().strip() or None
        return dict(TR=float(v['TR'].get()), PosiMax=float(v['PosiMax'].get()), THR=float(v['THR'].get()),
                    FIXED=1 if v['FIXED'].get() == 'fixed' else 0, Sm=float(v['Sm'].get()), reso=f(v['reso'].get()),
                    seed_mask=seed, mask_pct=float(v['mask_pct'].get()), lp_hz=f(v['lp_hz'].get()))

    # ---------- tabs
    def _build_pipeline(self):
        t = self.tab_pipe
        left = ttk.LabelFrame(t, text='Runs  (HCP layout: <run>/<run>.nii.gz + Movement_Regressors.txt + *SBRef)')
        left.grid(row=0, column=0, sticky='nsew', padx=6, pady=4)
        self.runs = tk.Listbox(left, height=7, width=75, selectmode='extended')
        self.runs.grid(row=0, column=0, columnspan=5, sticky='we')
        ttk.Button(left, text='Add files...', command=self._add_runs).grid(row=1, column=0, sticky='w')
        ttk.Button(left, text='Remove', command=lambda: [self.runs.delete(i) for i in reversed(self.runs.curselection())]).grid(row=1, column=1, sticky='w')
        self.subj = tk.StringVar(); self.pattern = tk.StringVar(value='_REST[12]'); self.nvols = tk.StringVar()
        self._entry(left, 'Subject dir', self.subj, 2, browse='dir', width=40)
        self._entry(left, 'run folder regex', self.pattern, 3, width=16, tip='e.g. _REST[12] (HCP) or BOLD_ (monkey)')
        self._entry(left, 'volumes per run', self.nvols, 4, width=8, tip='optional filter')
        ttk.Button(left, text='Find runs in subject dir', command=self._find_runs).grid(row=5, column=0, columnspan=2, sticky='w', pady=2)
        right = ttk.LabelFrame(t, text='Lag mapping'); right.grid(row=0, column=1, rowspan=2, sticky='nsew', padx=6, pady=4)
        self.pv = self._lag_params(right, 0)
        opt = ttk.LabelFrame(t, text='Pipeline options'); opt.grid(row=1, column=0, sticky='nsew', padx=6, pady=4)
        self.results_dir = tk.StringVar(); self.workname = tk.StringVar(value='Lag_concat_scrub')
        self.only_lag = tk.BooleanVar(); self.downsample = tk.BooleanVar(value=True); self.despike = tk.BooleanVar(value=True)
        self.spike_thr = tk.StringVar(value='1.5')
        self._entry(opt, 'Output folder', self.results_dir, 0, browse='dir', width=40, tip='empty = the Results folder of the runs (existing work folder is renamed old_*)')
        self._entry(opt, 'Work folder name', self.workname, 1, width=20)
        ttk.Checkbutton(opt, text='2x down-sampling before lag mapping (subsamp2offc)', variable=self.downsample).grid(row=2, column=0, columnspan=5, sticky='w')
        ttk.Checkbutton(opt, text='Spike regressors (DVARS > thr x median)', variable=self.despike).grid(row=3, column=0, columnspan=2, sticky='w')
        ttk.Entry(opt, textvariable=self.spike_thr, width=6).grid(row=3, column=2, sticky='w')
        ttk.Checkbutton(opt, text='Lag mapping only (skip deperfusioning)', variable=self.only_lag).grid(row=4, column=0, columnspan=5, sticky='w')
        self.jobs = tk.StringVar(value='1')
        self._entry(opt, 'Parallel runs (jobs)', self.jobs, 5, width=4, tip='scrub / deperf this many runs at once (memory ~3x one run per job)')
        t.columnconfigure(0, weight=3); t.columnconfigure(1, weight=2)

    def _build_lag(self):
        t = self.tab_lag
        f = ttk.LabelFrame(t, text='Input'); f.grid(row=0, column=0, sticky='nsew', padx=6, pady=4)
        self.lag_vols = tk.StringVar(); self.lag_name = tk.StringVar(value='run1'); self.lag_range = tk.StringVar(); self.lag_cwd = tk.StringVar()
        self._entry(f, '4D file', self.lag_vols, 0, browse='file', width=50, tip='motion corrected, normalised (e.g. REST4run.nii.gz)')
        self.lag_vols.trace_add('write', lambda *_: os.path.exists(self.lag_vols.get()) and self._tr_from_header(self.lag_vols.get(), self.lv['TR']))
        self._entry(f, 'Name tag', self.lag_name, 1, width=16, tip='appended to the result folder name')
        self._entry(f, 'Time range', self.lag_range, 2, width=16, tip="MATLAB style, e.g. 1:500 or 1:2:500 (empty = all)")
        self._entry(f, 'Work folder', self.lag_cwd, 3, browse='dir', width=50, tip='empty = folder of the 4D file')
        g = ttk.LabelFrame(t, text='Lag mapping'); g.grid(row=0, column=1, sticky='nsew', padx=6, pady=4)
        self.lv = self._lag_params(g, 0)

    def _build_dep(self):
        t = self.tab_dep
        f = ttk.LabelFrame(t, text='Deperfusion runs with an existing lag map (list the runs in the order used for lag mapping)')
        f.grid(row=0, column=0, sticky='nsew', padx=6, pady=4)
        self.dep_runs = tk.Listbox(f, height=5, width=80, selectmode='extended')
        self.dep_runs.grid(row=0, column=0, columnspan=4, sticky='we')
        ttk.Button(f, text='Add files...', command=lambda: [self.dep_runs.insert('end', p) for p in filedialog.askopenfilenames(filetypes=[('NIfTI', '*.nii *.nii.gz'), ('all', '*')])]).grid(row=1, column=0, sticky='w')
        ttk.Button(f, text='Remove', command=lambda: [self.dep_runs.delete(i) for i in reversed(self.dep_runs.curselection())]).grid(row=1, column=1, sticky='w')
        self.dep_lag = tk.StringVar(); self.dep_TR = tk.StringVar(value='0.72'); self.dep_lagdir = tk.StringVar()
        self.dep_reso = tk.StringVar(); self.dep_out = tk.StringVar()
        self._entry(f, 'rLagMap.nii', self.dep_lag, 2, browse='file', width=50, tip='lag map resliced on the grid of the runs')
        self._entry(f, 'Lag folder (Seeds.mat)', self.dep_lagdir, 3, browse='dir', width=50, tip='empty = folder of rLagMap.nii')
        self._entry(f, 'TR (s)', self.dep_TR, 4, width=8)
        self._entry(f, 'Tracking step (s)', self.dep_reso, 5, width=8, tip='as used for lag mapping; empty = TR')
        self._entry(f, 'Output folder', self.dep_out, 6, browse='dir', width=50, tip='empty = current folder')

    # ---------- run selection
    def _add_runs(self):
        for p in filedialog.askopenfilenames(filetypes=[('NIfTI', '*.nii *.nii.gz'), ('all', '*')]):
            self.runs.insert('end', p)
        self._tr_from_header(self.runs.get(0) if self.runs.size() else None, self.pv['TR'])

    def _tr_from_header(self, nii, var):
        """Set the TR field from the NIfTI header and report it in the log."""
        if not nii:
            return
        from .einsteining import header_tr
        h = header_tr(nii)
        if h:
            var.set(f'{h:g}')
            self.log.insert('end', f'TR read from the header of {os.path.basename(nii)}: {h:g} s\n')

    def _find_runs(self):
        from .einsteining import find_runs
        try:
            nv = int(self.nvols.get()) if self.nvols.get().strip() else None
            for r in find_runs(self.subj.get(), self.pattern.get(), nvols=nv):
                self.runs.insert('end', r)
            self._tr_from_header(self.runs.get(0) if self.runs.size() else None, self.pv['TR'])
        except Exception as e:
            messagebox.showerror('Find runs', str(e))

    # ---------- settings (JSON)
    def _vars(self):
        d = {}
        for grp, v in [('pipeline', self.pv), ('lag', self.lv)]:
            for k, var in v.items():
                d[f'{grp}.{k}'] = var
        for k in ['subj', 'pattern', 'nvols', 'results_dir', 'workname', 'only_lag', 'downsample', 'despike', 'spike_thr',
                  'jobs', 'lag_vols', 'lag_name', 'lag_range', 'lag_cwd', 'dep_lag', 'dep_TR', 'dep_lagdir', 'dep_reso', 'dep_out']:
            d[k] = getattr(self, k)
        return d

    def settings(self):
        d = {k: v.get() for k, v in self._vars().items()}
        d['runs'] = list(self.runs.get(0, 'end'))
        d['dep_runs'] = list(self.dep_runs.get(0, 'end'))
        d['tab'] = self.nb.index(self.nb.select())
        return d

    def apply_settings(self, d):
        for k, v in self._vars().items():
            if k in d:
                v.set(d[k])
        if 'runs' in d:
            self.runs.delete(0, 'end')
            for r in d['runs']:
                self.runs.insert('end', r)
        if 'dep_runs' in d:
            self.dep_runs.delete(0, 'end')
            for r in d['dep_runs']:
                self.dep_runs.insert('end', r)
        if 'tab' in d:
            self.nb.select(d['tab'])

    def save_settings(self, path=None):
        path = path or filedialog.asksaveasfilename(defaultextension='.json', filetypes=[('JSON', '*.json')])
        if path:
            with open(path, 'w') as f:
                json.dump(self.settings(), f, indent=1)
            self.log.insert('end', f'settings saved to {path}\n')

    def load_settings(self, path=None):
        path = path or filedialog.askopenfilename(filetypes=[('JSON', '*.json')])
        if path:
            with open(path) as f:
                self.apply_settings(json.load(f))
            self.log.insert('end', f'settings loaded from {path}\n')

    # ---------- execution
    def run(self):
        if self.worker and self.worker.is_alive():
            messagebox.showinfo('Busy', 'A job is already running.')
            return
        tab = self.nb.index(self.nb.select())
        try:
            job = [self._job_pipeline, self._job_lag, self._job_dep][tab]()
        except Exception as e:
            messagebox.showerror('Parameters', str(e))
            return
        self.result = None
        self.pvar.set(0.0); self.pstage.config(text='')
        self.status.config(text='running...'); self.run_btn.state(['disabled'])
        self.worker = threading.Thread(target=self._work, args=(job,), daemon=True)
        self.worker.start()

    def _job_pipeline(self):
        runs = list(self.runs.get(0, 'end'))
        if not runs:
            raise ValueError('add at least one run')
        kw = self._lag_kwargs(self.pv)
        from .einsteining import einsteining
        def job():
            d = einsteining(runs, kw['TR'], kw['PosiMax'], kw['THR'], kw['FIXED'], kw['Sm'], self.only_lag.get(),
                            self.downsample.get(), kw['reso'], kw['seed_mask'], kw['mask_pct'], kw['lp_hz'],
                            self.workname.get(), self.results_dir.get().strip() or None,
                            spike_thr=float(self.spike_thr.get()) if self.despike.get() else None, jobs=int(self.jobs.get() or 1))
            return dict(lagdir=d, underlay=os.path.join(os.path.dirname(d), 'Tmean.nii'))
        return job

    def _job_lag(self):
        vols = self.lag_vols.get().strip()
        if not os.path.exists(vols):
            raise ValueError('4D file not found')
        kw = self._lag_kwargs(self.lv)
        cwd = self.lag_cwd.get().strip() or os.path.dirname(os.path.abspath(vols))
        rng = self.lag_range.get().strip() or None
        from .lag4d import lag4d
        def job():
            d = lag4d(self.lag_name.get(), kw['TR'], vols, kw['PosiMax'], kw['THR'], kw['FIXED'], kw['Sm'], rng,
                      kw['reso'], kw['seed_mask'], kw['mask_pct'], kw['lp_hz'], cwd=cwd)
            return dict(lagdir=d, underlay=os.path.join(cwd, 'Tmean.nii'))
        return job

    def _job_dep(self):
        from .deperf import deperf
        runs, lag = list(self.dep_runs.get(0, 'end')), self.dep_lag.get().strip()
        if not runs or not os.path.exists(lag):
            raise ValueError('add the run files and the lag map')
        reso = float(self.dep_reso.get()) if self.dep_reso.get().strip() else None
        TR, lagdir, out = float(self.dep_TR.get()), self.dep_lagdir.get().strip() or None, self.dep_out.get().strip() or '.'
        def job():
            for i, run in enumerate(runs):
                print(f'Deperfusioning {os.path.basename(run)} (run {i + 1} of {len(runs)})', flush=True)
                f = deperf(run, lag, TR, i + 1, len(runs), lagdir, reso, out,
                           span=(i / len(runs), (i + 1) / len(runs)))
                print('written', f)
            return dict(lagdir=None)
        return job

    def _work(self, job):
        old = sys.stdout, sys.stderr
        sys.stdout = sys.stderr = _QueueWriter(self.q)
        try:
            self.result = job()
            print('\n*** done ***')
        except Exception:
            print(traceback.format_exc())
        finally:
            sys.stdout, sys.stderr = old
            self.q.put(None)

    def _poll(self):
        try:
            while True:
                s = self.q.get_nowait()
                if s is None:
                    self.status.config(text='finished'); self.run_btn.state(['!disabled'])
                    if self.result and self.result.get('lagdir'):
                        self.show_result()
                elif isinstance(s, tuple):
                    self.pvar.set(s[2]); self.pstage.config(text=f'{s[1]}  ({100 * s[2]:.0f} %)')
                else:
                    self.log.insert('end', s); self.log.see('end')
        except queue.Empty:
            pass
        self.after(200, self._poll)

    def show_result(self, lagdir=None):
        lagdir = lagdir or (self.result or {}).get('lagdir')
        if not lagdir:
            p = filedialog.askopenfilename(title='LagMap.nii', filetypes=[('NIfTI', '*.nii *.nii.gz')])
            if not p:
                return
            lagdir = os.path.dirname(p)
        under = (self.result or {}).get('underlay') or os.path.join(os.path.dirname(lagdir), 'Tmean.nii')
        try:
            ResultWindow(self, lagdir, under if os.path.exists(under) else None)
        except Exception as e:
            messagebox.showerror('Results', str(e))


class ResultWindow(tk.Toplevel):
    """Montage, slice viewer and lag-structure (sLFO) plot of one lag-map folder."""
    def __init__(self, master, lagdir, underlay=None):
        super().__init__(master)
        self.title(lagdir)
        self.lagdir, self.under = lagdir, underlay
        self.lagmap = os.path.join(lagdir, 'LagMap.nii')
        self.lim = tk.DoubleVar(value=4.0)
        top = ttk.Frame(self); top.pack(fill='x', padx=6, pady=4)
        ttk.Label(top, text='colour range ± s').pack(side='left')
        ttk.Spinbox(top, from_=0.5, to=20, increment=0.5, textvariable=self.lim, width=5, command=self.refresh).pack(side='left', padx=4)
        nb = ttk.Notebook(self); nb.pack(fill='both', expand=True)
        self.nb = nb
        # montage
        self.t_mont = ttk.Frame(nb); nb.add(self.t_mont, text=' Montage ')
        self.l_mont = ttk.Label(self.t_mont); self.l_mont.pack()
        # slice viewer
        self.t_slice = ttk.Frame(nb); nb.add(self.t_slice, text=' Slice viewer ')
        c = ttk.Frame(self.t_slice); c.pack(fill='x')
        self.axis = tk.IntVar(value=2)
        for k, nm in enumerate(['sagittal', 'coronal', 'axial']):
            ttk.Radiobutton(c, text=nm, variable=self.axis, value=k, command=self._axis_changed).pack(side='left', padx=4)
        self.idx = tk.IntVar(value=0)
        self.scale = ttk.Scale(c, from_=0, to=1, variable=self.idx, command=lambda _v: self.refresh_slice()); self.scale.pack(side='left', fill='x', expand=True, padx=8)
        self.l_slice = ttk.Label(self.t_slice); self.l_slice.pack()
        # lag structure
        self.t_lag = ttk.Frame(nb); nb.add(self.t_lag, text=' Lag structure (sLFO) ')
        c = ttk.Frame(self.t_lag); c.pack(fill='x')
        self.shifted = tk.BooleanVar(value=True)
        ttk.Checkbutton(c, text='time-shifted to lag (regressors)', variable=self.shifted, command=self.refresh_lag).pack(side='left', padx=4)
        ttk.Label(c, text='amplitude').pack(side='left')
        self.amp = tk.StringVar(value='normalised')
        cb = ttk.Combobox(c, textvariable=self.amp, values=['normalised', 'scaled', 'data'], width=10, state='readonly'); cb.pack(side='left', padx=4)
        cb.bind('<<ComboboxSelected>>', lambda _e: self.refresh_lag())
        ttk.Label(c, text='window start (s)').pack(side='left')
        self.t0 = tk.DoubleVar(value=0)
        self.t0scale = ttk.Scale(c, from_=0, to=1, variable=self.t0, command=lambda _v: self.refresh_lag()); self.t0scale.pack(side='left', fill='x', expand=True, padx=8)
        self.t0label = ttk.Label(c, text='0 s', width=7); self.t0label.pack(side='left')
        ttk.Label(c, text='length (min)').pack(side='left')
        self.nmin = tk.DoubleVar(value=10)
        ttk.Spinbox(c, from_=0.5, to=300, increment=0.5, textvariable=self.nmin, width=6, command=self.refresh_lag).pack(side='left', padx=4)
        self.l_lag = ttk.Label(self.t_lag); self.l_lag.pack()
        self._nsl = {}
        self.refresh()

    def _show(self, label, png):
        img = tk.PhotoImage(file=png)
        label.configure(image=img); label.image = img

    def _axis_changed(self):
        import nibabel as nib
        n = nib.load(self.lagmap).shape[self.axis.get()]
        self.scale.configure(to=n - 1); self.idx.set(n // 2); self.refresh_slice()

    def refresh(self):
        from .viewer import lagmap_montage
        self._show(self.l_mont, lagmap_montage(self.lagmap, self.under, lim=self.lim.get()))
        self._axis_changed()
        try:
            from .deperf import load_seeds
            from .viewer import _lag_step
            self.step = _lag_step(self.lagdir)
            self.t0scale.configure(to=max(0, (load_seeds(self.lagdir).shape[0] - 20) * self.step))
            self.refresh_lag()
        except Exception as e:
            self.l_lag.configure(text=f'no Seeds in {self.lagdir}: {e}')

    def refresh_slice(self):
        from .viewer import slice_image
        png, _ = slice_image(self.lagmap, self.under, self.axis.get(), int(self.idx.get()), self.lim.get())
        self._show(self.l_slice, png)

    def refresh_lag(self):
        from .viewer import lag_structure_plot
        self.t0label.configure(text=f'{self.t0.get():.0f} s')
        png = lag_structure_plot(self.lagdir, os.path.join(self.lagdir, '_lagstructure.png'), int(self.t0.get() / self.step),
                                 max(2, int(self.nmin.get() * 60 / self.step)), self.lim.get(), self.shifted.get(), amplitude=self.amp.get())
        self._show(self.l_lag, png)


def main(argv=None):
    """``boldlag-gui [settings.json]``"""
    argv = sys.argv[1:] if argv is None else argv
    app = App()
    if argv:
        app.load_settings(argv[0])
    app.mainloop()


if __name__ == '__main__':
    main()
