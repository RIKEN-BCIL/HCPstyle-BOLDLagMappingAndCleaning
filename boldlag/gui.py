"""Tkinter front end for boldlag (``boldlag-gui`` / ``python -m boldlag.gui``).

Three tabs -- the whole Einsteining pipeline, lag mapping of a single 4D file, and
deperfusioning of one run -- run in a background thread with their console output
shown in the log pane; the resulting lag map is displayed as a montage when done.
No dependencies beyond the package itself (matplotlib for the montage).
"""
import os, sys, glob, threading, queue, traceback
import tkinter as tk
from tkinter import ttk, filedialog, messagebox


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
        v = dict(TR=tk.StringVar(value='0.72'), PosiMax=tk.StringVar(value='9'), THR=tk.StringVar(value='0.2'),
                 FIXED=tk.StringVar(value='fixed'), Sm=tk.StringVar(value='8'), reso=tk.StringVar(value=''),
                 seed=tk.StringVar(value='hcp'), mask_pct=tk.StringVar(value='10'), lp_hz=tk.StringVar(value=''))
        self._entry(parent, 'TR (s)', v['TR'], row0, tip='repetition time')
        self._entry(parent, 'PosiMax', v['PosiMax'], row0 + 1, tip='tracking range +-PosiMax, in TR (or in s when a tracking step is given)')
        self._entry(parent, 'Min peak r (THR)', v['THR'], row0 + 2, tip='cross-correlogram peaks below this are ignored')
        ttk.Label(parent, text='Tracking').grid(row=row0 + 3, column=0, sticky='e', padx=4)
        ttk.Combobox(parent, textvariable=v['FIXED'], values=['fixed', 'recursive'], width=10, state='readonly').grid(row=row0 + 3, column=1, sticky='w')
        self._entry(parent, 'Smoothing FWHM (mm)', v['Sm'], row0 + 4, tip='8 for human 2 mm data, 4 for monkey; 0 = none')
        self._entry(parent, 'Tracking step (s)', v['reso'], row0 + 5, tip='empty = one TR (human); 0.5 = monkey (data resampled)')
        self._entry(parent, 'Seed mask', v['seed'], row0 + 6, browse='file', tip="'hcp' = bundled cerebral mask, a NIfTI file, or empty = whole brain")
        self._entry(parent, 'Brain mask %', v['mask_pct'], row0 + 7, tip='% of robust range of the mean image (10 human, 15 monkey)')
        self._entry(parent, 'Low-pass (Hz)', v['lp_hz'], row0 + 8, tip='empty = 0.9/(2*PosiMax s)')
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
        t.columnconfigure(0, weight=3); t.columnconfigure(1, weight=2)

    def _build_lag(self):
        t = self.tab_lag
        f = ttk.LabelFrame(t, text='Input'); f.grid(row=0, column=0, sticky='nsew', padx=6, pady=4)
        self.lag_vols = tk.StringVar(); self.lag_name = tk.StringVar(value='run1'); self.lag_range = tk.StringVar(); self.lag_cwd = tk.StringVar()
        self._entry(f, '4D file', self.lag_vols, 0, browse='file', width=50, tip='motion corrected, normalised (e.g. REST4run.nii.gz)')
        self._entry(f, 'Name tag', self.lag_name, 1, width=16, tip='appended to the result folder name')
        self._entry(f, 'Time range', self.lag_range, 2, width=16, tip="MATLAB style, e.g. 1:500 or 1:2:500 (empty = all)")
        self._entry(f, 'Work folder', self.lag_cwd, 3, browse='dir', width=50, tip='empty = folder of the 4D file')
        g = ttk.LabelFrame(t, text='Lag mapping'); g.grid(row=0, column=1, sticky='nsew', padx=6, pady=4)
        self.lv = self._lag_params(g, 0)

    def _build_dep(self):
        t = self.tab_dep
        f = ttk.LabelFrame(t, text='Deperfusion one run with an existing lag map'); f.grid(row=0, column=0, sticky='nsew', padx=6, pady=4)
        self.dep_vols = tk.StringVar(); self.dep_lag = tk.StringVar(); self.dep_TR = tk.StringVar(value='0.72')
        self.dep_sec = tk.StringVar(value='1'); self.dep_n = tk.StringVar(value='1'); self.dep_lagdir = tk.StringVar()
        self.dep_reso = tk.StringVar(); self.dep_out = tk.StringVar()
        self._entry(f, 'Original 4D run', self.dep_vols, 0, browse='file', width=50)
        self._entry(f, 'rLagMap.nii', self.dep_lag, 1, browse='file', width=50, tip='lag map on the grid of the run')
        self._entry(f, 'TR (s)', self.dep_TR, 2, width=8)
        self._entry(f, 'Run number', self.dep_sec, 3, width=8, tip='position of this run in the concatenation used for lag mapping')
        self._entry(f, 'Number of runs', self.dep_n, 4, width=8)
        self._entry(f, 'Lag folder (Seeds.mat)', self.dep_lagdir, 5, browse='dir', width=50, tip='empty = folder of rLagMap.nii')
        self._entry(f, 'Tracking step (s)', self.dep_reso, 6, width=8, tip='as used for lag mapping; empty = TR')
        self._entry(f, 'Output folder', self.dep_out, 7, browse='dir', width=50, tip='empty = current folder')

    # ---------- run selection
    def _add_runs(self):
        for p in filedialog.askopenfilenames(filetypes=[('NIfTI', '*.nii *.nii.gz'), ('all', '*')]):
            self.runs.insert('end', p)

    def _find_runs(self):
        from .einsteining import find_runs
        try:
            nv = int(self.nvols.get()) if self.nvols.get().strip() else None
            for r in find_runs(self.subj.get(), self.pattern.get(), nvols=nv):
                self.runs.insert('end', r)
        except Exception as e:
            messagebox.showerror('Find runs', str(e))

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
                            spike_thr=float(self.spike_thr.get()) if self.despike.get() else None)
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
        vols, lag = self.dep_vols.get().strip(), self.dep_lag.get().strip()
        if not (os.path.exists(vols) and os.path.exists(lag)):
            raise ValueError('run file or lag map not found')
        reso = float(self.dep_reso.get()) if self.dep_reso.get().strip() else None
        def job():
            out = deperf(vols, lag, float(self.dep_TR.get()), int(self.dep_sec.get()), int(self.dep_n.get()),
                         self.dep_lagdir.get().strip() or None, reso, self.dep_out.get().strip() or '.')
            print('written', out)
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
            lagmap = p
            under = os.path.join(os.path.dirname(os.path.dirname(p)), 'Tmean.nii')
        else:
            lagmap = os.path.join(lagdir, 'LagMap.nii')
            under = (self.result or {}).get('underlay') or os.path.join(os.path.dirname(lagdir), 'Tmean.nii')
        try:
            from .viewer import lagmap_montage
            png = lagmap_montage(lagmap, under if under and os.path.exists(under) else None, lim=4.0)
        except Exception as e:
            messagebox.showerror('Montage', str(e))
            return
        win = tk.Toplevel(self); win.title(lagmap)
        img = tk.PhotoImage(file=png)
        lbl = ttk.Label(win, image=img); lbl.image = img; lbl.pack()
        ttk.Label(win, text=png).pack()


def main():
    App().mainloop()


if __name__ == '__main__':
    main()
