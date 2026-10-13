#!/usr/bin/env python3
"""
stress_state_gui.py -- PILOT CRSS temperature-sweep workflow, from one window.

    python3 stress_state_gui.py

This is the graphical front end for the stress-state / CRSS half of PILOT.  It
drives the same engine (crss_core.py) as the command-line tool
(crss_search_alloy.py), so you can switch between them freely.

Everything for one search is kept in a study subfolder,
    studies/<study name>/
so the working directory stays uncluttered.  From here you can:

  * set the temperature range/step, stress bracket, number of rounds, chi/axis;
  * visualise and choose the stress state (launches the MRSSP explorer);
  * link the per-temperature alloy USFE files (from the PILOT lattice/USFE tool);
  * init / generate / update / extend / undo the bisection search;
  * grade moved / not-moved for each temperature in a dialog (after inspecting
    the dumps in OVITO);
  * fetch the dump files from a cluster, and delete the local copies afterwards;
  * plot CRSS vs temperature.

The heavy compute (PFDD) runs on the cluster via the SLURM array in slurm/; this
tool prepares its inputs and interprets its results.
"""

import json
import os
import subprocess
import sys
import threading

import crss_core as cc
import pilot_banner

HERE = os.path.dirname(os.path.abspath(__file__))
STUDIES = os.path.join(HERE, "studies")

# slip system -> (cases/jobs tag, dump-name pattern)
SLIP_DEFAULTS = {
    "{110}": ("Ta10W_110", "Ta10W_screw_110_{T}K_{tag}.dump"),
    "{112}": ("Ta10W_112", "Ta10W_screw_112_{T}K_{tag}.dump"),
}

# Remembered cluster fields (for the "remember" checkboxes) persist here so they
# survive between dialog opens and between sessions. A field is remembered iff
# its key is present in this file.
PREFS_PATH = os.path.expanduser("~/.pilot_stress_state.json")


def load_prefs():
    try:
        with open(PREFS_PATH) as f:
            return json.load(f)
    except Exception:
        return {}


def save_prefs(d):
    try:
        with open(PREFS_PATH, "w") as f:
            json.dump(d, f, indent=2)
    except Exception:
        pass


def run_gui():
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox, simpledialog
    from tkinter.scrolledtext import ScrolledText

    root = tk.Tk()
    root.title("PILOT  -  CRSS temperature sweep")
    root.geometry("900x860")
    root.lift(); root.focus_force()
    root.attributes("-topmost", True)
    root.after(200, lambda: root.attributes("-topmost", False))

    pad = {"padx": 6, "pady": 3}
    V = {}

    def sv(key, val):
        V[key] = tk.StringVar(value=str(val)); return V[key]

    # ---------- 1. Study & slip system ----------
    f1 = ttk.LabelFrame(root, text="1.  Study & slip system")
    f1.pack(fill="x", **pad)
    ttk.Label(f1, text="Study name").grid(row=0, column=0, sticky="w", **pad)
    ttk.Entry(f1, textvariable=sv("study", "study1"), width=22).grid(
        row=0, column=1, sticky="w", **pad)
    ttk.Label(f1, text="Slip system").grid(row=0, column=2, sticky="w", **pad)
    V["slip"] = tk.StringVar(value="{110}")
    slip_box = ttk.Combobox(f1, textvariable=V["slip"], values=list(SLIP_DEFAULTS),
                            state="readonly", width=8)
    slip_box.grid(row=0, column=3, sticky="w", **pad)
    ttk.Label(f1, text="Cases/jobs tag").grid(row=1, column=0, sticky="w", **pad)
    ttk.Entry(f1, textvariable=sv("tag", "Ta10W_110"), width=22).grid(
        row=1, column=1, sticky="w", **pad)
    ttk.Label(f1, text="(sbatch uses <tag>_jobs)", foreground="#555").grid(
        row=1, column=2, columnspan=2, sticky="w", **pad)

    # ---------- 2. Temperatures ----------
    f2 = ttk.LabelFrame(root, text="2.  Temperatures (K)")
    f2.pack(fill="x", **pad)
    for i, (key, lab, dflt) in enumerate([("t_min", "min", 0), ("t_max", "max", 1400),
                                          ("t_step", "step", 100)]):
        ttk.Label(f2, text=lab).grid(row=0, column=2 * i, sticky="e", **pad)
        ttk.Entry(f2, textvariable=sv(key, dflt), width=8).grid(
            row=0, column=2 * i + 1, sticky="w", **pad)

    # ---------- 3. Initial stress bracket ----------
    f3 = ttk.LabelFrame(root, text="3.  Initial stress bracket for binary search")
    f3.pack(fill="x", **pad)
    grid3 = [("tau_lo", "bracket lo (mu)", 0.0), ("tau_hi", "bracket hi (mu)", 0.08),
             ("n_rounds", "# of search rounds", 6),
             ("tau_extend_step", "extend step (mu)", 0.0),
             ("n_extend_rounds", "extend rounds", 0),
             ("seed_max", "max PFDD thermal seed", 1000)]
    for i, (key, lab, dflt) in enumerate(grid3):
        r, c = divmod(i, 2)
        ttk.Label(f3, text=lab).grid(row=r, column=2 * c, sticky="e", **pad)
        ttk.Entry(f3, textvariable=sv(key, dflt), width=10).grid(
            row=r, column=2 * c + 1, sticky="w", **pad)
    # auto-computed final resolution = initial bracket width * (1/2)^rounds
    ttk.Label(f3, text="final CRSS resolution (mu)").grid(
        row=3, column=0, sticky="e", **pad)
    V["crss_res"] = tk.StringVar(value="")
    ttk.Entry(f3, textvariable=V["crss_res"], width=12,
              state="readonly").grid(row=3, column=1, sticky="w", **pad)

    def _update_res(*_):
        try:
            w = abs(float(V["tau_hi"].get()) - float(V["tau_lo"].get()))
            n = int(float(V["n_rounds"].get()))
            V["crss_res"].set(f"{w * (0.5 ** n):.3e}")
        except (ValueError, ZeroDivisionError):
            V["crss_res"].set("")
    for _k in ("tau_lo", "tau_hi", "n_rounds"):
        V[_k].trace_add("write", _update_res)
    _update_res()
    ttk.Label(f3, text="resolution = initial bracket width x (1/2)^rounds",
              foreground="#555").grid(row=4, column=0, columnspan=4,
                                      sticky="w", **pad)

    # ---------- 4. Stress state ----------
    f4 = ttk.LabelFrame(root, text="4.  Stress state")
    f4.pack(fill="x", **pad)
    V["mode"] = tk.StringVar(value="shear")
    ttk.Radiobutton(f4, text="Shear couple (MRSSP)", value="shear",
                    variable=V["mode"]).grid(row=0, column=0, sticky="w", **pad)
    ttk.Radiobutton(f4, text="Uniaxial", value="uniaxial",
                    variable=V["mode"]).grid(row=0, column=1, sticky="w", **pad)
    ttk.Label(f4, text="chi (deg)").grid(row=1, column=0, sticky="e", **pad)
    ttk.Entry(f4, textvariable=sv("chi_deg", 0.0), width=8).grid(
        row=1, column=1, sticky="w", **pad)
    ttk.Label(f4, text="uniaxial axis").grid(row=1, column=2, sticky="e", **pad)
    ttk.Entry(f4, textvariable=sv("axis", "1, 0, 0"), width=12).grid(
        row=1, column=3, sticky="w", **pad)
    ttk.Button(f4, text="Visualize / set stress...",
               command=lambda: open_mrssp()).grid(row=0, column=2, columnspan=2,
                                                   sticky="w", **pad)

    # ---------- 5. Alloy USFE link ----------
    f5 = ttk.LabelFrame(root, text="5.  Alloy USFE link  (per-temperature USFE "
                                   "file written into each cases-file line)")
    f5.pack(fill="x", **pad)
    ttk.Label(f5, text="Names each temperature's USFE file (from the PILOT "
                       "lattice/USFE tool) so the cases file, and the SLURM job, "
                       "point PFDD at the right composition. Optionally pick the "
                       "folder to verify the files exist.",
              foreground="#555", wraplength=600).grid(
        row=0, column=0, columnspan=4, sticky="w", **pad)
    ttk.Label(f5, text="USFE name pattern").grid(row=1, column=0, sticky="w", **pad)
    ttk.Entry(f5, textvariable=sv("usfe_pattern", "USFE_{T}"), width=16).grid(
        row=1, column=1, sticky="w", **pad)
    ttk.Label(f5, text="({T} -> temperature)", foreground="#555").grid(
        row=1, column=2, sticky="w", **pad)
    ttk.Label(f5, text="USFE folder (opt.)").grid(row=2, column=0, sticky="w", **pad)
    ttk.Entry(f5, textvariable=sv("usfe_dir", ""), width=30).grid(
        row=2, column=1, columnspan=2, sticky="w", **pad)
    ttk.Button(f5, text="Browse...",
               command=lambda: _pick_dir(V["usfe_dir"])).grid(row=2, column=3, **pad)
    ttk.Button(f5, text="Check files exist",
               command=lambda: check_usfe()).grid(row=1, column=3, **pad)

    # ---------- cluster / dumps ----------
    f6 = ttk.LabelFrame(root, text="6.  Cluster output (for dump paths & fetch)")
    f6.pack(fill="x", **pad)
    ttk.Label(f6, text="cluster dump base").grid(row=0, column=0, sticky="w", **pad)
    ttk.Entry(f6, textvariable=sv("out_base", "OUTPUT/dumps"), width=40).grid(
        row=0, column=1, columnspan=3, sticky="w", **pad)

    V["dumps_name"] = tk.StringVar(value="dumps_local")

    # ---------- Actions ----------
    fa = ttk.LabelFrame(root, text="Actions")
    fa.pack(fill="x", **pad)
    btns = [("Initialize search", lambda: do(cc.init, force_confirm=True)),
            ("Generate round", lambda: do(cc.generate)),
            ("Update (grade)...", lambda: open_update()),
            ("Extend", lambda: do_extend()),
            ("Undo", lambda: do(cc.undo)),
            ("Report", lambda: show_report()),
            ("Send to cluster...", lambda: open_send()),
            ("Fetch dumps...", lambda: open_fetch()),
            ("Delete local dumps", lambda: delete_dumps()),
            ("Plot CRSS vs T", lambda: plot_crss())]
    for i, (lab, cmd) in enumerate(btns):
        ttk.Button(fa, text=lab, command=cmd, width=18).grid(
            row=i // 3, column=i % 3, **pad)

    # ---------- Log ----------
    log = ScrolledText(root, height=15, state="disabled", wrap="word",
                       font="TkFixedFont")
    log.pack(fill="both", expand=True, **pad)

    def append(msg):
        log.configure(state="normal")
        log.insert("end", msg + "\n")
        log.see("end")
        log.configure(state="disabled")

    # ---------- helpers ----------
    def _pick_dir(var):
        d = filedialog.askdirectory()
        if d:
            var.set(d)

    def study_dir():
        return os.path.join(STUDIES, V["study"].get().strip() or "study1")

    def build_cfg():
        def i(k):
            return int(float(V[k].get()))

        def f(k):
            return float(V[k].get())
        axis = (1.0, 0.0, 0.0)
        try:
            parts = [float(x) for x in V["axis"].get().replace(",", " ").split()]
            if len(parts) == 3:
                axis = tuple(parts)
        except ValueError:
            pass
        tag, dump = SLIP_DEFAULTS.get(V["slip"].get(), ("run", "dump_{T}K_{tag}.dump"))
        return cc.SearchConfig(
            work_dir=study_dir(),
            t_min=i("t_min"), t_max=i("t_max"), t_step=i("t_step"),
            tau_lo=f("tau_lo"), tau_hi=f("tau_hi"), n_rounds=i("n_rounds"),
            tau_extend_step=f("tau_extend_step"), n_extend_rounds=i("n_extend_rounds"),
            stress_mode=V["mode"].get(), chi_deg=f("chi_deg"), axis=axis,
            seed_min=1, seed_max=i("seed_max"),
            slip_system=V["slip"].get(), tag=V["tag"].get().strip() or tag,
            usfe_pattern=V["usfe_pattern"].get().strip() or "USFE_{T}",
            usfe_dir=V["usfe_dir"].get().strip(),
            out_base=V["out_base"].get().strip() or "OUTPUT/dumps",
            dump_name=dump)

    def do(fn, force_confirm=False):
        try:
            cfg = build_cfg()
        except ValueError as e:
            messagebox.showerror("Check inputs", str(e)); return
        try:
            if fn is cc.init:
                if os.path.exists(cfg.state_path()) and not messagebox.askyesno(
                        "Restart search?",
                        "A search already exists in this study. Overwrite it?"):
                    return
                fn(cfg, log=append, force=True)
            else:
                fn(cfg, log=append)
        except cc.SearchError as e:
            messagebox.showerror("Search", str(e)); append("! " + str(e))
        except Exception as e:            # noqa: BLE001
            messagebox.showerror("Error", str(e)); append("! " + str(e))

    def on_slip(*_):
        tag, _dump = SLIP_DEFAULTS[V["slip"].get()]
        V["tag"].set(tag)
    V["slip"].trace_add("write", on_slip)

    def open_mrssp():
        import mrssp_gui
        try:
            chi = float(V["chi_deg"].get())
        except ValueError:
            chi = 0.0
        axis = (1.0, 0.0, 0.0)
        try:
            parts = [float(x) for x in V["axis"].get().replace(",", " ").split()]
            if len(parts) == 3:
                axis = tuple(parts)
        except ValueError:
            pass

        def apply_stress(payload):
            V["mode"].set(payload["mode"])
            V["chi_deg"].set(f"{payload['chi_deg']:.0f}")
            V["axis"].set(", ".join(str(x) for x in payload["axis"]))
            append(f"Stress set from explorer: mode={payload['mode']}, "
                   f"chi={payload['chi_deg']:+.0f} deg, axis={payload['axis']}")
        mrssp_gui.run_mrssp(master=root, on_apply=apply_stress,
                            init={"mode": V["mode"].get(), "chi_deg": chi,
                                  "axis": axis})

    def check_usfe():
        cfg = build_cfg()
        if not cfg.usfe_dir:
            append("No USFE folder set -- cannot check existence."); return
        missing = []
        for T in cfg.temps():
            name = cfg.usfe_pattern.format(T=T) if "{T}" in cfg.usfe_pattern \
                else cfg.usfe_pattern
            if not os.path.exists(os.path.join(cfg.usfe_dir, name)):
                missing.append(name)
        if missing:
            append(f"MISSING {len(missing)} USFE file(s): {', '.join(missing)}")
        else:
            append(f"All {len(cfg.temps())} USFE files present in {cfg.usfe_dir}.")

    def show_report():
        try:
            append("\n" + cc.report_text(build_cfg()))
        except cc.SearchError as e:
            messagebox.showerror("Report", str(e))

    def open_update():
        try:
            cfg = build_cfg()
            pending = cc.list_pending(cfg)
        except cc.SearchError as e:
            messagebox.showerror("Update", str(e)); return
        if not pending:
            messagebox.showinfo("Update", "No pending results. Generate a round "
                                          "and run the jobs first.")
            return
        win = tk.Toplevel(root)
        win.title("Update -- did the dislocation move?")
        win.geometry("640x520")
        ttk.Label(win, text="Inspect each dump (e.g. in OVITO), then mark whether "
                            "the dislocation glided.", wraplength=600,
                  foreground="#333").pack(**pad)
        canv = tk.Canvas(win, borderwidth=0)
        frame = ttk.Frame(canv)
        vsb = ttk.Scrollbar(win, orient="vertical", command=canv.yview)
        canv.configure(yscrollcommand=vsb.set)
        vsb.pack(side="right", fill="y")
        canv.pack(side="left", fill="both", expand=True, **pad)
        canv.create_window((0, 0), window=frame, anchor="nw")
        frame.bind("<Configure>",
                   lambda _e: canv.configure(scrollregion=canv.bbox("all")))

        choices = {}
        for row in pending:
            rf = ttk.Frame(frame); rf.pack(fill="x", pady=2)
            ttk.Label(rf, text=f"T={row['T']}K  tau={cc.tau_str(row['tau'])}",
                      width=22).pack(side="left")
            var = tk.StringVar(value="")
            choices[row["T"]] = var
            ttk.Radiobutton(rf, text="moved", value="m", variable=var).pack(side="left")
            ttk.Radiobutton(rf, text="did not move", value="n",
                            variable=var).pack(side="left")
            ttk.Label(rf, text=os.path.basename(row["dump"]),
                      foreground="#777").pack(side="left", padx=6)

        def apply_all():
            results = {}
            for T, var in choices.items():
                if var.get() not in ("m", "n"):
                    messagebox.showwarning("Update", f"T={T}K is not marked yet.")
                    return
                results[T] = (var.get() == "m")
            try:
                cc.record_results(cfg, results, log=append)
            except cc.SearchError as e:
                messagebox.showerror("Update", str(e)); return
            win.destroy()
        ttk.Button(win, text="Apply results", command=apply_all).pack(**pad)

    def open_fetch():
        cfg = build_cfg()
        win = tk.Toplevel(root); win.title("Fetch dumps from cluster")
        win.geometry("660x330")
        prefs = load_prefs()
        fv = {}
        rem = {}
        rows = [("REMOTE_USER", "cluster username", ""),
                ("REMOTE_HOST", "cluster hostname", ""),
                ("REMOTE_ROOT", "remote dump root (contains <T>K/ dirs)", ""),
                ("FILE_PREFIX", "dump filename prefix",
                 cfg.dump_name.split("_{T}K")[0])]
        for i, (k, lab, dflt) in enumerate(rows):
            ttk.Label(win, text=lab).grid(row=i, column=0, sticky="w", **pad)
            fv[k] = tk.StringVar(value=prefs.get(k, dflt))
            ttk.Entry(win, textvariable=fv[k], width=32).grid(row=i, column=1, **pad)
            rem[k] = tk.BooleanVar(value=(k in prefs))
            ttk.Checkbutton(win, text="remember", variable=rem[k]).grid(
                row=i, column=2, sticky="w", **pad)
        ttk.Label(win, text="local dest folder name").grid(
            row=len(rows), column=0, sticky="w", **pad)
        ttk.Entry(win, textvariable=V["dumps_name"], width=40).grid(
            row=len(rows), column=1, **pad)
        ttk.Label(win, text=f"(created inside {cfg.work_dir})",
                  foreground="#555").grid(row=len(rows) + 1, column=0,
                                          columnspan=2, sticky="w", **pad)

        def run_fetch():
            pr = load_prefs()
            for k in fv:
                if rem[k].get():
                    pr[k] = fv[k].get().strip()
                else:
                    pr.pop(k, None)
            save_prefs(pr)
            local_dest = os.path.join(
                cfg.work_dir, V["dumps_name"].get().strip() or "dumps_local")
            env = dict(os.environ)
            for k in fv:
                env[k] = fv[k].get().strip()
            env["CASES"] = cfg.cases_path()
            env["LOCAL_DEST"] = local_dest
            script = os.path.join(HERE, "fetch_dumps.sh")
            append(f"Fetching dumps -> {local_dest} ...")
            win.destroy()

            def worker():
                try:
                    proc = subprocess.Popen(["bash", script], env=env,
                                            stdout=subprocess.PIPE,
                                            stderr=subprocess.STDOUT, text=True)
                    for line in proc.stdout:
                        append(line.rstrip())
                    proc.wait()
                    append(f"fetch_dumps.sh exited with code {proc.returncode}")
                except Exception as e:      # noqa: BLE001
                    append("! fetch failed: " + str(e))
            threading.Thread(target=worker, daemon=True).start()
        ttk.Button(win, text="Fetch", command=run_fetch).grid(
            row=len(rows) + 2, column=0, columnspan=2, **pad)

    def delete_dumps():
        cfg = build_cfg()
        initial = cfg.work_dir if os.path.isdir(cfg.work_dir) else STUDIES
        target = filedialog.askdirectory(
            initialdir=initial,
            title="Select the dumps folder to delete")
        if not target:
            return
        if messagebox.askyesno(
                "Delete folder",
                f"Permanently delete this folder and everything in it?\n\n{target}"):
            import shutil
            shutil.rmtree(target, ignore_errors=True)
            append(f"Deleted {target}.")

    def do_extend():
        try:
            step = float(V["tau_extend_step"].get())
            nr = int(float(V["n_extend_rounds"].get()))
        except ValueError:
            step = nr = 0
        if step <= 0 or nr <= 0:
            messagebox.showinfo(
                "Extend",
                "Set a nonzero 'extend step' and 'extend rounds' in section 3 "
                "first -- these control how far to widen temperatures whose CRSS "
                "fell outside the bracket, and how many extra rounds to grant.")
            return
        do(cc.extend)

    def open_send():
        cfg = build_cfg()
        win = tk.Toplevel(root)
        win.title("Send files to cluster (scp)")
        win.geometry("720x360")
        prefs = load_prefs()
        fv = {}
        rem = {}
        rows = [("REMOTE_USER", "cluster username", ""),
                ("REMOTE_HOST", "cluster hostname", ""),
                ("SUBMIT_DIR", "remote submit dir (cases / template / jobs)", ""),
                ("USFE_DIR", "remote USFE dir", "")]
        for i, (k, lab, dflt) in enumerate(rows):
            ttk.Label(win, text=lab).grid(row=i, column=0, sticky="w", **pad)
            fv[k] = tk.StringVar(value=prefs.get(k, dflt))
            ttk.Entry(win, textvariable=fv[k], width=34).grid(
                row=i, column=1, columnspan=2, **pad)
            rem[k] = tk.BooleanVar(value=(k in prefs))
            ttk.Checkbutton(win, text="remember", variable=rem[k]).grid(
                row=i, column=3, sticky="w", **pad)
        send_cases = tk.BooleanVar(value=True)
        send_tmpl = tk.BooleanVar(value=False)
        send_usfe = tk.BooleanVar(value=False)
        ttk.Checkbutton(win, text="cases file", variable=send_cases).grid(
            row=len(rows), column=0, sticky="w", **pad)
        ttk.Checkbutton(win, text="template + SLURM job", variable=send_tmpl).grid(
            row=len(rows), column=1, sticky="w", **pad)
        ttk.Checkbutton(win, text="USFE files (from folder)",
                        variable=send_usfe).grid(row=len(rows), column=2,
                                                 sticky="w", **pad)

        def run_send():
            pr = load_prefs()
            for k in fv:
                if rem[k].get():
                    pr[k] = fv[k].get().strip()
                else:
                    pr.pop(k, None)
            save_prefs(pr)
            user = fv["REMOTE_USER"].get().strip()
            host = fv["REMOTE_HOST"].get().strip()
            submit = fv["SUBMIT_DIR"].get().strip()
            usfe_remote = fv["USFE_DIR"].get().strip()
            if not user or not host:
                messagebox.showerror("Send", "Set cluster username and hostname.")
                return
            jobs = []
            if send_cases.get():
                if not os.path.exists(cfg.cases_path()):
                    append("! no cases file yet -- generate a round first.")
                elif not submit:
                    append("! set a remote submit dir for the cases file.")
                else:
                    jobs.append((["scp", cfg.cases_path(),
                                  f"{user}@{host}:{submit}/"], "cases file"))
            if send_tmpl.get():
                tmpl = os.path.join(HERE, "templates", f"in.{cfg.tag}_template")
                jobsh = os.path.join(HERE, "slurm", f"{cfg.tag}_jobs")
                files = [p for p in (tmpl, jobsh) if os.path.exists(p)]
                if files and submit:
                    jobs.append((["scp", *files, f"{user}@{host}:{submit}/"],
                                 "template + job"))
                else:
                    append("! template/job not found or no submit dir set.")
            if send_usfe.get():
                if not cfg.usfe_dir:
                    append("! set a USFE folder in section 5 to send USFE files.")
                elif not usfe_remote:
                    append("! set a remote USFE dir.")
                else:
                    files = []
                    for T in cfg.temps():
                        name = cfg.usfe_pattern.format(T=T) \
                            if "{T}" in cfg.usfe_pattern else cfg.usfe_pattern
                        pth = os.path.join(cfg.usfe_dir, name)
                        if os.path.exists(pth):
                            files.append(pth)
                    if files:
                        jobs.append((["scp", *files,
                                      f"{user}@{host}:{usfe_remote}/"],
                                     f"{len(files)} USFE files"))
                    else:
                        append("! no USFE files found in the folder.")
            if not jobs:
                append("Nothing to send.")
                return
            win.destroy()

            def worker():
                for cmd, desc in jobs:
                    append(f"scp {desc} -> {host} ...")
                    try:
                        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                                stderr=subprocess.STDOUT, text=True)
                        for line in proc.stdout:
                            append(line.rstrip())
                        proc.wait()
                        append(f"  {desc}: exit {proc.returncode}")
                    except Exception as e:      # noqa: BLE001
                        append("! scp failed: " + str(e))
            threading.Thread(target=worker, daemon=True).start()
        ttk.Button(win, text="Send", command=run_send).grid(
            row=len(rows) + 1, column=0, columnspan=3, **pad)

    def plot_crss():
        cfg = build_cfg()
        state = cfg.state_path()
        if not os.path.exists(state):
            messagebox.showerror("Plot", "No search state yet -- init and run "
                                         "at least one round first.")
            return
        default = f"{cfg.tag}  {cfg.slip_system}  {cfg.stress_mode}"
        label = simpledialog.askstring(
            "Plot legend",
            "Legend label for this study\n"
            "(e.g. material, slip plane of interest, applied stress state):",
            initialvalue=default, parent=root)
        script = os.path.join(HERE, "plot_crss.py")
        args = ["python3", script, state]
        if label:
            args += ["--labels", label]
        append(f"Launching plot for {state} ...")
        try:
            subprocess.Popen(args)
        except Exception as e:            # noqa: BLE001
            messagebox.showerror("Plot", str(e))

    # banner + welcome
    append(pilot_banner.BANNER)
    append("\nStress-state workflow ready. Study folder: " + STUDIES)
    append("Typical flow: set parameters -> Init -> Generate -> (run jobs) -> "
           "Fetch dumps -> Update -> repeat; Extend if needed; Plot CRSS vs T.")

    root.mainloop()


def main():
    try:
        import tkinter  # noqa: F401
    except Exception:
        print("Tkinter is unavailable, so the GUI can't open.\n"
              "Use the command-line driver instead:  python3 crss_search_alloy.py")
        sys.exit(1)
    run_gui()


if __name__ == "__main__":
    main()
