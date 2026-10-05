#!/usr/bin/env python3
"""
usfe_workflow_gui.py -- one guided workflow: build lattices (OTIS) and
interpolate USFE files, from a single window.

    python3 usfe_workflow_gui.py          # graphical interface
    python3 usfe_workflow_gui.py --cli     # text prompts (no Tk needed)

All the real work lives in usfe_core.py; this file is just the front-end.
Keep otis.py, the interp script(s) and your CSV reachable from where you run it
(same folder is simplest).
"""

from logging import root
import os
import queue
import sys
import threading
import traceback

import usfe_core as core

# Interpolation choices shown in the UI -> the script each one runs.
METHODS = {
    "{110} slip  (110_interp_usfe.py)": "110_interp_usfe.py",
    "{112} slip  (112_interp_usfe_nonbcc.py)": "112_interp_usfe_nonbcc.py",
}


# --------------------------- shared parsing helpers ----------------------

def _tokens(s):
    return [t for t in s.replace(",", " ").split() if t != ""]


def parse_config_from_strings(v):
    """Build a BuildConfig from a dict of raw strings (used by both front-ends).
    Raises ValueError with a readable message on bad input."""
    try:
        cell = tuple(int(x) for x in _tokens(v["cell"]))
        if len(cell) != 3:
            raise ValueError
    except ValueError:
        raise ValueError("Cell size must be three positive integers, e.g. 32 32 32")

    labels = _tokens(v["labels"])
    try:
        conc = [float(x) for x in _tokens(v["conc"])]
    except ValueError:
        raise ValueError("Concentrations must be numbers.")

    matrix = None
    input_type = "prob"
    if not v.get("random", True):
        input_type = v.get("input_type", "prob")
        rows = [r for r in v.get("matrix", "").splitlines() if r.strip()]
        try:
            matrix = [[float(x) for x in _tokens(r)] for r in rows]
        except ValueError:
            raise ValueError("SRO matrix entries must be numbers.")

    def _int(key, what):
        try:
            return int(v[key])
        except ValueError:
            raise ValueError(f"{what} must be an integer.")

    def _float(key, what):
        try:
            return float(v[key])
        except ValueError:
            raise ValueError(f"{what} must be a number.")

    return core.BuildConfig(
        lattice_type=v["lattice_type"].strip().upper(),
        cell=cell,
        lattice_param=_float("lattice_param", "Lattice parameter"),
        labels=labels,
        conc=conc,
        input_type=input_type,
        matrix=matrix,
        n_lat=_int("n_lat", "Number of lattices"),
        index_start=_int("index_start", "Index start"),
        index_step=_int("index_step", "Index step"),
        usfe_pattern=v["usfe_pattern"].strip(),
        lattice_pattern=v["lattice_pattern"].strip(),
        keep_lattices=bool(v.get("keep_lattices", False)),
        out_dir=v["out_dir"].strip() or ".",
        interp=v["interp"].strip(),
        csv=v["csv"].strip(),
        otis=v.get("otis", "otis.py").strip() or "otis.py",
        usfe_col=int(v.get("usfe_col", 4)),
        usfe_skip_header=int(v.get("usfe_skip_header", 1)),
    )


DEFAULTS = {
    "lattice_type": "BCC", "cell": "32 32 32", "lattice_param": "2.86",
    "labels": "Mo Nb Ti", "conc": "0.33 0.33 0.33",
    "random": True, "input_type": "prob", "matrix": "",
    "n_lat": "15", "index_start": "0", "index_step": "100",
    "usfe_pattern": "USFE_{i}", "lattice_pattern": "lattice_{i}",
    "keep_lattices": False, "out_dir": ".",
    "interp": "110_interp_usfe.py", "csv": "MoNbTi_sfe_data.csv",
    "otis": "otis.py", "usfe_col": "4", "usfe_skip_header": "1",
}


# ================================ GUI ====================================

def run_gui():
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox
    from tkinter.scrolledtext import ScrolledText

    root = tk.Tk()
    root.title("USFE workflow  -  OTIS lattice + interpolation")
    root.geometry("760x820")
    root.lift()
    root.focus_force()
    root.attributes("-topmost", True)
    root.after(200, lambda: root.attributes("-topmost", False))

    pad = {"padx": 6, "pady": 3}
    vars_ = {}

    def entry(parent, key, width=18):
        var = tk.StringVar(value=str(DEFAULTS[key]))
        vars_[key] = var
        return ttk.Entry(parent, textvariable=var, width=width)

    # ---------- 1. Lattice ----------
    f1 = ttk.LabelFrame(root, text="1.  Lattice  (OTIS)")
    f1.pack(fill="x", **pad)

    ttk.Label(f1, text="Lattice type").grid(row=0, column=0, sticky="w", **pad)
    vars_["lattice_type"] = tk.StringVar(value=DEFAULTS["lattice_type"])
    ttk.Combobox(f1, textvariable=vars_["lattice_type"], width=10,
                 values=["BCC", "FCC", "CUSTOM"], state="readonly"
                 ).grid(row=0, column=1, sticky="w", **pad)

    ttk.Label(f1, text="Cell  N1 N2 N3").grid(row=0, column=2, sticky="w", **pad)
    entry(f1, "cell", 14).grid(row=0, column=3, sticky="w", **pad)

    ttk.Label(f1, text="Lattice parameter").grid(row=1, column=0, sticky="w", **pad)
    entry(f1, "lattice_param", 10).grid(row=1, column=1, sticky="w", **pad)

    ttk.Label(f1, text="Elements").grid(row=2, column=0, sticky="w", **pad)
    entry(f1, "labels", 24).grid(row=2, column=1, columnspan=1, sticky="w", **pad)
    ttk.Label(f1, text="Concentrations").grid(row=2, column=2, sticky="w", **pad)
    entry(f1, "conc", 24).grid(row=2, column=3, sticky="w", **pad)

    vars_["random"] = tk.BooleanVar(value=DEFAULTS["random"])
    vars_["input_type"] = tk.StringVar(value=DEFAULTS["input_type"])
    matrix_txt = ScrolledText(f1, width=40, height=4)
    itype_box = ttk.Combobox(f1, textvariable=vars_["input_type"], width=6,
                             values=["prob", "wc"], state="readonly")

    def toggle_sro():
        on = not vars_["random"].get()
        state = "normal" if on else "disabled"
        matrix_txt.configure(state=state)
        itype_box.configure(state="readonly" if on else "disabled")

    ttk.Checkbutton(f1, text="Random solid solution (no SRO)",
                    variable=vars_["random"], command=toggle_sro
                    ).grid(row=3, column=0, columnspan=2, sticky="w", **pad)
    ttk.Label(f1, text="SRO type").grid(row=3, column=2, sticky="e", **pad)
    itype_box.grid(row=3, column=3, sticky="w", **pad)
    ttk.Label(f1, text="SRO matrix (prob or WC), one row per element:"
              ).grid(row=4, column=0, columnspan=4, sticky="w", **pad)
    matrix_txt.grid(row=5, column=0, columnspan=4, sticky="we", **pad)
    vars_["_matrix_widget"] = matrix_txt
    toggle_sro()

    # ---------- 2. How many / naming ----------
    f2 = ttk.LabelFrame(root, text="2.  How many lattices  /  naming")
    f2.pack(fill="x", **pad)

    ttk.Label(f2, text="Number of lattices").grid(row=0, column=0, sticky="w", **pad)
    entry(f2, "n_lat", 8).grid(row=0, column=1, sticky="w", **pad)
    ttk.Label(f2, text="Index start").grid(row=0, column=2, sticky="w", **pad)
    entry(f2, "index_start", 8).grid(row=0, column=3, sticky="w", **pad)
    ttk.Label(f2, text="Index step").grid(row=0, column=4, sticky="w", **pad)
    entry(f2, "index_step", 8).grid(row=0, column=5, sticky="w", **pad)

    ttk.Label(f2, text="USFE name pattern").grid(row=1, column=0, sticky="w", **pad)
    entry(f2, "usfe_pattern", 16).grid(row=1, column=1, sticky="w", **pad)
    ttk.Label(f2, text="( {i} -> index )").grid(row=1, column=2, sticky="w", **pad)

    vars_["keep_lattices"] = tk.BooleanVar(value=DEFAULTS["keep_lattices"])
    ttk.Checkbutton(f2, text="Keep lattice files", variable=vars_["keep_lattices"]
                    ).grid(row=2, column=0, sticky="w", **pad)
    ttk.Label(f2, text="Lattice name pattern").grid(row=2, column=1, sticky="e", **pad)
    entry(f2, "lattice_pattern", 16).grid(row=2, column=2, sticky="w", **pad)

    ttk.Label(f2, text="Output folder").grid(row=3, column=0, sticky="w", **pad)
    entry(f2, "out_dir", 34).grid(row=3, column=1, columnspan=3, sticky="w", **pad)
    ttk.Button(f2, text="Browse...",
               command=lambda: _pick_dir(vars_["out_dir"])
               ).grid(row=3, column=4, sticky="w", **pad)

    # ---------- 3. Interpolation ----------
    f3 = ttk.LabelFrame(root, text="3.  Interpolation")
    f3.pack(fill="x", **pad)

    ttk.Label(f3, text="Method").grid(row=0, column=0, sticky="w", **pad)
    method_var = tk.StringVar(value=list(METHODS)[0])
    interp_entry = entry(f3, "interp", 40)

    def on_method(*_):
        vars_["interp"].set(METHODS[method_var.get()])
    ttk.Combobox(f3, textvariable=method_var, width=44, state="readonly",
                 values=list(METHODS)).grid(row=0, column=1, columnspan=3,
                                            sticky="w", **pad)
    method_var.trace_add("write", on_method)

    ttk.Label(f3, text="Interp script").grid(row=1, column=0, sticky="w", **pad)
    interp_entry.grid(row=1, column=1, columnspan=2, sticky="w", **pad)
    ttk.Button(f3, text="Browse...",
               command=lambda: _pick_file(vars_["interp"])
               ).grid(row=1, column=3, sticky="w", **pad)

    ttk.Label(f3, text="CSV (SFE data)").grid(row=2, column=0, sticky="w", **pad)
    entry(f3, "csv", 40).grid(row=2, column=1, columnspan=2, sticky="w", **pad)
    ttk.Button(f3, text="Browse...",
               command=lambda: _pick_file(vars_["csv"])
               ).grid(row=2, column=3, sticky="w", **pad)

    # ---------- Advanced ----------
    f4 = ttk.LabelFrame(root, text="Advanced")
    f4.pack(fill="x", **pad)
    ttk.Label(f4, text="OTIS script").grid(row=0, column=0, sticky="w", **pad)
    entry(f4, "otis", 24).grid(row=0, column=1, sticky="w", **pad)
    ttk.Button(f4, text="Browse...",
               command=lambda: _pick_file(vars_["otis"])
               ).grid(row=0, column=2, sticky="w", **pad)
    ttk.Label(f4, text="USFE col").grid(row=0, column=3, sticky="e", **pad)
    entry(f4, "usfe_col", 5).grid(row=0, column=4, sticky="w", **pad)
    ttk.Label(f4, text="skip header").grid(row=0, column=5, sticky="e", **pad)
    entry(f4, "usfe_skip_header", 5).grid(row=0, column=6, sticky="w", **pad)

    # ---------- Controls + log ----------
    f5 = ttk.Frame(root)
    f5.pack(fill="both", expand=True, **pad)
    btn = ttk.Button(f5, text="Generate USFE files")
    btn.grid(row=0, column=0, sticky="w", **pad)
    stop_btn = ttk.Button(f5, text="Stop", state="disabled")
    stop_btn.grid(row=0, column=1, sticky="w", **pad)
    prog = ttk.Progressbar(f5, length=360, mode="determinate", maximum=1.0)
    prog.grid(row=0, column=2, sticky="we", **pad)
    status_var = tk.StringVar(value="Ready.")
    ttk.Label(f5, textvariable=status_var).grid(row=1, column=0, columnspan=3,
                                                sticky="w", **pad)
    log = ScrolledText(f5, height=16, state="disabled", wrap="none")
    log.grid(row=2, column=0, columnspan=3, sticky="nsew", **pad)
    f5.rowconfigure(2, weight=1)
    f5.columnconfigure(2, weight=1)

    def _pick_dir(var):
        d = filedialog.askdirectory()
        if d:
            var.set(d)

    def _pick_file(var):
        p = filedialog.askopenfilename()
        if p:
            var.set(p)

    def append_log(msg):
        log.configure(state="normal")
        log.insert("end", msg + "\n")
        log.see("end")
        log.configure(state="disabled")

    q = queue.Queue()
    stop_event = threading.Event()

    def worker(cfg):
        try:
            core.build_lattices(
                cfg,
                log=lambda m: q.put(("line", m)),
                status=lambda m: q.put(("status", m)),
                progress=lambda f: q.put(("progress", f)),
                should_stop=stop_event.is_set,
            )
            q.put(("done", None))
        except Exception as e:                      # noqa: BLE001
            q.put(("error", str(e)))
            q.put(("trace", traceback.format_exc()))

    def poll():
        try:
            while True:
                kind, payload = q.get_nowait()
                if kind == "line":
                    append_log(payload)
                elif kind == "status":
                    status_var.set(payload)
                elif kind == "progress":
                    prog["value"] = payload
                elif kind == "trace":
                    append_log(payload)
                elif kind == "error":
                    status_var.set("Error.")
                    append_log("ERROR: " + payload)
                    messagebox.showerror("USFE workflow", payload)
                    _finish()
                elif kind == "done":
                    status_var.set("Done.")
                    prog["value"] = 1.0
                    _finish()
        except queue.Empty:
            pass
        root.after(80, poll)

    def _finish():
        btn.configure(state="normal")
        stop_btn.configure(state="disabled")

    def on_generate():
        raw = {k: (v.get() if hasattr(v, "get") else v) for k, v in vars_.items()
               if not k.startswith("_")}
        raw["matrix"] = vars_["_matrix_widget"].get("1.0", "end")
        try:
            cfg = parse_config_from_strings(raw)
        except ValueError as e:
            messagebox.showerror("Check your inputs", str(e))
            return
        problems = core.validate(cfg)
        if problems:
            messagebox.showerror("Check your inputs", "\n".join(problems))
            return
        log.configure(state="normal")
        log.delete("1.0", "end")
        log.configure(state="disabled")
        prog["value"] = 0.0
        stop_event.clear()
        btn.configure(state="disabled")
        stop_btn.configure(state="normal")
        status_var.set("Running...")
        threading.Thread(target=worker, args=(cfg,), daemon=True).start()

    btn.configure(command=on_generate)
    stop_btn.configure(command=lambda: (stop_event.set(),
                                        status_var.set("Stopping after this step...")))
    root.after(80, poll)
    root.mainloop()


# ============================ CLI fallback ===============================

def cli_main():
    print("USFE workflow (command-line mode)\n"
          "Press Enter to accept the [default] shown.\n")

    def ask(key, prompt):
        d = DEFAULTS[key]
        r = input(f"{prompt} [{d}]: ").strip()
        return r if r else str(d)

    raw = {}
    raw["lattice_type"] = ask("lattice_type", "Lattice type (BCC/FCC/CUSTOM)")
    raw["cell"] = ask("cell", "Cell size N1 N2 N3")
    raw["lattice_param"] = ask("lattice_param", "Lattice parameter")
    raw["labels"] = ask("labels", "Elements")
    raw["conc"] = ask("conc", "Concentrations")
    raw["random"] = True                      # CLI keeps it simple: random SS
    raw["n_lat"] = ask("n_lat", "Number of lattices")
    raw["index_start"] = ask("index_start", "Index start")
    raw["index_step"] = ask("index_step", "Index step")
    raw["usfe_pattern"] = ask("usfe_pattern", "USFE name pattern ({i}=index)")
    kl = input("Keep lattice files? (y/N): ").strip().lower().startswith("y")
    raw["keep_lattices"] = kl
    raw["lattice_pattern"] = ask("lattice_pattern", "Lattice name pattern")
    raw["out_dir"] = ask("out_dir", "Output folder")
    print("\nInterpolation method:")
    for n, k in enumerate(METHODS):
        print(f"  {n+1}) {k}")
    sel = input(f"Choose 1-{len(METHODS)} [1]: ").strip() or "1"
    raw["interp"] = list(METHODS.values())[int(sel) - 1]
    raw["csv"] = ask("csv", "CSV (SFE data) file")
    raw["otis"] = ask("otis", "OTIS script")
    raw["usfe_col"] = DEFAULTS["usfe_col"]
    raw["usfe_skip_header"] = DEFAULTS["usfe_skip_header"]

    cfg = parse_config_from_strings(raw)
    problems = core.validate(cfg)
    if problems:
        print("\nCannot start:")
        for p in problems:
            print("  - " + p)
        return 1

    print()
    core.build_lattices(
        cfg,
        log=print,
        status=lambda m: print("\r" + m, end="", flush=True),
        progress=lambda f: None,
    )
    return 0


def main():
    if "--cli" in sys.argv:
        sys.exit(cli_main())
    try:
        import tkinter  # noqa: F401
    except Exception:
        print("Tkinter is not available, so the GUI can't open.\n"
              "Install it (e.g. 'sudo apt-get install python3-tk') for the window,\n"
              "or continue here in command-line mode.\n")
        sys.exit(cli_main())
    run_gui()


if __name__ == "__main__":
    main()
