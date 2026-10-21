#!/usr/bin/env python3
"""
usfe_workflow_gui.py -- PILOT lattice + USFE workflow, from a single window.

    python3 usfe_workflow_gui.py          # graphical interface
    python3 usfe_workflow_gui.py --cli    # text prompts (no Tk needed)

Pick a slip system, describe the alloy, and the tool builds random
solid-solution lattices and their per-site USFE files.  The slip system fixes
both the lattice type and the interpolation, so they can never be mismatched:

    {110}  ->  BCC lattice              ->  110_interp_usfe.py
    {112}  ->  CUSTOM_NON_ORTHO lattice ->  112_interp_usfe_nonbcc.py

The real work lives in usfe_core.py; this file is just the front end.  A
"View generated files..." button opens viewer.py to inspect any result.
Keep lattice_generator.py, the interp scripts, viewer.py and your CSV reachable
from where you run this (same folder is simplest).
"""

import os
import queue
import sys
import threading
import traceback

import usfe_core as core
import pilot_banner

# Slip system shown in the UI -> (lattice type, interpolation script).
# The interpolation is also derivable in core; this table drives the UI choice.
SLIP_SYSTEMS = {
    "{110}   (BCC lattice)":              ("BCC", "110_interp_usfe.py"),
    "{112}   (CUSTOM_NON_ORTHO lattice)": ("CUSTOM_NON_ORTHO",
                                           "112_interp_usfe_nonbcc.py"),
}


# --------------------------- shared parsing helpers ----------------------

def _tokens(s):
    return [t for t in s.replace(",", " ").split() if t != ""]


def parse_config_from_strings(vals):
    """Build a BuildConfig from a dict of raw strings (used by both front-ends).
    Raises ValueError with a readable message on bad input."""
    try:
        cell = tuple(int(x) for x in _tokens(vals["cell"]))
        if len(cell) != 3:
            raise ValueError
    except ValueError:
        raise ValueError("Cell size must be three positive integers, e.g. 32 32 32")

    labels = _tokens(vals["labels"])
    try:
        conc = [float(x) for x in _tokens(vals["conc"])]
    except ValueError:
        raise ValueError("Concentrations must be numbers.")

    def _int(key, what):
        try:
            return int(vals[key])
        except ValueError:
            raise ValueError(f"{what} must be an integer.")

    def _float(key, what):
        try:
            return float(vals[key])
        except ValueError:
            raise ValueError(f"{what} must be a number.")

    return core.BuildConfig(
        lattice_type=vals["lattice_type"].strip().upper(),
        cell=cell,
        lattice_param=_float("lattice_param", "Lattice parameter"),
        labels=labels,
        conc=conc,
        n_lat=_int("n_lat", "Number of lattices"),
        index_start=_int("index_start", "Index start"),
        index_step=_int("index_step", "Index step"),
        usfe_pattern=vals["usfe_pattern"].strip(),
        lattice_pattern=vals["lattice_pattern"].strip(),
        comp_pattern=vals["comp_pattern"].strip(),
        keep_lattices=bool(vals.get("keep_lattices", True)),
        keep_comp=bool(vals.get("keep_comp", True)),
        out_dir=vals["out_dir"].strip() or ".",
        csv=vals["csv"].strip(),
        generator=vals.get("generator", "lattice_generator.py").strip()
        or "lattice_generator.py",
        usfe_col=int(vals.get("usfe_col", 4)),
        usfe_skip_header=int(vals.get("usfe_skip_header", 1)),
    )


DEFAULTS = {
    "slip": list(SLIP_SYSTEMS)[0],
    "lattice_type": "BCC", "cell": "32 32 32", "lattice_param": "3.16",
    "labels": "Ta W", "conc": "0.5 0.5",
    "n_lat": "1", "index_start": "0", "index_step": "100",
    "usfe_pattern": "USFE_{i}", "lattice_pattern": "lattice_{i}",
    "comp_pattern": "comp_{i}",
    "keep_lattices": True, "keep_comp": True, "out_dir": ".",
    "csv": "TaW_110_SFE_NequIP.csv",
    "generator": "lattice_generator.py", "usfe_col": "4", "usfe_skip_header": "1",
}


# ================================ GUI ====================================

def run_gui():
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox
    from tkinter.scrolledtext import ScrolledText

    root = tk.Tk()
    root.title("PILOT  -  lattice + USFE workflow")
    root.geometry("780x760")
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

    # ---------- 1. Slip system + alloy ----------
    f1 = ttk.LabelFrame(root, text="1.  Slip system  &  alloy")
    f1.pack(fill="x", **pad)

    ttk.Label(f1, text="Slip system").grid(row=0, column=0, sticky="w", **pad)
    vars_["slip"] = tk.StringVar(value=DEFAULTS["slip"])
    vars_["lattice_type"] = tk.StringVar(value=DEFAULTS["lattice_type"])
    slip_box = ttk.Combobox(f1, textvariable=vars_["slip"], width=34,
                            values=list(SLIP_SYSTEMS), state="readonly")
    slip_box.grid(row=0, column=1, columnspan=3, sticky="w", **pad)

    ttk.Label(f1, text="Cell  N1 N2 N3").grid(row=1, column=0, sticky="w", **pad)
    entry(f1, "cell", 14).grid(row=1, column=1, sticky="w", **pad)
    ttk.Label(f1, text="Lattice parameter").grid(row=1, column=2, sticky="w", **pad)
    entry(f1, "lattice_param", 10).grid(row=1, column=3, sticky="w", **pad)

    ttk.Label(f1, text="# of alloying elements").grid(row=2, column=0,
                                                      sticky="w", **pad)
    n_elem_var = tk.StringVar(value="2")
    ttk.Spinbox(f1, from_=1, to=8, width=5, textvariable=n_elem_var).grid(
        row=2, column=1, sticky="w", **pad)
    ttk.Label(f1, text="(last concentration auto-fills to 1 - the rest)",
              foreground="#555").grid(row=2, column=2, columnspan=2,
                                      sticky="w", **pad)

    # dynamic element / concentration grid (row of elements over a row of concs)
    alloy_frame = ttk.Frame(f1)
    alloy_frame.grid(row=3, column=0, columnspan=4, sticky="w", **pad)
    elem_vars = []
    conc_vars = []

    def _fmt_conc(x):
        txt = f"{x:.6f}".rstrip("0").rstrip(".")
        return txt if txt else "0"

    def recompute_last(*_):
        n = len(conc_vars)
        if n == 0:
            return
        if n == 1:
            conc_vars[0].set("1")
            return
        total = 0.0
        for c in range(n - 1):
            try:
                total += float(conc_vars[c].get())
            except ValueError:
                return
        conc_vars[n - 1].set(_fmt_conc(1.0 - total))

    def rebuild_alloy(*_):
        try:
            n = int(float(n_elem_var.get()))
        except ValueError:
            return
        n = max(1, min(n, 8))
        old_e = [v.get() for v in elem_vars]
        old_c = [v.get() for v in conc_vars]
        for w in alloy_frame.winfo_children():
            w.destroy()
        elem_vars.clear()
        conc_vars.clear()
        seed_e = old_e if any(old_e) else ["Ta", "W", "Mo", "Nb", "Ti", "V",
                                           "Cr", "Fe"]
        seed_c = old_c if any(old_c) else ["0.5"]
        ttk.Label(alloy_frame, text="element type").grid(
            row=0, column=0, sticky="e", **pad)
        ttk.Label(alloy_frame, text="concentration").grid(
            row=1, column=0, sticky="e", **pad)
        for c in range(n):
            ev = tk.StringVar(value=seed_e[c] if c < len(seed_e) else "")
            cv = tk.StringVar(value=seed_c[c] if c < len(seed_c) else "")
            elem_vars.append(ev)
            conc_vars.append(cv)
            ttk.Entry(alloy_frame, textvariable=ev, width=8).grid(
                row=0, column=c + 1, **pad)
            is_last = (c == n - 1 and n >= 2)
            ttk.Entry(alloy_frame, textvariable=cv, width=8,
                      state="readonly" if is_last else "normal").grid(
                row=1, column=c + 1, **pad)
            if not is_last:
                cv.trace_add("write", recompute_last)
        recompute_last()

    n_elem_var.trace_add("write", rebuild_alloy)
    rebuild_alloy()

    ttk.Label(f1, text="Random solid solution (no short-range order).",
              foreground="#555").grid(row=4, column=0, columnspan=4,
                                      sticky="w", **pad)

    # ---------- 2. How many / naming ----------
    f2 = ttk.LabelFrame(root, text="2.  How many lattices  /  naming")
    f2.pack(fill="x", **pad)

    ttk.Label(f2, text="Number of lattices").grid(row=0, column=0, sticky="w", **pad)
    entry(f2, "n_lat", 8).grid(row=0, column=1, sticky="w", **pad)
    ttk.Label(f2, text="Index start").grid(row=0, column=2, sticky="w", **pad)
    entry(f2, "index_start", 8).grid(row=0, column=3, sticky="w", **pad)
    ttk.Label(f2, text="Index step").grid(row=0, column=4, sticky="w", **pad)
    entry(f2, "index_step", 8).grid(row=0, column=5, sticky="w", **pad)

    ttk.Label(f2, text="USFE name").grid(row=1, column=0, sticky="w", **pad)
    entry(f2, "usfe_pattern", 14).grid(row=1, column=1, sticky="w", **pad)
    ttk.Label(f2, text="Lattice name").grid(row=1, column=2, sticky="w", **pad)
    entry(f2, "lattice_pattern", 14).grid(row=1, column=3, sticky="w", **pad)
    ttk.Label(f2, text="Comp name").grid(row=1, column=4, sticky="w", **pad)
    entry(f2, "comp_pattern", 14).grid(row=1, column=5, sticky="w", **pad)
    ttk.Label(f2, text="( {i} -> index in every name )", foreground="#555"
              ).grid(row=2, column=0, columnspan=3, sticky="w", **pad)

    vars_["keep_lattices"] = tk.BooleanVar(value=DEFAULTS["keep_lattices"])
    ttk.Checkbutton(f2, text="Keep lattice files", variable=vars_["keep_lattices"]
                    ).grid(row=2, column=3, sticky="w", **pad)
    vars_["keep_comp"] = tk.BooleanVar(value=DEFAULTS["keep_comp"])
    ttk.Checkbutton(f2, text="Keep composition files",
                    variable=vars_["keep_comp"]
                    ).grid(row=2, column=4, columnspan=2, sticky="w", **pad)

    ttk.Label(f2, text="Output folder").grid(row=3, column=0, sticky="w", **pad)
    entry(f2, "out_dir", 40).grid(row=3, column=1, columnspan=4, sticky="w", **pad)
    ttk.Button(f2, text="Browse...",
               command=lambda: _pick_dir(vars_["out_dir"])
               ).grid(row=3, column=5, sticky="w", **pad)

    # ---------- 3. SFE data ----------
    f3 = ttk.LabelFrame(root, text="3.  SFE data (CSV)")
    f3.pack(fill="x", **pad)
    ttk.Label(f3, text="CSV").grid(row=0, column=0, sticky="w", **pad)
    entry(f3, "csv", 46).grid(row=0, column=1, columnspan=2, sticky="w", **pad)
    ttk.Button(f3, text="Browse...",
               command=lambda: _pick_file(vars_["csv"])
               ).grid(row=0, column=3, sticky="w", **pad)
    ttk.Label(f3, text="Header names the elements; 2 columns = binary, "
                       "3 = ternary.", foreground="#555"
              ).grid(row=1, column=0, columnspan=4, sticky="w", **pad)

    # keep the generator path configurable but out of the way
    vars_["generator"] = tk.StringVar(value=DEFAULTS["generator"])
    vars_["usfe_col"] = tk.StringVar(value=DEFAULTS["usfe_col"])
    vars_["usfe_skip_header"] = tk.StringVar(value=DEFAULTS["usfe_skip_header"])

    # keep lattice_type in sync with the slip system choice
    def on_slip(*_):
        lt, _interp = SLIP_SYSTEMS[vars_["slip"].get()]
        vars_["lattice_type"].set(lt)
        # nudge the default CSV name to match the slip system on first switches
    vars_["slip"].trace_add("write", on_slip)
    on_slip()

    # ---------- Controls + log ----------
    f5 = ttk.Frame(root)
    f5.pack(fill="both", expand=True, **pad)
    btn = ttk.Button(f5, text="Generate USFE files")
    btn.grid(row=0, column=0, sticky="w", **pad)
    stop_btn = ttk.Button(f5, text="Stop", state="disabled")
    stop_btn.grid(row=0, column=1, sticky="w", **pad)
    prog = ttk.Progressbar(f5, length=260, mode="determinate", maximum=1.0)
    prog.grid(row=0, column=2, sticky="we", **pad)

    # View button: its own full-width row, larger (extra vertical padding)
    view_btn = ttk.Button(f5, text="View generated files...")
    view_btn.grid(row=1, column=0, columnspan=3, sticky="we", ipady=10, **pad)

    status_var = tk.StringVar(value="Ready.")
    ttk.Label(f5, textvariable=status_var).grid(row=2, column=0, columnspan=3,
                                                sticky="w", **pad)
    log = ScrolledText(f5, height=14, state="disabled", wrap="none",
                   font="TkFixedFont")
    log.grid(row=3, column=0, columnspan=3, sticky="nsew", **pad)
    f5.rowconfigure(3, weight=1)
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
                    messagebox.showerror("PILOT workflow", payload)
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

    def _collect_cfg():
        raw = {k: (val.get() if hasattr(val, "get") else val)
               for k, val in vars_.items() if not k.startswith("_")}
        raw["labels"] = " ".join(v.get().strip() for v in elem_vars
                                 if v.get().strip())
        raw["conc"] = " ".join(v.get().strip() for v in conc_vars
                               if v.get().strip())
        return parse_config_from_strings(raw)

    def on_generate():
        try:
            cfg = _collect_cfg()
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
        append_log(pilot_banner.BANNER)
        prog["value"] = 0.0
        stop_event.clear()
        btn.configure(state="disabled")
        stop_btn.configure(state="normal")
        status_var.set("Running...")
        threading.Thread(target=worker, args=(cfg,), daemon=True).start()

    def on_view():
        import viewer
        slip = "{112}" if vars_["lattice_type"].get() == "CUSTOM_NON_ORTHO" else "{110}"
        folder = vars_["out_dir"].get().strip() or "."
        viewer.run_viewer(master=root, folder=folder, slip=slip,
                          csv=vars_["csv"].get().strip())

    btn.configure(command=on_generate)
    stop_btn.configure(command=lambda: (stop_event.set(),
                                        status_var.set("Stopping after this step...")))
    view_btn.configure(command=on_view)
    root.after(80, poll)
    root.mainloop()


# ============================ CLI fallback ===============================

def cli_main():
    print("PILOT lattice + USFE workflow (command-line mode)\n"
          "Press Enter to accept the [default] shown.\n")

    def ask(key, prompt):
        d = DEFAULTS[key]
        r = input(f"{prompt} [{d}]: ").strip()
        return r if r else str(d)

    print("Slip system:")
    for n, k in enumerate(SLIP_SYSTEMS):
        print(f"  {n+1}) {k}")
    sel = input(f"Choose 1-{len(SLIP_SYSTEMS)} [1]: ").strip() or "1"
    slip = list(SLIP_SYSTEMS)[int(sel) - 1]
    lattice_type, _interp = SLIP_SYSTEMS[slip]

    raw = {}
    raw["slip"] = slip
    raw["lattice_type"] = lattice_type
    raw["cell"] = ask("cell", "Cell size N1 N2 N3")
    raw["lattice_param"] = ask("lattice_param", "Lattice parameter")
    raw["labels"] = ask("labels", "Elements")
    raw["conc"] = ask("conc", "Concentrations")
    raw["n_lat"] = ask("n_lat", "Number of lattices")
    raw["index_start"] = ask("index_start", "Index start")
    raw["index_step"] = ask("index_step", "Index step")
    raw["usfe_pattern"] = ask("usfe_pattern", "USFE name pattern ({i}=index)")
    raw["lattice_pattern"] = ask("lattice_pattern", "Lattice name pattern")
    raw["comp_pattern"] = ask("comp_pattern", "Composition name pattern")
    kl = input("Keep lattice files? (Y/n): ").strip().lower()
    raw["keep_lattices"] = not kl.startswith("n")
    kc = input("Keep composition files? (Y/n): ").strip().lower()
    raw["keep_comp"] = not kc.startswith("n")
    raw["out_dir"] = ask("out_dir", "Output folder")
    raw["csv"] = ask("csv", "CSV (SFE data) file")
    raw["generator"] = DEFAULTS["generator"]
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
    print(pilot_banner.BANNER)
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
