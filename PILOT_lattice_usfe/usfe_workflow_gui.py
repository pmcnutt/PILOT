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

The PILOT banner is built into this file (no separate pilot_banner.py needed).
In the window, drag the bar above the terminal up to enlarge the terminal, or
down to shrink it; the controls above it scroll when they don't all fit.
"""

import os
import queue
import sys
import threading
import traceback

import usfe_core as core

# --------------------------- PILOT banner --------------------------------

BANNER = r"""----------------------------------------------------------------------

                                 /||\
                                / -- \
                               /      \
                               |  ::  |
                               |  __  |
                               | /  \ |
     ____________________ _|_  | |  | |  _|_ ____________________
    /--------------------/   \----------/   \--------------------\
    |                    |   |          |   |                    |
    |     :              |   |    ::    |   |            :       |
    |_______________ ___ \   /          \   / ___ _______________|
     \_____________|___|__| |____________| |__|___|_____________/
                          | |  |  ||  |  | |
                               |  ||  |
                               |      |
                               / \  / \
                             -(   )(   )-
                               \ /  \ /
                                \    /
                           _____/ || \_____
                          / :    \||/    : \
                          | ____  ||  ____ |
                           \___|__||__|___/
                                  ||

     ____           ___          _               ___           _____
    |  _ \         |_ _|        | |             / _ \         |_   _|
    | |_) |         | |         | |            | | | |          | |
    |  __/          | |         | |            | | | |          | |
    | |             | |         | |___         | |_| |          | |
    |_|            |___|        |_____|         \___/           |_|

Phase-field dislocation dynamics Input, Lattice & Optimization Toolkit
                  
                  3D Lattice and USFE Cell Builder 
----------------------------------------------------------------------"""


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
    "csv": "SFE_data/TaW_110_SFE_NequIP.csv",
    "generator": "lattice_generator.py", "usfe_col": "4", "usfe_skip_header": "1",
}


# ================================ GUI ====================================

# ------------------- resizable terminal + scrollable panel -----------------

def _split_layout(root, tk, ttk):
    """Split the window into a scrollable control panel (top) and the
    terminal (bottom), with a draggable bar between them.

        drag the bar UP    -> bigger terminal
        drag the bar DOWN  -> smaller terminal

    The control panel has its own vertical scrollbar (and mouse-wheel
    scrolling), so nothing is lost when the terminal takes up more room.

    Returns (form, term_frame, finish):
        form        -- put the controls in here
        term_frame  -- put the terminal (log) in here
        finish()    -- call once every widget exists; sets the starting bar
                       position and the mouse-wheel routing.
    """
    paned = tk.PanedWindow(root, orient="vertical", sashwidth=10,
                           sashrelief="raised", sashpad=1, borderwidth=0,
                           showhandle=False, handlesize=10, handlepad=40,
                           opaqueresize=True, sashcursor="sb_v_double_arrow")
    paned.pack(fill="both", expand=True)

    # ---- top pane: scrollable control panel ----
    top = ttk.Frame(paned)
    try:
        bg = ttk.Style().lookup("TFrame", "background") or root.cget("bg")
    except tk.TclError:
        bg = root.cget("bg")
    canvas = tk.Canvas(top, highlightthickness=0, borderwidth=0,
                       background=bg, yscrollincrement=20)
    vsb = ttk.Scrollbar(top, orient="vertical", command=canvas.yview)
    canvas.configure(yscrollcommand=vsb.set)
    vsb.pack(side="right", fill="y")
    canvas.pack(side="left", fill="both", expand=True)
    form = ttk.Frame(canvas)
    form_id = canvas.create_window((0, 0), window=form, anchor="nw")

    def _sync(_e=None):
        # form spans the full canvas width; scroll region = form's full height
        w = max(canvas.winfo_width(), form.winfo_reqwidth())
        h = max(canvas.winfo_height(), form.winfo_reqheight())
        canvas.itemconfigure(form_id, width=w)
        canvas.configure(scrollregion=(0, 0, w, h))
    form.bind("<Configure>", _sync)
    canvas.bind("<Configure>", _sync)

    # ---- bottom pane: terminal ----
    term_frame = ttk.Frame(paned)
    ttk.Label(term_frame,
              text="Terminal   (drag the bar above up / down to resize)",
              foreground="#555").pack(side="top", anchor="w", padx=6)

    paned.add(top, minsize=60)
    paned.add(term_frame, minsize=80)    # terminal gets extra room on resize

    # ---- mouse wheel scrolls the control panel when the pointer is over it ----
    def _scroll(e):
        if form.winfo_reqheight() <= canvas.winfo_height():
            return                                   # everything fits already
        if e.num == 4:
            step = -1
        elif e.num == 5:
            step = 1
        elif e.delta:
            step = -1 if e.delta > 0 else 1
        else:
            return
        canvas.yview_scroll(step, "units")

    def _over_panel(e):
        try:
            w = root.winfo_containing(e.x_root, e.y_root)
        except (KeyError, tk.TclError):
            return False
        if w is None:
            return False
        p, t = str(w), str(top)
        return p == t or p.startswith(t + ".")

    def _wheel(e):
        if _over_panel(e):
            _scroll(e)

    for seq in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
        root.bind_all(seq, _wheel, add="+")

    def _guard_wheel(widget):
        # Comboboxes / spinboxes change their value on mouse-wheel.  Inside a
        # scrolling panel that would silently flip settings (e.g. the slip
        # system) while the user is just scrolling, so scroll the panel instead.
        for child in widget.winfo_children():
            if child.winfo_class() in ("TCombobox", "TSpinbox", "Spinbox"):
                for seq in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
                    child.bind(seq, lambda e: (_scroll(e), "break")[1])
            _guard_wheel(child)

    def finish():
        _guard_wheel(form)
        # Start the bar right under the controls, like the original layout:
        # the controls get their natural height, the terminal gets the rest.
        form.update_idletasks()
        paned.paneconfigure(top, height=form.winfo_reqheight())

    return form, term_frame, finish


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

    # controls go in a scrollable panel; the terminal sits below a draggable bar
    form, term_frame, finish_layout = _split_layout(root, tk, ttk)

    pad = {"padx": 6, "pady": 3}
    vars_ = {}

    def entry(parent, key, width=18):
        var = tk.StringVar(value=str(DEFAULTS[key]))
        vars_[key] = var
        return ttk.Entry(parent, textvariable=var, width=width)

    # ---------- 1. Slip system + alloy ----------
    f1 = ttk.LabelFrame(form, text="1.  Slip system  &  alloy")
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
    f2 = ttk.LabelFrame(form, text="2.  How many lattices  /  naming")
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
    f3 = ttk.LabelFrame(form, text="3.  SFE data (CSV)")
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

    # ---------- Controls ----------
    f5 = ttk.Frame(form)
    f5.pack(fill="x", **pad)
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
    f5.columnconfigure(2, weight=1)

    # ---------- Terminal (log), below the draggable bar ----------
    log = ScrolledText(term_frame, height=14, state="disabled", wrap="none",
                       font="TkFixedFont")
    log.pack(fill="both", expand=True, **pad)

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
        append_log(BANNER)
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
    finish_layout()
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
    print(BANNER)
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
