#!/usr/bin/env python3
"""
viewer.py -- interactive viewer for the lattice / composition / USFE files that
the PILOT USFE workflow produces.  It reproduces the two Jupyter notebooks
(one for {110}, one for {112}) inside a single Tk window.

For a chosen (lattice, composition, USFE) triple it shows:
  * a 2-D atomic plane colored by local composition,
  * the same plane colored by USFE,
  * a full-cell 3-D view of the atom types,
  * a full-cell 3-D view of the USFE field.

Slip planes
-----------
Each slip system has three distinct slip planes, and the plane selector names
them explicitly:

    {110}:  0 -> (110)     1 -> (01-1)     2 -> (101)
    {112}:  0 -> (211)     1 -> (-1-21)    2 -> (-11-2)

Selecting a plane colors the USFE views with that plane's USFE column and, for
the 2-D plots, slices the cell along that plane's normal (all three planes
contain the Burgers vector b = [-111], which stays the in-plane horizontal
axis).  The 3-D views use a fixed viewing frame per slip system, so switching
the USFE index only recolors the cell -- it never reorients the camera.

Composition coloring
---------------------
Coloring (Coloring)
  * binary (2 components):     color by a chosen element fraction (1-D bar),
  * ternary (3 components):    ternary RGB triangle (each element = one corner),
  * 4+ components:             color by a chosen element fraction (1-D bar).

Color scale
------------
Every 1-D color bar auto-scales to the data by default; the "Color scale"
boxes let you fix a min/max (leave blank for auto) so the bar stays consistent
as you step through planes or USFE indices.

The 3-D views are live: drag to rotate, or use the elevation / azimuth sliders.
Axes are shown by default.

Open it from the workflow GUI ("View generated files..."), or run standalone:

    python3 viewer.py [folder]

The data-loading and projection helpers below have no Tk dependency and can be
imported and tested on their own.
"""

import os
import re
import sys

import numpy as np


# Burgers vector: lies in every slip plane listed below, so it is the natural
# in-plane horizontal axis for all of them.
B_DIR = np.array([-1, 1, 1], float)

# For each slip system: (plane normal, ascii label, TeX label), indexed 0/1/2 to
# match the l column in the USFE files.
SLIP_PLANES = {
    "{110}": [
        ((1, 1, 0),   "(110)",   "(110)"),
        ((0, 1, -1),  "(01-1)",  r"(01$\bar{1}$)"),
        ((1, 0, 1),   "(101)",   "(101)"),
    ],
    "{112}": [
        ((2, 1, 1),   "(211)",   "(211)"),
        ((-1, -2, 1), "(-1-21)", r"($\bar{1}\bar{2}$1)"),
        ((-1, 1, -2), "(-11-2)", r"($\bar{1}$1$\bar{2}$)"),
    ],
}
LATTICE_TO_SLIP = {"BCC": "{110}", "CUSTOM_NON_ORTHO": "{112}"}

# Plane-0 frames, kept for backward compatibility / standalone use.
SLIP_FRAMES = {
    "{110}": (B_DIR.copy(), np.cross(np.array([1, 1, 0], float), B_DIR)),
    "{112}": (B_DIR.copy(), np.cross(np.array([2, 1, 1], float), B_DIR)),
}


# ------------------------------ pure loaders -----------------------------

def read_lattice(path):
    """element + real-space xyz. Returns (coords[M,3], elems[M])."""
    coords = np.genfromtxt(path, usecols=[1, 2, 3], skip_header=2)
    elems = np.genfromtxt(path, usecols=[0], dtype=str, skip_header=2)
    return coords, elems


def read_comp(path, ncols):
    """Parse '(idx, j): [c0 c1 ...]' lines into comp[M, ncols], ordered by idx."""
    rows = {}
    pat = re.compile(r"\((\d+),\s*\d+\):\s*\[([^\]]+)\]")
    with open(path) as f:
        for line in f:
            m = pat.search(line)
            if m:
                rows[int(m.group(1))] = [float(x) for x in m.group(2).split()]
    if not rows:
        raise ValueError(f"No '(idx, j): [...]' rows found in {path}")
    comp = np.full((max(rows) + 1, ncols), np.nan)
    for i, val in rows.items():
        comp[i] = val
    return comp


def comp_ncols(path):
    """How many composition columns a comp file has (from its first data row)."""
    pat = re.compile(r"\((\d+),\s*\d+\):\s*\[([^\]]+)\]")
    with open(path) as f:
        for line in f:
            m = pat.search(line)
            if m:
                return len(m.group(2).split())
    raise ValueError(f"No composition rows found in {path}")


def read_usfe(path, slip_plane):
    """Keep one slip plane; return (usfe[N1*N2*N3], (N1,N2,N3)) in idx order."""
    with open(path) as f:
        N1, N2, N3, _ = (int(x) for x in f.readline().split())
    d = np.genfromtxt(path, skip_header=1)
    d = d[d[:, 3] == slip_plane]
    idx = (N1 * N2 * d[:, 2] + N1 * d[:, 1] + d[:, 0]).astype(int)
    usfe = np.full(N1 * N2 * N3, np.nan)
    usfe[idx] = d[:, 4]
    return usfe, (N1, N2, N3)


def read_usfe_all(path):
    """Return ({l: usfe[N]} for l in 0,1,2, (N1,N2,N3)) so the viewer can switch
    slip planes without re-reading the file."""
    with open(path) as f:
        N1, N2, N3, _ = (int(x) for x in f.readline().split())
    d = np.genfromtxt(path, skip_header=1)
    out = {}
    for l in (0, 1, 2):
        dl = d[d[:, 3] == l]
        arr = np.full(N1 * N2 * N3, np.nan)
        if len(dl):
            idx = (N1 * N2 * dl[:, 2] + N1 * dl[:, 1] + dl[:, 0]).astype(int)
            arr[idx] = dl[:, 4]
        out[l] = arr
    return out, (N1, N2, N3)


def comp_columns_from_csv(csv_path):
    """Composition-column labels = reversed SFE-CSV header (matches the interp)."""
    with open(csv_path, "r", encoding="utf-8-sig") as f:
        header = f.readline()
    return [c.strip() for c in header.strip().split(",")][1:][::-1]


# ------------------------------ projection -------------------------------

def _miller(v, tol=1e-6):
    v = np.array(v, float)
    v[np.abs(v) < tol * np.max(np.abs(v))] = 0.0
    nz = np.abs(v[v != 0])
    if len(nz):
        v = v / nz.min()
    return np.round(v).astype(int)


def _dir_label(v):
    s = "".join((r"$\bar{%d}$" % abs(int(c)) if c < 0 else "%d" % int(c)) for c in v)
    return f"[{s}]"


def project(coords, x_dir, y_dir):
    """Project coords onto the viewing frame.  Returns U, V, Wp, planes and the
    three axis labels.  The out-of-plane (viewing) direction is X x Y."""
    xhat = x_dir / np.linalg.norm(x_dir)
    yhat = y_dir / np.linalg.norm(y_dir)
    what = np.cross(xhat, yhat)
    what /= np.linalg.norm(what)
    U = coords @ xhat
    V = coords @ yhat
    Wp = coords @ what
    planes = np.unique(np.round(Wp, 3))
    labels = (_dir_label(_miller(x_dir)),
              _dir_label(_miller(y_dir)),
              _dir_label(_miller(what)))
    return U, V, Wp, planes, labels


def frame_for(system, idx):
    """(x_dir, y_dir) for slip plane `idx` of `system`: x = b = [-111],
    y = normal x b, so the viewing direction is the plane normal."""
    normal = np.array(SLIP_PLANES[system][idx][0], float)
    return B_DIR.copy(), np.cross(normal, B_DIR)


# ================================ Tk GUI =================================

def run_viewer(master=None, folder=".", slip="{110}",
               lattice="", comp="", usfe="", csv=""):
    """Open the viewer.  If master is given it opens as a Toplevel child of the
    workflow GUI; otherwise it creates its own root window."""
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox
    import matplotlib
    matplotlib.use("TkAgg")
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import (
        FigureCanvasTkAgg, NavigationToolbar2Tk)
    from matplotlib.lines import Line2D
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (registers 3d)

    win = tk.Toplevel(master) if master is not None else tk.Tk()
    win.title("PILOT  -  lattice / composition / USFE viewer")
    win.geometry("1140x880")

    if slip not in SLIP_PLANES:
        slip = "{110}"

    state = {"coords": None, "elems": None, "comp": None, "comp_cols": None,
             "usfe_by_l": None, "cell": None, "2d": None, "3d": None,
             "loaded": False}

    pad = {"padx": 5, "pady": 2}
    left = ttk.Frame(win)
    left.pack(side="left", fill="y", **pad)
    right = ttk.Frame(win)
    right.pack(side="right", fill="both", expand=True, **pad)

    def plane_choices(system):
        return [f"{i}: {asc}" for i, (_n, asc, _tex)
                in enumerate(SLIP_PLANES[system])]

    v = {
        "lattice": tk.StringVar(value=lattice),
        "comp": tk.StringVar(value=comp),
        "usfe": tk.StringVar(value=usfe),
        "csv": tk.StringVar(value=csv),
        "slip": tk.StringVar(value=slip),
        "plane": tk.StringVar(value=plane_choices(slip)[0]),
        "mode": tk.StringVar(value="Composition (2-D plane)"),
        "color_elem": tk.StringVar(value=""),
        "plane_idx": tk.IntVar(value=0),
        "elev": tk.IntVar(value=22),
        "azim": tk.IntVar(value=180),
        "subsample": tk.DoubleVar(value=1.0),
        "comp_min": tk.StringVar(value=""),
        "comp_max": tk.StringVar(value=""),
        "usfe_min": tk.StringVar(value=""),
        "usfe_max": tk.StringVar(value=""),
    }

    def cur_l():
        try:
            return int(v["plane"].get().split(":")[0])
        except (ValueError, IndexError):
            return 0

    def cur_plane_tex():
        return SLIP_PLANES[v["slip"].get()][cur_l()][2]

    def _lim(minvar, maxvar):
        """Parse the min/max entry boxes; blank or invalid -> None (auto)."""
        def parse(s):
            s = s.strip()
            try:
                return float(s) if s else None
            except ValueError:
                return None
        return parse(minvar.get()), parse(maxvar.get())

    # ---- file pickers ----
    ff = ttk.LabelFrame(left, text="Files")
    ff.pack(fill="x", **pad)

    def _row(parent, label, var, r):
        ttk.Label(parent, text=label).grid(row=r, column=0, sticky="w", **pad)
        ttk.Entry(parent, textvariable=var, width=26).grid(row=r, column=1, **pad)

        def pick():
            p = filedialog.askopenfilename(initialdir=folder)
            if p:
                var.set(p)
        ttk.Button(parent, text="...", width=3, command=pick).grid(
            row=r, column=2, **pad)

    _row(ff, "Lattice", v["lattice"], 0)
    _row(ff, "Composition", v["comp"], 1)
    _row(ff, "USFE", v["usfe"], 2)
    _row(ff, "SFE CSV (opt.)", v["csv"], 3)

    idxf = ttk.Frame(ff)
    idxf.grid(row=4, column=0, columnspan=3, sticky="w", **pad)
    ttk.Label(idxf, text="or fill index").pack(side="left")
    idx_var = tk.StringVar(value="0")
    ttk.Entry(idxf, textvariable=idx_var, width=6).pack(side="left", padx=3)

    def fill_by_index():
        i = idx_var.get().strip()
        v["lattice"].set(os.path.join(folder, f"lattice_{i}"))
        v["comp"].set(os.path.join(folder, f"comp_{i}"))
        v["usfe"].set(os.path.join(folder, f"USFE_{i}"))
    ttk.Button(idxf, text="fill", width=5, command=fill_by_index).pack(side="left")

    # ---- slip system + plane ----
    sf = ttk.LabelFrame(left, text="Slip system & plane")
    sf.pack(fill="x", **pad)
    ttk.Label(sf, text="System").grid(row=0, column=0, sticky="w", **pad)
    slip_box = ttk.Combobox(sf, textvariable=v["slip"], values=list(SLIP_PLANES),
                            state="readonly", width=8)
    slip_box.grid(row=0, column=1, sticky="w", **pad)
    ttk.Label(sf, text="Slip plane").grid(row=1, column=0, sticky="w", **pad)
    plane_box = ttk.Combobox(sf, textvariable=v["plane"],
                             values=plane_choices(slip), state="readonly", width=12)
    plane_box.grid(row=1, column=1, sticky="w", **pad)
    ttk.Label(sf, text="(index -> slip plane; sets USFE column & 2-D slice)",
              foreground="#555", wraplength=250).grid(
        row=2, column=0, columnspan=2, sticky="w", **pad)

    ttk.Button(left, text="Load files", command=lambda: load()).pack(
        fill="x", **pad)

    # ---- view mode ----
    mf = ttk.LabelFrame(left, text="View")
    mf.pack(fill="x", **pad)
    ttk.Combobox(mf, textvariable=v["mode"], state="readonly", width=26,
                 values=["Composition (2-D plane)", "USFE (2-D plane)",
                         "Atom types (3-D cell)", "USFE (3-D cell)"]).pack(**pad)
    ttk.Label(mf, text="Color by composition").pack(anchor="w", **pad)
    color_box = ttk.Combobox(mf, textvariable=v["color_elem"], state="readonly",
                             width=12)
    color_box.pack(anchor="w", **pad)
    color_note = ttk.Label(mf, text="", foreground="#555", wraplength=250)
    color_note.pack(anchor="w", **pad)

    # ---- color scale (manual limits) ----
    cs = ttk.LabelFrame(left, text="Color scale  (blank = auto)")
    cs.pack(fill="x", **pad)
    ttk.Label(cs, text="Composition").grid(row=0, column=0, sticky="w", **pad)
    ttk.Label(cs, text="min").grid(row=0, column=1, sticky="e")
    ttk.Entry(cs, textvariable=v["comp_min"], width=8).grid(row=0, column=2, **pad)
    ttk.Label(cs, text="max").grid(row=0, column=3, sticky="e")
    ttk.Entry(cs, textvariable=v["comp_max"], width=8).grid(row=0, column=4, **pad)
    ttk.Label(cs, text="USFE").grid(row=1, column=0, sticky="w", **pad)
    ttk.Label(cs, text="min").grid(row=1, column=1, sticky="e")
    ttk.Entry(cs, textvariable=v["usfe_min"], width=8).grid(row=1, column=2, **pad)
    ttk.Label(cs, text="max").grid(row=1, column=3, sticky="e")
    ttk.Entry(cs, textvariable=v["usfe_max"], width=8).grid(row=1, column=4, **pad)
    ttk.Button(cs, text="Apply scale", command=lambda: draw(full=True)).grid(
        row=2, column=0, columnspan=5, sticky="w", **pad)
    ttk.Label(cs, text="(ternary RGB triangle ignores these)",
              foreground="#555").grid(row=3, column=0, columnspan=5,
                                      sticky="w", **pad)

    # ---- 2-D controls ----
    twod = ttk.LabelFrame(left, text="2-D plane")
    twod.pack(fill="x", **pad)
    ttk.Label(twod, text="Atomic plane index (along the slip-plane normal)").pack(
        anchor="w", **pad)
    plane_scale = ttk.Scale(twod, from_=0, to=1, orient="horizontal",
                            variable=v["plane_idx"])
    plane_scale.pack(fill="x", **pad)
    plane_lbl = ttk.Label(twod, text="plane 0")
    plane_lbl.pack(anchor="w", **pad)

    # ---- 3-D controls ----
    threed = ttk.LabelFrame(left, text="3-D view")
    threed.pack(fill="x", **pad)
    ttk.Label(threed, text="elevation").pack(anchor="w", **pad)
    ttk.Scale(threed, from_=-90, to=90, orient="horizontal",
              variable=v["elev"]).pack(fill="x", **pad)
    ttk.Label(threed, text="azimuth").pack(anchor="w", **pad)
    ttk.Scale(threed, from_=-180, to=180, orient="horizontal",
              variable=v["azim"]).pack(fill="x", **pad)
    ttk.Label(threed, text="subsample (0-1, lower=faster)").pack(anchor="w", **pad)
    ttk.Scale(threed, from_=0.05, to=1.0, orient="horizontal",
              variable=v["subsample"]).pack(fill="x", **pad)
    ttk.Label(threed, text="drag the plot to rotate too").pack(anchor="w", **pad)

    status = tk.StringVar(value="Choose files and click Load.")
    ttk.Label(left, textvariable=status, wraplength=270,
              foreground="#555").pack(fill="x", **pad)

    # ---- matplotlib canvas ----
    fig = Figure(figsize=(7.6, 7.2))
    canvas = FigureCanvasTkAgg(fig, master=right)
    canvas.get_tk_widget().pack(fill="both", expand=True)
    NavigationToolbar2Tk(canvas, right)

    rng = np.random.default_rng(0)

    def _project_for(idx):
        xd, yd = frame_for(v["slip"].get(), idx)
        return project(state["coords"], xd, yd)

    def reproject():
        """2-D uses the selected plane's frame; 3-D uses a fixed per-system frame
        (plane 0), so switching the USFE index never reorients the 3-D cell."""
        state["2d"] = _project_for(cur_l())
        state["3d"] = _project_for(0)
        planes = state["2d"][3]
        plane_scale.configure(to=max(len(planes) - 1, 1))
        if v["plane_idx"].get() > len(planes) - 1:
            v["plane_idx"].set(len(planes) // 2)

    def _update_color_controls():
        cols = state["comp_cols"] or []
        ncomp = len(cols)
        if ncomp == 3:
            color_box.configure(state="disabled")
            color_note.configure(text="Ternary: colored by RGB triangle "
                                      "(one element per corner).")
        else:
            color_box.configure(state="readonly")
            color_note.configure(
                text="Binary/4+: colored by one element fraction."
                if ncomp else "")

    def load():
        try:
            latp, compp, usfep = v["lattice"].get(), v["comp"].get(), v["usfe"].get()
            for lbl, p in (("lattice", latp), ("composition", compp),
                           ("USFE", usfep)):
                if not p or not os.path.exists(p):
                    raise ValueError(f"{lbl} file not found: {p}")
            coords, elems = read_lattice(latp)
            ncols = comp_ncols(compp)
            comp = read_comp(compp, ncols)
            usfe_by_l, cell = read_usfe_all(usfep)
            n = cell[0] * cell[1] * cell[2]
            if not (len(coords) == len(comp) == len(usfe_by_l[0]) == n):
                raise ValueError(
                    f"file sizes disagree: coords={len(coords)}, comp={len(comp)}, "
                    f"usfe={len(usfe_by_l[0])}, N1*N2*N3={n}")

            if v["csv"].get() and os.path.exists(v["csv"].get()):
                cols = comp_columns_from_csv(v["csv"].get())
                if len(cols) != ncols:
                    cols = [f"col{c}" for c in range(ncols)]
            else:
                cols = [f"col{c}" for c in range(ncols)]

            state.update(coords=coords, elems=elems, comp=comp, comp_cols=cols,
                         usfe_by_l=usfe_by_l, cell=cell, loaded=True)

            color_box["values"] = cols
            if v["color_elem"].get() not in cols:
                v["color_elem"].set(cols[-1] if cols else "")
            _update_color_controls()
            reproject()
            status.set(f"Loaded {cell[0]}x{cell[1]}x{cell[2]} = {n} sites, "
                       f"{len(cols)}-component. Composition is the slip-plane-0 "
                       f"local average.")
            draw(full=True)
        except Exception as e:            # noqa: BLE001
            messagebox.showerror("Viewer", str(e))
            status.set("Error: " + str(e))

    def _draw_triangle_legend(triax, labels3):
        """Fill a dedicated (external) axes with the ternary RGB triangle."""
        n = 44
        xs, ys, cs = [], [], []
        h = np.sqrt(3) / 2
        for a in range(n + 1):
            for b in range(n + 1 - a):
                c0 = a / n
                c1 = b / n
                c2 = 1.0 - c0 - c1
                xs.append(c1 * 0.0 + c2 * 1.0 + c0 * 0.5)
                ys.append(c0 * h)
                cs.append((c0, c1, c2))
        triax.scatter(xs, ys, c=np.clip(cs, 0, 1), s=10, marker="o", linewidths=0)
        triax.text(0.5, h + 0.06, labels3[0], ha="center", va="bottom", fontsize=10)
        triax.text(-0.04, -0.06, labels3[1], ha="right", va="top", fontsize=10)
        triax.text(1.04, -0.06, labels3[2], ha="left", va="top", fontsize=10)
        triax.set_xlim(-0.35, 1.35)
        triax.set_ylim(-0.2, h + 0.22)
        triax.set_aspect("equal")
        triax.set_title("composition", fontsize=10)
        triax.axis("off")

    def _style_3d(ax, title, labels, U, V, Wp):
        ax.set_xlabel(labels[0], labelpad=10)
        ax.set_ylabel(labels[1], labelpad=10)
        ax.set_zlabel(labels[2], labelpad=10)
        ax.set_title(title, pad=10)
        try:
            ax.set_box_aspect((np.ptp(U), np.ptp(V), np.ptp(Wp)))
        except Exception:
            pass

    def draw(full=True):
        if not state["loaded"]:
            return
        mode = v["mode"].get()
        plane_tex = cur_plane_tex()
        l = cur_l()
        is3d = "3-D" in mode

        # fast path: 3-D angle change only -> just move the camera
        if is3d and not full and getattr(fig, "_is3d", False) and fig.axes:
            fig.axes[0].view_init(elev=v["elev"].get(), azim=v["azim"].get())
            canvas.draw_idle()
            return

        fig.clear()
        fig._is3d = is3d

        if not is3d:
            U, V, Wp, planes, labels = state["2d"]
            p = int(v["plane_idx"].get())
            p = max(0, min(p, len(planes) - 1))
            mask = np.abs(Wp - planes[p]) <= 0.05
            plane_lbl.configure(text=f"plane {p}/{len(planes)-1} "
                                     f"({int(mask.sum())} atoms)")

            if mode.startswith("Composition"):
                cols = state["comp_cols"]
                ncomp = len(cols)
                if ncomp == 3:
                    # main scatter + a separate triangle panel (outside the box)
                    gs = fig.add_gridspec(1, 2, width_ratios=[5, 1.5])
                    ax = fig.add_subplot(gs[0])
                    triax = fig.add_subplot(gs[1])
                    rgb = np.clip(state["comp"][:, :3], 0, 1)
                    ax.scatter(U[mask], V[mask], c=rgb[mask], s=80, marker="H",
                               edgecolors="k", linewidths=0.3)
                    _draw_triangle_legend(triax, cols)
                    ax.set_title(f"Local composition (ternary RGB)  -  "
                                 f"{plane_tex} slice  -  plane {p}/{len(planes)-1}")
                else:
                    ax = fig.add_subplot(111)
                    vmin, vmax = _lim(v["comp_min"], v["comp_max"])
                    col = v["color_elem"].get()
                    ci = cols.index(col) if col in cols else len(cols) - 1
                    field = state["comp"][:, ci]
                    sc = ax.scatter(U[mask], V[mask], c=field[mask], cmap="rainbow",
                                    vmin=vmin, vmax=vmax, s=80, marker="H",
                                    edgecolors="k", linewidths=0.3)
                    cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
                    cb.set_label(f"{cols[ci]} fraction")
                    ax.set_title(f"Local {cols[ci]} fraction  -  {plane_tex} slice "
                                 f" -  plane {p}/{len(planes)-1}")
            else:
                ax = fig.add_subplot(111)
                vmin, vmax = _lim(v["usfe_min"], v["usfe_max"])
                field = state["usfe_by_l"][l]
                sc = ax.scatter(U[mask], V[mask], c=field[mask], cmap="viridis",
                                vmin=vmin, vmax=vmax, s=80, marker="H",
                                edgecolors="k", linewidths=0.3)
                cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
                cb.set_label(r"USFE_%d (mJ/m$^2$)" % l)
                ax.set_title(f"USFE_{l} on {plane_tex}  -  plane {p}/{len(planes)-1}")
            ax.set_aspect("equal")
            ax.set_xlabel(labels[0])
            ax.set_ylabel(labels[1])
        else:
            U, V, Wp, planes, labels = state["3d"]     # fixed per-system frame
            f = float(v["subsample"].get())
            keep = (np.ones(len(U), bool) if f >= 1.0 else rng.random(len(U)) < f)
            ax = fig.add_subplot(111, projection="3d")
            if mode.startswith("Atom"):
                elems = state["elems"]
                uniq = sorted(set(elems.tolist()))
                cmap = plt.get_cmap("tab10")
                cmap_e = {e: cmap(i % 10) for i, e in enumerate(uniq)}
                cvec = np.array([cmap_e[e] for e in elems[keep]])
                ax.scatter(U[keep], V[keep], Wp[keep], c=cvec, s=25,
                           depthshade=True, edgecolors="k", linewidths=0.3)
                handles = [Line2D([], [], ls="", marker="o", markersize=7,
                                  markerfacecolor=cmap_e[e], markeredgecolor="k",
                                  label=f"{e} ({int((elems==e).sum())})")
                           for e in uniq]
                ax.legend(handles=handles, loc="upper left", fontsize=9)
                _style_3d(ax, "3-D cell - atom types", labels, U, V, Wp)
            else:
                vmin, vmax = _lim(v["usfe_min"], v["usfe_max"])
                usfe = state["usfe_by_l"][l]
                m = keep & ~np.isnan(usfe)
                sc = ax.scatter(U[m], V[m], Wp[m], c=usfe[m], cmap="viridis",
                                vmin=vmin, vmax=vmax, s=25, depthshade=True,
                                edgecolors="k", linewidths=0.3)
                cb = fig.colorbar(sc, ax=ax, fraction=0.03, pad=0.10)
                cb.set_label(r"USFE_%d (mJ/m$^2$)" % l)
                _style_3d(ax, f"3-D cell - USFE_{l} on {plane_tex}",
                          labels, U, V, Wp)
            ax.view_init(elev=v["elev"].get(), azim=v["azim"].get())
        canvas.draw_idle()

    # ---- wire controls ----
    def on_system(*_):
        plane_box.configure(values=plane_choices(v["slip"].get()))
        v["plane"].set(plane_choices(v["slip"].get())[0])   # triggers on_plane
        if state["loaded"]:
            reproject()
            draw(full=True)

    def on_plane(*_):
        # only the 2-D slice frame depends on the plane; the 3-D frame is fixed,
        # so the 3-D camera stays put and only the USFE coloring changes.
        if state["loaded"]:
            reproject()
            draw(full=True)

    v["slip"].trace_add("write", on_system)
    v["plane"].trace_add("write", on_plane)
    v["mode"].trace_add("write", lambda *_: draw(full=True))
    v["color_elem"].trace_add("write", lambda *_: draw(full=True))
    v["plane_idx"].trace_add("write", lambda *_: draw(full=True))
    v["subsample"].trace_add("write", lambda *_: draw(full=True))
    v["elev"].trace_add("write", lambda *_: draw(full=False))
    v["azim"].trace_add("write", lambda *_: draw(full=False))
    for key in ("comp_min", "comp_max", "usfe_min", "usfe_max"):
        # apply on Enter as well as via the Apply button
        pass

    # bind Enter in the scale boxes to redraw
    for child in cs.winfo_children():
        if isinstance(child, ttk.Entry):
            child.bind("<Return>", lambda _e: draw(full=True))

    if lattice and comp and usfe:
        win.after(200, load)

    if master is None:
        win.mainloop()
    return win


def main():
    folder = sys.argv[1] if len(sys.argv) > 1 else "."
    run_viewer(master=None, folder=folder)


if __name__ == "__main__":
    main()
