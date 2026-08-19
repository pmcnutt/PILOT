#!/usr/bin/env python3
"""
usfe_core.py -- the actual work behind the USFE workflow (no GUI dependencies).

One pass = OTIS builds a random SRO lattice, then an interp script maps local
composition to USFE and writes <lattice>_usfe.  build_lattices() repeats that
for however many lattices you ask for and records mean/std of each.

This module is imported by usfe_workflow_gui.py; it has no tkinter import so it
can be tested and reused on its own.
"""

import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from typing import Callable, List, Optional, Tuple

import numpy as np


# ------------------------------ configuration ----------------------------

@dataclass
class BuildConfig:
    # --- lattice (OTIS) ---
    lattice_type: str                      # BCC | FCC | CUSTOM
    cell: Tuple[int, int, int]             # N1, N2, N3
    lattice_param: float                   # nearest-neighbour distance
    labels: List[str]                      # element names
    conc: List[float]                      # concentrations (need not pre-sum to 1)
    input_type: str = "prob"               # prob | wc
    matrix: Optional[List[List[float]]] = None   # components x components; None -> random

    # --- how many / naming ---
    n_lat: int = 1                         # number of lattices to build
    index_start: int = 0                   # first index value
    index_step: int = 1                    # step between index values
    usfe_pattern: str = "USFE_{i}"         # final USFE filename; must contain {i}
    lattice_pattern: str = "lattice_{i}"   # OTIS lattice filename; must contain {i}
    keep_lattices: bool = False
    out_dir: str = "."

    # --- interpolation ---
    interp: str = "110_interp_usfe.py"
    csv: str = "MoNbTi_sfe_data.csv"

    # --- tooling / misc ---
    otis: str = "otis.py"
    usfe_col: int = 4
    usfe_skip_header: int = 1
    stats_file: str = "usfe_stats.csv"
    interp_side_files: List[str] = field(default_factory=lambda: [
        "monbti_110_comp_values.txt",
        "monbti_112_1323K_comp_values.txt",
    ])

    def indices(self) -> List[int]:
        return [self.index_start + self.index_step * n for n in range(self.n_lat)]


# ------------------------------ small helpers ----------------------------

def build_otis_input_text(cfg: BuildConfig, out_name: str) -> str:
    """Return the text of an OTIS input file for this config + output name."""
    nc = len(cfg.labels)
    if cfg.matrix is None:                       # random solid solution
        s = sum(cfg.conc)
        row = [c / s for c in cfg.conc]
        matrix = [row[:] for _ in range(nc)]
        input_type = "prob"
    else:
        matrix = cfg.matrix
        input_type = cfg.input_type

    N1, N2, N3 = cfg.cell
    lines = [
        f"{cfg.lattice_type} #lattice type (BCC, FCC, or CUSTOM)",
        f"{N1} {N2} {N3} #number of primitive cells in each direction",
        f"{cfg.lattice_param} #distance between nearest neighbors",
        f"{nc} #number of elements",
        " ".join(cfg.labels) + " #names of elements",
        " ".join(f"{c:.8f}" for c in cfg.conc) + " #concentration of each element",
        f"{input_type} #prob (pij) or wc (Warren-Cowley)",
    ]
    for r in matrix:
        lines.append(" ".join(f"{v:.8f}" for v in r))
    lines.append(f"{out_name} #output file name")
    return "\n".join(lines) + "\n"


def usfe_stats(path: str, col: int, skip: int) -> Tuple[float, float, int]:
    """Mean / std / N of the USFE column (same approach as avg_usfe.py)."""
    data = np.genfromtxt(path, usecols=[col], skip_header=skip)
    return float(np.mean(data)), float(np.std(data)), int(data.size)


def stream_process(cmd: List[str], emit: Callable[[str, bool], None]) -> None:
    """Run cmd, forwarding output as it appears.  emit(text, is_status): is_status
    is True for carriage-return progress updates (e.g. OTIS's bar), False for
    completed lines.  Raises RuntimeError on nonzero exit."""
    env = {**os.environ, "PYTHONUNBUFFERED": "1"}
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, env=env)
    buf = ""
    assert proc.stdout is not None
    while True:
        ch = proc.stdout.read(1)
        if ch == "":
            break
        if ch == "\n":
            emit(buf, False)
            buf = ""
        elif ch == "\r":
            emit(buf, True)
            buf = ""
        else:
            buf += ch
    if buf:
        emit(buf, False)
    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError(f"command failed ({proc.returncode}): {' '.join(cmd)}")


# ------------------------------ validation -------------------------------

def validate(cfg: BuildConfig) -> List[str]:
    """Return a list of human-readable problems; empty list means OK."""
    errs = []
    if cfg.lattice_type not in ("BCC", "FCC", "CUSTOM"):
        errs.append("Lattice type must be BCC, FCC, or CUSTOM.")
    if any(int(x) <= 0 for x in cfg.cell):
        errs.append("Cell dimensions N1 N2 N3 must all be positive integers.")
    if len(cfg.labels) == 0:
        errs.append("At least one element is required.")
    if len(cfg.labels) != len(cfg.conc):
        errs.append(f"Got {len(cfg.labels)} elements but {len(cfg.conc)} "
                    f"concentrations -- counts must match.")
    if cfg.conc and abs(sum(cfg.conc) - 1.0) > 0.1:
        errs.append(f"Concentrations sum to {sum(cfg.conc):.3f}; must be ~1 "
                    f"(within 0.1). They are renormalised automatically.")
    if cfg.matrix is not None:
        nc = len(cfg.labels)
        if len(cfg.matrix) != nc or any(len(r) != nc for r in cfg.matrix):
            errs.append(f"SRO matrix must be {nc}x{nc} to match the elements.")
        if cfg.input_type not in ("prob", "wc"):
            errs.append("SRO input type must be 'prob' or 'wc'.")
    if cfg.n_lat <= 0:
        errs.append("Number of lattices must be at least 1.")
    if "{i}" not in cfg.usfe_pattern:
        errs.append("USFE naming pattern must contain '{i}'.")
    if cfg.keep_lattices and "{i}" not in cfg.lattice_pattern:
        errs.append("Lattice naming pattern must contain '{i}'.")
    for label, p in (("OTIS script", cfg.otis), ("interp script", cfg.interp),
                     ("CSV file", cfg.csv)):
        if not os.path.exists(p):
            errs.append(f"{label} not found: {p}")
    # naming collisions across the run
    idxs = cfg.indices()
    names = [cfg.usfe_pattern.format(i=i) for i in idxs]
    if len(set(names)) != len(names):
        errs.append("USFE naming pattern + start/step produce duplicate names.")
    return errs


# ------------------------------ main routine -----------------------------

def _noop(*_a, **_k):
    pass


def build_lattices(cfg: BuildConfig,
                   log: Callable[[str], None] = _noop,
                   status: Callable[[str], None] = _noop,
                   progress: Callable[[float], None] = _noop,
                   should_stop: Callable[[], bool] = lambda: False):
    """Run the full workflow.  Callbacks:
        log(msg)      -> a completed log line
        status(msg)   -> a transient one-line status (OTIS progress etc.)
        progress(f)   -> overall fraction in [0, 1]
        should_stop() -> return True to abort cleanly between steps
    Returns the list of (index, usfe_file, lattice_file, mean, std, n)."""
    problems = validate(cfg)
    if problems:
        raise ValueError("Cannot start:\n  - " + "\n  - ".join(problems))

    os.makedirs(cfg.out_dir, exist_ok=True)
    idxs = cfg.indices()
    N1, N2, N3 = cfg.cell
    total_steps = len(idxs) * 2
    step = 0
    stats = []

    def emit(text, is_status):
        (status if is_status else log)(text)

    for n, i in enumerate(idxs):
        if should_stop():
            log("Stopped by user.")
            break

        usfe_name = cfg.usfe_pattern.format(i=i)
        lattice_name = cfg.lattice_pattern.format(i=i)
        log(f"[{n+1}/{len(idxs)}] OTIS -> {lattice_name}")

        # per-lattice OTIS input (temp file in cwd, removed afterwards)
        with tempfile.NamedTemporaryFile("w", suffix=f"_otis_{i}",
                                         delete=False, dir=".") as tf:
            otis_in = tf.name
            tf.write(build_otis_input_text(cfg, lattice_name))
        try:
            stream_process(["python3", cfg.otis, otis_in], emit)
            if not os.path.exists(lattice_name):
                raise RuntimeError(f"OTIS did not produce {lattice_name}")
            step += 1
            progress(step / total_steps)

            log(f"[{n+1}/{len(idxs)}] interp ({os.path.basename(cfg.interp)}) "
                f"-> {usfe_name}")
            stream_process(["python3", cfg.interp, lattice_name,
                            str(N1), str(N2), str(N3), cfg.csv], emit)
            interp_out = f"{lattice_name}_usfe"
            if not os.path.exists(interp_out):
                raise RuntimeError(f"interp did not produce {interp_out}")

            usfe_final = os.path.join(cfg.out_dir, usfe_name)
            shutil.move(interp_out, usfe_final)

            mean, std, npts = usfe_stats(usfe_final, cfg.usfe_col,
                                         cfg.usfe_skip_header)
            log(f"       {usfe_name}: mean={mean:.4f}  std={std:.4f}  (N={npts})")
            stats.append((i, usfe_name, lattice_name, mean, std, npts))

            if cfg.keep_lattices:
                dest = os.path.join(cfg.out_dir, lattice_name)
                if os.path.abspath(dest) != os.path.abspath(lattice_name):
                    shutil.move(lattice_name, dest)
            elif os.path.exists(lattice_name):
                os.remove(lattice_name)

            step += 1
            progress(step / total_steps)
        finally:
            for junk in [otis_in, *cfg.interp_side_files]:
                if os.path.exists(junk):
                    os.remove(junk)

    # summary file
    if stats:
        stats_path = os.path.join(cfg.out_dir, cfg.stats_file)
        with open(stats_path, "w") as f:
            f.write("index,usfe_file,lattice_file,mean_usfe,std_usfe,n_points\n")
            for i, uf, lat, mean, std, npts in stats:
                f.write(f"{i},{uf},{lat},{mean:.6f},{std:.6f},{npts}\n")
        log(f"Wrote {len(stats)} USFE files + summary -> {stats_path}")
    return stats
