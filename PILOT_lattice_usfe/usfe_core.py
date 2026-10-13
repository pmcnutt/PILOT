#!/usr/bin/env python3
"""
usfe_core.py -- the work behind the PILOT USFE workflow (no GUI dependencies).

One pass = lattice_generator builds a random solid-solution (RSS) lattice, then
the matching interpolation script maps each site's local composition to a USFE
value and writes <lattice>_usfe.  build_lattices() repeats that for as many
lattices as you ask for and records the mean/std of each.

Two things are now automatic, so the user cannot mismatch them:

    lattice_type BCC              -> {110} interpolation (110_interp_usfe.py)
    lattice_type CUSTOM_NON_ORTHO -> {112} interpolation (112_interp_usfe_nonbcc.py)

Only random solid solutions are produced; there is no short-range-order option.

This module has no tkinter import, so it can be tested and reused on its own.
"""

import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from typing import Callable, List, Optional, Tuple

import numpy as np


# lattice type -> (interpolation script, human-readable slip system)
LATTICE_INTERP = {
    "BCC":              ("110_interp_usfe.py",        "{110}"),
    "CUSTOM_NON_ORTHO": ("112_interp_usfe_nonbcc.py", "{112}"),
}


# ------------------------------ configuration ----------------------------

@dataclass
class BuildConfig:
    # --- lattice ---
    lattice_type: str                      # BCC | CUSTOM_NON_ORTHO
    cell: Tuple[int, int, int]             # N1, N2, N3
    lattice_param: float                   # nearest-neighbour distance
    labels: List[str]                      # element names
    conc: List[float]                      # concentrations (need not pre-sum to 1)

    # --- how many / naming ---
    n_lat: int = 1                         # number of lattices to build
    index_start: int = 0                   # first index value
    index_step: int = 1                    # step between index values
    usfe_pattern: str = "USFE_{i}"         # final USFE filename; must contain {i}
    lattice_pattern: str = "lattice_{i}"   # lattice filename; must contain {i}
    comp_pattern: str = "comp_{i}"         # composition filename; must contain {i}
    keep_lattices: bool = True             # keep the RSS lattice files
    keep_comp: bool = True                 # keep the per-site composition files
    out_dir: str = "."

    # --- interpolation ---
    csv: str = "MoNbTi_sfe_data.csv"       # SFE data; header names the elements
    interp: str = ""                       # auto-set from lattice_type (see resolve)

    # --- tooling / misc ---
    generator: str = "lattice_generator.py"
    usfe_col: int = 4
    usfe_skip_header: int = 1
    stats_file: str = "usfe_stats.csv"

    def indices(self) -> List[int]:
        return [self.index_start + self.index_step * n for n in range(self.n_lat)]

    def slip_system(self) -> str:
        return LATTICE_INTERP.get(self.lattice_type, ("", "?"))[1]

    def resolve_interp(self) -> str:
        """The interpolation script implied by the lattice type."""
        if self.lattice_type in LATTICE_INTERP:
            return LATTICE_INTERP[self.lattice_type][0]
        return self.interp


# ------------------------------ small helpers ----------------------------

def build_generator_input_text(cfg: BuildConfig, out_name: str) -> str:
    """Text of a lattice_generator input file for this config + output name."""
    N1, N2, N3 = cfg.cell
    nc = len(cfg.labels)
    lines = [
        f"{cfg.lattice_type} #lattice type (BCC or CUSTOM_NON_ORTHO)",
        f"{N1} {N2} {N3} #number of primitive cells in each direction",
        f"{cfg.lattice_param} #distance between nearest neighbors",
        f"{nc} #number of elements",
        " ".join(cfg.labels) + " #names of elements",
        " ".join(f"{c:.8f}" for c in cfg.conc) + " #concentration of each element",
        f"{out_name} #output file name",
    ]
    return "\n".join(lines) + "\n"


def usfe_stats(path: str, col: int, skip: int) -> Tuple[float, float, int]:
    """Mean / std / N of the USFE column."""
    data = np.genfromtxt(path, usecols=[col], skip_header=skip)
    return float(np.mean(data)), float(np.std(data)), int(data.size)


def stream_process(cmd: List[str], emit: Callable[[str, bool], None]) -> None:
    """Run cmd, forwarding output as it appears.  emit(text, is_status):
    is_status True for carriage-return progress updates, False for finished
    lines.  Raises RuntimeError on nonzero exit."""
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
    if cfg.lattice_type not in LATTICE_INTERP:
        errs.append("Lattice type must be BCC or CUSTOM_NON_ORTHO.")
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
    if cfg.n_lat <= 0:
        errs.append("Number of lattices must be at least 1.")
    for label, pat in (("USFE", cfg.usfe_pattern),
                       ("Lattice", cfg.lattice_pattern),
                       ("Composition", cfg.comp_pattern)):
        if "{i}" not in pat:
            errs.append(f"{label} naming pattern must contain '{{i}}'.")
    interp = cfg.resolve_interp()
    for label, p in (("Lattice generator", cfg.generator),
                     ("Interpolation script", interp),
                     ("CSV file", cfg.csv)):
        if not p or not os.path.exists(p):
            errs.append(f"{label} not found: {p}")
    # naming collisions across the run
    idxs = cfg.indices()
    for label, pat in (("USFE", cfg.usfe_pattern),
                       ("lattice", cfg.lattice_pattern),
                       ("composition", cfg.comp_pattern)):
        names = [pat.format(i=i) for i in idxs]
        if len(set(names)) != len(names):
            errs.append(f"{label} naming pattern + start/step produce duplicate "
                        f"names.")
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
        status(msg)   -> a transient one-line status (generator progress etc.)
        progress(f)   -> overall fraction in [0, 1]
        should_stop() -> return True to abort cleanly between steps
    Returns a list of (index, usfe_file, lattice_file, comp_file, mean, std, n)."""
    problems = validate(cfg)
    if problems:
        raise ValueError("Cannot start:\n  - " + "\n  - ".join(problems))

    interp = cfg.resolve_interp()
    log(f"Slip system {cfg.slip_system()}  ->  {interp}")

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
        comp_name = cfg.comp_pattern.format(i=i)
        log(f"[{n+1}/{len(idxs)}] lattice_generator -> {lattice_name}")

        # per-lattice generator input (temp file in cwd, removed afterwards)
        with tempfile.NamedTemporaryFile("w", suffix=f"_geninput_{i}",
                                         delete=False, dir=".") as tf:
            gen_in = tf.name
            tf.write(build_generator_input_text(cfg, lattice_name))

        # where the interp writes its per-site composition file.  If we're
        # keeping it, write straight into out_dir with the comp_{i} name; if
        # not, write to a temp name in cwd and delete it afterwards.
        if cfg.keep_comp:
            comp_out = os.path.join(cfg.out_dir, comp_name)
        else:
            comp_out = f".__comp_tmp_{i}"

        try:
            stream_process(["python3", cfg.generator, gen_in], emit)
            if not os.path.exists(lattice_name):
                raise RuntimeError(f"generator did not produce {lattice_name}")
            step += 1
            progress(step / total_steps)

            log(f"[{n+1}/{len(idxs)}] interp ({os.path.basename(interp)}) "
                f"-> {usfe_name}")
            stream_process(["python3", interp, lattice_name,
                            str(N1), str(N2), str(N3), cfg.csv, comp_out], emit)
            interp_out = f"{lattice_name}_usfe"
            if not os.path.exists(interp_out):
                raise RuntimeError(f"interp did not produce {interp_out}")

            usfe_final = os.path.join(cfg.out_dir, usfe_name)
            shutil.move(interp_out, usfe_final)

            mean, std, npts = usfe_stats(usfe_final, cfg.usfe_col,
                                         cfg.usfe_skip_header)
            comp_final = comp_out if cfg.keep_comp else ""
            log(f"       {usfe_name}: mean={mean:.4f}  std={std:.4f}  (N={npts})")
            stats.append((i, usfe_name, lattice_name, comp_name,
                          mean, std, npts))

            # lattice file: keep (move to out_dir) or discard
            if cfg.keep_lattices:
                dest = os.path.join(cfg.out_dir, lattice_name)
                if os.path.abspath(dest) != os.path.abspath(lattice_name):
                    shutil.move(lattice_name, dest)
            elif os.path.exists(lattice_name):
                os.remove(lattice_name)

            step += 1
            progress(step / total_steps)
        finally:
            if os.path.exists(gen_in):
                os.remove(gen_in)
            if not cfg.keep_comp and os.path.exists(comp_out):
                os.remove(comp_out)

    # summary file
    if stats:
        stats_path = os.path.join(cfg.out_dir, cfg.stats_file)
        with open(stats_path, "w") as f:
            f.write("index,usfe_file,lattice_file,comp_file,"
                    "mean_usfe,std_usfe,n_points\n")
            for i, uf, lat, cmp_, mean, std, npts in stats:
                f.write(f"{i},{uf},{lat},{cmp_},{mean:.6f},{std:.6f},{npts}\n")
        log(f"Wrote {len(stats)} USFE files + summary -> {stats_path}")
    return stats
