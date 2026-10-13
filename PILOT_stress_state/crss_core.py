#!/usr/bin/env python3
"""
crss_core.py -- the engine behind the PILOT CRSS temperature sweep (no GUI deps).

A bisection search for the critical resolved shear stress (CRSS) of a screw
dislocation as a function of temperature, driving PFDD SLURM job arrays.  Both
the command-line driver (crss_search_alloy.py) and the graphical workflow
(stress_state_gui.py) call the functions here, so the search logic lives in one
place and behaves identically either way.

Everything for one search lives in a work directory (a study subfolder):
    <work_dir>/search_state.json     the search state (brackets, history)
    <work_dir>/search_state.undo.json snapshot stack powering undo
    <work_dir>/cases_<tag>.txt        the per-round cases file fed to SLURM

Key ideas carried over from the original driver:
  * The initial bracket [lo, hi] may be either sign but must not straddle zero
    (0 may sit at one end).  The regime (positive/negative) is inferred from the
    sign and the move/not-move logic flips to match.
  * Alloys: every case line carries the USFE file for that temperature, so the
    SLURM script can point PFDD at the right composition-resolved USFE field.
  * Two loadings: a pure shear couple on the MRSSP (chi about b), or a uniaxial
    stress along a chosen axis.  The scalar being bisected ("tau") is the load
    magnitude in both cases.
"""

import json
import math
import os
import random
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple


# ------------------------------ configuration ----------------------------

@dataclass
class SearchConfig:
    # --- where everything is written ---
    work_dir: str = "."

    # --- temperatures ---
    t_min: int = 0
    t_max: int = 1400
    t_step: int = 100

    # --- initial bracket (units of mu); the two ends in any order, must not
    #     straddle zero (0 may sit at one end) ---
    tau_lo: float = 0.0
    tau_hi: float = 0.08
    n_rounds: int = 6

    # --- extension when CRSS lies outside the initial range ---
    tau_extend_step: float = 0.02
    n_extend_rounds: int = 4

    # --- stress state ---
    stress_mode: str = "shear"          # "shear" | "uniaxial"
    chi_deg: float = 0.0                 # shear: MRSSP angle about b (0 -> (110))
    axis: Tuple[float, float, float] = (1.0, 0.0, 0.0)   # uniaxial loading axis

    # --- randomisation ---
    seed_min: int = 1
    seed_max: int = 1000

    # --- slip system / naming ---
    slip_system: str = "{110}"          # "{110}" | "{112}"
    tag: str = "Ta10W_110"              # used in cases/dump file names

    # --- alloy USFE linkage (per-temperature file the SLURM job feeds PFDD) ---
    usfe_pattern: str = "USFE_{T}"      # {T} -> temperature; or a constant name
    usfe_dir: str = ""                  # optional: local dir to check files exist

    # --- cluster-side output naming (recorded in the cases/history) ---
    out_base: str = "OUTPUT/dumps"      # cluster path where dumps will land
    dump_name: str = "Ta10W_screw_110_{T}K_{tag}.dump"

    # --- file names inside work_dir ---
    state_file: str = "search_state.json"
    undo_file: str = "search_state.undo.json"
    undo_depth: int = 12

    def temps(self) -> List[int]:
        return list(range(self.t_min, self.t_max + 1, self.t_step))

    def cases_file(self) -> str:
        return f"cases_{self.tag}.txt"

    # absolute paths within the work dir
    def state_path(self) -> str:
        return os.path.join(self.work_dir, self.state_file)

    def undo_path(self) -> str:
        return os.path.join(self.work_dir, self.undo_file)

    def cases_path(self) -> str:
        return os.path.join(self.work_dir, self.cases_file())


# ------------------------------ stress math ------------------------------

def _unit(v):
    n = math.sqrt(sum(c * c for c in v))
    return [c / n for c in v]


def stress_components_shear(tau0, chi_deg):
    """sigma = tau0 (b (x) m + m (x) b), b || [-1 1 1], m = (110) normal rotated
    by chi about b (positive toward b x n110).  Returns the 6 cube-axis
    components (sxx, syy, szz, sxy, sxz, syz)."""
    b = _unit((-1.0, 1.0, 1.0))
    n = _unit((1.0, 1.0, 0.0))                       # (110) normal, chi = 0
    bxn = (b[1] * n[2] - b[2] * n[1],
           b[2] * n[0] - b[0] * n[2],
           b[0] * n[1] - b[1] * n[0])
    c, s = math.cos(math.radians(chi_deg)), math.sin(math.radians(chi_deg))
    m = [c * n[i] + s * bxn[i] for i in range(3)]
    S = [[tau0 * (b[i] * m[j] + m[i] * b[j]) for j in range(3)] for i in range(3)]
    return (S[0][0], S[1][1], S[2][2], S[0][1], S[0][2], S[1][2])


def stress_components_uniaxial(mag, axis):
    """sigma = mag * (a (x) a), a = axis/|axis|.  Returns the 6 cube-axis
    components (sxx, syy, szz, sxy, sxz, syz)."""
    a = _unit(axis)
    S = [[mag * a[i] * a[j] for j in range(3)] for i in range(3)]
    return (S[0][0], S[1][1], S[2][2], S[0][1], S[0][2], S[1][2])


def tau_str(tau):
    return f"{tau:.6f}".rstrip("0").rstrip(".") or "0"


def tau_tag(tau):
    return tau_str(tau).replace(".", "p").replace("-", "m")


def sigma_line(cfg: SearchConfig, tau):
    """The full PFDD 'sigma initial ... delta ... #tau ...' line for this load."""
    if cfg.stress_mode == "uniaxial":
        c = stress_components_uniaxial(tau, cfg.axis)
        tail = ("#tau %s uniaxial axis %g %g %g"
                % (tau_str(tau), cfg.axis[0], cfg.axis[1], cfg.axis[2]))
    else:
        c = stress_components_shear(tau, cfg.chi_deg)
        tail = "#tau %s chi %+g deg" % (tau_str(tau), cfg.chi_deg)
    comps = " ".join(f"{x:.6f}" for x in c)
    return (f"sigma     initial {comps} "
            f"delta 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 "
            f"{tail}")


# ------------------------------ state I/O --------------------------------

class SearchError(Exception):
    """Raised for user-facing problems (no state, unrecorded rounds, etc.)."""


def load_state(cfg: SearchConfig):
    p = cfg.state_path()
    if not os.path.exists(p):
        raise SearchError(f"No search state found in {cfg.work_dir}. "
                          f"Initialise the search first.")
    with open(p) as f:
        return json.load(f)


def _push_undo(cfg: SearchConfig, prev_state):
    stack = []
    up = cfg.undo_path()
    if os.path.exists(up):
        try:
            with open(up) as f:
                stack = json.load(f)
        except (ValueError, OSError):
            stack = []
    stack.append(prev_state)
    stack = stack[-cfg.undo_depth:]
    with open(up, "w") as f:
        json.dump(stack, f)


def save_state(cfg: SearchConfig, state, snapshot=True):
    os.makedirs(cfg.work_dir, exist_ok=True)
    if snapshot and os.path.exists(cfg.state_path()):
        try:
            with open(cfg.state_path()) as f:
                _push_undo(cfg, json.load(f))
        except (ValueError, OSError):
            pass
    with open(cfg.state_path(), "w") as f:
        json.dump(state, f, indent=2)


# --------------------------- bracket helpers -----------------------------

def normalized_bracket(a, b):
    """Sort ends into (lo, hi), reject a straddle of zero, infer the regime."""
    lo, hi = (a, b) if a <= b else (b, a)
    if lo == hi:
        raise SearchError(f"Invalid bracket: the two ends must differ "
                          f"(got {a}, {b}).")
    if lo < 0 < hi:
        raise SearchError(f"Invalid bracket [{lo}, {hi}]: 0 must not lie strictly "
                          f"inside the bracket. Put 0 at an end, or keep both "
                          f"ends the same sign.")
    return lo, hi, ("positive" if hi > 0 else "negative")


def state_regime(state):
    r = state.get("regime")
    if r:
        return r
    any_d = next(iter(state["temps"].values()))
    return "negative" if any_d["hi"] <= 0 else "positive"


def active_temps(state):
    return {T: d for T, d in state["temps"].items() if d.get("active", True)}


def is_unbounded_above(d):
    h = d["history"]
    return bool(h) and all(x.get("moved") is False for x in h)


def is_unbounded_below(d):
    h = d["history"]
    return bool(h) and all(x.get("moved") is True for x in h)


def apply_result(d, tau, moved, regime):
    """Update one temperature's bracket from a move/not-move result."""
    if regime == "positive":
        if moved:
            d["hi"] = tau
        else:
            d["lo"] = tau
    else:  # negative
        if moved:
            d["lo"] = tau
        else:
            d["hi"] = tau


# ------------------------------- commands --------------------------------

def _noop(*_a, **_k):
    pass


def init(cfg: SearchConfig, log: Callable[[str], None] = _noop, force=False):
    """Create a fresh search_state.json in the work dir."""
    os.makedirs(cfg.work_dir, exist_ok=True)
    if os.path.exists(cfg.state_path()) and not force:
        raise SearchError(f"{cfg.state_file} already exists in {cfg.work_dir} -- "
                          f"delete it (or use force) to restart the search.")
    if os.path.exists(cfg.undo_path()):
        os.remove(cfg.undo_path())
    lo, hi, regime = normalized_bracket(cfg.tau_lo, cfg.tau_hi)
    state = {
        "round": 0,
        "chi_deg": cfg.chi_deg,
        "regime": regime,
        "n_rounds": cfg.n_rounds,
        "stress_mode": cfg.stress_mode,
        "slip_system": cfg.slip_system,
        "temps": {str(T): {"lo": lo, "hi": hi, "history": [], "active": True}
                  for T in cfg.temps()},
    }
    save_state(cfg, state)
    log(f"Initialized {regime} search: {len(cfg.temps())} temperatures, "
        f"bracket [{lo}, {hi}] mu, {cfg.n_rounds} rounds planned.")
    return state


def generate(cfg: SearchConfig, log: Callable[[str], None] = _noop):
    """Write the cases file for the next round.  Returns a dict with the count,
    cases path, and the exact sbatch command."""
    state = load_state(cfg)
    n_rounds = state.get("n_rounds", cfg.n_rounds)
    act = active_temps(state)
    if not act:
        raise SearchError("No active temperatures left. See the report.")
    for T, d in act.items():
        if d["history"] and "moved" not in d["history"][-1]:
            raise SearchError(f"T={T}K has an unrecorded result from round "
                              f"{state['round']} -- run update first.")
    if state["round"] >= n_rounds:
        raise SearchError(f"Already completed {n_rounds} rounds. See the report, "
                          f"or extend to push out-of-range temps further.")

    state["round"] += 1
    lines = [f"# Round {state['round']} of {n_rounds} -- generated by PILOT",
             "# <T_K>  <seed>  <USFE_file>  <sigma line>"]
    for T, d in act.items():
        tau = 0.5 * (d["lo"] + d["hi"])
        seed = random.randint(cfg.seed_min, cfg.seed_max)
        tag = tau_tag(tau)
        usfe = cfg.usfe_pattern.format(T=T) if "{T}" in cfg.usfe_pattern \
            else cfg.usfe_pattern
        dump = os.path.join(cfg.out_base, f"{T}K",
                            cfg.dump_name.format(T=T, tag=tag))
        d["history"].append({"round": state["round"], "tau": tau, "seed": seed,
                             "usfe": usfe, "dump": dump})
        lines.append(f"{T}  {seed}  {usfe}  {sigma_line(cfg, tau)}")

    os.makedirs(cfg.work_dir, exist_ok=True)
    with open(cfg.cases_path(), "w") as f:
        f.write("\n".join(lines) + "\n")
    save_state(cfg, state)

    n = len(act)
    sbatch = f"sbatch --array=0-{n-1} {cfg.tag}_jobs"
    log(f"Round {state['round']}: wrote {n} cases to {cfg.cases_path()}")
    log(f"Submit with:  {sbatch}")
    return {"n": n, "cases_path": cfg.cases_path(), "sbatch": sbatch,
            "round": state["round"], "n_rounds": n_rounds}


def list_pending(cfg: SearchConfig):
    """Rows awaiting a moved/not-moved grade: dicts with T, tau, dump, usfe."""
    state = load_state(cfg)
    out = []
    for T, d in sorted(state["temps"].items(), key=lambda x: int(x[0])):
        if d["history"] and "moved" not in d["history"][-1]:
            run = d["history"][-1]
            out.append({"T": T, "tau": run["tau"], "dump": run["dump"],
                        "usfe": run.get("usfe", "")})
    return out


def record_results(cfg: SearchConfig, results: Dict[str, bool],
                   log: Callable[[str], None] = _noop):
    """Apply a {T: moved} mapping to the pending round and save.

    results keys are temperature strings; values are True (moved) / False."""
    state = load_state(cfg)
    regime = state_regime(state)
    n_rounds = state.get("n_rounds", cfg.n_rounds)
    pending = {T: d for T, d in state["temps"].items()
               if d["history"] and "moved" not in d["history"][-1]}
    if not pending:
        raise SearchError("No pending results. Generate the next round first.")
    missing = [T for T in pending if T not in results]
    if missing:
        raise SearchError("Missing results for T = "
                          + ", ".join(sorted(missing, key=int)) + " K.")

    for T, moved in results.items():
        if T not in pending:
            continue
        run = pending[T]["history"][-1]
        run["moved"] = bool(moved)
        apply_result(pending[T], run["tau"], bool(moved), regime)
        log(f"T={T}K  tau={tau_str(run['tau'])}  -> "
            f"{'MOVED' if moved else 'did not move'}")

    save_state(cfg, state)
    if state["round"] >= n_rounds:
        above = [T for T, d in active_temps(state).items() if is_unbounded_above(d)]
        below = [T for T, d in active_temps(state).items() if is_unbounded_below(d)]
        log(f"Search complete ({n_rounds} rounds).")
        if above:
            log("CRSS above the probed range for T = "
                + ", ".join(sorted(above, key=int)) + " K -- use extend.")
        if below:
            log("CRSS below the probed range for T = "
                + ", ".join(sorted(below, key=int)) + " K -- use extend.")
    else:
        log("Brackets updated. Generate the next round when ready.")
    return state


def update_cli(cfg: SearchConfig, answer_fn: Callable[[dict], bool],
               log: Callable[[str], None] = _noop):
    """CLI update: answer_fn(pending_row) -> bool for each pending temperature."""
    pend = list_pending(cfg)
    if not pend:
        raise SearchError("No pending results. Generate the next round first.")
    results = {row["T"]: bool(answer_fn(row)) for row in pend}
    return record_results(cfg, results, log)


def extend(cfg: SearchConfig, log: Callable[[str], None] = _noop):
    """Widen out-of-range temps, granting extra rounds; retire finished temps."""
    state = load_state(cfg)
    regime = state_regime(state)
    pending = [T for T, d in state["temps"].items()
               if d["history"] and "moved" not in d["history"][-1]]
    if pending:
        raise SearchError("Unrecorded results for T = "
                          + ", ".join(sorted(pending, key=int))
                          + " K -- run update first.")

    above, below, bracketed = [], [], []
    for T, d in active_temps(state).items():
        if is_unbounded_above(d):
            above.append((T, d))
        elif is_unbounded_below(d):
            below.append((T, d))
        else:
            bracketed.append((T, d))

    if not above and not below:
        log("No active temperature has CRSS outside its bracket -- "
            "nothing to extend.")
        return state

    for _, d in bracketed:
        d["active"] = False

    extended, degenerate = [], []
    step = cfg.tau_extend_step

    for T, d in above:                      # never moved: push strong end out
        if regime == "positive":
            d["lo"] = d["hi"]
            d["hi"] = round(d["hi"] + step, 12)
        else:
            d["hi"] = d["lo"]
            d["lo"] = round(d["lo"] - step, 12)
        d["active"] = True
        extended.append((T, d))

    for T, d in below:                      # always moved: pull weak end to zero
        if regime == "positive":
            new_hi = d["lo"]
            new_lo = max(0.0, round(d["lo"] - step, 12))
        else:
            new_lo = d["hi"]
            new_hi = min(0.0, round(d["hi"] + step, 12))
        if new_lo >= new_hi:
            d["active"] = False
            degenerate.append((T, d))
            continue
        d["lo"], d["hi"] = new_lo, new_hi
        d["active"] = True
        extended.append((T, d))

    if extended:
        state["n_rounds"] = state.get("n_rounds", cfg.n_rounds) + cfg.n_extend_rounds
    save_state(cfg, state)

    if extended:
        ts = ", ".join(T for T, _ in sorted(extended, key=lambda x: int(x[0])))
        log(f"Extended {len(extended)} temp(s): T = {ts} K "
            f"(+{cfg.n_extend_rounds} rounds, budget now {state['n_rounds']}).")
    if bracketed:
        rt = ", ".join(T for T, _ in sorted(bracketed, key=lambda x: int(x[0])))
        log(f"Retired {len(bracketed)} already-bracketed temp(s): T = {rt} K.")
    if degenerate:
        dg = ", ".join(T for T, _ in sorted(degenerate, key=lambda x: int(x[0])))
        log(f"Cannot extend T = {dg} K: bracket already reaches zero.")
    return state


def report_rows(cfg: SearchConfig):
    """Per-temperature summary rows for tables/plots."""
    state = load_state(cfg)
    regime = state_regime(state)
    rows = []
    for T in sorted(state["temps"], key=int):
        d = state["temps"][T]
        status = "active" if d.get("active", True) else "done"
        if is_unbounded_above(d):
            regime_flag = "above range"
            est = d["hi"] if regime == "positive" else d["lo"]
        elif is_unbounded_below(d):
            regime_flag = "below range"
            est = d["lo"] if regime == "positive" else d["hi"]
        else:
            regime_flag = "bracketed"
            est = 0.5 * (d["lo"] + d["hi"])
        rows.append({"T": int(T), "lo": d["lo"], "hi": d["hi"],
                     "mid": 0.5 * (d["lo"] + d["hi"]),
                     "half": 0.5 * (d["hi"] - d["lo"]),
                     "est": est, "regime": regime_flag, "status": status})
    return {"round": state["round"], "n_rounds": state.get("n_rounds", cfg.n_rounds),
            "chi_deg": state.get("chi_deg", cfg.chi_deg), "regime": regime,
            "rows": rows}


def report_text(cfg: SearchConfig):
    r = report_rows(cfg)
    lines = [f"Round {r['round']} / {r['n_rounds']}   chi = {r['chi_deg']:+g} deg"
             f"   regime = {r['regime']}",
             f"{'T (K)':>7} {'lo (mu)':>10} {'hi (mu)':>10} "
             f"{'CRSS est.':>14} {'+/-':>12} {'status':>8}"]
    for row in r["rows"]:
        if row["regime"] == "bracketed":
            lines.append(f"{row['T']:>7} {row['lo']:>10.6f} {row['hi']:>10.6f} "
                         f"{row['mid']:>14.6f} {row['half']:>12.6f} "
                         f"{row['status']:>8}")
        else:
            note = ("> " if row["regime"] == "above range" else "< ") \
                + f"{row['est']:.4f}"
            lines.append(f"{row['T']:>7} {row['lo']:>10.6f} {row['hi']:>10.6f} "
                         f"{note:>14} {row['regime']:>12} {row['status']:>8}")
    return "\n".join(lines)


def undo(cfg: SearchConfig, log: Callable[[str], None] = _noop):
    """Step back one command (revert the most recent generate/update/extend)."""
    up = cfg.undo_path()
    if not os.path.exists(up):
        raise SearchError("Nothing to undo -- no snapshots recorded yet.")
    with open(up) as f:
        stack = json.load(f)
    if not stack:
        raise SearchError("Nothing to undo -- the undo history is empty.")
    prev = stack.pop()
    with open(up, "w") as f:
        json.dump(stack, f)
    save_state(cfg, prev, snapshot=False)
    pend = sum(1 for d in prev["temps"].values()
               if d["history"] and "moved" not in d["history"][-1])
    log(f"Reverted to round {prev['round']} "
        f"({len(stack)} further undo step(s) available).")
    if pend:
        log(f"{pend} temperature(s) now have an ungraded round -- "
            f"run update to (re-)record them, then generate.")
    log(f"Note: {cfg.cases_file()} on disk may be stale until you regenerate.")
    return prev
