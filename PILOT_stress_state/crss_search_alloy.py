#!/usr/bin/env python3
"""
crss_search_alloy.py -- command-line driver for the PILOT CRSS temperature sweep.

This is the text-mode front end; the graphical workflow (stress_state_gui.py)
does the same thing with buttons.  Both call crss_core.py, so the search logic
is identical either way.

Edit the constants below, then run one command at a time:

    python3 crss_search_alloy.py init        # create search_state.json
    python3 crss_search_alloy.py generate    # write the cases file for this round
    sbatch --array=0-14 Ta10W_110_jobs       # (generate prints the exact command)
    ... wait for the array to finish ...
    python3 crss_search_alloy.py update      # record moved/not-moved per temperature
    python3 crss_search_alloy.py report      # show current brackets / CRSS estimates

If a temperature's CRSS lies outside the initial range:

    python3 crss_search_alloy.py extend      # widen out-of-range temps, retire finished
    python3 crss_search_alloy.py generate    # then continue the generate/update loop

Made a mistake?  `undo` reverts the most recent generate/update/extend:

    python3 crss_search_alloy.py undo

Marking results: implement did_dislocation_move() below to parse your dump files
automatically; until then, `update` prompts you for m/n per temperature after you
have inspected the dumps (e.g. in OVITO -- see fetch_dumps.sh).
"""

import sys

import crss_core as cc

# ----------------------------- configuration -----------------------------
# Everything for one search is written into WORK_DIR.  Leave it "." to work in
# the current directory (classic behaviour), or point it at a study subfolder.
WORK_DIR = "."

CFG = cc.SearchConfig(
    work_dir=WORK_DIR,

    # temperatures (K)
    t_min=0, t_max=1400, t_step=100,

    # initial bracket (units of mu); the two ends in any order, must not
    # straddle zero (0 may sit at one end)
    tau_lo=0.0, tau_hi=0.08, n_rounds=6,

    # extension when CRSS lies outside the initial range
    tau_extend_step=0.02, n_extend_rounds=4,

    # stress state: "shear" (MRSSP couple, chi about b) or "uniaxial" (along axis)
    stress_mode="shear", chi_deg=0.0, axis=(1.0, 0.0, 0.0),

    # slip system + file naming
    slip_system="{110}", tag="Ta10W_110",

    # alloy USFE linkage: per-temperature USFE file name ({T} -> temperature),
    # produced by the PILOT lattice/USFE tool and numbered to match the temps
    usfe_pattern="USFE_{T}",

    # cluster-side dump naming recorded in the cases/history
    out_base="OUTPUT/dumps",
    dump_name="Ta10W_screw_110_{T}K_{tag}.dump",

    seed_min=1, seed_max=1000,
)


# --------------------------- result detection ----------------------------

def did_dislocation_move(dumpfile):
    """Return True if the dislocation moved in this run, False if not.

    TODO: implement using your dump format (e.g. locate the dislocation from the
    xir disregistry fields and compare with the initial position).  Until then,
    `update` falls back to asking you interactively after you inspect the dumps.
    """
    raise NotImplementedError


# ------------------------------- commands --------------------------------

def _answer(row):
    """Interactive grader for one pending temperature."""
    try:
        return bool(did_dislocation_move(row["dump"]))
    except NotImplementedError:
        pass
    ans = ""
    while ans not in ("m", "n"):
        ans = input(f"T={row['T']}K  tau={cc.tau_str(row['tau'])}  "
                    f"[{row['dump']}]\n  moved? (m/n): ").strip().lower()
    return ans == "m"


def main():
    cmds = {
        "init":     lambda: cc.init(CFG, log=print),
        "generate": lambda: cc.generate(CFG, log=print),
        "update":   lambda: cc.update_cli(CFG, _answer, log=print),
        "extend":   lambda: cc.extend(CFG, log=print),
        "report":   lambda: print(cc.report_text(CFG)),
        "undo":     lambda: cc.undo(CFG, log=print),
        "redo":     lambda: cc.undo(CFG, log=print),   # alias for undo
    }
    if len(sys.argv) != 2 or sys.argv[1] not in cmds:
        sys.exit(f"Usage: python3 {sys.argv[0]} "
                 f"{{init|generate|update|extend|report|undo}}")
    try:
        cmds[sys.argv[1]]()
    except cc.SearchError as e:
        sys.exit(str(e))


if __name__ == "__main__":
    main()
