#!/usr/bin/env python3
r"""
plot_crss.py -- plot critical resolved shear stress (CRSS) vs temperature from
one or more search_state.json files (e.g. different slip planes) on shared axes.

CRSS at each temperature is the midpoint of the current bisection bracket
[lo, hi]; the error bar is the half-width (hi - lo)/2 -- the uncertainty from
the stress regime tested so far. crss_search.py only tightens lo/hi on
`update`, so this reflects whatever rounds have been recorded: run it now for
the current picture, re-run after the next update for smaller error bars.

Usage:
    # CURRENT LINE: python3 plot_crss.py \(211\)T.json \(-11-2\)T.json \(-1-21\)T.json --abs --no-dodge --out crss_112T.png
    # single plane
    python3 plot_crss.py state_110.json

    # three planes on one plot, with labels
    python3 plot_crss.py state_110.json state_112AT.json state_112TW.json \\
        --labels "{110}" "{112} anti-twinning" "{112} twinning"

    # label inline instead (path=Label), any order
    python3 plot_crss.py "s110.json={110}" "s112.json={112} TW"

    # save, and/or convert mu -> MPa
    python3 plot_crss.py *.json --out crss.pdf --mu 161000 --unit MPa

Options:
    --labels ...   one label per file, in the same order (overrides path=Label)
    --out FILE     save the figure (.png/.pdf/.eps/.svg); implies no window
    --mu VALUE     scale all stresses by VALUE (e.g. shear modulus)
    --unit NAME    y-axis stress unit label (default: "mu")
    --title TEXT   override the auto title
    --no-dodge     do not offset series horizontally
    --no-show      do not open an interactive window
"""

import argparse
import json
import os
import sys

from matplotlib.ticker import AutoMinorLocator

# distinct, print-friendly styles cycled across files
COLORS = ["#1f5fa8", "#c05621", "#2f855a", "#6b46c1", "#b7791f", "#a02c50"]
MARKERS = ["o", "s", "^", "D", "v", "P"]


def load_state(path):
    if not os.path.exists(path):
        sys.exit(f"File not found: {path}")
    with open(path) as f:
        return json.load(f)


def summarize(state):
    """Per-temperature rows sorted by T, with bracket midpoint/half-width and
    a regime flag: 'bracketed', 'upper_bounded' (moved at every tau), or
    'lower_bounded' (never moved)."""
    rows = []
    for T_str, d in state["temps"].items():
        lo, hi = float(d["lo"]), float(d["hi"])
        flags = [h["moved"] for h in d.get("history", []) if "moved" in h]
        pending = bool(d.get("history")) and "moved" not in d["history"][-1]
        if flags and all(flags):
            regime = "upper_bounded"
        elif flags and not any(flags):
            regime = "lower_bounded"
        else:
            regime = "bracketed"
        rows.append({"T": int(T_str), "lo": lo, "hi": hi,
                     "mid": 0.5 * (lo + hi), "half": 0.5 * (hi - lo),
                     "n_done": len(flags), "pending": pending, "regime": regime})
    rows.sort(key=lambda r: r["T"])
    return rows


def default_label(state, path):
    """Fall back to a readable label when the user gives none."""
    stem = os.path.splitext(os.path.basename(path))[0]
    if stem and stem != "search_state":
        return stem
    chi = state.get("chi_deg")
    return fr"$\chi={chi:+g}^\circ$" if chi is not None else stem or path


def parse_inputs(tokens):
    """Split positional args of the form 'path' or 'path=Label'."""
    files, inline = [], []
    for tok in tokens:
        if "=" in tok:
            path, lab = tok.split("=", 1)
            files.append(path)
            inline.append(lab)
        else:
            files.append(tok)
            inline.append(None)
    return files, inline


def main():
    ap = argparse.ArgumentParser(description="Plot CRSS vs temperature for one or more slip planes")
    ap.add_argument("states", nargs="+", help="search_state.json file(s), optionally path=Label")
    ap.add_argument("--labels", nargs="+", default=None, help="one label per file, in order")
    ap.add_argument("--out", default=None, help="save figure to this file")
    ap.add_argument("--mu", type=float, default=1.0, help="scale factor for stresses")
    ap.add_argument("--unit", default="mu", help='y-axis stress unit (default: "mu")')
    ap.add_argument("--title", default=None, help="override the auto title")
    ap.add_argument("--abs", action="store_true",
                    help="plot |CRSS| (magnitude), e.g. for negative-shear conventions")
    ap.add_argument("--no-dodge", action="store_true", help="do not offset series in T")
    ap.add_argument("--no-show", action="store_true", help="do not open a window")
    args = ap.parse_args()

    import matplotlib
    if args.no_show or args.out:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    files, inline = parse_inputs(args.states)
    if args.labels and len(args.labels) != len(files):
        sys.exit(f"--labels has {len(args.labels)} entries but {len(files)} files were given.")

    series = []
    for i, path in enumerate(files):
        st = load_state(path)
        if args.labels:
            label = args.labels[i]
        elif inline[i] is not None:
            label = inline[i]
        else:
            label = default_label(st, path)
        series.append({"path": path, "label": label, "state": st,
                       "rows": summarize(st)})

    s, unit = args.mu, args.unit
    mag = (lambda x: abs(x)) if args.abs else (lambda x: x)   # magnitude toggle
    ylabel = (r"$\tau_{\mathrm{CRSS}}$" + f"  ({unit})") if unit else r"$\tau_{\mathrm{CRSS}}$"

    # horizontal dodge so error bars at the same T don't overlap
    all_T = sorted({r["T"] for ser in series for r in ser["rows"]})
    if len(all_T) > 1:
        dt = min(b - a for a, b in zip(all_T, all_T[1:]))
    else:
        dt = 100.0
    n = len(series)
    step = 0.0 if (args.no_dodge or n == 1) else dt * 0.11
    offsets = [(i - (n - 1) / 2.0) * step for i in range(n)]

    fig, ax = plt.subplots(figsize=(15.2, 10.4))
    any_onesided = False

    for i, ser in enumerate(series):
        color = COLORS[i % len(COLORS)]
        marker = MARKERS[i % len(MARKERS)]
        dx = offsets[i]
        rows = ser["rows"]
        normal = [r for r in rows if r["regime"] == "bracketed"]
        upper = [r for r in rows if r["regime"] == "upper_bounded"]
        lower = [r for r in rows if r["regime"] == "lower_bounded"]

        if normal:
            ax.errorbar([r["T"] + dx for r in normal], [mag(r["mid"]) * s for r in normal],
                        yerr=[r["half"] * s for r in normal],
                        fmt=marker, ms=6, capsize=3.5, lw=1.2, elinewidth=1.1,
                        color=color, ecolor=color, mfc=color, mec=color,
                        label=ser["label"], zorder=3)
        # one-sided points: known bound + caret toward the unknown side
        if upper:
            any_onesided = True
            ax.errorbar([r["T"] + dx for r in upper], [mag(r["hi"]) * s for r in upper],
                        yerr=[abs(r["hi"]) * s * 0.06 for r in upper], uplims=True,
                        fmt=marker, ms=6, lw=1.1, color=color, ecolor=color,
                        mfc="none", mec=color, zorder=3,
                        label=None if normal else ser["label"])
        if lower:
            any_onesided = True
            ax.errorbar([r["T"] + dx for r in lower], [mag(r["lo"]) * s for r in lower],
                        yerr=[max(abs(r["lo"]), 1e-9) * s * 0.06 for r in lower], lolims=True,
                        fmt=marker, ms=6, lw=1.1, color=color, ecolor=color,
                        mfc="none", mec=color, zorder=3,
                        label=None if normal else ser["label"])

    ax.set_xlabel("Temperature (K)", fontsize=32)
    ax.set_ylabel(r'$|\sigma^{crit}|/\mu$' if args.abs else r'$\sigma^{crit}/\mu$', fontsize=32)
    ax.grid(True, alpha=0.5, linewidth=.7)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))          
    ax.grid(True, which="minor", axis="x", alpha=0.5, linewidth=.7)
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))          
    ax.grid(True, which="minor", axis="y", alpha=0.5, linewidth=.7)
    ax.tick_params(axis='both', labelsize=28)
    ax.margins(x=0.04)
    if args.abs or all(r["lo"] >= 0 for ser in series for r in ser["rows"]):
        ax.set_ylim((0, 0.14))

    # if args.title:
    #     ax.set_title(args.title)
    # else:
    #     chis = {ser["state"].get("chi_deg") for ser in series}
    #     title = "CRSS vs temperature"
    #     if len(series) == 1 and None not in chis:
    #         title += fr"  --  $\chi = {chis.pop():+g}^\circ$"
    #     ax.set_title(title)

    # completeness footnote across all series
    max_done = max((r["n_done"] for ser in series for r in ser["rows"]), default=0)
    n_pending = sum(r["pending"] for ser in series for r in ser["rows"])
    note = f"{len(series)} plane(s); up to {max_done} recorded round(s)"
    if n_pending:
        note += f"; {n_pending} point(s) awaiting the current round"
    if any_onesided:
        note += "; open carets = one-sided bound (widen initial range)"
    fig.text(0.99, 0.01, note, ha="right", va="bottom", fontsize=8, color="#666666")

    ax.legend(frameon=False, fontsize=28, loc="best", title="slip plane" if n > 4 else None)
    fig.tight_layout()

    # console tables
    for ser in series:
        print(f"\n# {ser['label']}   ({ser['path']})")
        print(f"{'T (K)':>7} {'lo':>10} {'hi':>10} {'CRSS':>11} {'+/-':>10}  regime")
        for r in ser["rows"]:
            print(f"{r['T']:>7} {r['lo']*s:>10.5g} {r['hi']*s:>10.5g} "
                  f"{r['mid']*s:>11.5g} {r['half']*s:>10.4g}  {r['regime']}"
                  + ("  [pending]" if r["pending"] else ""))

    if args.out:
        fig.savefig(args.out, dpi=200, bbox_inches="tight")
        print(f"\nsaved figure -> {args.out}")
    if not args.no_show and not args.out:
        plt.show()


if __name__ == "__main__":
    main()
