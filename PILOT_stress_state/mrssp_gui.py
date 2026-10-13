#!/usr/bin/env python3
"""
mrssp_gui.py -- interactive stress-state explorer for a BCC screw dislocation.

View down the Burgers vector b = [-1 1 1]: every slip plane of the [-111] zone
appears as a trace through the origin.  Two loadings are offered:

  * Shear couple (MRSSP):  sigma = tau0 (b (x) m + m (x) b), with the MRSSP
    orientation chi set by the slider.  Zero Escaig shear, zero normal stress and
    zero pressure on every plane of the zone -- the clean loading for isolating
    the chi-dependence of glide.

  * Uniaxial:  sigma = s0 (a (x) a) for a chosen loading axis a.  Generic
    non-glide components appear; the effective MRSSP angle is computed and drawn.

For each plane the bar chart shows either the resolved shear tau/tau0 (shear
mode) or the glide / Escaig / normal traction decomposition (uniaxial mode), so
you can see what the {110} and {112} planes would feel during a PFDD run.

Run standalone:

    python3 mrssp_gui.py

Or open it from the PILOT stress-state workflow; when launched with a callback it
shows an "Apply to search" button that hands the chosen stress state (mode, chi,
axis) back to the search configuration.

The geometry helpers have no Tk dependency and can be imported on their own.
"""

import subprocess
import sys

import numpy as np


# ------------------------------ geometry ---------------------------------

b = np.array([-1.0, 1.0, 1.0]); b /= np.linalg.norm(b)      # Burgers vector
n1 = np.array([1.0, 1.0, 0.0]); n1 /= np.linalg.norm(n1)    # (110) normal, chi=0
e2 = np.cross(b, n1); e2 /= np.linalg.norm(e2)              # completes the frame

# Six distinct planes of the [-111] zone; theta measured about b from (110).
PLANES = [
    dict(theta=0.0,    label=r"$(110)$",             ascii="(110)",   family="110"),
    dict(theta=60.0,   label=r"$(01\bar{1})$",       ascii="(011̅)",  family="110"),
    dict(theta=-60.0,  label=r"$(101)$",             ascii="(101)",   family="110"),
    dict(theta=-30.0,  label=r"$(211)$",             ascii="(211)",   family="112"),
    dict(theta=-150.0, label=r"$(\bar{1}\bar{2}1)$", ascii="(1̅2̅1)", family="112", flip_label=True),
    dict(theta=90.0,   label=r"$(\bar{1}1\bar{2})$", ascii="(1̅12̅)", family="112"),
]

COL = {"110": "#1D9E75", "112": "#2246bd", "mrssp": "#EF9F27"}
COMP_COL = {"escaig": "#C4553D", "normal": "#7B3FA0"}
R = 1.0


def plane_normal(theta_deg):
    c, s = np.cos(np.radians(theta_deg)), np.sin(np.radians(theta_deg))
    return c * n1 + s * e2


def mrssp_normal(chi_deg):
    return plane_normal(chi_deg)


def stress_tensor_shear(tau0, chi_deg):
    m = mrssp_normal(chi_deg)
    return tau0 * (np.outer(b, m) + np.outer(m, b))


def stress_tensor_uniaxial(mag, axis):
    a = np.asarray(axis, float)
    a = a / np.linalg.norm(a)
    return mag * np.outer(a, a)


def resolved_shear(chi_deg, theta_deg):
    return np.cos(np.radians(theta_deg - chi_deg))


def decompose(S, theta_deg):
    """(glide, Escaig, normal) traction components on the plane at theta."""
    n = plane_normal(theta_deg)
    g = np.cross(n, b)
    return b @ S @ n, g @ S @ n, n @ S @ n


def mrssp_angle(S):
    """Effective MRSSP angle chi (deg) and max resolved shear for any tensor."""
    A, B = b @ S @ n1, b @ S @ e2
    if abs(A) < 1e-12 and abs(B) < 1e-12:
        return 0.0, 0.0
    return float(np.degrees(np.arctan2(B, A))), float(np.hypot(A, B))


def sigma_line(mode, mag, chi_deg, axis):
    """PFDD 'sigma initial ... delta ... #tau ...' line for the current state."""
    if mode == "uniaxial":
        S = stress_tensor_uniaxial(mag, axis)
        tail = "#tau %g uniaxial axis %g %g %g" % (mag, axis[0], axis[1], axis[2])
    else:
        S = stress_tensor_shear(mag, chi_deg)
        tail = "#tau %g chi %+d deg" % (mag, int(round(chi_deg)))
    comps = "%f %f %f %f %f %f" % (S[0, 0], S[1, 1], S[2, 2],
                                   S[0, 1], S[0, 2], S[1, 2])
    return ("sigma     initial " + comps
            + " delta 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 "
            + tail)


def _to_clipboard(text, tk_widget=None):
    """Best-effort clipboard copy; returns True on success."""
    try:
        import pyperclip
        pyperclip.copy(text)
        return True
    except Exception:
        pass
    if sys.platform == "darwin":
        try:
            subprocess.run(["pbcopy"], input=text.encode("utf-8"), check=True)
            return True
        except Exception:
            pass
    if tk_widget is not None:
        try:
            tk_widget.clipboard_clear()
            tk_widget.clipboard_append(text)
            return True
        except Exception:
            pass
    return False


# ================================ Tk GUI =================================

def run_mrssp(master=None, on_apply=None, init=None):
    """Open the explorer.  master given -> Toplevel child; else its own root.
    on_apply(dict) is called by the "Apply to search" button with keys
    mode, chi_deg, magnitude, axis.  init seeds the initial controls."""
    import tkinter as tk
    from tkinter import ttk, messagebox
    import matplotlib
    matplotlib.use("TkAgg")
    from matplotlib.figure import Figure
    from matplotlib.patches import Circle, Patch
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

    init = init or {}
    win = tk.Toplevel(master) if master is not None else tk.Tk()
    win.title("PILOT  -  stress state / MRSSP explorer")
    win.geometry("1080x720")

    v = {
        "mode": tk.StringVar(value=init.get("mode", "shear")),
        "chi": tk.DoubleVar(value=float(init.get("chi_deg", 0.0))),
        "mag": tk.StringVar(value=str(init.get("magnitude", 0.055))),
        "axis": tk.StringVar(value=", ".join(
            str(x) for x in init.get("axis", (10.0, -2.0, 5.0)))),
    }

    pad = {"padx": 5, "pady": 3}
    left = ttk.Frame(win); left.pack(side="left", fill="y", **pad)
    right = ttk.Frame(win); right.pack(side="right", fill="both", expand=True, **pad)

    # ---- mode ----
    mf = ttk.LabelFrame(left, text="Loading mode")
    mf.pack(fill="x", **pad)
    ttk.Radiobutton(mf, text="Shear couple (MRSSP)", value="shear",
                    variable=v["mode"]).pack(anchor="w", **pad)
    ttk.Radiobutton(mf, text="Uniaxial (along axis)", value="uniaxial",
                    variable=v["mode"]).pack(anchor="w", **pad)

    # ---- shear controls ----
    sf = ttk.LabelFrame(left, text="Shear: MRSSP angle chi (deg)")
    sf.pack(fill="x", **pad)
    chi_scale = ttk.Scale(sf, from_=-180, to=180, orient="horizontal",
                          variable=v["chi"])
    chi_scale.pack(fill="x", **pad)
    chi_lbl = ttk.Label(sf, text="chi = 0 deg")
    chi_lbl.pack(anchor="w", **pad)
    presets = ttk.Frame(sf); presets.pack(fill="x", **pad)
    for p in PLANES:
        tk.Button(presets, text=f"{int(p['theta'])}\n{p['ascii']}",
                  width=8, height=2, justify="center",
                  command=lambda t=p["theta"]: v["chi"].set(t)).pack(
            side="left", padx=1)
    ttk.Label(sf, text="presets: chi that puts the MRSSP on each plane",
              foreground="#555").pack(anchor="w", **pad)

    # ---- uniaxial controls ----
    uf = ttk.LabelFrame(left, text="Uniaxial: loading axis")
    uf.pack(fill="x", **pad)
    ttk.Entry(uf, textvariable=v["axis"], width=18).pack(side="left", **pad)
    ttk.Label(uf, text="e.g. 1, 0, 0").pack(side="left", **pad)

    # ---- magnitude ----
    gf = ttk.LabelFrame(left, text="Magnitude (tau0 / s0, units of mu)")
    gf.pack(fill="x", **pad)
    ttk.Entry(gf, textvariable=v["mag"], width=12).pack(side="left", **pad)

    # ---- readout ----
    rf = ttk.LabelFrame(left, text="Stress tensor (cube axes) & sigma line")
    rf.pack(fill="both", expand=True, **pad)
    readout = tk.Text(rf, height=11, width=40, wrap="word", font="TkFixedFont")
    readout.pack(fill="both", expand=True, **pad)

    # ---- buttons ----
    bf = ttk.Frame(left); bf.pack(fill="x", **pad)

    def do_copy():
        line = state["line"]
        print(line)
        ok = _to_clipboard(line, win)
        status.set("copied to clipboard (and printed)" if ok
                   else "clipboard unavailable -- printed to terminal")
    ttk.Button(bf, text="Copy sigma line", command=do_copy).pack(
        side="left", **pad)

    def do_apply():
        try:
            mag = float(v["mag"].get())
        except ValueError:
            mag = 0.0
        payload = {"mode": v["mode"].get(), "chi_deg": float(v["chi"].get()),
                   "magnitude": mag, "axis": _parse_axis()}
        if on_apply:
            on_apply(payload)
            status.set("applied to search")
        else:
            status.set("no search to apply to (standalone)")
    apply_btn = ttk.Button(bf, text="Apply to search", command=do_apply)
    if on_apply is not None:
        apply_btn.pack(side="left", **pad)

    status = tk.StringVar(value="")
    ttk.Label(left, textvariable=status, foreground="#555",
              wraplength=260).pack(fill="x", **pad)

    # ---- figure ----
    fig = Figure(figsize=(6.6, 6.6))
    canvas = FigureCanvasTkAgg(fig, master=right)
    canvas.get_tk_widget().pack(fill="both", expand=True)

    state = {"line": ""}

    def _parse_axis():
        try:
            vals = [float(t) for t in v["axis"].get().replace(",", " ").split()]
            if len(vals) == 3 and np.linalg.norm(vals) > 0:
                return tuple(vals)
        except ValueError:
            pass
        return (1.0, 0.0, 0.0)

    def _mag():
        try:
            return float(v["mag"].get())
        except ValueError:
            return 0.0

    def draw(*_):
        mode = v["mode"].get()
        mag = _mag()
        axis = _parse_axis()
        if mode == "uniaxial":
            S = stress_tensor_uniaxial(mag, axis)
            chi_eff, tau_max = mrssp_angle(S)
        else:
            chi = float(v["chi"].get())
            chi_lbl.configure(text=f"chi = {chi:+.0f} deg")
            S = stress_tensor_shear(mag, chi)
            chi_eff = chi

        fig.clear()
        ax_p = fig.add_axes([0.02, 0.30, 0.46, 0.66])   # polar-style view
        ax_b = fig.add_axes([0.58, 0.30, 0.40, 0.60])   # bars

        # ---- view down b ----
        ax_p.set_aspect("equal"); ax_p.axis("off")
        ax_p.set_xlim(-1.5, 1.5); ax_p.set_ylim(-1.5, 1.5)
        ax_p.add_patch(Circle((0, 0), R * 1.08, fill=False,
                              color="0.8", lw=0.8))
        ax_p.plot(0, 0, "o", ms=5, color="0.2")
        ax_p.set_title("view down b = [-111]", fontsize=10)
        for p in PLANES:
            a = np.radians(p["theta"] + 90.0)
            x, y = R * np.cos(a), R * np.sin(a)
            if mode == "uniaxial":
                g, _, _ = decompose(S, p["theta"])
                w = abs(g) / max(1e-12, max(abs(decompose(S, q["theta"])[0])
                                            for q in PLANES))
            else:
                w = abs(resolved_shear(v["chi"].get(), p["theta"]))
            ax_p.plot([-x, x], [-y, y], color=COL[p["family"]],
                      lw=0.8 + 4.5 * w, alpha=0.30 + 0.70 * w,
                      solid_capstyle="round")
            la = a + (np.pi if p.get("flip_label") else 0.0)
            ax_p.text(1.30 * np.cos(la), 1.30 * np.sin(la), p["label"],
                      ha="center", va="center", fontsize=10)
        am = np.radians(chi_eff + 90.0)
        ax_p.plot([-R * np.cos(am), R * np.cos(am)],
                  [-R * np.sin(am), R * np.sin(am)],
                  color=COL["mrssp"], lw=2, ls=(0, (5, 4)), zorder=5,
                  label="MRSSP")
        ax_p.legend(loc="lower left", frameon=False, fontsize=9)

        # ---- bar chart ----
        ypos = np.arange(len(PLANES))[::-1]
        ax_b.set_yticks(ypos)
        ax_b.set_yticklabels([p["label"] for p in PLANES], fontsize=10)
        ax_b.axvline(0, color="0.6", lw=0.8)
        ax_b.spines[["top", "right", "left"]].set_visible(False)
        if mode == "uniaxial":
            vals = {"shear": [], "escaig": [], "normal": []}
            for p in PLANES:
                g, e, nrm = decompose(S, p["theta"])
                vals["shear"].append(g); vals["escaig"].append(e)
                vals["normal"].append(nrm)
            off = {"shear": 0.22, "escaig": 0.0, "normal": -0.22}
            shear_cols = [COL[p["family"]] for p in PLANES]   # green {110}/blue {112}
            ax_b.barh(ypos + off["shear"], vals["shear"], height=0.16,
                      color=shear_cols)
            ax_b.barh(ypos + off["escaig"], vals["escaig"], height=0.16,
                      color=COMP_COL["escaig"])
            ax_b.barh(ypos + off["normal"], vals["normal"], height=0.16,
                      color=COMP_COL["normal"])
            for y in ypos[:-1]:                      # separators between planes
                ax_b.axhline(y - 0.5, color="0.9", lw=0.6)
            vmax = max(1e-9, *(abs(x) for c in vals for x in vals[c]))
            ax_b.set_xlim(-1.2 * vmax, 1.2 * vmax)
            ax_b.set_xlabel("traction on plane (units of mu)")
            handles = [Patch(color=COL["110"], label="shear {110}"),
                       Patch(color=COL["112"], label="shear {112}"),
                       Patch(color=COMP_COL["escaig"], label=r"Escaig $\perp$b"),
                       Patch(color=COMP_COL["normal"], label=r"normal $\perp$plane")]
            ax_b.legend(handles=handles, loc="lower center",
                        bbox_to_anchor=(0.5, 1.01), ncols=4, frameon=False,
                        fontsize=7)
        else:
            chi = float(v["chi"].get())
            ts = [resolved_shear(chi, p["theta"]) for p in PLANES]
            ax_b.barh(ypos, ts, color=[COL[p["family"]] for p in PLANES],
                      height=0.55)
            ax_b.set_xlim(-1.05, 1.05)
            ax_b.set_xlabel(r"resolved shear  $\tau/\tau_0$")
            ax_b.set_title("resolved (glide) shear on each plane", fontsize=9)
            for y, t in zip(ypos, ts):
                ax_b.text(t + (0.04 if t >= 0 else -0.04), y, f"{t:+.2f}",
                          va="center", ha="left" if t >= 0 else "right",
                          fontsize=9)

        canvas.draw_idle()

        # ---- text readout + sigma line ----
        state["line"] = sigma_line(mode, mag, float(v["chi"].get()), axis)
        txt = (f"mode: {mode}\n"
               f"magnitude = {mag:g}   MRSSP chi = {chi_eff:+.1f} deg\n"
               f"mean stress tr(S)/3 = {np.trace(S)/3:+.5f}\n"
               f"line stress b.S.b = {b @ S @ b:+.5f}\n\n"
               f"sxx {S[0,0]:+.5f}  syy {S[1,1]:+.5f}  szz {S[2,2]:+.5f}\n"
               f"sxy {S[0,1]:+.5f}  sxz {S[0,2]:+.5f}  syz {S[1,2]:+.5f}\n\n"
               + state["line"])
        readout.delete("1.0", "end")
        readout.insert("1.0", txt)

    for key in ("mode", "chi", "mag", "axis"):
        v[key].trace_add("write", lambda *_: draw())
    draw()

    if master is None:
        win.mainloop()
    return win


def main():
    run_mrssp(master=None)


if __name__ == "__main__":
    main()
