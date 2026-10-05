"""
Interactive slip-plane / MRSSP explorer for a BCC screw dislocation.

View down the Burgers vector b = [-1 1 1]: every plane of the zone appears
as a trace (a line through the origin). Drag the chi slider to rotate the
maximum resolved shear stress plane (MRSSP) and watch:

  * the resolved shear tau = tau0*cos(theta - chi) on each plane
    (trace thickness/opacity and the bar chart),
  * the applied stress tensor sigma = tau0*(b (x) m + m (x) b) expressed
    in cube axes -- the same six numbers the original script prints as
    its "sigma initial" line.

Requires: numpy, matplotlib.  Run:  python slip_plane_mrssp_explorer.py
"""

import subprocess
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, TextBox

# ----------------------------------------------------------------------
# Geometry (identical to the original stress-generation script)
# ----------------------------------------------------------------------
# Applied shear magnitude (sigma in the generation script); editable in
# the GUI via the "sigma" text box.
state = {"tau0": 0.055, "line": ""}

b  = np.array([-1.0, 1.0, 1.0]); b  /= np.linalg.norm(b)    # Burgers vector
n1 = np.array([ 1.0, 1.0, 0.0]); n1 /= np.linalg.norm(n1)   # (110) normal, chi = 0 reference
e2 = np.cross(b, n1);            e2 /= np.linalg.norm(e2)   # completes the frame about b

# All six distinct planes of the [-111] zone:
#   the three {110} planes of the original script, plus the three {112}
#   planes (211), (-1-21), (-11-2) of the "new 112 slip" set.
#   angle theta is measured about b from the (110) normal.
#   NOTE: the old [1,2,-1] normal is the same plane as (-1-21) with the
#   normal flipped (theta -150 vs +30), and [2,1,1]+[-1,-2,1] = [1,-1,2]
#   is (-11-2) flipped -- so these six cover every plane in the zone
#   with no duplicates.
PLANES = [
    dict(theta=0.0,    label=r"$(110)$",             family="110"),
    dict(theta=60.0,   label=r"$(01\bar{1})$",       family="110"),
    dict(theta=-60.0,  label=r"$(101)$",             family="110"),
    dict(theta=-30.0,  label=r"$(211)$",             family="112"),
    dict(theta=-150.0, label=r"$(\bar{1}\bar{2}1)$", family="112",
         flip_label=True),   # label at the upper trace end, with the others
    dict(theta=90.0,   label=r"$(\bar{1}1\bar{2})$", family="112"),
]

COL = {"110": "#1D9E75", "112": "#2246bd", "mrssp": "#EF9F27"}
R = 1.0          # trace radius in the polar view


def mrssp_normal(chi_deg):
    """Unit normal m of the MRSSP at angle chi from (110), in cube axes."""
    c = np.radians(chi_deg)
    return np.cos(c) * n1 + np.sin(c) * e2


def stress_tensor(chi_deg):
    """Pure shear couple sigma = tau0*(b (x) m + m (x) b) in cube axes."""
    m = mrssp_normal(chi_deg)
    return state["tau0"] * (np.outer(b, m) + np.outer(m, b))


def resolved_shear(chi_deg, theta_deg):
    """tau/tau0 on the plane at angle theta for MRSSP at angle chi."""
    return np.cos(np.radians(theta_deg - chi_deg))


# ----------------------------------------------------------------------
# Figure layout
# ----------------------------------------------------------------------
def default_figsize(frac=0.72, aspect=13.0 / 7.5, dpi=100.0):
    """Figure size (inches) covering `frac` of the screen width, keeping
    the layout's aspect ratio and never exceeding `frac` of the height.
    Falls back to (13, 7.5) if the screen size can't be determined."""
    w_px = h_px = None
    try:                        # cross-platform, if installed
        from screeninfo import get_monitors
        mon = get_monitors()[0]
        w_px, h_px = mon.width, mon.height
    except Exception:
        pass
    if w_px is None and sys.platform == "darwin":
        try:                    # macOS native via CoreGraphics (no GUI toolkit)
            import ctypes
            cg = ctypes.CDLL("/System/Library/Frameworks/"
                             "CoreGraphics.framework/CoreGraphics")
            cg.CGMainDisplayID.restype = ctypes.c_uint32
            did = cg.CGMainDisplayID()
            w_px = cg.CGDisplayPixelsWide(did)
            h_px = cg.CGDisplayPixelsHigh(did)
        except Exception:
            pass
    if not w_px or not h_px:
        return (13.0, 7.5)
    w = frac * w_px / dpi
    h = w / aspect
    if h > frac * h_px / dpi:       # height-limited screen: shrink to fit
        h = frac * h_px / dpi
        w = h * aspect
    return (w, h)


fig = plt.figure(figsize=default_figsize())
fig.canvas.manager.set_window_title("BCC screw dislocation: slip planes and MRSSP")

ax_polar = fig.add_axes([0.04, 0.22, 0.46, 0.72])   # view down b
ax_bars  = fig.add_axes([0.58, 0.40, 0.38, 0.52])   # resolved shear bars
ax_chi   = fig.add_axes([0.10, 0.08, 0.55, 0.035])  # slider

# ----- polar view (view down the dislocation line) --------------------
ax_polar.set_aspect("equal")
ax_polar.set_xlim(-1.45, 1.45)
ax_polar.set_ylim(-1.45, 1.45)
ax_polar.axis("off")
ax_polar.add_patch(plt.Circle((0, 0), R * 1.08, fill=False,
                              color="0.8", lw=0.8))
ax_polar.plot(0, 0, "o", ms=5, color="0.2")
ax_polar.text(0.06, 0.06, r"$\mathbf{b}\ \odot\ [\bar{1}11]$",
              fontsize=10, color="0.25")
ax_polar.set_title("View down the dislocation line", fontsize=11)

trace_lines = []
for p in PLANES:
    a = np.radians(p["theta"] + 90.0)            # trace is perpendicular to n
    x, y = R * np.cos(a), R * np.sin(a)
    (ln,) = ax_polar.plot([-x, x], [-y, y], color=COL[p["family"]],
                          lw=4, solid_capstyle="round")
    trace_lines.append(ln)
    la = a + (np.pi if p.get("flip_label") else 0.0)
    ax_polar.text(1.28 * np.cos(la), 1.28 * np.sin(la), p["label"],
                  ha="center", va="center", fontsize=11)

(mrssp_line,) = ax_polar.plot([0, 0], [-R, R], color=COL["mrssp"],
                              lw=2, ls=(0, (5, 4)), zorder=5,
                              label="MRSSP")
ax_polar.legend(loc="lower left", frameon=False, fontsize=9)

# ----- bar chart -------------------------------------------------------
ypos = np.arange(len(PLANES))[::-1]
bars = ax_bars.barh(ypos, [1] * len(PLANES),
                    color=[COL[p["family"]] for p in PLANES], height=0.55)
ax_bars.set_yticks(ypos)
ax_bars.set_yticklabels([p["label"] for p in PLANES], fontsize=11)
ax_bars.set_xlim(-1.05, 1.05)
ax_bars.axvline(0, color="0.6", lw=0.8)
ax_bars.set_xlabel(r"resolved shear  $\tau/\tau_0$")
ax_bars.spines[["top", "right", "left"]].set_visible(False)
bar_texts = [ax_bars.text(0, y, "", va="center", fontsize=10) for y in ypos]

# ----- tensor readout --------------------------------------------------
tensor_text = fig.text(0.58, 0.24, "", fontsize=10, family="monospace",
                       va="top")
fig.text(0.58, 0.30,
         r"$\sigma = \tau_0(\mathbf{b}\otimes\mathbf{m}"
         r" + \mathbf{m}\otimes\mathbf{b})$ in cube axes",
         fontsize=10, va="top")

# ----- slider and preset buttons ---------------------------------------
chi_slider = Slider(ax_chi, r"MRSSP angle $\chi$ (deg)", -180, 180,
                    valinit=0, valstep=1,
                    color=COL["mrssp"],          # filled part of the track
                    track_color="0.9",           # unfilled part of the track
                    initcolor=COL["110"],            # the tick marking valinit (or a color)
                    handle_style={"facecolor": "white",
                                  "edgecolor": COL["mrssp"],
                                  "size": 12})   # the circular grab handle
presets = [(f"{int(p['theta'])} {p['label']}", p["theta"]) for p in PLANES]

buttons = []
for i, (lab, val) in enumerate(presets):
    ax_b = fig.add_axes([0.68 + 0.052 * i, 0.065, 0.048, 0.05])
    btn = Button(ax_b, lab)
    btn.label.set_fontsize(10)
    btn.on_clicked(lambda _e, v=val: chi_slider.set_val(v))
    buttons.append(btn)

# ----- sigma input and copy button -------------------------------------
ax_sig = fig.add_axes([0.16, 0.02, 0.09, 0.045])
sigma_box = TextBox(ax_sig, r"$\sigma$ (= $\tau_0$)  ",
                    initial=str(state["tau0"]))

ax_copy = fig.add_axes([0.28, 0.02, 0.16, 0.045])
copy_btn = Button(ax_copy, "copy sigma line for PFDD input")
copy_btn.label.set_fontsize(9)
status_text = fig.text(0.46, 0.043, "", fontsize=9, color="0.35",
                       va="center")


def to_clipboard(text):
    """Best-effort clipboard copy; returns True on success.

    Never creates a new Tk root: instantiating Tk inside a non-Tk GUI
    event loop (e.g. matplotlib's native 'macosx' backend) corrupts the
    process and segfaults on the next mouse event.
    """
    try:                                   # pyperclip, if installed
        import pyperclip
        pyperclip.copy(text)
        return True
    except Exception:
        pass
    if sys.platform == "darwin":           # macOS native, no GUI toolkit
        try:
            subprocess.run(["pbcopy"], input=text.encode("utf-8"),
                           check=True)
            return True
        except Exception:
            pass
    try:                                   # only if the canvas IS Tk-based
        w = fig.canvas.get_tk_widget()
        w.clipboard_clear()
        w.clipboard_append(text)
        return True
    except Exception:
        return False


def on_sigma(text):
    try:
        state["tau0"] = float(text)
    except ValueError:
        sigma_box.set_val(str(state["tau0"]))   # revert bad input
        return
    update()


def on_copy(_event):
    line = state["line"]
    print(line)                                # always echo to the terminal
    if to_clipboard(line):
        status_text.set_text("copied to clipboard (also printed to terminal)")
    else:
        status_text.set_text("clipboard unavailable -- printed to terminal")
    fig.canvas.draw_idle()


sigma_box.on_submit(on_sigma)
copy_btn.on_clicked(on_copy)


# ----------------------------------------------------------------------
# Update loop
# ----------------------------------------------------------------------
def update(_val=None):
    chi = chi_slider.val

    # MRSSP trace
    a = np.radians(chi + 90.0)
    mrssp_line.set_data([-R * np.cos(a), R * np.cos(a)],
                        [-R * np.sin(a), R * np.sin(a)])

    # plane traces + bars
    for p, ln, bar, txt in zip(PLANES, trace_lines, bars, bar_texts):
        t = resolved_shear(chi, p["theta"])
        ln.set_linewidth(0.8 + 4.5 * abs(t))
        ln.set_alpha(0.30 + 0.70 * abs(t))
        bar.set_width(t)
        txt.set_x(t + (0.04 if t >= 0 else -0.04))
        txt.set_ha("left" if t >= 0 else "right")
        txt.set_text(f"{t:+.2f}")

    # cube-frame tensor (the "sigma initial" line of the original script)
    s = stress_tensor(chi)
    init_part = ("sigma     initial %f %f %f %f %f %f"
                 % (s[0, 0], s[1, 1], s[2, 2], s[0, 1], s[0, 2], s[1, 2]))
    delta_part = ("delta %f %f %f %f %f %f #tau %g chi %+d deg"
              % (0, 0, 0, 0, 0, 0, state["tau0"], chi))
    state["line"] = init_part + " " + delta_part
    tensor_text.set_text(
        f"tau0 = {state['tau0']:g}   chi = {chi:+.0f} deg\n"
        f"sxx {s[0,0]:+.4f}   syy {s[1,1]:+.4f}   szz {s[2,2]:+.4f}\n"
        f"sxy {s[0,1]:+.4f}   sxz {s[0,2]:+.4f}   syz {s[1,2]:+.4f}\n\n"
        + init_part + "\n  " + delta_part
    )
    status_text.set_text("")
    fig.canvas.draw_idle()


chi_slider.on_changed(update)
update()

if __name__ == "__main__":
    plt.show()
