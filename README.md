# USFE Workflow

Generate **unstable stacking fault energy (USFE) lattice files** for random BCC
alloys, end to end, from a single window. You describe the alloy, the tool builds
the atomic lattices with [OTIS](#otis-the-lattice-generator), maps local
composition to USFE by interpolating your SFE data, and writes one `USFE_<i>`
file per lattice — ready for downstream PFDD runs.

It works for **2- and 3-component alloys** with no code changes: everything is
driven by what you type into the GUI and by the SFE CSV you point it at.

---

## What it does

For each lattice it runs a two-step pipeline:

```
          you describe the alloy in the GUI
                        │
                        ▼
   ┌─────────────────────────────────────────┐
   │  1. OTIS                                 │   random solid-solution
   │     otis.py  →  lattice_<i>              │   BCC / CUSTOM lattice
   └─────────────────────────────────────────┘
                        │
                        ▼
   ┌─────────────────────────────────────────┐
   │  2. Interpolation                        │   local composition → USFE
   │     110_interp_usfe.py       ({110})     │   from your SFE CSV
   │     112_interp_usfe_nonbcc.py ({112})    │
   │        →  USFE_<i>                        │
   └─────────────────────────────────────────┘
```

It repeats that for however many lattices you ask for, records the mean and
standard deviation of each USFE file, cleans up the intermediates, and writes a
`usfe_stats.csv` summary.

---

## Repository contents

| File | Role |
|------|------|
| `usfe_workflow_gui.py` | The application. Launches the GUI (or a text-prompt fallback) and drives the whole workflow. **This is what you run.** |
| `usfe_core.py` | All the real logic (builds OTIS inputs, runs OTIS + interpolation, computes stats). Imported by the GUI; no GUI dependencies of its own. |
| `otis.py` | Order-Through-Informed-Swapping lattice generator. Builds the random solid-solution atomic lattice. |
| `110_interp_usfe.py` | Interpolates USFE on the **{110}** slip plane. |
| `112_interp_usfe_nonbcc.py` | Interpolates USFE on the **{112}** slip plane (non-orthogonal grid). |
| `MoNbTi_sfe_data.csv` | Example SFE data — 3-component alloy. |
| `TaW_110_SFE.csv`, `TaW_112_SFE.csv` | Example SFE data — 2-component alloy, {110} and {112}. |

Keep all of these in the same folder, and run everything from that folder.

---

## Requirements

- **Python 3**
- **numpy**, **scipy**, **matplotlib** (used by the interpolation scripts)
- **tkinter** for the graphical interface (optional — see below)

```bash
pip install numpy scipy matplotlib
```

`tkinter` ships with most Python installs. If yours doesn't have it, install it
from your OS package manager, e.g. on Debian/Ubuntu:

```bash
sudo apt-get install python3-tk
```

On macOS the python.org build or `brew install python-tk` provides it. If
tkinter is unavailable, the tool automatically falls back to
[command-line mode](#command-line-mode) — you're never blocked.

---

## Quick start (GUI)

```bash
python3 usfe_workflow_gui.py
```

A window opens with three sections, top to bottom:

1. **Lattice (OTIS)** — lattice type, cell size, lattice parameter, elements,
   concentrations.
2. **How many / naming** — number of lattices, naming pattern, output folder.
3. **Interpolation** — pick {110} or {112}, and choose your SFE CSV.

Then click **Generate USFE files**. A live log and progress bar track each
lattice; per-file mean/std are printed as they finish, and a **Stop** button
halts cleanly between lattices.

---

## Command-line mode

Same workflow, text prompts, no tkinter required:

```bash
python3 usfe_workflow_gui.py --cli
```

The tool also switches to this automatically if it can't open a window.

---

## GUI reference

### 1. Lattice (OTIS)

| Field | Meaning | Example |
|-------|---------|---------|
| Lattice type | `BCC`, `FCC`, or `CUSTOM` — **see the note below** | `BCC` |
| Cell N1 N2 N3 | Number of primitive cells in each direction | `32 32 32` |
| Lattice parameter | Nearest-neighbour distance | `2.86` |
| Elements | Element names, space/comma separated | `Mo Nb Ti` |
| Concentrations | Fractions, same count as elements (renormalised if needed) | `0.33 0.33 0.33` |
| Random solid solution | Leave checked for a random alloy; uncheck to enter a Warren–Cowley / probability matrix | ✔ |

### 2. How many / naming

| Field | Meaning | Example |
|-------|---------|---------|
| Number of lattices | How many distinct lattices/USFE files to build | `15` |
| Index start / step | The `{i}` values walk `start, start+step, …` | `0`, `100` |
| USFE name pattern | Output filename; must contain `{i}` | `USFE_{i}` |
| Keep lattice files | Keep the intermediate OTIS lattices | off |
| Output folder | Where `USFE_<i>` and `usfe_stats.csv` are written | `.` |

With start `0` and step `100`, 15 lattices produce `USFE_0, USFE_100, … USFE_1400`.
For plain indexing use start `0`, step `1`.

### 3. Interpolation

| Field | Meaning |
|-------|---------|
| Method | `{110}` → `110_interp_usfe.py`, `{112}` → `112_interp_usfe_nonbcc.py` |
| CSV | Your SFE data file ([format below](#sfe-csv-format)) |

---

## ⚠️ Lattice type must match the slip system

The OTIS lattice type and the interpolation method go together:

| Slip system | OTIS **Lattice type** | Interpolation method |
|-------------|-----------------------|----------------------|
| **{110}**   | **`BCC`**             | `110_interp_usfe.py` |
| **{112}**   | **`CUSTOM`**          | `112_interp_usfe_nonbcc.py` |

- Doing **{110}** → always select **`BCC`** in OTIS.
- Doing **{112}** → always select **`CUSTOM`** in OTIS (this builds the
  non-orthogonal, non-BCC grid the {112} interpolation expects).

Mixing these (e.g. BCC with the {112} interpolation) will not give correct USFE
values.

---

## SFE CSV format

The CSV is the only alloy-specific input to the interpolation. Its **header names
the composition columns**, and that's what makes the workflow generic:

```
SFE,<element1>,<element2>[,<element3>]
<sfe>,<frac1>,<frac2>[,<frac3>]
...
```

- **2 element columns → binary alloy** (1-D interpolation of SFE vs composition).
- **3 element columns → ternary alloy** (2-D composition-triangle interpolation).
- A UTF-8 BOM and Windows (CRLF) line endings are both fine.
- Blank lines and lines starting with `#` are ignored.

**Element names must match** the element labels you type into the GUI
(case-sensitive). The interpolation maps each atom to a composition column by
name; a mismatch produces a clear error such as
`Lattice has elements not in the CSV header: ['Mo', 'Nb', 'Ti']`.

### Include the pure-element points

The interpolation should span the full composition range, so include the
pure-element rows:

- **Ternary** — add the three corners (fractions in header column order):

  ```
  # for a header of  sfe,Ti,Nb,Mo
  1213.605,0,0,1     # pure Mo
  775.399,0,1,0      # pure Nb
  465.379,1,0,0      # pure Ti
  ```

- **Binary** — include the two endpoints (pure element rows). The example
  `TaW_*.csv` files already do.

If a local composition falls **outside** the range your CSV covers, the tool
does **not** fail — it falls back to the nearest known value and prints a
warning telling you how many sites were affected. Adding the pure-element rows
removes the warning and gives exact interpolation near the edges.

---

## Outputs

Written to the output folder you choose:

| Output | Description |
|--------|-------------|
| `USFE_<i>` | One USFE file per lattice. First line is `N1 N2 N3 3`; each subsequent line is `i j k l usfe`. |
| `usfe_stats.csv` | Summary: `index, usfe_file, lattice_file, mean_usfe, std_usfe, n_points`. |
| `lattice_<i>` | The OTIS lattices — only if **Keep lattice files** is on. |

Intermediate files (per-lattice OTIS inputs, the interpolation's
`*_comp_values.txt`) are removed automatically.

---

## How the interpolation works

Each interpolation script:

1. Reads your lattice and the SFE CSV (element names come from the CSV header).
2. Computes each site's **local composition** by averaging over its slip-plane
   neighbours (the {110} script uses 30 neighbours; the {112} script uses 16).
3. Maps local composition → USFE:
   - **3 components:** embeds the composition simplex in 2-D and does a linear
     interpolation over your CSV points.
   - **2 components:** 1-D interpolation of SFE vs the composition fraction.
4. Any site outside the CSV's data range uses nearest-neighbour extrapolation
   (with a warning) rather than producing `NaN`.

Run standalone if you like:

```bash
python3 110_interp_usfe.py <lattice_file> N1 N2 N3 [SFE_csv]
```

---

## OTIS (the lattice generator)

`otis.py` builds a random solid-solution lattice by informed atom swapping. The
GUI/core generate its input automatically, but if you run it by hand the input
format is:

```
BCC                 # lattice type: BCC, FCC, or CUSTOM
32 32 32            # N1 N2 N3
2.86                # nearest-neighbour distance
3                   # number of components
Mo Nb Ti            # element names
0.33 0.33 0.33      # concentrations
prob                # 'prob' (pij) or 'wc' (Warren-Cowley)
0.33 0.33 0.33      # one matrix row per component
0.33 0.33 0.33
0.33 0.33 0.33
my_lattice          # output filename
```

```bash
python3 otis.py <input_file>
```

Each run is a fresh random realisation, so building N lattices gives N distinct
configurations.

---

## Troubleshooting

**Buttons need two clicks (macOS).** The window isn't focused; the first click
just focuses it. Bring the window to the front first, or add this after the
window is created so it grabs focus on launch:

```python
root.lift()
root.focus_force()
root.attributes("-topmost", True)
root.after(200, lambda: root.attributes("-topmost", False))
```

**`ModuleNotFoundError: tkinter`.** Install `python3-tk`, or use `--cli`.

**USFE files are full of `NaN`, or you see “nearest-neighbour extrapolation”
warnings.** Your CSV doesn't span the compositions being queried — add the
pure-element rows (corners for a ternary, endpoints for a binary). See
[SFE CSV format](#sfe-csv-format).

**`QhullError` / “flat simplex” during interpolation.** This comes from feeding
a binary CSV to the old ternary-only interpolation. Make sure you're using the
current `110_interp_usfe.py` / `112_interp_usfe_nonbcc.py`, which handle 2- and
3-component alloys.

**`Lattice has elements not in the CSV header`.** The element names in the GUI
don't match the CSV header names (case-sensitive). Make them identical.
