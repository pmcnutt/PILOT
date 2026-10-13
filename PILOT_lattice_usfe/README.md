# PILOT — Lattice & USFE generation

This is the **lattice / USFE half of PILOT** (Phase field dislocation dynamics
Input, Lattice & Optimization Toolkit). It builds **random solid-solution (RSS)
lattices** for a bcc alloy and turns each one into a **per-site unstable
stacking fault energy (USFE) field** ready for a Phase Field Dislocation
Dynamics (PFDD) run — end to end, from a single window.

You describe the alloy, pick a slip system, and the tool:

1. generates one or more RSS lattices,
2. computes each site's local composition over its slip-plane neighbours,
3. interpolates your stacking fault energy (SFE) data to that composition, and
4. writes one `USFE_<i>` file per lattice, plus the lattice and composition
   files if you want to keep them.

A built-in **viewer** then lets you inspect any generated cell — 2-D planes and
rotatable 3-D views of atom types, composition, and USFE — without leaving the
workflow.

It works for **2- and 3-component alloys** with no code changes: everything is
driven by what you enter and by the SFE CSV you point it at.

---

## Slip system fixes everything

PILOT models **random solid solutions only** — there is no atom swapping and no
short-range-order option. Choosing the slip system automatically selects the
lattice type *and* the interpolation, so they can never be mismatched:

| Slip system | Lattice type       | Interpolation script          | Neighbours | Slip planes (l = 0, 1, 2) |
|-------------|--------------------|-------------------------------|------------|---------------------------|
| **{110}**   | `BCC`              | `110_interp_usfe.py`          | 30         | (110), (01-1), (101)      |
| **{112}**   | `CUSTOM_NON_ORTHO` | `112_interp_usfe_nonbcc.py`   | 16         | (211), (-1-21), (-11-2)   |

The `{112}` case needs a **non-orthogonal grid** (`CUSTOM_NON_ORTHO`), which is
built automatically for that slip system. You never pair the lattice and the
interpolation by hand.

---

## What it does

For each lattice it runs a two-step pipeline:

```
             you describe the alloy + slip system
                          │
                          ▼
   ┌──────────────────────────────────────────────┐
   │  1. lattice_generator.py                      │  random solid-solution
   │       → lattice_<i>                            │  BCC or CUSTOM_NON_ORTHO
   └──────────────────────────────────────────────┘
                          │
                          ▼
   ┌──────────────────────────────────────────────┐
   │  2. interpolation  (auto-selected)            │  local composition → USFE
   │       110_interp_usfe.py        ({110})       │  from your SFE CSV
   │       112_interp_usfe_nonbcc.py ({112})       │  → USFE_<i>, comp_<i>
   └──────────────────────────────────────────────┘
```

It repeats that for however many lattices you ask for, records the mean and
standard deviation of each USFE file, and writes a `usfe_stats.csv` summary.

Generating **many** lattices is the point: each is an independent random
realisation, so the set samples the statistical variation in the local
compositions a dislocation encounters in the bulk.

---

## Repository contents

| File | Role |
|------|------|
| `usfe_workflow_gui.py` | The application. Launches the GUI (or a text-prompt fallback) and drives the whole workflow. **This is what you run.** |
| `usfe_core.py` | All the logic (builds generator inputs, runs the generator + interpolation, names/keeps files, computes stats). Imported by the GUI; no GUI dependencies of its own. |
| `viewer.py` | Interactive viewer (Tk + matplotlib) for lattice / composition / USFE files. Opened from the GUI, or run standalone. |
| `lattice_generator.py` | Random solid-solution lattice generator (replaces the old OTIS generator). BCC and CUSTOM_NON_ORTHO only; no swapping. |
| `110_interp_usfe.py` | Interpolates USFE on the **{110}** slip plane (BCC grid, 30 neighbours). |
| `112_interp_usfe_nonbcc.py` | Interpolates USFE on the **{112}** slip plane (non-orthogonal grid, 16 neighbours). |
| `TaW_110_SFE.csv`, `TaW_112_SFE.csv` | Example SFE data — 2-component Ta–W, {110} and {112}. |
| `AB_110_SFE.csv` | Example SFE data — generic 2-component A–B. |
| `MoNbTi_sfe_data.csv` | Example SFE data — 3-component alloy. |

Keep these in the same folder, and run everything from that folder.

---

## Requirements

- **Python 3**
- **numpy**, **scipy**, **matplotlib**
- **tkinter** for the GUI and the viewer (optional — see below)

```bash
pip install numpy scipy matplotlib
```

`tkinter` ships with most Python installs. If yours doesn't have it:

```bash
sudo apt-get install python3-tk      # Debian/Ubuntu
brew install python-tk               # macOS (or use the python.org build)
```

If tkinter is unavailable, the workflow automatically falls back to
[command-line mode](#command-line-mode) — you're never blocked. (The viewer
does need tkinter.)

---

## Quick start (GUI)

```bash
python3 usfe_workflow_gui.py
```

A window opens with three sections:

1. **Slip system & alloy** — pick `{110}` or `{112}`, then set cell size,
   lattice parameter, elements, and concentrations.
2. **How many / naming** — number of lattices, index start/step, the name
   patterns for the USFE, lattice, and composition files, and whether to keep
   the lattice and composition files.
3. **SFE data (CSV)** — choose your SFE CSV.

Click **Generate USFE files**. A live log and progress bar track each lattice.
When it's done, click **View generated files...** to open the viewer.

---

## Command-line mode

Same workflow, text prompts, no tkinter required:

```bash
python3 usfe_workflow_gui.py --cli
```

The tool also switches to this automatically if it can't open a window.

---

## Naming and keeping files

Every output is indexed by the same `{i}`, so a set of lattices lines up
one-to-one across the three file types (and, if you like, with a downstream
temperature series):

| Field | Default pattern | Example (start 0, step 100) |
|-------|-----------------|-----------------------------|
| USFE name | `USFE_{i}` | `USFE_0, USFE_100, … USFE_1400` |
| Lattice name | `lattice_{i}` | `lattice_0, lattice_100, …` |
| Comp name | `comp_{i}` | `comp_0, comp_100, …` |

"Keep lattice files" and "Keep composition files" are on by default so the
viewer has everything it needs. Turn them off to write only the `USFE_<i>`
files and the summary.

With index start `0` and step `1` you get plain `USFE_0, USFE_1, …`.

---

## The viewer

Open it from the GUI ("View generated files...") or run it standalone:

```bash
python3 viewer.py [folder]
```

In the viewer:

- **Files** — point it at a lattice, composition, and USFE file. The
  *fill index* box fills all three from `lattice_{i}` / `comp_{i}` /
  `USFE_{i}` in the output folder in one click. The optional **SFE CSV**
  field lets it label the composition columns by element name.
- **Slip system & plane** — pick the system (`{110}` or `{112}`), then the
  **slip plane**, which is named explicitly:
  - {110}: `0: (110)`, `1: (01-1)`, `2: (101)`
  - {112}: `0: (211)`, `1: (-1-21)`, `2: (-11-2)`

  Selecting a plane colors the USFE views with that plane's USFE column *and*
  sets the 2-D slice to that plane (all three contain the Burgers vector
  `[-111]`, which stays the in-plane horizontal axis).
- **View** — choose one of:
  - *Composition (2-D plane)* — one atomic plane colored by local
    composition,
  - *USFE (2-D plane)* — the same plane colored by USFE,
  - *Atom types (3-D cell)* — the whole cell, atoms colored by element,
  - *USFE (3-D cell)* — the whole cell colored by USFE.

  USFE color bars are labelled `USFE_0` / `USFE_1` / `USFE_2` for the
  selected slip plane. The 3-D views use a fixed viewing frame per slip
  system, so switching the USFE index only recolors the cell — it does not
  reorient the camera.
- **Color by composition** — how the composition views are colored:
  - *binary*: by one element fraction (1-D color bar),
  - *ternary*: a **ternary RGB triangle** shown in its own panel beside the
    plot (so it never overlaps the atoms) — each element sits at a corner
    (red / green / blue), so every atom's color is its full composition,
  - *4+ components*: by one chosen element fraction (1-D color bar).
- **Color scale** — the 1-D color bars auto-scale to the data by default;
  the *Composition* and *USFE* min/max boxes fix the range instead (leave a
  box blank for auto). Fixing it keeps the bar consistent as you step
  through planes or USFE indices. The ternary RGB triangle ignores these.
- **2-D plane** — slide through the atomic planes *along the selected slip
  plane's normal* (index 0 = one edge of the cell, up to the far edge).
- **3-D view** — **axes are on**, the default camera matches the notebooks
  (elevation 22°, azimuth 180°), and you can adjust it live with the
  **elevation / azimuth sliders** or by **dragging the plot**. Lower the
  *subsample* for snappier rotation on large cells.

The viewer reproduces the two Jupyter notebooks (`{110}` and `{112}`) in one
window; you no longer need to open a notebook to look at a result.

---

## SFE CSV format

The CSV is the only alloy-specific input to the interpolation. Its **header
names the composition columns**, which is what makes the workflow generic:

```
SFE,<element1>,<element2>[,<element3>]
<sfe>,<frac1>,<frac2>[,<frac3>]
...
```

- **2 element columns → binary alloy** (1-D interpolation of SFE vs composition).
- **3 element columns → ternary alloy** (2-D composition-triangle interpolation).
- **4+ element columns → nearest-neighbour lookup**: each site is assigned the
  USFE of the closest composition in the CSV (a lookup table).
- A UTF-8 BOM and Windows (CRLF) line endings are both fine.
- Blank lines and lines starting with `#` are ignored.
- **Element names must match** the labels you enter (case-sensitive). A mismatch
  gives a clear error, e.g. `Lattice has elements not in the CSV header:
  ['Mo', 'Nb', 'Ti']`.

**Include the pure-element points** so the interpolation spans the full range:

- **Binary** — include the two endpoints (pure-element rows). The example
  `TaW_*.csv` files already do.
- **Ternary** — add the three corners (fractions in header column order):

  ```
  # for a header of  sfe,Ti,Nb,Mo
  1213.605,0,0,1     # pure Mo
  775.399,0,1,0      # pure Nb
  465.379,1,0,0      # pure Ti
  ```

If a local composition falls **outside** the CSV range, the tool does **not**
fail — it falls back to the nearest known value and prints a warning saying how
many sites were affected. Adding the pure-element rows removes the warning.

---

## Outputs

Written to the output folder you choose:

| Output | Description |
|--------|-------------|
| `USFE_<i>` | One USFE file per lattice. First line is `N1 N2 N3 3`; each subsequent line is `i j k l usfe`. Consumed directly by PFDD. |
| `lattice_<i>` | The RSS lattice (extended-XYZ style). Written if **Keep lattice files** is on. |
| `comp_<i>` | Per-site local composition, one `(idx, 0): [c0 c1 …]` line per site. Written if **Keep composition files** is on. |
| `usfe_stats.csv` | Summary: `index, usfe_file, lattice_file, comp_file, mean_usfe, std_usfe, n_points`. |

---

## How the interpolation works

Each interpolation script:

1. Reads the lattice and the SFE CSV (element names come from the CSV header).
2. Computes each site's **local composition** by averaging over its slip-plane
   neighbours (**{110}** uses 30 neighbours; **{112}** uses 16).
3. Maps local composition → USFE:
   - **3 components:** embeds the composition simplex in 2-D and does a linear
     interpolation over your CSV points.
   - **2 components:** 1-D interpolation of SFE vs the composition fraction.
   - **4+ components:** nearest-neighbour lookup — each site takes the USFE of
     the closest CSV composition. Add more CSV rows for finer resolution.
4. For 2- and 3-component alloys, any site outside the CSV's range uses
   nearest-neighbour extrapolation (with a warning) rather than producing `NaN`.

Run one standalone if you like (6th argument is the composition-output path):

```bash
python3 110_interp_usfe.py <lattice_file> N1 N2 N3 <SFE_csv> [comp_out]
```

---

## The lattice generator

`lattice_generator.py` builds a **random solid solution**: atoms are placed on
the grid in the exact target proportions and then shuffled. There is no atom
swapping and no short-range order. The GUI/core write its input automatically,
but if you run it by hand the input format is:

```
BCC                 # lattice type: BCC or CUSTOM_NON_ORTHO
32 32 32            # N1 N2 N3
2.86                # nearest-neighbour distance
2                   # number of components
Ta W                # element names
0.5 0.5             # concentrations
my_lattice          # output filename
```

```bash
python3 lattice_generator.py <input_file>
```

Each run is a fresh random realisation, so building N lattices gives N distinct
configurations.

---

## How it fits into PILOT

This half produces the composition-resolved **USFE lattice files** that PFDD
consumes. The other half of PILOT (the CRSS temperature sweep) takes those
files, applies a stress state, and searches for the critical resolved shear
stress as a function of temperature. Together they form the end-to-end PILOT
workflow, with PFDD at the center.

---

## Troubleshooting

**`ModuleNotFoundError: tkinter`.** Install `python3-tk`, or use `--cli` for the
workflow. (The viewer needs tkinter.)

**USFE files full of `NaN`, or "nearest-neighbour extrapolation" warnings.**
The CSV doesn't span the compositions being queried — add the pure-element rows
(endpoints for a binary, corners for a ternary).

**`Lattice has elements not in the CSV header`.** The element names you entered
don't match the CSV header (case-sensitive). Make them identical.

**Viewer says "file sizes disagree".** The lattice, composition, and USFE files
are from different runs or different cell sizes. Use a matched `{i}` triple (the
*fill index* box does this for you).

**Buttons need two clicks (macOS).** The window isn't focused yet; the first
click just focuses it. Bring the window to the front first.
