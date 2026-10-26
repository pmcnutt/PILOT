# PILOT — Phase field dislocation dynamics Input, Lattice & Optimization Toolkit

PILOT is an end-to-end workflow for modeling plasticity in body-centered cubic
(bcc) metals and alloys with **Phase Field Dislocation Dynamics (PFDD)**. It
automates the two tasks that otherwise make PFDD alloy studies slow and manual:
building the composition-resolved **unstable stacking fault energy (USFE)**
inputs PFDD needs, and locating the **critical resolved shear stress (CRSS)** as
a function of temperature.

The toolkit comes in two halves that meet at the PFDD solver:

```
   alloy + slip system                         stress state + temperatures
          │                                              │
          ▼                                              ▼
  ┌───────────────────────┐   USFE_<i> files   ┌───────────────────────┐
  │  PILOT_lattice_usfe   │ ─────────────────► │ PILOT_stress_state    │
  │  RSS lattices +       │   (one per temp)   │ bisection CRSS sweep  │
  │  per-site USFE fields │                    │ driving PFDD on SLURM │
  └───────────────────────┘                    └───────────────────────┘
          │                                              │
          ▼                                              ▼
     3-D viewer                                   CRSS vs temperature
```

* **`PILOT_lattice_usfe/`** builds random solid-solution (RSS) lattices for a
  bcc alloy and maps each site's local composition to a USFE value, producing
  the per-site `USFE_<i>` files PFDD reads. It includes an interactive 3-D
  viewer for the generated lattices, compositions, and USFE fields.

* **`PILOT_stress_state/`** takes those USFE files, applies a stress state, and
  runs a **bisection search** across temperature to locate the CRSS, driving
  PFDD through a SLURM job array. It includes a stress-state explorer, a
  GUI-based result grader, and CRSS-vs-temperature plotting.

Each half runs on its own and has its **own detailed README** in its
subdirectory. This file is the overview and the map between them.

---

## Repository layout

```
PILOT/
├── README.md                 ← you are here (overview + how the halves connect)
├── PILOT_lattice_usfe/       ← lattice + USFE generation (and 3-D viewer)
│   └── README.md             ← full details for this half
└── PILOT_stress_state/       ← CRSS temperature sweep (and stress explorer)
    └── README.md             ← full details for this half
```

---

## Requirements

* **Python 3** with **numpy**, **scipy**, and **matplotlib**.
* **tkinter** for the graphical interfaces.
  **For the GUIs to load and render properly, Tk 8.6 or newer is required.**
  You can check your version with:

  ```bash
  python3 -c "import tkinter; print(tkinter.TkVersion)"
  ```

  If this prints `8.5` (common with the macOS *system* Python), the GUI windows
  may open blank. Use a Python that ships Tk 8.6+ — for example the installer
  from python.org, or `brew install python-tk` on macOS, or
  `sudo apt-get install python3-tk` on Debian/Ubuntu. The command-line
  interfaces do not need tkinter.
* The CRSS half additionally requires a **SLURM** scheduler and a user-built
  **PFDD executable** for the actual simulations; the rest of the toolkit runs
  on a laptop.

Install the Python dependencies with:

```bash
pip install numpy scipy matplotlib
```

---

## Quick start

### 1. Build the USFE inputs (`PILOT_lattice_usfe/`)

```bash
cd PILOT_lattice_usfe
python3 usfe_workflow_gui.py          # GUI  (or: --cli for command-line mode)
```

Pick a slip system (`{110}` or `{112}`), set the cell size and lattice
parameter, choose the number of alloying elements and their concentrations, and
point the tool at your stacking fault energy CSV. It writes one `USFE_<i>` file
per lattice (plus the lattice and composition files), and the built-in viewer
("View generated files...") lets you inspect any result in 2-D and 3-D.

Number the USFE files to line up with your temperature series (for example
`USFE_0, USFE_100, …`) so the stress half can pick them up per temperature.

See `PILOT_lattice_usfe/README.md` for the full description.

### 2. Run the CRSS temperature sweep (`PILOT_stress_state/`)

```bash
cd PILOT_stress_state
python3 stress_state_gui.py           # GUI  (or: python3 crss_search_alloy.py for CLI)
```

Name a study, set the temperature range and the initial stress bracket, choose
the stress state (visualize it with the MRSSP explorer), and link the
per-temperature USFE files from step 1. Then initialize the search and iterate:
generate a round, submit the SLURM array, fetch the dumps, grade which
temperatures glided, and repeat. The tool plots CRSS versus temperature when
you are done.

See `PILOT_stress_state/README.md` for the full description.

---

## How the two halves connect

The link between the halves is the **USFE file per temperature**. The lattice
half writes `USFE_<i>` files; the stress half's cases file carries, on each line,
the USFE file name for that temperature, so the SLURM job points PFDD at the
right composition-resolved USFE field. Numbering the `USFE_<i>` files to match
your temperature list (via the index start/step in the lattice GUI) is what makes
them line up one-to-one.

Both halves share the same look and feel: a graphical interface backed by a
Tk-free engine, with an equivalent command-line path, and the PILOT banner at
launch.

---

## Command-line use

Everything the GUIs do is also available without a display:

* Lattice/USFE: `python3 usfe_workflow_gui.py --cli`
* CRSS sweep: `python3 crss_search_alloy.py {init|generate|update|extend|report|undo}`

This is convenient on headless clusters or when Tk is unavailable.

---

## Citing and further documentation

PILOT is described in an accompanying manuscript (in preparation). Each
subdirectory's README documents that half's inputs, outputs, file formats, and
options in full, including the SFE CSV format, the slip-system/lattice pairing,
the viewer controls, the bisection search, and the SLURM setup.
