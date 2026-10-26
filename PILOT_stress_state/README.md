# PILOT — Stress state & CRSS temperature sweep

This is the **stress-state half of PILOT** (Phase field dislocation dynamics
Input, Lattice & Optimization Toolkit). It finds the **critical resolved shear
stress (CRSS)** of a screw dislocation as a function of temperature by driving
Phase Field Dislocation Dynamics (PFDD) SLURM job arrays through a **bisection
search**, for **alloys** — each temperature is run with its own
composition-resolved USFE file produced by the PILOT lattice/USFE tool.

It pairs with the other half (`PILOT_lattice_usfe`): that tool builds the
`USFE_<i>` files; this tool applies a stress state to them and searches for
CRSS(T). Together, with PFDD at the center, they form the end-to-end workflow.

You can run everything from a **GUI** or from the **command line** — both call
the same engine (`crss_core.py`), so they behave identically.

---

## What it does

For each temperature it maintains a stress bracket `[lo, hi]` known to contain
the CRSS. Each round tests the midpoint: if the dislocation glides, the CRSS is
at or below it (bring the strong end in); if not, it is beyond it (push the weak
end out). Bracketing halves the interval each round, so the resolution improves
as `2^-n` in the number of rounds `n` — a handful of rounds pins CRSS to a tight
range instead of an exhaustive stress scan.

```
   set temps, bracket, rounds, stress state, USFE link
                     │
     init ──► generate ──► (submit SLURM array) ──► fetch dumps
                     ▲                                   │
                     │                                   ▼
                  extend ◄──────── update ◄──── inspect in OVITO
                     │           (moved / not-moved per T)
                     ▼
                  report ──► plot CRSS vs T
```

`extend` widens temperatures whose CRSS fell outside the initial bracket; `undo`
reverts the most recent generate/update/extend if you make a mistake.

---

## Files

| File | Role |
|------|------|
| `stress_state_gui.py` | The graphical workflow. **Run this for the GUI.** |
| `crss_search_alloy.py` | The command-line driver (init/generate/update/extend/report/undo). |
| `crss_core.py` | The search engine shared by the GUI and CLI. No GUI dependencies. |
| `mrssp_gui.py` | Stress-state explorer (shear couple **and** uniaxial). Launched from the GUI to choose/visualize the stress state, or run standalone. |
| `plot_crss.py` | Plot CRSS vs temperature from one or more `search_state.json` files. |
| `fetch_dumps.sh` | Copy the PFDD dump files for a study from a cluster to a local folder (for OVITO). |
| `templates/in.Ta10W_110_template`, `templates/in.Ta10W_112_template` | PFDD input templates (with `temperature`, `seed`, `sigma`, `DUMPFILE_PLACEHOLDER`, `USFEFILE_PLACEHOLDER` lines the SLURM job fills in). |
| `slurm/Ta10W_110_jobs`, `slurm/Ta10W_112_jobs` | SLURM array scripts. One array task per case line; builds a run-specific input from the template and launches PFDD. |
| `pilot_banner.py` | The PILOT banner (shown when the GUI opens). |

Generated files for each search are written under `studies/<study name>/` so the
top level stays clean:

```
studies/<name>/
    search_state.json         the search state (brackets, per-round history)
    search_state.undo.json    snapshot stack powering undo
    cases_<tag>.txt           the per-round cases file fed to SLURM
    dumps_local/              dumps fetched from the cluster (delete when done)
```

---

## Requirements

- **Python 3** with **numpy** and **matplotlib** (matplotlib only for the MRSSP
  explorer and `plot_crss.py`; the core search needs only the standard library).
- **tkinter** for the GUIs (`sudo apt-get install python3-tk`, or the python.org
  build / `brew install python-tk` on macOS). The CLI needs no tkinter.
- A cluster with **SLURM** and a built **PFDD** executable for the actual runs.

---

## Quick start (GUI)

```bash
python3 stress_state_gui.py
```

1. **Study & slip system** — name the study (creates `studies/<name>/`) and pick
   `{110}` or `{112}`. The cases/jobs tag follows the slip system.
2. **Temperatures** — min, max, step.
3. **Initial stress bracket for binary search** — the initial bracket `[lo, hi]` (same sign;
   0 may sit at one end, but the bracket must not straddle 0), the number of
   the number of search rounds, and the extend settings. The **final CRSS resolution** box updates automatically as bracket width x (1/2)^rounds.
4. **Stress state** — shear couple or uniaxial; click **Visualize / set
   stress...** to open the MRSSP explorer, choose an orientation while watching
   what each {110}/{112} plane feels, and **Apply to search**.
5. **Alloy USFE files** — the per-temperature USFE name pattern (`USFE_{T}`),
   produced by the PILOT lattice/USFE tool and numbered to match your
   temperatures. Point it at the USFE folder and **Check files exist** to verify.
6. **Actions** — Init, Generate, Update (grade moved/not-moved after inspecting
   dumps), Extend, Undo, Report, **Send to cluster** (scp the cases file, and
   optionally the template/job and the USFE files, to remote directories),
   Fetch dumps, Delete local dumps, and Plot CRSS (prompts for a legend label).

The heavy PFDD compute runs on the cluster; this tool prepares each round's
inputs and interprets the results.

---

## Command-line workflow

Edit the constants at the top of `crss_search_alloy.py`, then:

```bash
python3 crss_search_alloy.py init         # create search_state.json
python3 crss_search_alloy.py generate     # write the cases file for this round
sbatch --array=0-14 Ta10W_110_jobs        # (generate prints the exact command)
#   ... wait for the array to finish, then fetch + inspect the dumps ...
python3 crss_search_alloy.py update       # grade moved/not-moved per temperature
python3 crss_search_alloy.py report       # show brackets / CRSS estimates
python3 crss_search_alloy.py extend       # widen out-of-range temps (if needed)
python3 crss_search_alloy.py undo         # revert the last generate/update/extend
```

---

## The stress-state explorer (shear vs uniaxial)

Open it with **Visualize / set stress...**.  The bar chart shows what each
{110} and {112} plane of the [-111] zone would feel:

- **Shear couple mode** applies sigma = tau0 (b (x) m + m (x) b) with the MRSSP
  set by chi.  The bars are the **resolved (glide) shear** tau/tau0 on each
  plane -- the shear stress driving the screw dislocation along its Burgers
  vector b on that plane.  {110} bars are green, {112} bars are blue.  The chi
  preset buttons are labeled with the plane each angle puts the MRSSP on.

- **Uniaxial mode** applies sigma = s0 (a (x) a) along an axis a.  A general
  load resolves onto each plane as three tractions, shown as three bars per plane:
  - **shear** -- the shear along b (this *is* the "resolved shear"; it is what
    drives glide, and equals the shear-couple bar when the couple is applied),
  - **Escaig** -- the in-plane shear perpendicular to b (does not drive glide,
    but can affect the core / cross-slip),
  - **normal** -- the stress acting perpendicular to the plane (opens/closes it;
    a non-glide component absent from the pure shear couple).

"Escaig" and "normal" are the extra non-glide components a uniaxial load introduces 
that a pure shear couple does not.

## Cases file format (the alloy link)

Each non-comment line fully specifies one PFDD run:

```
# <T_K>  <seed>  <USFE_file>  <sigma line>
0    734  USFE_0    sigma     initial ... delta ... #tau 0.08 chi +0 deg
100  512  USFE_100  sigma     initial ... delta ... #tau 0.08 chi +0 deg
```

The third column is the **USFE file for that temperature** — this is what makes
the sweep alloy-aware. The SLURM script reads it and points PFDD at
`USFE_DIR/<USFE_file>`, so make sure the USFE files you generated with the
lattice/USFE tool are in the `USFE_DIR` set in the SLURM script, named to match
your `usfe_pattern`.

The `sigma` line ends with `#tau <value>`, which both tags the dump filename and
records the load magnitude; a shear couple also records `chi`, a uniaxial load
records its `axis`.

---

## Setting up the SLURM run

The templates and SLURM scripts are provided and their array logic does not need
changing, but a few paths are cluster-specific. In each `slurm/*_jobs` script,
edit the clearly-marked `# EDIT` items:

- the `#SBATCH --account`, `--partition`, and `-o` / `-e` log paths,
- `SCRATCH`, `CODE_SRC` (your PFDD build), `USFE_DIR` (where your USFE files are),
  `DUMP_BASE`, `LOG_BASE`, and the `module load` lines.

The script parses each case line into `T`, `SEED`, `USFE_FILE`, and the `sigma`
line, then `sed`-substitutes the template's `temperature`, `seed`, `sigma`,
`DUMPFILE_PLACEHOLDER`, and `USFEFILE_PLACEHOLDER`. Personal cluster details from
the original have been replaced with placeholders.

---

## Fetching dumps and grading

After a round's array finishes, copy its dumps down to inspect in OVITO:

- **GUI:** *Fetch dumps...* — fill in the cluster user/host/root; files land in
  `studies/<name>/dumps_local/`. Then *Update (grade)...* lists each temperature
  with **moved / did not move** radio buttons.
  Each cluster field has a **remember** checkbox: tick it to save that value and
  auto-fill it next time (shared across the Fetch and Send dialogs for the username
  and hostname); untick to forget it. Remembered values live in `~/.pilot_stress_state.json`.
- **CLI:** edit and run `fetch_dumps.sh` (it reads `REMOTE_USER`, `REMOTE_HOST`,
  `REMOTE_ROOT`, `FILE_PREFIX`, `CASES`, `LOCAL_DEST` from its top or the
  environment), then `python3 crss_search_alloy.py update`.

When you no longer need the local dumps, **Delete local dumps** (GUI) or
`rm -rf studies/<name>/dumps_local` removes the large files.

To grade automatically instead of by eye, implement `did_dislocation_move()` in
`crss_search_alloy.py` to parse your dump format; `update` will then use it.

---

## Plotting

```bash
python3 plot_crss.py studies/<name>/search_state.json
# several planes on shared axes, with labels and unit conversion:
python3 plot_crss.py s110.json s112.json --labels "{110}" "{112}" --mu 161000 --unit MPa
```

CRSS at each temperature is the bracket midpoint; the error bar is the
half-width, so it shrinks as you record more rounds. In the GUI, **Plot CRSS vs
T** opens this for the current study.

---

## How it fits into PILOT

`PILOT_lattice_usfe` builds the composition-resolved `USFE_<i>` files; this tool
links one per temperature, applies a stress state, and drives PFDD across
temperature to locate CRSS(T). The two halves meet at the PFDD solver and share
the same banner and look-and-feel.
