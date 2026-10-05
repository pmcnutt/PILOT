# PFDD CRSS Temperature Sweep

Automation for computing the **critical resolved shear stress (CRSS)** of a
screw dislocation as a function of temperature using
[PFDD](https://github.com/lanl/Phase-Field-Dislocation-Dynamics-PFDD-2.x)
phase-field dislocation dynamics simulations on a SLURM cluster.

Instead of hand-editing one input file per run, this workflow drives a **SLURM
job array**: each temperature/stress condition runs as an independent task
(one per node allocation), and a Python driver implements a **bisection search**
that locates each temperature's CRSS in a handful of rounds rather than an
exhaustive stress scan.

The example files target a ⟨111⟩ screw dislocation in bcc tungsten loaded on the
(110) plane, but the approach generalizes to any material, glide plane, or
dislocation character.

---

## Contents

| File | Purpose |
|------|---------|
| `crss_search.py` | Bisection-search driver. Generates the cases file each round, tracks per-temperature stress brackets, and updates them from run results. |
| `run_W_110_TSarray.sh` | SLURM job-array script. One array task per line of the cases file; builds a run-specific input from a template and launches the simulation. |
| `cases_W_110.txt` | The per-run manifest: one line per simulation, listing temperature, random seed, and the full stress (`sigma`) line. |
| `in.W_110_template` | The PFDD input template. Its `temperature`, `seed`, `sigma`, and dump-path lines are overwritten per run; all other settings pass through unchanged. |

The only thing **not** included is the **PFDD executable** itself — you supply
that (see [Setup](#setup)). PFDD is available from LANL:
[Phase-Field-Dislocation-Dynamics-PFDD-2.x](https://github.com/lanl/Phase-Field-Dislocation-Dynamics-PFDD-2.x).

---

## How it works

Each line of `cases_W_110.txt` fully specifies one simulation:

```
# <T_K>  <seed>  <sigma line, including a trailing "#tau <value>" comment>
300  356  sigma     initial -0.032660 0.032660 0.000000 0.000000 0.016330 0.016330 delta 0.0 0.0 0.0 0.0 0.0 0.0 #tau 0.040 chi +0 deg
```

When submitted, the SLURM array launches one task per line. Each task:

1. reads **its own line** (indexed by `$SLURM_ARRAY_TASK_ID`);
2. parses the temperature, seed, stress line, and the `#tau` value;
3. copies the input **template** and uses `sed` to substitute the
   `temperature`, `seed`, `sigma`, and dump-path lines;
4. runs PFDD, writing a uniquely named dump, e.g.
   `.../W/300K/W_screw_110_300K_0p040.dump`.

Because the tasks are independent, the cluster starts as many as node
availability allows and runs the rest as space frees up — no task waits for any
other, and a crash in one does not affect the others.

The **bisection driver** wraps this in a search. For each temperature it holds a
bracket `[τ_lo, τ_hi]` known to contain the CRSS. Each round tests the bracket
midpoint; if the dislocation glides, the CRSS is at or below that stress (upper
bound moves down), otherwise it is above (lower bound moves up). Six rounds
narrow an initial `[0, 0.08 μ]` bracket to a width of
0.08 μ / 2⁶ ≈ 1.25 × 10⁻³ μ — 6 simulations per temperature instead of the
64 an equivalent-resolution scan would need.

---

## Setup

1. **PFDD executable.** Build PFDD from the LANL source
   ([Phase-Field-Dislocation-Dynamics-PFDD-2.x](https://github.com/lanl/Phase-Field-Dislocation-Dynamics-PFDD-2.x))
   and note the path to the resulting binary — this is the one piece not
   included here. Point `EXE` in `run_W_110_TSarray.sh` at it.

2. **Input template.** `in.W_110_template` is included and ready to use. It is a
   standard PFDD input file with four lines set up to be overwritten each run:
   the `temperature`, `seed`, and `sigma` lines (whose values act as
   placeholders) and the `dump` line, whose path is the literal token
   `DUMPFILE_PLACEHOLDER`. Every other setting — cell size, lattice, solver,
   material preset, mobilities — passes through unchanged, so adapt those to
   your own system by editing the template directly. If you rename it, update
   `IN_TEMPLATE` in the SLURM script to match.

   > If you build your own template from scratch, keep exactly one line each
   > beginning with `temperature`, `seed`, and `sigma`, and put
   > `DUMPFILE_PLACEHOLDER` where the dump path goes — these are the anchors the
   > `sed` substitutions rely on.

3. **Edit the paths** near the top of both scripts:
   - `run_W_110_TSarray.sh`: `EXE`, `IN_TEMPLATE`, `CASES_FILE`, `OUT_BASE`,
     and the `module load` / SLURM header lines for your cluster.
   - `crss_search.py`: `OUT_BASE`, `DUMP_NAME`, and the search parameters
     `TEMPS`, `TAU_LO`, `TAU_HI`, `N_ROUNDS`, `CHI_DEG`, `SEED_MIN/MAX`.

   > The `DUMP_NAME` in the Python driver **must match** the dump naming in the
   > SLURM script so the driver can find each run's output.

4. **Stress lines.** The `sigma` lines encode a pure shear couple on the glide
   plane. Generate them with the companion GUI
   ([pfdd-mrssp-gui](https://github.com/pmcnutt/pfdd-mrssp-gui)); the driver can
   also import that module directly to compute them, and falls back to
   equivalent built-in math if it is not on the `PYTHONPATH`.

---

## Usage

### Option A — automated bisection search (recommended)

```bash
# One time: initialize the search state
python3 crss_search.py init

# --- each round ---
python3 crss_search.py generate          # writes cases_W_110.txt for this round
sbatch --array=0-14 run_W_110_TSarray.sh # submit (adjust range to your # of temps)
#   ...wait for ALL array tasks to finish (check: squeue -u $USER)...
python3 crss_search.py update            # record moved/not-moved per temperature
# --- repeat generate/submit/update for N_ROUNDS rounds ---

python3 crss_search.py report            # CRSS estimate ± uncertainty per temperature
```

At the `update` step the driver shows each run's dump path and asks `m` (moved)
or `n` (didn't). To automate this, implement `did_dislocation_move(dumpfile)`
in `crss_search.py` to parse your dump format and return `True`/`False`; the
`update` step then records all results with no prompting.

All state (brackets, seeds, dump paths, per-round history) is stored in
`search_state.json`, so every run remains reproducible after the fact. Delete
that file to start a fresh search.

### Option B — run a fixed set of conditions

If you just want to run a specific list rather than search, write the conditions
directly into `cases_W_110.txt` (one line per run) and submit:

```bash
sbatch --array=0-$(( $(grep -cvE '^[[:space:]]*(#|$)' cases_W_110.txt) - 1 )) run_W_110_TSarray.sh
```

The `--array` range on the command line overrides the one in the script header
and is computed from the number of (non-comment) lines in the cases file, so you
never have to keep the two in sync by hand.

---

## Cases file format

- One run per line: `<temperature_K>  <seed>  <sigma line>`.
- Lines beginning with `#` and blank lines are ignored.
- The line index (ignoring comments/blanks) is the SLURM array index: first
  line → task 0, second → task 1, and so on.
- The `sigma` line **must** end with a `#tau <value>` comment — this value is
  parsed to tag the dump filename (`0.055` → `0p055`).
- Nothing requires temperatures to be unique. Multiple lines at the same
  temperature (e.g. several seeds or stresses) are fine; each produces its own
  input, log, and dump.

---

## Output

Per array task, in `OUT_BASE/<T>K/`:

- `W_screw_110_<T>K_<tag>.dump` — the simulation dump.

In the submission directory:

- `in.W_110_T<T>_<tag>` — the generated input file for that run.
- `out_W_110_T<T>_<tag>.txt` — PFDD stdout.
- `slurm_W_110.<jobid>_<taskid>.{out,err}` — SLURM logs (the applied
  temperature, seed, and stress are echoed here for the record).

---

## Notes and caveats

- **Keep `update` for after a full round.** It expects a result for every
  temperature in the round; run it only once all array tasks have finished.
  `sbatch --wait` can be used to block until an array completes if you script
  the rounds.
- **Random seeds and finite-temperature scatter.** With a random seed per run,
  the moved/not-moved outcome near the CRSS is stochastic (thermal activation is
  probabilistic). The converged value is effectively one sample; for publication
  you may want to repeat the final round with several seeds per temperature and
  report the spread. The cases format supports this with no code changes.
- **Define "moved" up front.** Near the athermal transition the classification
  can depend on run length. Fix a quantitative criterion (e.g. core displacement
  beyond one Peierls valley by the end of the run) and apply it consistently.
- **Array range vs. temperature count.** In Option A, remember to set
  `--array=0-(N-1)` to match your number of temperatures (15 temps → `0-14`).
  Use `--array=0-14%K` to cap how many run at once if your cluster limits
  concurrent jobs.

---

## Citation / related

Stress states are generated with the companion tool
[pfdd-mrssp-gui](https://github.com/pmcnutt/pfdd-mrssp-gui).
