#!/bin/bash
#SBATCH --job-name=W_110_TSsweep
#SBATCH --nodes=1 --ntasks-per-node=8
#SBATCH --time=500:00:00
#SBATCH --array=0-3
#SBATCH -o slurm_W_110.%A_%a.out
#SBATCH -e slurm_W_110.%A_%a.err

cd $SLURM_SUBMIT_DIR
module load intel/20
module load openmpi
ulimit -s unlimited

EXE=pfdd_openmpi
IN_TEMPLATE=in.W_110_template
CASES_FILE=cases_W_110.txt
OUT_BASE=/home/OUTPUT/W

# ---- Read this task's case line (skipping comments/blank lines) ----
mapfile -t CASES < <(grep -vE '^[[:space:]]*(#|$)' "$CASES_FILE")

NCASES=${#CASES[@]}
if [ "$SLURM_ARRAY_TASK_ID" -ge "$NCASES" ]; then
    echo "ERROR: array index $SLURM_ARRAY_TASK_ID >= number of cases ($NCASES)" >&2
    exit 1
fi

CASE="${CASES[$SLURM_ARRAY_TASK_ID]}"

# Field 1 = temperature; field 2 = seed; the rest = the full sigma line
T=$(echo "$CASE" | awk '{print $1}')
SEED=$(echo "$CASE" | awk '{print $2}')
SIGMA_LINE=$(echo "$CASE" | awk '{$1=""; $2=""; sub(/^ +/,""); print}')

# Sanity-check the seed (must be a positive integer)
if ! [[ "$SEED" =~ ^[0-9]+$ ]]; then
    echo "ERROR: second field of case line is not an integer seed:" >&2
    echo "  $CASE" >&2
    exit 1
fi

# Extract tau from the trailing "#tau <value>" comment, e.g. 0.055
TAU=$(echo "$SIGMA_LINE" | sed -n 's/.*#tau[[:space:]]*\([0-9.eE+-]*\).*/\1/p')
if [ -z "$TAU" ]; then
    echo "ERROR: could not parse '#tau <value>' from sigma line:" >&2
    echo "  $SIGMA_LINE" >&2
    exit 1
fi
TAU_TAG="${TAU//./p}"        # 0.055 -> 0p055

echo "Array task $SLURM_ARRAY_TASK_ID / job $SLURM_ARRAY_JOB_ID"
echo "T = ${T} K   tau = ${TAU}   seed = ${SEED}   sigma: $SIGMA_LINE"
echo "START $(date) on $(hostname)"

DUMP_DIR="${OUT_BASE}/${T}K"
mkdir -p "$DUMP_DIR"

INFILE="in.W_110_T${T}_${TAU_TAG}"
DUMPFILE="${DUMP_DIR}/W_screw_110_${T}K_${TAU_TAG}.dump"

# Build the run-specific input from the template:
#   1) replace the temperature line
#   2) replace the seed line with this case's seed
#   3) replace the entire sigma line with this case's sigma line
#   4) replace the dump placeholder with the real path
sed \
    -e "s/^temperature[[:space:]].*/temperature     ${T}/" \
    -e "s/^seed[[:space:]].*/seed            ${SEED}/" \
    -e "s|^sigma[[:space:]].*|${SIGMA_LINE}|" \
    -e "s#DUMPFILE_PLACEHOLDER#${DUMPFILE}#" \
    "$IN_TEMPLATE" > "$INFILE"
    # optional per-T material preset:
    # -e "s/^material[[:space:]].*/material         preset W110_${T}K dAdc 0.0 0.0 0.0 aflag 2/" \

echo "Input file: $INFILE"
echo "Dump file:  $DUMPFILE"

mpirun -np $SLURM_NTASKS "$EXE" -in "$INFILE" > "out_W_110_T${T}_${TAU_TAG}.txt"

echo "END $(date)"
