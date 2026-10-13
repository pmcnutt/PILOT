#!/bin/bash
set -euo pipefail

# Fetch PFDD dump files for one study from a cluster into a local folder, so you
# can inspect them (e.g. in OVITO) and grade moved/not-moved for `update`.
#
# Configure by editing the values below OR by exporting the same names as
# environment variables before running (the PILOT stress-state GUI does the
# latter). Anything left blank must be filled in before this will run.
#
#   REMOTE_USER   your username on the cluster
#   REMOTE_HOST   the cluster hostname
#   REMOTE_ROOT   absolute path on the cluster that contains the <T>K/ dump dirs
#   FILE_PREFIX   dump filename prefix, before _<T>K_<tau>.dump
#   CASES         local cases file listing the runs to fetch
#   LOCAL_DEST    local folder the dumps are copied into

# ============ EDIT THESE (or export as env vars) ============
REMOTE_USER="${REMOTE_USER:-}"                 # e.g. your cluster username
REMOTE_HOST="${REMOTE_HOST:-}"                 # e.g. login.mycluster.edu
REMOTE_ROOT="${REMOTE_ROOT:-}"                 # e.g. /scratch/$USER/OUTPUT/dumps/<study>
FILE_PREFIX="${FILE_PREFIX:-Ta10W_screw_110}"  # filename before _<T>K_<tau>.dump
CASES="${CASES:-cases_Ta10W_110.txt}"          # local copy of the cases file
LOCAL_DEST="${LOCAL_DEST:-./dumps_local}"      # where files land locally
# ============================================================

if [ -z "$REMOTE_USER" ] || [ -z "$REMOTE_HOST" ] || [ -z "$REMOTE_ROOT" ]; then
    echo "ERROR: set REMOTE_USER, REMOTE_HOST and REMOTE_ROOT (edit this file" >&2
    echo "       or export them / use the GUI) before running." >&2
    exit 1
fi
if [ ! -f "$CASES" ]; then
    echo "ERROR: cases file not found: $CASES" >&2
    exit 1
fi

LIST="filelist.txt"

# Build the list of files to pull, one relative path per line
: > "$LIST"
while read -r T SEED USFE REST; do
    tau=$(sed -n 's/.*#tau[[:space:]]*\([0-9.eE+-]*\).*/\1/p' <<< "$REST")
    if [ -z "$tau" ]; then
        echo "WARN: no #tau found for T=${T}, skipping" >&2
        continue
    fi
    tau_tag="${tau//./p}"                       # 0.04 -> 0p04
    echo "${T}K/${FILE_PREFIX}_${T}K_${tau_tag}.dump" >> "$LIST"
done < <(grep -vE '^[[:space:]]*(#|$)' "$CASES")

echo "Fetching $(wc -l < "$LIST") files from ${REMOTE_HOST}:${REMOTE_ROOT}"

mkdir -p "$LOCAL_DEST"
rsync -avz --partial --progress \
    --files-from="$LIST" \
    "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_ROOT}/" \
    "$LOCAL_DEST/"

echo "Done -> $LOCAL_DEST"
