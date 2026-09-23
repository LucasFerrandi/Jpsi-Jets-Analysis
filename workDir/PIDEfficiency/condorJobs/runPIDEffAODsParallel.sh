#!/usr/bin/env bash

# This substitutes AO2D and EFFMAP_FILE in a json and runs PID efficiency task with it

set -euo pipefail

AO2D="$1"
EFFMAP_FILE="$2"
BASE_JSON="$3"
SUB_OUTDIR="$4"
CLUSTER="$5"
PROCESS="$6"

printf "\n===== PID efficiency job =====\n"
printf "AO2D file   : %s\n" "$AO2D"
printf "Eff map     : %s\n" "$EFFMAP_FILE"
printf "Cluster/Proc: %s / %s\n" "$CLUSTER" "$PROCESS"
printf "================================\n\n"

TMP_JSON="$(mktemp --suffix=.json)"
trap 'rm -f "$TMP_JSON"' EXIT

cp "$BASE_JSON" "$TMP_JSON"

# substitute effMap
sed -i "s|\"effMapPath\":.*|\"effMapPath\": \"${EFFMAP_FILE}\",|" "$TMP_JSON"

# substitute AO2D
sed -i "s|\"aod-file-private\":.*|\"aod-file-private\": \"${AO2D}\",|" "$TMP_JSON"

JOBDIR="${SUB_OUTDIR}/cluster${CLUSTER}_process${PROCESS}"
mkdir -p "$JOBDIR"
cd "$JOBDIR"
# cd $SUB_OUTDIR

time o2-analysis-pid-efficiency \
    --configuration json://"$TMP_JSON" \
    --shm-segment-size 8000000000 \
    --aod-memory-rate-limit 500000000 -b