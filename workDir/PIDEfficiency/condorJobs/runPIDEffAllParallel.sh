#!/usr/bin/env bash
'''
This creates a json with 1 AOD and 1 electron PID-efficiency map and runs pid-efficiency taks
'''
set -euo pipefail

AOD="$1"
EFFMAP_FILE="$2"
BASE_JSON="$3"

TMP_JSON="$(mktemp --suffix=.json)"
trap 'rm -f "$TMP_JSON"' EXIT

cp "$BASE_JSON" "$TMP_JSON"

python3 - "$TMP_JSON" "$EFFMAP_FILE" "$AOD" <<'PY'
import json, sys

path, effmap, aod = sys.argv[1:]

with open(path) as f:
    data = json.load(f)

data["effMapPath"] = effmap
data["aod-file-private"] = aod

with open(path, "w") as f:
    json.dump(data, f, indent=2)
PY

time o2-analysis-pid-efficiency \
    --configuration json://"$TMP_JSON" \
    --shm-segment-size 8000000000 \
    --aod-memory-rate-limit 500000000 -b