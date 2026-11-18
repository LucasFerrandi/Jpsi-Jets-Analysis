#!/bin/bash
JSON="$1"
BATCH="$2" #"-b" has to be the 2nd argument!
time \

o2-analysis-event-selection-service "$BATCH" --configuration json://"$JSON" --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000 --aod-parent-access-level 2 | \
o2-analysis-tracks-extra-v002-converter "$BATCH" --configuration json://"$JSON" --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000 --aod-parent-access-level 2