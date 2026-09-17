#!/bin/bash
JSON="$1"
BATCH="$2" #"-b" has to be the 2nd argument!
time \

o2-analysis-dq-table-reader-with-assoc "$BATCH" --configuration json://"$JSON" --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000
#  | \
# o2-analysis-dq-model-converter-event-extended "$BATCH" --configuration json://"$JSON" --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000