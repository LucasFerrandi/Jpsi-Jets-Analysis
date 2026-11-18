#!/bin/bash
JSON="$1"
BATCH="$2" #"-b" has to be the 2nd argument!
# time \

o2-analysis-je-dq-fragmentation-function "$BATCH" --configuration json://"$JSON" --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000 | \
o2-analysis-je-jet-finder-dielectron-data-charged "$BATCH" --configuration json://"$JSON" --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000 | \
o2-analysis-je-jet-substructure-dielectron "$BATCH" --configuration json://"$JSON" --aod-file @input_data.txt --aod-writer-json OutputDirector.json --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000
