#!/bin/bash
JSON="$1"
BATCH="$2"

o2-analysis-dq-efficiency-with-assoc "$BATCH" --configuration json://"$JSON" --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000 --aod-writer-keep AOD/RTDIELECTRONALL/0| \
o2-analysis-dq-model-converter-mc-reduced-event "$BATCH" --configuration json://"$JSON" --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000