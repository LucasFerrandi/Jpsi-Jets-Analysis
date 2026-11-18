#!/bin/bash
JSON="$1"
BATCH="$2" #"-b" has to be the 2nd argument!
time \

o2-analysis-je-dq-fragmentation-function "$BATCH" --configuration json://"$JSON" --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000 --aod-parent-access-level 1 --aod-writer-keep AOD/JDIELALLPLUSJET/0 | \
o2-analysis-mccollision-converter "$BATCH" --configuration json://"$JSON" --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000 --aod-parent-access-level 1 | \
o2-analysis-tracks-extra-v002-converter "$BATCH" --configuration json://"$JSON" --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000 --aod-parent-access-level 1 | \
o2-analysis-je-jet-deriveddata-producer "$BATCH" --configuration json://"$JSON" --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000 --aod-parent-access-level 1 | \
o2-analysis-lf-mc-centrality "$BATCH" --configuration json://"$JSON" --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000 --aod-parent-access-level 1 | \
o2-analysis-je-jet-finder-dielectron-mcd-charged "$BATCH" --configuration json://"$JSON" --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000 --aod-parent-access-level 1 | \
o2-analysis-je-jet-finder-dielectron-mcp-charged "$BATCH" --configuration json://"$JSON" --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000 --aod-parent-access-level 1 | \
o2-analysis-je-jet-matching-mc-dielectron-ch "$BATCH" --configuration json://"$JSON"  --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000 --aod-parent-access-level 1 | \
o2-analysis-je-estimator-rho "$BATCH" --configuration json://"$JSON"  --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000 --aod-parent-access-level 1 | \
o2-analysis-je-jet-luminosity-calculator "$BATCH" --configuration json://"$JSON" --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000 --aod-parent-access-level 1 | \
o2-analysis-je-jet-luminosity-producer "$BATCH" --configuration json://"$JSON" --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000 --aod-parent-access-level 1 | \
o2-analysis-ft0-corrected-table "$BATCH" --configuration json://"$JSON" --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000 --aod-parent-access-level 1 | \
o2-analysis-multcenttable "$BATCH" --configuration json://"$JSON" --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000 --aod-parent-access-level 1 | \
o2-analysis-event-selection-service "$BATCH" --configuration json://"$JSON" --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000 --aod-parent-access-level 1 | \
o2-analysis-trackselection "$BATCH" --configuration json://"$JSON" --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000 --aod-parent-access-level 1 | \
o2-analysis-propagationservice "$BATCH" --configuration json://"$JSON" --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000 --aod-parent-access-level 1
# --aod-parent-base-path-replacement "alien://;/home/alitrain/train-workdir/testdata/LFN"
