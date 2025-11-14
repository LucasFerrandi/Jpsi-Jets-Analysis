#!/bin/bash
JSON="$1"
o2-analysis-bc-converter -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-timestamp -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-event-selection -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-ft0-corrected-table -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-track-propagation -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-trackselection -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-track-to-collision-associator -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-tracks-extra-converter -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-centrality-table -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-collision-converter -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-zdc-converter -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-pid-tpc-base -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-pid-tpc-full -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-pid-tof-base -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-pid-tof-full -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-hf-track-index-skim-creator  -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-hf-candidate-creator-2prong -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-hf-candidate-selector-d0 -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-hf-derived-data-creator-d0-to-k-pi -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-hf-pid-creator -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-multiplicity-table -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-je-jet-deriveddata-producer -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-je-jet-deriveddata-trigger-producer -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-je-jet-finder-d0-data-charged -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error | \
o2-analysis-je-jet-hf-fragmentation -b --aod-writer-keep AOD/JETDISTTABLE/0 --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error
