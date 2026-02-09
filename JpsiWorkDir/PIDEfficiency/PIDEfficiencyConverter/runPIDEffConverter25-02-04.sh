#!/bin/bash
JSON="$1"
BATCH="$2" #"-b" has to be the 2nd argument!
# FileIn="$2" #Example: ~/alice/JPsiDerivedDatasetSemiMerged/alice/cern.ch/user/a/alihyperloop/jobs/0122/hy_1220648/AOD/001/AO2D.root  --aod-file $2
#OutputDirectory="$3"
#echo "FileIn: $FileIn"
#--aod-writer-json test.json
time \
# o2-analysis-dq-efficiency-with-assoc $BATCH --configuration json://$JSON --aod-writer-keep AOD/RTDIELEEXTRA/0,AOD/RTDIELECTRON/0,AOD/RTDIELECTRONALL/0 --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000 

o2-analysis-pid-efficiency --configuration json://"$JSON" "$BATCH" --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000 | \
o2-analysis-dq-model-converter-mc-reduced-event --configuration json://"$JSON" "$BATCH" --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000