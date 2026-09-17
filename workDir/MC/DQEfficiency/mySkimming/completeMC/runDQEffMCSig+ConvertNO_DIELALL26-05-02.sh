#!/bin/bash
JSON="$1"
BATCH="$2" #"-b" has to be the 2nd argument!
echo "Used JSON: $JSON"
#--aod-writer-json test.json

time \
o2-analysis-dq-efficiency-with-assoc $BATCH --configuration json://$JSON --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000 | \
#--aod-writer-keep AOD/RTDIELECTRONALL/0 
o2-analysis-dq-model-converter-mc-reduced-event $BATCH --configuration json://$JSON --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000

# o2-analysis-dq-efficiency-with-assoc $BATCH --configuration json://$JSON --aod-writer-keep AOD/RTDIELEEXTRA/0,AOD/RTDIELECTRON/0,AOD/RTDIELECTRONALL/0 --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000 

# o2-analysis-dq-efficiency-with-assoc $BATCH --configuration json://$JSON  ‐‐aod-writer-resfile="~/alice/O2Physics/PWGJE/Tasks/JPsiWorkDir/JPsiMC/DQEfficiency/AOD.root" #--severity error --shm-segment-size 12000000000 --aod-writer-json aodWriterTempConfig.json -b | \