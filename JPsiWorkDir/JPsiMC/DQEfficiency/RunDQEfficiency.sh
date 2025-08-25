#!/bin/bash
JSON="$1"
BATCH="$2" #"-b" has to be the 2nd argument!
# FileIn="$2" #Example: ~/alice/JPsiDerivedDatasetSemiMerged/alice/cern.ch/user/a/alihyperloop/jobs/0122/hy_1220648/AOD/001/AO2D.root  --aod-file $2
#OutputDirectory="$3"
echo "JSON: $JSON"
#echo "FileIn: $FileIn"
#--aod-writer-json test.json
time \
o2-analysis-dq-efficiency-with-assoc $BATCH --configuration json://$JSON --aod-writer-keep AOD/RTDIELEEXTRA/0,AOD/RTDIELECTRON/0 --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000

# o2-analysis-dq-efficiency $BATCH --configuration json://$JSON --aod-writer-keep AOD/RTDIELECTRONALL/0
# o2-analysis-dq-efficiency-with-assoc $BATCH --configuration json://$JSON  ‐‐aod-writer-resfile="~/alice/O2Physics/PWGJE/Tasks/JPsiWorkDir/JPsiMC/DQEfficiency/AOD.root" \
#--severity error --shm-segment-size 12000000000 --aod-writer-json aodWriterTempConfig.json -b | \
# o2-analysis-timestamp --configuration json://$JSON -b | \
# o2-analysis-event-selection --configuration json://$JSON -b | \
# o2-analysis-multiplicity-table --configuration json://$JSON -b | \
# o2-analysis-trackselection --configuration json://$JSON -b | \
# o2-analysis-pid-tof-base --configuration json://$JSON -b | \
# o2-analysis-pid-tof --configuration json://$JSON -b | \
# o2-analysis-pid-tof-full --configuration json://$JSON -b | \
# o2-analysis-pid-tof-beta --configuration json://$JSON -b | \
# o2-analysis-track-propagation --configuration json://$JSON -b | \
# o2-analysis-tracks-extra-v002-converter --configuration json://$JSON -b | \
# o2-analysis-mccollision-converter --configuration json://$JSON -b | \
# o2-analysis-track-to-collision-associator --configuration json://$JSON -b | \
# o2-analysis-pid-tpc --configuration json://$JSON -b | \
# o2-analysis-pid-tpc-base --configuration json://$JSON -b | \
# o2-analysis-dq-table-maker-mc-with-assoc --configuration json://$JSON -b 

# o2-analysis-trackextension --configuration json://$JSON -b | \
# o2-analysis-pid-tpc-full --configuration json://$JSON -b | \





# Giacomo Alocco 2022:
# o2-analysis-dq-table-maker-mc --configuration json://configLHC21i3f3.json --severity error --shm-segment-size 12000000000 --aod-writer-json aodWriterTempConfig.json -b | o2-analysis-timestamp --configuration json://configLHC21i3f3.json -b | o2-analysis-event-selection --configuration json://configLHC21i3f3.json -b | o2-analysis-multiplicity-table --configuration json://configLHC21i3f3.json -b | o2-analysis-trackselection --configuration json://configLHC21i3f3.json -b | o2-analysis-trackextension --configuration json://configLHC21i3f3.json -b | o2-analysis-pid-tof-base --configuration json://configLHC21i3f3.json -b | o2-analysis-pid-tof --configuration json://configLHC21i3f3.json -b | o2-analysis-pid-tof-full --configuration json://configLHC21i3f3.json -b | o2-analysis-pid-tof-beta --configuration json://configLHC21i3f3.json -b | o2-analysis-pid-tpc-full --configuration json://configLHC21i3f3.json -b | o2-analysis-track-propagation --configuration json://configLHC21i3f3.json -b
