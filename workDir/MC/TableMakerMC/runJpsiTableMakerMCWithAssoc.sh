#!/bin/bash
JSON="$1"
BATCH="$2" #"-b" has to be the 2nd argument!
# FileIn="$2" #Example: ~/alice/JPsiDerivedDatasetSemiMerged/alice/cern.ch/user/a/alihyperloop/jobs/0122/hy_1220648/AOD/001/AO2D.root  --aod-file $2
#OutputDirectory="$3"
echo "JSON: $JSON"
#echo "FileIn: $FileIn"
#--aod-writer-json test.json
time \
# o2-analysis-dq-efficiency-with-assoc $BATCH --configuration  json://$JSON  | \
#--severity error --shm-segment-size 12000000000 --aod-writer-json aodWriterTempConfig.json -b  | \
o2-analysis-dq-table-maker-mc-with-assoc $BATCH --configuration  json://$JSON  | \
o2-analysis-timestamp $BATCH --configuration  json://$JSON  | \
o2-analysis-event-selection $BATCH --configuration  json://$JSON  | \
o2-analysis-multiplicity-table $BATCH --configuration  json://$JSON  | \
o2-analysis-trackselection $BATCH --configuration  json://$JSON  | \
o2-analysis-track-propagation $BATCH --configuration  json://$JSON  | \
o2-analysis-fwdtrackextension $BATCH --configuration  json://$JSON  | \
# o2-analysis-track-propagation-tester $BATCH --configuration  json://$JSON  | \
# o2-analysis-onthefly-tracker $BATCH --configuration  json://$JSON  | \
o2-analysis-tracks-extra-v002-converter $BATCH --configuration  json://$JSON  | \
o2-analysis-mccollision-converter $BATCH --configuration  json://$JSON  | \
o2-analysis-track-to-collision-associator $BATCH --configuration  json://$JSON | \
#Bloco de pid
o2-analysis-pid-tpc-base $BATCH --configuration  json://$JSON  | \
o2-analysis-ft0-corrected-table $BATCH --configuration  json://$JSON  | \
o2-analysis-pid-tof-full $BATCH --configuration  json://$JSON  | \
o2-analysis-pid-tof-base $BATCH --configuration  json://$JSON  | \
o2-analysis-pid-tof-beta $BATCH --configuration  json://$JSON  | \
o2-analysis-pid-tpc $BATCH --configuration  json://$JSON



# o2-analysis-pid-tpase $BATCH --configuration  json://$JSON
# o2-analysis-trackextension $BATCH --configuration  json://$JSON | \
# o2-analysis-pid-toase $BATCH --configuration  json://$JSON | \
# o2-analysis-pid-tof $BATCH --configuration  json://$JSON | \
# o2-analysis-pid-tof-full $BATCH --configuration  json://$JSON | \
# o2-analysis-pid-toeta $BATCH --configuration  json://$JSON | \
# o2-analysis-pid-tpc-full $BATCH --configuration  json://$JSON | \
# o2-analysis-hf-mc-pid-tof $BATCH --configuration  json://$JSON | \
# o2-analysis-pid-tof-merge $BATCH --configuration  json://$JSON | \


# Giacomo Alocco 2022:
# o2-analysis-dq-table-maker-mc $BATCH --configuration  json://configLHC21i3f3.json --severity error --shm-segment-size 12000000000 --aod-writer-json aodWriterTempConfig.json -b | o2-analysis-timestamp $BATCH --configuration  json://configLHC21i3f3.json -b | o2-analysis-event-selection $BATCH --configuration  json://configLHC21i3f3.json -b | o2-analysis-multiplicity-table $BATCH --configuration  json://configLHC21i3f3.json -b | o2-analysis-trackselection $BATCH --configuration  json://configLHC21i3f3.json -b | o2-analysis-trackextension $BATCH --configuration  json://configLHC21i3f3.json -b | o2-analysis-pid-tof-base $BATCH --configuration  json://configLHC21i3f3.json -b | o2-analysis-pid-tof $BATCH --configuration  json://configLHC21i3f3.json -b | o2-analysis-pid-tof-full $BATCH --configuration  json://configLHC21i3f3.json -b | o2-analysis-pid-tof-beta $BATCH --configuration  json://configLHC21i3f3.json -b | o2-analysis-pid-tpc-full $BATCH --configuration  json://configLHC21i3f3.json -b | o2-analysis-track-propagation $BATCH --configuration  json://configLHC21i3f3.json -b


