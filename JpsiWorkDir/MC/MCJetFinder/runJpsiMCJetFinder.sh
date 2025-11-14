#!/bin/bash
JSON="$1"
BATCH="$2" #"-b" has to be the 2nd argument!
# FileIn="$2" #Example: ~/alice/JPsiDerivedDatasetSemiMerged/alice/cern.ch/user/a/alihyperloop/jobs/0122/hy_1220648/AOD/001/AO2D.root  --aod-file $2
#OutputDirectory="$3"
echo "JSON: $JSON"
#echo "FileIn: $FileIn"
#--aod-writer-json test.json

time \
o2-analysis-je-jpsi-fragmentation $BATCH --configuration json://$JSON | \
o2-analysis-je-jet-finder-dielectron-mcp-charged $BATCH --configuration json://$JSON | \
o2-analysis-je-jet-deriveddata-producer $BATCH --configuration json://$JSON | \
o2-analysis-mccollision-converter $BATCH --configuration json://$JSON | \  #--aod-parent-access-level 1 | \
o2-analysis-dq-efficiency-with-assoc $BATCH --configuration json://$JSON #| \
# o2-analysis-mc-converter $BATCH --configuration json://$JSON --aod-parent-access-level 1 #| \
# o2-analysis-je-jet-finder-dielectron-data-charged $BATCH --configuration json://$JSON
# o2-analysis-je-jet-deriveddata-writer  $BATCH --configuration json://$JSON | \
# o2-analysis-je-jet-finder-dielectron-data-charged  $BATCH --configuration json://$JSON

# o2-analysis-je-jet-deriveddata-producer  $BATCH --configuration json://$JSON | \
# o2-analysis-track-propagation  $BATCH --configuration json://$JSON | \
# o2-analysis-tracks-extra-v002-converter $BATCH --configuration json://$JSON | \
# o2-analysis-timestamp  $BATCH --configuration json://$JSON | \
# o2-analysis-trackselection $BATCH --configuration json://$JSON | \
# o2-analysis-ft0-corrected-table $BATCH --configuration json://$JSON | \
# o2-analysis-event-selection $BATCH --configuration json://$JSON | \
# o2-analysis-multiplicity-table $BATCH --configuration json://$JSON #| \
# o2-analysis-centrality-table $BATCH --configuration json://$JSON | \
# o2-analysis-run2bcinfos-converter $BATCH --configuration json://$JSON






# o2-analysis-je-jet-finder-dielectron-data-charged $BATCH --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# #--aod-writer-resdir $OutputDirectory | \
# o2-analysis-je-jet-deriveddata-producer $BATCH --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-je-jet-deriveddata-producer-dummy-dielectron $BATCH --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-je-jet-deriveddata-writer $BATCH --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-dq-table-reader-with-assoc $BATCH --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-centrality-table $BATCH --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-ft0-corrected-table $BATCH --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-event-selection $BATCH --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-multiplicity-table $BATCH --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-timestamp $BATCH --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-track-propagation $BATCH --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-trackselection $BATCH --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-je-jet-deriveddata-trigger-producer $BATCH --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \

#TRAIN 240251:
# o2-analysis-dq-filter-pp-with-association $BATCH --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-dq-table-maker-with-assoc $BATCH --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-fwdtrack-to-collision-associator $BATCH --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-fwdtrackextension $BATCH --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| 
# o2-analysis-pid-tof-base $BATCH --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-pid-tof-beta $BATCH --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-pid-tof-full $BATCH --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-pid-tpc $BATCH --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-pid-tpc-base $BATCH --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1

#ExtraTasks(maybe only for old AODs):


#o2-analysis-mft-tracks-converter $BATCH --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error
#o2-analysis-mft-tracks-converter $BATCH --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error





















##o2-analysis-track-to-collision-associator $BATCH --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
##o2-analysis-centrality-table $BATCH --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
##o2-analysis-bc-converter $BATCH --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
#o2-analysis-timestamp $BATCH --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
#o2-analysis-tracks-extra-converter $BATCH --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
#o2-analysis-ft0-corrected-table  $BATCH --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
#o2-analysis-trackselection $BATCH --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
#o2-analysis-track-propagation $BATCH --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
#o2-analysis-event-selection $BATCH --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
#o2-analysis-multiplicity-table $BATCH --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
#o2-analysis-je-jet-deriveddata-producer $BATCH --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
#o2-analysis-je-jet-deriveddata-producer-dummy $BATCH --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
#o2-analysis-je-jet-deriveddata-trigger-producer $BATCH --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
#o2-analysis-je-jet-finder-data-charged $BATCH --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
#o2-analysis-je-jet-finder-d0-data-charged $BATCH --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
#o2-analysis-je-jet-hf-fragmentation $BATCH --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error
##o2-analysis-je-jet-finder-charged-qa $BATCH --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
##o2-analysis-je-jet-substructure $BATCH --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error
