#!/bin/bash
JSON="$1"
# FileIn="$2" #Example: ~/alice/JPsiDerivedDatasetSemiMerged/alice/cern.ch/user/a/alihyperloop/jobs/0122/hy_1220648/AOD/001/AO2D.root  --aod-file $2
#OutputDirectory="$3"
echo "JSON: $JSON"
#echo "FileIn: $FileIn"
#--aod-writer-json test.json

time \
o2-analysis-je-jet-jpsi-fragmentation -b --configuration json://$JSON | \
o2-analysis-je-jet-finder-dielectron-data-charged -b --configuration json://$JSON
# o2-analysis-je-jet-deriveddata-writer  --configuration json://$JSON --aod-parent-access-level 2 | \
# o2-analysis-je-jet-finder-dielectron-data-charged  --configuration json://$JSON --aod-parent-access-level 2

# o2-analysis-je-jet-deriveddata-producer  --configuration json://$JSON --aod-parent-access-level 2 | \
# o2-analysis-track-propagation  --configuration json://$JSON --aod-parent-access-level 2 | \
# o2-analysis-tracks-extra-v002-converter --configuration json://$JSON --aod-parent-access-level 2 | \
# o2-analysis-timestamp  --configuration json://$JSON --aod-parent-access-level 2 | \
# o2-analysis-trackselection --configuration json://$JSON --aod-parent-access-level 2 | \
# o2-analysis-ft0-corrected-table --configuration json://$JSON --aod-parent-access-level 2 | \
# o2-analysis-event-selection --configuration json://$JSON --aod-parent-access-level 2 | \
# o2-analysis-multiplicity-table --configuration json://$JSON --aod-parent-access-level 2 | \
# o2-analysis-centrality-table --configuration json://$JSON --aod-parent-access-level 2 | \
# o2-analysis-run2bcinfos-converter --configuration json://$JSON --aod-parent-access-level 2






# o2-analysis-je-jet-finder-dielectron-data-charged -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# #--aod-writer-resdir $OutputDirectory | \
# o2-analysis-je-jet-deriveddata-producer -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-je-jet-deriveddata-producer-dummy-dielectron -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-je-jet-deriveddata-writer -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-dq-table-reader-with-assoc -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-centrality-table -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-ft0-corrected-table -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-event-selection -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-multiplicity-table -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-timestamp -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-track-propagation -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-trackselection -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-je-jet-deriveddata-trigger-producer -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \

#TRAIN 240251:
# o2-analysis-dq-filter-pp-with-association -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-dq-table-maker-with-assoc -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-fwdtrack-to-collision-associator -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-fwdtrackextension -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| 
# o2-analysis-pid-tof-base -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-pid-tof-beta -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-pid-tof-full -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-pid-tpc -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
# o2-analysis-pid-tpc-base -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1

#ExtraTasks(maybe only for old AODs):


#o2-analysis-mft-tracks-converter -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error
#o2-analysis-mft-tracks-converter -b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error





















##o2-analysis-track-to-collision-associator -b --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
##o2-analysis-centrality-table -b --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
##o2-analysis-bc-converter -b --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
#o2-analysis-timestamp -b --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
#o2-analysis-tracks-extra-converter -b --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
#o2-analysis-ft0-corrected-table  -b --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
#o2-analysis-trackselection -b --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
#o2-analysis-track-propagation -b --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
#o2-analysis-event-selection -b --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
#o2-analysis-multiplicity-table -b --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
#o2-analysis-je-jet-deriveddata-producer -b --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
#o2-analysis-je-jet-deriveddata-producer-dummy -b --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
#o2-analysis-je-jet-deriveddata-trigger-producer -b --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
#o2-analysis-je-jet-finder-data-charged -b --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
#o2-analysis-je-jet-finder-d0-data-charged -b --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
#o2-analysis-je-jet-hf-fragmentation -b --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error
##o2-analysis-je-jet-finder-charged-qa -b --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error --aod-parent-access-level 1| \
##o2-analysis-je-jet-substructure -b --configuration json://$JSON --aod-writer-json $OutputDirectory --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --min-failure-level error
