#!/bin/bash

# Runs o2-analysis-pid-efficiency for all Electron Maps. Ouput saved in output/IdasMapsPromptNonPNevts<...>/LHC25b<...>/AnalysisResults_SUBDIRNAME.root

MAPS_DIR="$HOME/alice/datasets/IdasElectronMaps_DQ_LHC24_pass1_skimmed_V0candidates"
CONFIG_TEMPLATE="configPIDEff.json"
N="$1" #Number of maps (for test)
OUTPUT_DIR="output/IdasMapsPromptNonP20Runs" #! ATENTION: Must be changed by hand!


MC_LIST=(LHC25b14 LHC25b15 LHC25b16 LHC25b17)

DATASETS_DIR="/home/ferrandi/alice/datasets/MC/derivedMC/myTableMakerMCNoPID25-11-19/higherStats/20Runs/"
AO2D_LIST=()

for i in LHC25b14 LHC25b15 LHC25b16 LHC25b17; do
    # AO2D_LIST+=("$DATASETS_DIR/$i/jobs/034${i}/hy_343${i}000/0001/AO2D.root")
    AO2D_LIST+=("@$DATASETS_DIR$i/input.txt")
done #! ATENTION: Must be changed by hand!


for i in "${!MC_LIST[@]}"; do

    MC="${MC_LIST[$i]}"
    AO2D="${AO2D_LIST[$i]}"

    echo "===== Processing $MC ====="
    echo "AO2D: $AO2D"
    echo "=========================="

    COUNT=0

    for SUBDIR in "$MAPS_DIR"/*/; do # Loop over maps

        ((COUNT++))
        if [[ -n "$N" && "$COUNT" -gt "$N" ]]; then
            break
        fi

        echo "Processing $SUBDIR"

        EFFMAP_FILE=$(find "$SUBDIR" -maxdepth 1 -type f -name "effMap*.root" | head -n 1)

        if [[ -z "$EFFMAP_FILE" ]]; then
            continue
        fi

        TMP_JSON="tmp_config.json"
        cp "$CONFIG_TEMPLATE" "$TMP_JSON"

        # substitute effMap
        sed -i "s|\"effMapPath\":.*|\"effMapPath\": \"${EFFMAP_FILE}\",|" "$TMP_JSON"

        # substitute AO2D
        sed -i "s|\"aod-file-private\":.*|\"aod-file-private\": \"${AO2D}\",|" "$TMP_JSON"

        time \
        o2-analysis-pid-efficiency \
            --configuration json://"$TMP_JSON" \
            --shm-segment-size 8000000000 \
            --aod-memory-rate-limit 500000000 -b | \
        o2-analysis-dq-model-converter-mc-reduced-event \
            --configuration json://"$TMP_JSON" \
            --shm-segment-size 8000000000 \
            --aod-memory-rate-limit 500000000 -b

        SUBDIR_NAME=$(basename "$SUBDIR")

        mkdir -p "$OUTPUT_DIR/$MC"
        mv AnalysisResults.root \
        "$OUTPUT_DIR/$MC/AnalysisResults_${SUBDIR_NAME}.root"

        rm "$TMP_JSON"

    done
done