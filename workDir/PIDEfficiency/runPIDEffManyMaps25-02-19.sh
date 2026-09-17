#!/bin/bash

# Runs o2-analysis-pid-efficiency for all Electron Maps. Ouput saved in output/IdasMapsPromptNonPNevts<...>/LHC25b<...>/AnalysisResults_SUBDIRNAME.root

BASE_DIR="$HOME/alice/datasets/IdasElectronMaps_DQ_LHC24_pass1_skimmed_V0candidates"
CONFIG_TEMPLATE="configPIDEff.json"
OUTPUT_DIR="output/IdasMapsPromptNonPNevts2"
MC="$1" # e.g LHC25b14
N="$2" # numero de subdirs a processar, se quiser processar todos, deixe em branco


# mkdir -p "$OUTPUT_DIR"

COUNT=0

for SUBDIR in "$BASE_DIR"/*/; do # Loop over maps

    ((COUNT++))

    if [[ -n "$N" && "$COUNT" -gt "$N" ]]; then
        break
    fi

    echo "Processing $SUBDIR"


    # pega o primeiro arquivo que começa com effMap
    EFFMAP_FILE=$(find "$SUBDIR" -maxdepth 1 -type f -name "effMap*.root" | head -n 1)

    if [ -z "$EFFMAP_FILE" ]; then
        echo "No effMap found in $SUBDIR"
        continue
    fi

    echo "Found map: $EFFMAP_FILE"

    # cria json temporário
    TMP_JSON="tmp_config.json"
    cp "$CONFIG_TEMPLATE" "$TMP_JSON"

    # substitui linha effMapPath
    sed -i "s|\"effMapPath\":.*|\"effMapPath\": \"${EFFMAP_FILE}\",|" "$TMP_JSON"

    # roda
    time \
    o2-analysis-pid-efficiency --configuration json://"$TMP_JSON" --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000 -b | \
    o2-analysis-dq-model-converter-mc-reduced-event --configuration json://"$TMP_JSON" --shm-segment-size 8000000000 --aod-memory-rate-limit 500000000 -b

    # nome baseado na subdir
    SUBDIR_NAME=$(basename "$SUBDIR")

    mv AnalysisResults.root "$OUTPUT_DIR/$MC/AnalysisResults_${SUBDIR_NAME}.root"

    rm "$TMP_JSON"

done
