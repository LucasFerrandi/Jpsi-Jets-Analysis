#!/bin/bash

BASE_DIR="$HOME/alice/Jpsi-Jets-Analysis/JpsiWorkDir/PIDEfficiency/PIDEfficiencyConverter/IdasElectronMaps_DQ_LHC24_pass1_skimmed_V0candidates"
CONFIG_TEMPLATE="configPIDEffConverter.json"
OUTPUT_DIR="output/IdasMaps"
N="$1" # numero de subdirs a processar, se quiser processar todos, deixe em branco


# mkdir -p "$OUTPUT_DIR"

COUNT=0

for SUBDIR in "$BASE_DIR"/*/; do

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

    mv AnalysisResults.root "$OUTPUT_DIR/AnalysisResults_${SUBDIR_NAME}.root"

    rm "$TMP_JSON"

done
