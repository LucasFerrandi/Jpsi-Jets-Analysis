#!/usr/bin/env bash

#This creates .sub files and runs condor_submit for each tuple (MC Dataset, ElectronEffMap)


# Arguments
NMaps="${1:-}" # If none, all maps will be analyzed
N_AODS="${2:-}" # If none, all AODs will be analyzed


# Script stops when errors
set -euo pipefail

MC_LIST=(LHC25b16)

MAPS_DIR="/storage/ferrandi/datasets/IdasElectronMaps_DQ_LHC24_pass1_skimmed_V0candidates"
DATASETS_DIR="/storage/ferrandi/datasets/MC/derivedMC/myTableMakerMCNoPID25-11-19/completeDatasets"
BASE_JSON="/sampa/ferrandi/alice/Jpsi-Jets-Analysis/workDir/PIDEfficiency/condorJobs/configPIDEff.json"

OUTDIR="$PWD/outputCondor2"
mkdir -p "$OUTDIR/joblists"
# Loads ALICE env
export WORK_DIR="$HOME/alice/sw"

set +u # Does not interrupt when finding unbound variables in init.sh
source $HOME/alice/sw/slc9_x86-64/O2Physics/master-local1/etc/profile.d/init.sh
set -u

# cd $JpsiPIDDir

for MC in "${MC_LIST[@]}"; do
    AOD_LIST_FILE="$DATASETS_DIR/$MC/inputCorrect.txt"

    if [[ ! -f "$AOD_LIST_FILE" ]]; then
        echo "AOD list não encontrada: $AOD_LIST_FILE" >&2
        continue
    fi

    # Reads AOD file, ignoring comments and empty lines
    mapfile -t AODS < <(grep -vE '^[[:space:]]*($|#)' "$AOD_LIST_FILE")

    # Se o arquivo tiver só um .txt com outra lista dentro, faz flatten
    # Maybe useful for inputing chunks of aod
    # if (( ${#AODS[@]} == 1 )) && [[ "${AODS[0]}" == *.txt ]] && [[ -f "${AODS[0]}" ]]; then
    #     mapfile -t AODS < <(grep -vE '^[[:space:]]*($|#)' "${AODS[0]}")
    # fi

    if (( ${#AODS[@]} == 0 )); then
        echo "Nenhuma AOD válida em $AOD_LIST_FILE" >&2
        continue
    fi

    JOBLIST="$OUTDIR/joblists/AODs${MC}.txt"
    printf '%s\n' "${AODS[@]}" > "$JOBLIST"

    # Limits number of AODs
    if [[ -n "$N_AODS" ]]; then
        printf '%s\n' "${AODS[@]:0:$N_AODS}" > "$JOBLIST"
    else
        printf '%s\n' "${AODS[@]}" > "$JOBLIST"
    fi

    COUNT_MAPS=0
    for SUBDIR in "$MAPS_DIR"/*/; do
        # removes last "/"
        SUBDIR=${SUBDIR%/}
        ((++COUNT_MAPS))
        # Runs until Nmaps > COUNT_MAPS
        if [[ -n "$NMaps" && "$COUNT_MAPS" -gt "$NMaps" ]]; then
            echo "Final map reached"
            break
        fi

        echo "Processing MC=$MC map=$SUBDIR"
        #Gets the name after "_"
        MAP_NAME=${SUBDIR##*_}
        SUB_OUTDIR="$OUTDIR/${MC}_${MAP_NAME}"
        mkdir -p "$SUB_OUTDIR/logs"

        EFFMAP_FILE=$(find "$SUBDIR" -maxdepth 1 -type f -name "effMap*.root" | head -n 1)
        if [[ -z "$EFFMAP_FILE" ]]; then
            echo "  nenhum effMap*.root em $SUBDIR" >&2
            continue
        fi


        SUBFILE="$SUB_OUTDIR/${MC}_map${MAP_NAME}.sub"

        # Create 1 sub file for each (MC, Map), with all AODs as arguments
        cat > "$SUBFILE" <<EOF
universe   = vanilla
executable = runPIDEffAODsParallel.sh

should_transfer_files   = YES
when_to_transfer_output  = ON_EXIT
getenv = True

request_cpus   = 1
request_memory = 4GB

log    = $SUB_OUTDIR/logs/${MC}_map${MAP_NAME}.log
output = $SUB_OUTDIR/logs/${MC}_map${MAP_NAME}.\$(Cluster).\$(Process).out
error  = $SUB_OUTDIR/logs/${MC}_map${MAP_NAME}.\$(Cluster).\$(Process).err

arguments = \$(AOD) $EFFMAP_FILE $BASE_JSON $SUB_OUTDIR \$(Cluster) \$(Process)
queue AOD from $JOBLIST
EOF

        echo "Submitting: MC=$MC map=$MAP_NAME jobs=$(wc -l < "$JOBLIST")"
        condor_submit "$SUBFILE"
    done
done