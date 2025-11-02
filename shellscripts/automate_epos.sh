#!/bin/bash

# Define paths
BASE_DIR="/star/u/momatyas/work/epos/simulation_runs/sequential_epos_runs/200GeV_5nfull5nfreeze"
TEMP_DIR="$BASE_DIR/tmp"
OUTPUT_DIR="$BASE_DIR/output"
OPTNSNM="prueba_EOS2"
OPTNS_TEMPLATE="$BASE_DIR/${OPTNSNM}.optns"
EPOS_EXEC="/gpfs/mnt/gpfs01/star/pwg_tasks/cf01/EPOS4/epos4-mod.img"

# Number of sequential runs
N_RUNS=400

# Create necessary directories
mkdir -p "$TEMP_DIR"
mkdir -p "$OUTPUT_DIR"

# Start the loop
for JOBNUM in $(seq 1 $N_RUNS); do
    echo " UUU* --- Running job #$JOBNUM --- "

    # Define run directory
    RUN_DIR="$TEMP_DIR/run_${JOBNUM}"
    mkdir -p "$RUN_DIR"
    cd "$RUN_DIR"

    # Copy the original configuration file
    cp "$OPTNS_TEMPLATE" ./"${OPTNSNM}.optns"

    # Generate unique seed values
    SEEDJ=$(( JOBNUM * 12 + 34 ))
    SEEDI=$(( JOBNUM * 56 + 78 ))
    echo " UUU* SEEDJ: $SEEDJ, SEEDI: $SEEDI"

    # Modify the .optns file with the new seeds (in-place modification)
    sed -i -E "s/set seedj [0-9]+/set seedj $SEEDJ/" "$OPTNSNM.optns"
    sed -i -E "s/set seedi [0-9]+/set seedi $SEEDI/" "$OPTNSNM.optns"

    # Run the containerized simulation
    "$EPOS_EXEC" "$OPTNSNM"

    # Rename and move the output .root file
    if [[ -e "z-${OPTNSNM}.root" ]]; then
        ROOT_OUT="$OUTPUT_DIR/z-${OPTNSNM}_${JOBNUM}.root"
        mv "z-${OPTNSNM}.root" "$ROOT_OUT"
        echo " UUU* Saved output: $ROOT_OUT"
    else
        echo " UUU* WARNING: No output .root file found for JOBNUM = $JOBNUM"
    fi

    # Cleanup temp files (optional)
    rm -rf "$RUN_DIR"
done

echo " UUU* All runs completed!"

