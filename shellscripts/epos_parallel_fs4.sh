#!/bin/bash

# Check if FROM is provided as a command line argument
if [ -z "$1" ]; then
    echo "Usage: $0 <FROM>"
    exit 1
fi

FROM=$1
#FROM=1 # always change to largest (#+1) in resultfile_#

N_RUNS=10 # 20
OPTNSNM="prueba_EOS2"
BASE_DIR=$PWD
OPTNS="${BASE_DIR}/${OPTNSNM}.optns"
PIDS=()

mkdir -p "${BASE_DIR}/logfiles"
mkdir -p "${BASE_DIR}/output"
OUTPUT_DIR="${BASE_DIR}/output"

for JOBNUM in $(seq 0 $(($N_RUNS-1))); do
    OUTPUT_NUM=$(($FROM + $JOBNUM))
    echo "UUU* --- Starting job #$JOBNUM with output #$OUTPUT_NUM --- "
    echo "UUU* With random seed values near seedj = "
    echo $(date '+%N')
    echo "UUU* and seedi = "
    echo $(date '+%N')

    # Run directory creation
    RUN_DIR="${BASE_DIR}/run_${JOBNUM}"
    mkdir -p "${RUN_DIR}"
    cd "$RUN_DIR"

    # Renaming & copying options file
    NEW_OPTNS="${OPTNSNM}_${OUTPUT_NUM}"
    cp "$OPTNS" ./"${NEW_OPTNS}.optns"

    # START THE SIMULATION
    #${EPO}scripts/epos -root "${NEW_OPTNS}" &> "${BASE_DIR}/logfiles/run_${JOBNUM}_out${OUTPUT_NUM}.log" &
    $MYDIR/epos_install/bin/epos -root "${NEW_OPTNS}" &> "${BASE_DIR}/logfiles/run_${JOBNUM}_out${OUTPUT_NUM}.log" &
    PIDS[$JOBNUM]=$!

    echo "UUU* Job with PID $! started. Now I will wait 2 seconds to make sure that different random seed is generated."
    sleep 2s
    cd "$BASE_DIR"
done

echo "UUU* All jobs sent. Now waiting for execution."
for pid in "${PIDS[@]}"; do
    wait $pid
    echo "Job with PID ${pid} done."
done

echo "UUU* All jobs done. Collecting output root files & deleting temporary run directories."

for JOBNUM in $(seq 0 $(($N_RUNS-1))); do
    OUTPUT_NUM=$(($FROM + $JOBNUM))
    RUN_DIR="${BASE_DIR}/run_${JOBNUM}"
    mv "${RUN_DIR}/z-${OPTNSNM}_${OUTPUT_NUM}.root" "${OUTPUT_DIR}/"
    rm -r "$RUN_DIR"
done

echo "KONIEC."
