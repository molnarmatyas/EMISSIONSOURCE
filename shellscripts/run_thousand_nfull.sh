#!/bin/bash

for BATCH in $(seq 2 19); do
    echo "Doing job batch $(($BATCH * 10 + 1))"
    ./epos_parallel_fs4.sh $(($BATCH * 10 + 1))
done
echo "All job batches done."
