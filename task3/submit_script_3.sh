#!/bin/bash
# Number of cores
#SBATCH -c 80
# Runtime of this job is less than 10 minutes
#             (hh:mm:ss)
#SBATCH --time=00:10:00
# Clear Environment
module purge > /dev/null 2>&1

# Specify parameters
XMIN="-1.570796326794897"  # -pi/2
XMAX="1.570796326794897"   # pi/2
SAMPLES="10000000"

# Array of functions to run
FUNCTIONS=("SINX" "COS2XINV" "X4M5")

# Array of thread counts to use
THREAD_COUNTS=(1 5 10 20 40 80)


for FUNCTION in "${FUNCTIONS[@]}"; do
    rm "${FUNCTION}.txt"
    OUTPUT_FILE="${FUNCTION}.txt"

    for THREAD in "${THREAD_COUNTS[@]}"; do
        export OMP_NUM_THREADS=$THREAD
        ./mcint "$FUNCTION" "$XMIN" "$XMAX" "$SAMPLES" >> "$OUTPUT_FILE"
    done
done