#!/bin/bash

# Define the remote server address
REMOTE_SERVER="nssc_e11801674@tcad30.iue.tuwien.ac.at"
REMOTE_DIR="."  # Assumes the files are in the current directory on the remote machine

# Define the output file pattern
THREAD_COUNTS=(1 2 4 8 10 20 40)
SCHEDULE_POLICIES=("static" "static,1" "dynamic")

# Loop through each combination of thread counts and scheduling policies
for THREADS in "${THREAD_COUNTS[@]}"; do
  for SCHEDULE in "${SCHEDULE_POLICIES[@]}"; do
    # File pattern on remote server
    FILE_PATTERN="solver_output_${THREADS}_${SCHEDULE}.txt"

    # Use scp to transfer the files to the local machine
    echo "Transferring files for ${THREADS} threads with ${SCHEDULE} scheduling policy..."
    scp ${REMOTE_SERVER}:${REMOTE_DIR}/${FILE_PATTERN} ~/Desktop  # Transfers to current local directory
  done
done

echo "File transfer completed."
