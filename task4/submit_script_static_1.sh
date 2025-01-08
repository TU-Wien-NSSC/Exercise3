#!/bin/bash
#SBATCH --time=00:35:00
#SBATCH --job-name=solver_test
#SBATCH --output=solver_output_%x.txt  # Output file based on the scheduling policy name
#SBATCH -c 40  # Reserve 40 CPUs for the entire job

THREAD_COUNTS=(1 2 4 8 10 20 40)
SCHEDULE="static,1"

# Output file
OUTPUT_FILE="solver_output_static,1.txt"


# Write headers to each output file
echo "Results for scheduling policy: static,1" > $OUTPUT_FILE
echo "Threads,Total Runtime (s)" >> $OUTPUT_FILE

# Loop through each thread count
for THREADS in "${THREAD_COUNTS[@]}"; do
  # Set OpenMP environment variables
  export OMP_NUM_THREADS=${THREADS}
  export OMP_SCHEDULE="${SCHEDULE}"

  # Display debug information
  echo "Running with ${THREADS} threads and ${SCHEDULE} scheduling policy"

  # Run the solver program and capture runtime
  RUNTIME=$(./solver reference 2048 100000 | grep "Total runtime" | awk '{print $3}')

  echo "${THREADS},${RUNTIME}" >> $OUTPUT_FILE_STATIC

done
