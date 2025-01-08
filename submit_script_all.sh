#!/bin/bash

# Array of thread counts
THREAD_COUNTS=(1 2 4 8 10 20 40)

# Array of scheduling policies
SCHEDULE_POLICIES=("static" "static,1" "dynamic")

# Loop through each combination and create an individual batch script
for THREADS in "${THREAD_COUNTS[@]}"; do
  for SCHEDULE in "${SCHEDULE_POLICIES[@]}"; do
    # Create a unique batch script for this configuration
    cat << EOF > solver_batch_${THREADS}_${SCHEDULE}.sh
#!/bin/bash
#SBATCH --time=00:25:00
#SBATCH --job-name=solver_test
#SBATCH --output=solver_output_${THREADS}_${SCHEDULE}.txt
#SBATCH -c ${THREADS}  # Set the number of cores to the current thread count

# Clear the environment
module purge > /dev/null 2>&1

# Set OpenMP environment variables
export OMP_NUM_THREADS=${THREADS}
export OMP_SCHEDULE="${SCHEDULE}"

# Display debug information
echo "Running with ${OMP_NUM_THREADS} threads and ${OMP_SCHEDULE} scheduling policy"

# Run the solver program
./solver reference 2048 100000
EOF

    # Submit the batch script
    sbatch solver_batch_${THREADS}_${SCHEDULE}.sh
  done
done
