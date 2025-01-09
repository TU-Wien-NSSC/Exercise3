#! /bin/bash
# Number of cores (up to 40)
#SBATCH -c 40
# Time
#SBATCH --time=00:59:00
# Output file for results
#SBATCH --output=task4_results_%x.txt  # Output file based on the scheduling policy name
# Clear the environment
module purge > /dev/null 2>&1

# Parameters
RESOLUTION=2048
ITERATIONS=100000
NAME=reference
SCHEDULE_NAME="static"

# Loop through thread counts
for THREADS in 1 2 4 8 10 20 40; do
    export OMP_NUM_THREADS=$THREADS
    export OMP_SCHEDULE="${SCHEDULE_NAME}"

    # Append thread and policy information to the result file
    echo "Threads: $THREADS" >> task4_results_${SCHEDULE_NAME}.txt
    echo "Policy: $SCHEDULE_NAME" >> task4_results_${SCHEDULE_NAME}.txt

    # Run the solver and redirect its output to the result file
    ./solver $NAME $RESOLUTION $ITERATIONS >> task4_results_${SCHEDULE_NAME}.txt
done
