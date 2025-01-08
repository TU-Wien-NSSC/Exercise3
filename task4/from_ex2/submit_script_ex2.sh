#! /bin/bash
# Number of cores (up to 40)
#SBATCH -c 40
# Output file for results
#SBATCH --output=task4_results.txt
# Clear the environment
module purge > /dev/null 2>&1
# Set OMP_NUM_THREADS to match SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

RESOLUTION=2048
ITERATIONS=100000
NAME=reference

for THREADS in 1 2 4 8 10 20 40; do
    export OMP_NUM_THREADS=$THREADS
    echo "Running with $THREADS threads..." >> task4_results.txt
    
    for SCHEDULE in 1 2 3; do
        case $SCHEDULE in
            1) SCHEDULE_NAME="static";;
            2) SCHEDULE_NAME="static,1";;
            3) SCHEDULE_NAME="dynamic";;
        esac
        export OMP_SCHEDULE="${SCHEDULE_NAME}"
        echo "  Testing schedule policy: $SCHEDULE_NAME" >> task4_results.txt
        echo "./solver $NAME $RESOLUTION $ITERATIONS" >> task4_results.txt
        ./solver $NAME $RESOLUTION $ITERATIONS 
        
    done
done