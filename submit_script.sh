#!/bin/bash
# Number of cores
#SBATCH -c 4
# Runtime of this jobs is less then 10 minutes
#            (hh:mm:ss)
#SBATCH --time=00:10:00
# Clear the environment
module purge > /dev/null 2>&1
# Set OMP_NUM_THREADS to the same value as -c
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_SCHEDULe="dynamic"

echo "Number of threads: ${OMP_NUM_THREADS}"
echo "Scheduling policy: ${OMP_SCHEDULE}"

# You can start several programs with one script file/submission
./solver num1 256 10000
./solver num2 512 10000
#./solver num3 1024 10000
