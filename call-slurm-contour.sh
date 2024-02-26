#!/bin/bash

# Check if the argument is provided
if [ $# -ne 3 ]; then
    echo "Usage: $0 <job_name> <delay> <mode>"
    exit 1
fi

# Set the job name variable
job_name="$1"
delay=$2
mode=$3

# Replace the comment in the script
sed -i "s/^#SBATCH --job-name=.*/#SBATCH --job-name=$job_name/g" ./autorun_scripts/slurm-contour.sh
echo "Job name updated to: $job_name"

sbatch ./autorun_scripts/slurm-contour.sh "$delay" "$mode"
