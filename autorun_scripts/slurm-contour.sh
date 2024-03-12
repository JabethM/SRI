#!/bin/bash
#
#
#SBATCH --job-name=Low_2
#SBATCH --partition=long
#SBATCH --time=3-00:00:00
#SBATCH --array=0-9
#SBATCH --mem=16G
#
########################################################################
 
# Decide which input file this Slurm task will process.
# We use the special $SLURM_ARRAY_TASK_ID variable, which tells the
# task which one it is among the whole job array.
directories=("data/data-ContourMaps/RelationStrength.vs.TimeDelay" "data/data-ContourMaps/RelationStrength.vs.R0" "data/data-ContourMaps/R0.vs.TimeDelay")
delay=$1
mode=$2
# Now do our "processing" on the input file
python3 -m src.parallel.Contours-repeats "${directories[$(( mode % 3 ))]}/${delay}Delay/" $mode ${SLURM_ARRAY_TASK_ID} "10"
 
# Finally let's sleep for 60 seconds so that we can watch the tasks
# get processed for the purposes of this tutorial.
echo "$delay $mode ${SLURM_ARRAY_TASK_ID}"
echo "Done"
