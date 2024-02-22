#!/bin/bash

PYTHON_SCRIPT="parallel-Contours.py"
directories=("data/data-ContourMaps/R0.vs.RelationStrength" "data/data-ContourMaps/R0.vs.TimeDelay" "data/data-ContourMaps/RelationStrength.vs.TimeDelay")

run_cases(){
    for ((i=1; i<4; i++)); do
    	screen -S "contourscreen.mode$i"
    	
    	screen -S "contourscreen.mode$i" -X stuff "python3 $PYTHON_SCRIPT ${directories[$i]} $i\n"
    done
}

run_cases

echo "end"

