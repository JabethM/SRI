#!/bin/bash
PYTHON_SCRIPT="parallel-Contours.py"
directories=("data/data-ContourMaps/R0.vs.RelationStrength" "data/data-ContourMaps/R0.vs.TimeDelay" "data/data-ContourMaps/RelationStrength.vs.TimeDelay")

run_cases(){
    for ((i=1; i<4; i++)); do
	echo "$PYTHON_SCRIPT ${directories[$(( i % 3 ))]} $i" > "contours/screenrunner.$i.txt"

	screen -S "$i.contourscreen.mode" -d -m
    	#screen -S "$i.contourscreen.mode" -X stuff "longjob -28day -c ./parse_contours.sh\n"
    done
}

run_cases

echo "end"

