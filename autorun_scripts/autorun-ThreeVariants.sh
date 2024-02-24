#!/bin/bash

folder="../data/data-TripleVariants"
PYTHON_SCRIPT="../src/parallel-ThreeDelays.py"

count=0
run_cases(){
    for dir in "$folder"/*/; do
	echo "$PYTHON_SCRIPT $dir" > "screen_loaders/screenrunner.$count.txt"
    	screen -S "$count.triplescreen" -d -m
    	
    	screen -S "$count.triplescreen" -X stuff "longjob -28day -c ./parse_screen.sh\n"
    	((count++))
    done
}

run_cases

echo "end"
