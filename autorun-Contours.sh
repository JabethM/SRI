#!/bin/bash

folder="data/data-TripleVariants"
PYTHON_SCRIPT="parallel-ThreeDelays.py"

count=0
run_cases(){
    for dir in "$folder"/*/; do
    	screen -S "triplescreen.$count"
    	
    	screen -S "triplescreen.$count" -X stuff "python3 $PYTHON_SCRIPT $dir\n"
    	((count++))
    done
}

run_cases

echo "end"

