#!/bin/bash
screen_num="$STY"

IFS='.' read -r -a parts <<< "$screen_num"


second_section="${parts[1]}"
integer_value=$((second_section))

directory=$(<./screen_loaders/screenrunner.$integer_value.txt)

PYTHON_SCRIPT="src/parallel.py"
for ((i=0; i<20; i++)); do
    echo "script $i:"
    time_delay_path="$directory/Time_delay_$i"
    python3 "$PYTHON_SCRIPT" "$time_delay_path"
   
done
