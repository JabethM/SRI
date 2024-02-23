#!/bin/bash
screen_num="$STY"

IFS='.' read -r -a parts <<< "$screen_num"


second_section="${parts[1]}"
integer_value=$((second_section))

directory=$(<./contours/screenrunner.$integer_value.txt)
echo "$integer_value"
echo "$directory"
read -r -a args_array < $directory
#python3 "${args_array[@]}"
