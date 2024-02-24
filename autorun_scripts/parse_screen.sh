#!/bin/bash
screen_num="$STY"

IFS='.' read -r -a parts <<< "$screen_num"


second_section="${parts[1]}"
integer_value=$((second_section))

directory="screen_loaders/screenrunner.$integer_value.txt"


read -r -a args_array < $directory
echo ${args_array[@]}
python3 "${args_array[@]}"

