#!/bin/bash

for ((i=0;i<5;i++)); do 
    file_contents=$(<./auto_run_files/load_directory.$i.txt)
    echo "this is contents: $file_contents"
done
