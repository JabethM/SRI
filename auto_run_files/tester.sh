#!/bin/bash

for ((i=0;i<20;i++)); do
    echo $i > "load_directory.$i.txt"
done
