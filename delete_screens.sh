#!/bin/bash

for ((i=0;i<20;i++)); do
    screen -S "$i.contourscreen" -X quit
done
