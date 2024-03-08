#!/bin/bash


for file in *_*.*_*.json; do
	# Get the filename from the argument
	

	# Extract the parts of the filename using parameter expansion
	name="${file%.*}"   # Remove the extension
	extension="${file##*.}"  # Extract the extension

	numbers="${name#*.}"
	numbers="${numbers%.*}"

	numbers="${numbers//_/}"


	prefix="${name%.*}"    # Extract the prefix before the last underscore

	# Replace underscore with hyphen
	new_name="${prefix}-${numbers//_/}.json"

	# Rename the file
	mv "$file" "${prefix}-${numbers}.json"
done
