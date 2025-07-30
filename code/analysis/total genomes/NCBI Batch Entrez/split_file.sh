#!/bin/bash

# Prompt the user to enter the input file
read -p "Enter the input file: " input_file

# Check if the input file exists
if [ ! -f "$input_file" ]; then
    echo "Input file not found!"
    exit 1
fi

# Create a directory to store the output files
output_dir="output_metascope"
mkdir -p "$output_dir"

# Counter for naming output files
file_count=1

# Loop through each line in the input file
while IFS= read -r line; do
    # Check if the line is "//"
    if [ "$line" == "//" ]; then
        # Write "//" to corresponding output file
        echo "//" >> "$output_dir/section$file_count.txt"
        # Increment file count
        ((file_count++))
    else
        # Write line to corresponding output file
        echo "$line" >> "$output_dir/section$file_count.txt"
    fi
done < "$input_file"

echo "Separation complete. Output files are stored in the '$output_dir' directory."
