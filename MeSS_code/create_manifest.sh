#!/bin/bash

# Define the output file
output_file="mess_manifest.txt"

file_path="/restricted/projectnb/pathoscope/work/aubrey/newMB_02_25_meta_benchmark/MeSS/out_samples"

# Find all .fastq.gz files recursively and write to output file
find $file_path -type f -path "*/fastq/*_R1.fq.gz" > "$output_file"

sort $output_file -o $output_file

# Optional: Print a message
echo "File paths written to $output_file"
