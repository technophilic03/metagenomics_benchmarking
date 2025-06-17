#!/bin/bash

# Define the output file
output_file="mess_manifest_test.txt"

file_path="/restricted/projectnb/pathoscope/work/aubrey/newMB_02_25_meta_benchmark/MeSS/out_samples/mess_sample_100k_err/fastq"

# Find all .fastq.gz files recursively and write to output file
find "$file_path" -type f -name "*_R1.fq.gz" > "$output_file"
sort $output_file -o $output_file

echo "File paths written to $output_file"
