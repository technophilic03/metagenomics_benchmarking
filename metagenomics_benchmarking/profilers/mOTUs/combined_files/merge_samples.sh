#!/bin/bash

# Directory containing the mOTUs profile files

stem=/restricted/projectnb/pathoscope/work/aubrey/newMB_02_25_meta_benchmark

PROFILE_DIR=${stem}/profilers/mOTUs/outputs  # Change this to your actual directory
OUTPUT_FILE=${stem}/profilers/mOTUs/combined_files/merged_mOTUs_table.txt

# Find all .txt files in the directory (assumed to be mOTUs profile outputs)
PROFILE_FILES=$(find "$PROFILE_DIR" -type f -name "*.profile" | paste -sd "," -)

# Load dependencies
module load samtools
module load bwa/0.7.17
module load htslib
module load boost
module load miniconda

# Activate conda environment
loc=${stem}/my_conda_env
mamba activate $loc

# Run the motus merge command
motus merge -i "$PROFILE_FILES" > "$OUTPUT_FILE"

echo "Merged profile written to $OUTPUT_FILE"
