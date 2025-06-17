#!bin/bash/sh

module load python3

stem=/restricted/projectnb/pathoscope/work/aubrey/newMB_02_25_meta_benchmark/profilers/Bracken

files=$(ls ${stem}/outputs/*.bracken)

${stem}/Bracken-2.9/analysis_scripts/combine_bracken_outputs.py \
    --files $files \
    --output ${stem}/combined_files/test_lluch.txt
