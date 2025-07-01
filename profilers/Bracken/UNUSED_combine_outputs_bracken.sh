#!bin/bash/sh

module load python3

stem=/restricted/projectnb/pathoscope/work/aubrey/newMB_02_25_meta_benchmark/profilers/Bracken
code_stem=/restricted/projectnb/pathoscope/work/aubrey/newMB_02_25_meta_benchmark/code/Bracken-2.9
## files=$(ls ${stem}/outputs/*.bracken)

${code_stem}/analysis_scripts/combine_bracken_outputs.py \
    --files $files \
    --output ${stem}/combined_files/test_lluch.txt

for subfolder in ${stem}/outputs/*/; do
    files=$(ls "${subfolder}"/*.bracken)
        ${code_stem}/analysis_scripts/combine_bracken_outputs.py \
            --files $files \
            --output "${stem}/combined_files/$(basename "$subfolder")_combined.txt"
done
