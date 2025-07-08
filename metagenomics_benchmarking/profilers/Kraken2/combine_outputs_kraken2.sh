#!bin/bash/sh

module load python3

stem=/restricted/projectnb/pathoscope/work/aubrey/newMB_02_25_meta_benchmark

for this_batch in ${stem}/profilers/Kraken2/outputs/*/; do

    batch_name=$(basename "$this_batch")

    files=$(ls "${this_batch}"/*_taxa.txt)

    echo $files > "${stem}/profilers/Kraken2/combined_files/${batch_name}_samplenames.txt"

    ${stem}/code/KrakenTools-1.2/combine_mpa.py \
        --input $files \
        --output "${stem}/profilers/Kraken2/combined_files/${batch_name}_combined.txt"
done
