#!bin/bash/sh

module load python3

stem=/restricted/projectnb/pathoscope/work/aubrey/newMB_02_25_meta_benchmark

files=$(ls ${stem}/profilers/Kraken2/outputs/*_taxa.txt)

echo $files > ${stem}/profilers/Kraken2/combined_files/test_lluch_samplenames.txt

${stem}/code/KrakenTools-1.2/combine_mpa.py \
    --input $files \
    --output ${stem}/profilers/Kraken2/combined_files/test_lluch.txt
