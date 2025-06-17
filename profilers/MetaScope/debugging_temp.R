# TEMPORARY SCRIPT
readPath1 <- "/restricted/projectnb/pathoscope/data/Mock_Lluch_2015/run_mock__5EB-95GR-2_R1.fastq.gz"
readPath2 <- "/restricted/projectnb/pathoscope/data/Mock_Lluch_2015/run_mock__5EB-95GR-2_R2.fastq.gz"
indexDir <- "/restricted/projectnb/pathoscope/reflib/2020_index_bowtie/"
expTag <- "run_mock__5EB-95GR-2"
outDir <- "/restricted/projectnb/pathoscope/work/aubrey/newMB_02_25_meta_benchmark/profilers/MetaScope/outputs"
tmpDir <- "/scratch/run_mock__5EB-95GR-2_tmp/"
threads <- 1
targets <- stringr::str_split("bacteria,fungi,viral", ",")[[1]]
filters <- stringr::str_split("human_mouse,phix174", ",")[[1]]
taxDB <- "/restricted/projectnb/pathoscope/data/blastdb/2024_accession_taxa/accessionTaxa.sql"


# ID step issues

# Time this!
now <- Sys.time()

suppressPackageStartupMessages({
  library(MetaScope)
  library(tidyverse)
})

# Using bt2 params
bt2_params <- "--local -R 2 -N 0 -L 25 -i S,1,0.75 -k 9 --score-min L,0,1.7"

# Align to targets
target_map <- align_target_bowtie(read1 = readPath1,
                                  read2 = readPath2,
                                  lib_dir = indexDir,
                                  libs =  targets,
                                  align_dir = tmpDir,
                                  align_file = expTag,
                                  overwrite = TRUE,
                                  threads = threads,
                                  quiet = FALSE,
                                  bowtie2_options = bt2_params)
message("TARGET STEP COMPLETE")
# Align to filters
output <- paste(file.path(tmpDir, expTag), "filtered", sep = ".")
final_map <- filter_host_bowtie(reads_bam = target_map,
                                lib_dir = indexDir,
                                libs = filters,
                                make_bam = FALSE,
                                output = output,
                                bowtie2_options = bt2_params,
                                threads = threads,
                                overwrite = TRUE,
                                quiet = FALSE)
message("FILTER STEP COMPLETE")

# MetaScope ID
message("running id step")
metascope_id(final_map, input_type = "csv.gz", aligner = "bowtie2",
             db = "ncbi",
             num_species_plot = 1,
             maxitsEM = 100,
             accession_path = taxDB,
             tmp_dir = tmpDir,
             update_bam = FALSE,
             quiet = FALSE,
             out_dir = outDir)

# input_file <- final_map
input_file <- final_map
input_type = "csv.gz"
aligner = "bowtie2"
db = "ncbi"
db_feature_table = NULL
accession_path = taxDB
tmp_dir = tmpDir
out_dir = outDir
convEM = 1 / 10000
maxitsEM = 25
update_bam = FALSE
num_species_plot = 1
out_fastas = FALSE
num_genomes = 100
num_reads = 50
quiet = FALSE

message(capture.output(Sys.time() - now))