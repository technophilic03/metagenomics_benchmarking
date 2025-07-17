# Take in arguments from bash script
args <- commandArgs(trailingOnly = TRUE)

readPath1 <- args[1]
readPath2 <- args[2]
indexDir <- args[3]
expTag <- args[4] 
outDir <- args[5]
tmpDir <- args[6]
threads <- args[7]
targets <- stringr::str_split(args[8], ",")[[1]]
filters <- stringr::str_split(args[9], ",")[[1]]
taxDB <- args[10]
priors_df <- args[11]
blast_db <- args[12]

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
output_id_csv <- metascope_id(final_map, input_type = "csv.gz", aligner = "bowtie2",
                              db = "ncbi",
                              num_species_plot = 0,
                              maxitsEM = 100,
                              accession_path = taxDB,
                              tmp_dir = tmpDir,
                              update_bam = TRUE,
                              quiet = FALSE,
                              out_dir = outDir,
                              priors_df = priors_df)

# MetaBlast
message("running MetaBlast")
metascope_id_path <- output_id_csv
metascope_blast(metascope_id_path = metascope_id_path, 
                bam_file_path = list.files(tmpDir, ".updated.bam$",full.names = TRUE)[1], # THIS ASSUMES EACH SAMPLE HAS THEIR OWN TMP DIRECTORY
                tmp_dir = tmpDir, 
                out_dir = outDir, 
                sample_name = expTag, 
                fasta_dir = NULL,
                num_results = 100, 
                num_reads = 50, 
                hit_list = 10,
                num_threads = threads, 
                db_path = blast_db, 
                quiet = FALSE,
                db = "ncbi", 
                accession_path = taxDB)
metascope_blast_path <- gsub("metascope_id.csv","metascope_blast.csv", metascope_id_path) 
message("MetaBlast COMPLETE")

message("Running blast reassignment")
blast_reassignment(metascope_blast_path, 
                   species_threshold = 0.5, 
                   num_hits = 100, 
                   blast_tmp_dir = paste0(tmpDir, "/blast"), 
                   out_dir = outDir, 
                   sample_name = expTag)
message("BLAST REASSIGNMENT COMPLETE")

message(capture.output(Sys.time() - now))
