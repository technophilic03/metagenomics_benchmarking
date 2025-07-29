library(rentrez)

na_accession <- "code/analysis/total genomes/metascope_na_accession.txt"
success_file <- "code/analysis/total genomes/metascope_success_conted.txt"
failed_file <- "code/analysis/total genomes/metascope_failed_conted.txt"

accession <- readLines(na_accession)

if (file.exists(success_file)) {
  existing_success <- readLines(success_file)
  result <- existing_success
} else {
  result <- c()
}

if (file.exists(failed_file)) {
  existing_failed <- readLines(failed_file)
  failed <- existing_failed
} else {
  failed <- c()
}

processed_count <- length(result) + length(failed)
start_index <- processed_count + 1

message("Starting from accession ", start_index, " out of ", length(accession))

if (start_index <= length(accession)) {
  for (i in start_index:length(accession)) {
    search <- entrez_search(db = "nucleotide",
                            term = accession[i],
                            use_history = TRUE)
    summary <- tryCatch({
      entrez_summary(db = "nucleotide", web_history = search$web_history)
    }, error = function(e){
      message("Error: ", e$message)
      return(NULL)
    })
    
    if (!is.null(summary) && !is.null(summary$taxid)) {
      result <- c(result, summary$taxid)
    } else {
      failed <- c(failed, accession[i])
    }
    
    message(i, " out of ", length(accession))
    
    write.table(result, file = success_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(failed, file = failed_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}