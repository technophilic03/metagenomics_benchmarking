suppressPackageStartupMessages({
  library(tidyverse)
  library(rentrez)
})
na_accession <- read.table("code/analysis/total genomes/centrifuge_na_accession.txt")

library(rentrez)

retry_with_backoff <- function(expr, max_attempts = 6, base_delay = 1) {
  for (attempt in 1:max_attempts) {
    result <- tryCatch({
      return(expr)
    }, error = function(e) {
      if (attempt == max_attempts) {
        message("Failed after ", max_attempts, " attempts: ", e$message)
        return(NULL)
      }
      delay <- base_delay * (2 ^ (attempt - 1))
      message("Attempt ", attempt, " failed: ", e$message, ". Retrying in ", delay, " seconds...")
      Sys.sleep(delay)
      return(NULL)
    })
    if (!is.null(result)) return(result)
  }
  return(NULL)
}

result <- c()
failed_accessions <- c()
for (i in seq_along(na_accession$V1)) {
  current_accession <- na_accession$V1[i]
  message("Processing ", current_accession, " (", i, " out of ", nrow(na_accession), " )")
  
  # Retry NCBI search
  ncbi_id <- retry_with_backoff({
    search_result <- entrez_search(db = "nucleotide", term = current_accession)
    if (length(search_result$ids) == 0) stop("No IDs found")
    search_result$ids
  })
  
  if (is.null(ncbi_id) || length(ncbi_id) == 0) {
    message("No NCBI ID found for ", current_accession)
    failed_accessions <- c(failed_accessions, current_accession)
    next
  }
  
  Sys.sleep(2)
  
  # Retry entrez_link
  our_id <- retry_with_backoff({
    link_result <- entrez_link(dbfrom = "nuccore", id = ncbi_id, db = "all")
    if (is.null(link_result$links$nuccore_nuccore_rsgb)) stop("No links found")
    link_result
  })
  
  if (is.null(our_id) || length(our_id) == 0) {
    message("No links found for ", current_accession)
    failed_accessions <- c(failed_accessions, current_accession)
    next
  }
  
  Sys.sleep(2)
  
  # Retry entrez_summary
  extract_taxid <- retry_with_backoff({
    summary_result <- entrez_summary(db="nuccore", id=our_id$links$nuccore_nuccore_rsgb)
    taxid <- extract_from_esummary(summary_result, "taxid")
    if (is.null(taxid) || length(taxid) == 0) stop("No taxID found")
    taxid
  })
  
  if (is.null(extract_taxid)) {
    message("No taxID found for ", current_accession)
    failed_accessions <- c(failed_accessions, current_accession)
    next
  } else {
    message("Found taxID for ", current_accession)
  }
  
  result <- c(result, extract_taxid)
  Sys.sleep(2)
}

write.table(result, "na_result.txt", row.names = FALSE, col.names = FALSE)
if (length(failed_accessions) > 0) {
  write.table(failed_accessions, "failed_accessions.txt", row.names = FALSE, col.names = FALSE)
  message("Failed to process ", length(failed_accessions), " accessions. See failed_accessions.txt")
}
message("Successfully processed ", length(result), " accessions")

