#' Calculates exclusive counts form Kraken results and adds taxonomy tables
#'
#' @param kraken_res A dataframe of Kraken results
#' @returns A dataframe containing exclusive kracken counts and taxonomy levels
#' 
#' 
kraken_exclusive_counts <- function(kraken_res){
  colnames(kraken_res) <- c("taxon", "counts")
  
  # Add depth
  kraken_res$depth <- stringr::str_count(kraken_res$taxon, "\\|")
  
  # Preprocess: split each taxon into its parent string
  kraken_res$parent <- sapply(kraken_res$taxon, function(x) {
    parts <- strsplit(x, "\\|")[[1]]
    if (length(parts) == 1) return(NA)
    paste(parts[-length(parts)], collapse = "|")
  })
  
  # Create exclusive_count column
  kraken_res$exclusive_count <- kraken_res$counts
  
  # Loop from deepest to shallowest, subtracting only **direct child** counts
  kraken_res <- kraken_res[order(-kraken_res$depth), ]  # deepest first
  for (i in 1:nrow(kraken_res)) {
    parent_taxon <- kraken_res$parent[i]
    if (!is.na(parent_taxon)) {
      parent_row <- which(kraken_res$taxon == parent_taxon)
      if (length(parent_row) == 1) {
        kraken_res$exclusive_count[parent_row] <- kraken_res$exclusive_count[parent_row] - kraken_res$counts[i]
      }
    }
  }
  
  # Sort for output
  kraken_res <- kraken_res[order(kraken_res$depth, kraken_res$taxon), ]
  kraken_res <- kraken_res |>
    dplyr::select("taxon", "counts", "exclusive_count")
  
  kraken_res <- kraken_res |>
    dplyr::mutate(domain = stringr::str_extract(taxon, "(?<=d__)[^|]+"), 
                  phylum = stringr::str_extract(taxon, "(?<=p__)[^|]+"),
                  class = stringr::str_extract(taxon, "(?<=c__)[^|]+"),
                  order = stringr::str_extract(taxon, "(?<=o__)[^|]+"),
                  family = stringr::str_extract(taxon, "(?<=f__)[^|]+"),
                  genus = stringr::str_extract(taxon, "(?<=g__)[^|]+"),
                  species = stringr::str_extract(taxon, "(?<=s__)[^|]+")) |>
    dplyr::relocate(taxon, domain, phylum, class, order, family, genus, species)
  
  return(kraken_res)
}


#' Adds Kraken Taxonomy to Bracken Results
#'
#' @param kraken_res A dataframe of Kraken results that has been cleaned through kraken_exclusive_counts
#' @param bracken_res A dataframe of Bracken results 
#' @returns A dataframe of Bracken Results with Taxonomy added from the Kraken results

add_taxonomy_bracken <- function(kraken_res, bracken_res) {
  bracken_res <- dplyr::left_join(bracken_res, kraken_res, 
                                  by=dplyr::join_by(name == species)) |>
    dplyr::select(-counts, -exclusive_count, -taxon) |>
    dplyr::relocate(domain, phylum, class, order, family, genus) |>
    dplyr::rename(species = name)
  
  return(bracken_res)
}