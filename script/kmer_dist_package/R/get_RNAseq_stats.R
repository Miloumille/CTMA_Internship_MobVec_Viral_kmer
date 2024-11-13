# R/spe_kmers_list.R

#' Title of the Function
#'
#' Description of what the function does.
#' @param x Description of parameter x
#' @return Description of what the function returns
#' @export

library(data.table)
library(stringr)
library(ShortRead)

count_kmer_in_fastq <- function(fq, specific_kmer) {
  kmer_count <- sum(str_count(fq, specific_kmer))
  total_size <- sum(nchar(fq))
  freq_m <- if (total_size > 0) kmer_count / total_size else 0
  return(list(kmer_count = kmer_count, freq_m = freq_m))
}

get_kmer_freq <- function(v_kmers_stats, fastq_path, fastq_id_list) {
  all_results_df <- list()
  for (id in fastq_id_list) {
    fq1 <- readFastq(paste0(fastq_path, id, "_1.fastq"))
    fq2 <- readFastq(paste0(fastq_path, id, "_2.fastq"))
    fq <- c(as.character(sread(fq1)), as.character(sread(fq2)))
    temp_results <- v_kmers_stats
    temp_results$kmer_count <- integer(length(temp_results$Kmer))
    temp_results$freq_m <- numeric(length(temp_results$Kmer))
    for (i in seq_along(temp_results$Kmer)) {
      kmer <- temp_results$Kmer[i]
      fastq_count <- count_kmer_in_fastq(fq, kmer)
      temp_results$kmer_count[i] <- fastq_count$kmer_count
      temp_results$freq_m[i] <- fastq_count$freq_m
    }
    temp_results$accession_id <- id
    all_results_df[[length(all_results_df) + 1]] <- temp_results
  }
  all_results_df <- rbindlist(all_results_df)
  return(all_results_df)
}
