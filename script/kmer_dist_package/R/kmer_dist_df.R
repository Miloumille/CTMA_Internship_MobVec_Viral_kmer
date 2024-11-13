# R/spe_kmers_list.R

#' Title of the Function
#'
#' Description of what the function does.
#' @param x Description of parameter x
#' @return Description of what the function returns
#' @export


gen_kmers_pos_df <- function(fasta_file, kmers_list, virus_type) {  
  find_kmer_positions <- function(sequence, kmer) {  
    positions <- gregexpr(kmer, sequence)[[1]]  
    positions <- positions[positions > 0]  
    return(positions)  
  }  
  
  genomes <- readDNAStringSet(fasta_file)  
  kmers_df <- data.frame(SequenceName = character(), Kmer = character(), Position = integer(), Type = character())  
  
  for (i in seq_along(genomes)) {  
    sequence_name <- names(genomes)[i]  
    sequence <- as.character(genomes[i])  
    sequence_length <- nchar(sequence)  
    
    for (kmer in kmers_list) {  
      positions <- find_kmer_positions(sequence, kmer)  
      if (length(positions) > 0) {  
        temp_df <- data.frame(SequenceName = sequence_name, Kmer = kmer, Position = positions, Type = virus_type)  
        kmers_df <- rbind(kmers_df, temp_df)  
      }  
      final_position_df <- data.frame(SequenceName = sequence_name, Kmer = kmer, Position = sequence_length, Type = virus_type)  
      kmers_df <- rbind(kmers_df, final_position_df)  
    }  
  }  
  kmers_df <- kmers_df %>%  
    filter(Kmer %in% kmers_list) %>%  
    arrange(SequenceName, Kmer, Position) %>%  
    group_by(SequenceName, Kmer) %>%  
    mutate(  
      segm_size = Position - lag(Position),  
      log10_segm_size = ifelse(is.na(segm_size), NA, log10(segm_size + 1))
    ) %>%  
    ungroup()  
  
  return(kmers_df)  
}
