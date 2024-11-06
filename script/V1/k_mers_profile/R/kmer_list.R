# R/spe_kmers_list.R

#' Title of the Function
#'
#' Description of what the function does.
#' @param x Description of parameter x
#' @return Description of what the function returns
#' @export

gen_obs_kmers <- function(fasta_file,k){
  extract_kmers <- function(sequence, k) {
    kmers <- sapply(1:(nchar(sequence) - k + 1), function(i) {
      substring(sequence, i, i + k - 1)
    })
    return(kmers)
  }
  
  genomes <- readDNAStringSet(fasta_file)
  kmers_list <- list()
  for (i in seq_along(genomes)) {
    sequence <- as.character(genomes[i])
    kmers <- extract_kmers(sequence, k)
    kmers_list[[names(genomes)[i]]] <- kmers
  }
  common_kmers <- Reduce(intersect, kmers_list)
  return(common_kmers)
}

gen_spe_first_obs_kmers <- function(fasta_file,k,x){
  s_obs_kmers <- gen_obs_kmers(fasta_file,k)
  v_kmer_first_pos_df <- gen_kmers_first_pos_df(fasta_file,s_obs_kmers)
  kmer_max_position_df <- v_kmer_first_pos_df %>%
    group_by(SequenceName, Kmer) %>%
    summarize(Position = min(Position)) %>%
    ungroup() %>%
    group_by(Kmer) %>%
    summarize(Position = max(Position), .groups = 'drop') %>%
    slice_min(order_by = Position, n = x)
  best_kmers_list <- kmer_max_position_df$Kmer 
  return(best_kmers_list)
}


gen_gen_first_obs_kmers <- function(fasta_folder,k,x){
  
  v_fasta_files <- tools::file_path_sans_ext(list.files(path = fasta_folder ,pattern = "\\.(fasta)$",full.names = FALSE))
  g_obs_kmers <- NULL 
  for (v_fasta_file in v_fasta_files) {  
    obs_kmers <- gen_obs_kmers(paste("../data/virus/reference/",v_fasta_file,".fasta",sep = ""), k) 
    if (is.null(g_obs_kmers)) {
      g_obs_kmers <- obs_kmers
    } else {
      g_obs_kmers <- intersect(g_obs_kmers, obs_kmers)
    }
  }
  all_v_kmer_first_pos_df <- list() 
  for (v_fasta_file in v_fasta_files) {  
    v_kmer_pos_df <- gen_kmers_first_pos_df(paste("../data/virus/reference/",v_fasta_file,".fasta",sep = ""), g_obs_kmers) 
    all_v_kmer_first_pos_df[[length(all_v_kmer_first_pos_df) + 1]] <- v_kmer_pos_df
  }
  all_v_kmer_first_pos_df <- do.call(rbind, all_v_kmer_first_pos_df)
  
  kmer_max_position_df <- all_v_kmer_first_pos_df %>%
    group_by(SequenceName, Kmer) %>%
    summarize(Position = min(Position)) %>%
    ungroup() %>%
    group_by(Kmer) %>%
    summarize(Position = max(Position), .groups = 'drop') %>%
    slice_min(order_by = Position, n = x)
  
  best_kmers_list <- kmer_max_position_df$Kmer
  return(best_kmers_list)
}



