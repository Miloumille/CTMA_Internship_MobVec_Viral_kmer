library(Biostrings)
library(stringr)
library(dplyr)
library(ggplot2)
library(purrr)
library(kmer.dist)

v_fasta_folder <- "../data/virus/all_dengue/"

best_kmers_list <- gen_gen_first_obs_kmers(v_fasta_folder,k,x)
comb_v_kmer_pos_df <- list() 
v_fasta_files <- tools::file_path_sans_ext(list.files(path = v_fasta_folder ,pattern = "\\.(fasta)$",full.names = FALSE))
for (v_fasta_file in v_fasta_files) {    
  v_kmer_pos_df <- gen_kmers_pos_df(paste(v_fasta_folder,v_fasta_file,".fasta",sep = ""), best_kmers_list,v_fasta_file)
  comb_v_kmer_pos_df[[length(comb_v_kmer_pos_df) + 1]] <- v_kmer_pos_df
}
comb_v_kmer_pos_df <- do.call(rbind, comb_v_kmer_pos_df)