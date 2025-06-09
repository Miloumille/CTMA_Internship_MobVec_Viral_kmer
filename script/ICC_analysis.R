library(KmerEnrich)
library(dplyr)
library(psych)
library(ggplot2)

######### ICC #########

virus_path <- "data_processed/virus/dengue/dengue_type_1/NC_001477.1.fasta"
# Comparison conditions
fastq_path_list <- c("data_processed/vector/aedes/aegypty/SRR29420386",
                     "data_processed/vector/aedes/aegypty/SRR29420375",
                     "data_processed/vector/aedes/aegypty/SRR29420334",
                     "data_processed/vector/aedes/aegypty/SRR29420368",
                     "data_processed/vector/aedes/aegypty/SRR29420325")

# Comparison experiments
fastq_path_list <- c("data_processed/vector/aedes/aegypty/SRR29420334",
                     "data_processed/vector/aedes/aegypty/SRR7499952",
                     "data_processed/vector/aedes/aegypty/ERR15075903",
                     "data_processed/vector/aedes/aegypty/SRR23079314",
                     "data_processed/vector/aedes/aegypty/SRR17786679",
                     "data_processed/vector/aedes/aegypty/SRR14655772")

k <- 7




generate_kmers <- function(k) {
  bases <- c("A", "C", "T", "G")
  all_kmers <- as.vector(apply(expand.grid(rep(list(bases), k)), 1, paste, collapse = ""))
  return(all_kmers)
}

# Example usage
kmers <- generate_kmers(k)


v_kmer_pos_df <- kmers_pos_df(virus_path, kmers)
v_kmers_stats <- get_virus_stats_kmers(v_kmer_pos_df)


complete_kmer_stats <- get_vector_stats_kmers(v_kmers_stats, fastq_path_list)

df_selected <- complete_kmer_stats %>%
  select(accession_id, Enrichment_score)

data_for_icc <- df_selected %>%
  group_by(accession_id) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = accession_id, values_from = Enrichment_score) %>%
  select(-row)


icc_result <- ICC(data_for_icc)
print(icc_result)


########  GRAPH ICC  ######

random_kmer_sample <-  sample(kmers, 100)
random_kmer_stats <- complete_kmer_stats %>% filter(Kmer %in% random_kmer_sample)

ggplot(random_kmer_stats, aes(x = factor(Kmer), y = Enrichment_score, color = accession_id)) + 
  geom_line(aes(group = Kmer), color = "gray", position = position_dodge(width = 0.2)) +  # Lines are gray
  geom_point() + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(x = "Randomly selected K-mers", y = "freq virus / freq mos", title = "ICC freq virus / freq mosquito of kmer") + 
  theme( 
    axis.text.x = element_text(angle = 90, hjust = 1), 
    plot.margin = margin(10, 30, 10, 10), 
    plot.title = element_text(size = 10)
  )




