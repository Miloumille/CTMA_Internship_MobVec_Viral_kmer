library(kmer.dist)
library(ShortRead)
library(psych) 
library(dplyr)
library(ggplot2)
library(tidyr)



virus_name <- "dengue_virus_type_1"
virus_path <- "data_processed/virus/reference/dengue_virus_type_1.fasta"
fastq_path <- "data_processed/vector/RNA_seq/"
fastq_id_list <- c("SRR23079314","SRR8482202","SRR29420355","SRR8482204")

k <- 6
x <- 100


generate_kmers <- function(k,sample_size) {
  bases <- c("A", "C", "T", "G")
  all_kmers <- as.vector(apply(expand.grid(rep(list(bases), k)), 1, paste, collapse = ""))
  kmers <- sample(all_kmers, sample_size)
  return(kmers)
}
kmers <- generate_kmers(6,40)


v_kmer_pos_df <- gen_kmers_pos_df(virus_path, kmers, virus_name)
v_kmers_stats <- get_stats_kmers(v_kmer_pos_df)

kmer_freq_list <- list()


complete_kmer_stats <- get_kmer_freq(v_kmers_stats, fastq_path, fastq_id_list)
complete_kmer_stats$y_value <- log10(complete_kmer_stats$mean_freq_virus) - log10(complete_kmer_stats$freq_m)

df_selected <- complete_kmer_stats %>%
  select(accession_id, y_value)

data_for_icc <- df_selected %>%
  group_by(accession_id) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = accession_id, values_from = y_value) %>%
  select(-row)


icc_result <- ICC(data_for_icc)
print(icc_result)

####    GRAPH     ####

complete_kmer_stats <- complete_kmer_stats[, !colnames(complete_kmer_stats) %in% "mean_count"]

ggplot(complete_kmer_stats, aes(x = factor(Kmer), y = y_value, color = accession_id)) + 
  geom_line(aes(group = Kmer), color = "gray", position = position_dodge(width = 0.2)) +  # Lines are gray
  geom_point() + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(x = "K-mer Name", y = "log10(freq virus) - log10(freq mos)", title = "ICC log10(freq virus) - log10(freq mosquito) of kmer appearing more than 5 times/virus") + 
  theme( 
    axis.text.x = element_text(angle = 90, hjust = 1), 
    plot.margin = margin(10, 30, 10, 10), 
    plot.title = element_text(size = 10)
  )



