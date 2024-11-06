library(Biostrings)
library(stringr)
library(dplyr)
library(ggplot2)
library(purrr)
library(kmer.dist)
library(ShortRead)

virus_name <- "dengue_virus_type_1"
k <- 6
x <- 100
size <- "1000000"

v_fasta_file <- paste("data_processed//virus/reference/", virus_name, ".fasta", sep = "")

best_kmers_list <- gen_spe_first_obs_kmers(v_fasta_file,k,x)

v_kmer_pos_df <- gen_kmers_pos_df(v_fasta_file,best_kmers_list,virus_name)

# ---
get_stats_segments <- function(df_virus) {  
  # Filter out rows with NA in segm_size
  df_virus <- df_virus %>% 
    filter(!is.na(segm_size)) 
  
  # Group by SequenceName and Kmer to calculate mean and sd on the original segm_size
  virus_values <- df_virus %>%
    group_by(SequenceName, Kmer) %>%
    summarize(mean_segm_size = mean(segm_size, na.rm = TRUE),
              sd_segm_size = ifelse(n() > 1, sd(segm_size, na.rm = TRUE), NA),
              count = n(),
              .groups = 'drop')
  
  # Now group by Kmer and calculate the mean and sd of stats across different SequenceName values
  results <- virus_values %>%
    group_by(Kmer) %>%
    summarize(mean_virus = log10(mean(mean_segm_size, na.rm = TRUE)),  # Take log10 of the mean of segm_size
              sd_virus = sd(mean_segm_size, na.rm = TRUE),             # sd remains in original scale
              mean_count = mean(count, na.rm = TRUE),
              sd_count = ifelse(n() > 1, sd(count, na.rm = TRUE), NA),
              .groups = 'drop')
  
  return(results)
}





results <- get_stats_segments(v_kmer_pos_df)


count_specific_kmer_in_fastq <- function(fq, specific_kmer) { 
  sequences <- as.character(sread(fq)) 
  kmer_count <- 0 
  total_size <- 0 
  for (seq in sequences) { 
    total_size <- total_size + nchar(seq)
    kmer_count <- kmer_count + (gregexpr(specific_kmer, seq, fixed = TRUE)[[1]] != -1) %>% sum() 
  } 
  ratio <- if (total_size > 0) { 
    kmer_count / total_size 
  } else { 
    0 
  } 
  return(list(kmer_count = kmer_count, ratio = ratio))
}


calculate_kmer_frequencies <- function(results, fq) {
  kmers <- results$Kmer
  
  results$kmer_count <- numeric(length(kmers))
  results$ratio <- numeric(length(kmers))
  
  for (i in seq_along(kmers)) {
    kmer <- kmers[i]
    result <- count_specific_kmer_in_fastq(fq, kmer)  
    results$kmer_count[i] <- result$kmer_count     
    results$ratio[i] <- result$ratio                 
  }
  
  return(results)
}

fq <- readFastq("data_processed//testing_data/SRR23079314_1.fastq")
results <- calculate_kmer_frequencies(results, fq)

results$frequency_virus <- 1 / (10 ^ results$mean_virus)

library(ggplot2)

# Create new columns for the X and Y axes
results$x_value <- ifelse(results$sd_virus == 0,
                          sqrt(1 / min(results$sd_virus[results$sd_virus != 0]) + 0.01),
                          sqrt(1 / results$sd_virus))
results$y_value <- results$frequency_virus / results$ratio  # Y-axis values

# This code adds a dark red horizontal line at y = 1.
ggplot(results, aes(x = x_value, y = y_value)) +
  geom_text(aes(label = paste0(Kmer, "\n", mean_count)), 
            vjust = -0.5, size = 2) +
  labs(title = "frequency_virus / frequency_mos",
       x = "sqrt(1/sd_virus)",
       y = "frequency virus / frequency mos") +
  geom_hline(yintercept = 1, color = "darkred", linetype = "solid") +  # Dark red horizontal line at y = 1
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 1)) +
  expand_limits(x = max(results$x_value) * 1.05, 
                y = max(results$y_value) * 1.1)

