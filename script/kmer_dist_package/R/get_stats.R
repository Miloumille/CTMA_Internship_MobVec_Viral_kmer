# R/spe_kmers_list.R

#' Title of the Function
#'
#' Description of what the function does.
#' @param x Description of parameter x
#' @return Description of what the function returns
#' @export

get_stats_kmers <- function(df_virus) {  
  # Filter out rows with NA in segm_size
  df_virus <- df_virus %>% 
    filter(!is.na(segm_size)) 
  
  # Group by SequenceName and Kmer to calculate mean and sd on the original segm_size
  virus_values <- df_virus %>%
    group_by(SequenceName, Kmer) %>%
    summarize(mean_segm_size = mean(log10(segm_size), na.rm = TRUE),
              sd_segm_size = ifelse(n() > 1, sd(log10(segm_size), na.rm = TRUE), 0),
              count = n(),
              freq = n()/max(Position),
              .groups = 'drop')
  
  # Now group by Kmer and calculate the mean and sd of stats across different SequenceName values
  results <- virus_values %>%
    group_by(Kmer) %>%
    summarize(
      mean_virus = mean(mean_segm_size, na.rm = TRUE), 
      mean_sd_virus = mean(sd_segm_size, na.rm = TRUE), 
      mean_count = mean(count, na.rm = TRUE),
      sd_count = ifelse(n() > 1, sd(count, na.rm = TRUE), NA), 
      mean_freq_virus = mean(freq, na.rm = TRUE), 
      .groups = 'drop'
    )
  return(results)
}


