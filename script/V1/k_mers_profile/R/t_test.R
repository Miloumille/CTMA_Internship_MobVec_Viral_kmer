# R/spe_kmers_list.R

#' Title of the Function
#'
#' Description of what the function does.
#' @param x Description of parameter x
#' @return Description of what the function returns
#' @export

get_t_test <- function(df_virus, df_moust) {
  virus_values <- df_virus %>%
    group_by(Kmer) %>%
    summarize(log10_segm_size = list(na.omit(log10_segm_size)), .groups = 'drop')
  
  moustique_values <- df_moust %>%
    group_by(Kmer) %>%
    summarize(log10_segm_size = list(na.omit(log10_segm_size)), .groups = 'drop')
  
  results <- data.frame(Kmer = character(),
                        mean_virus = numeric(),
                        sd_virus = numeric(),
                        mean_moustique = numeric(),
                        sd_moustique = numeric(),
                        difference = numeric(),
                        t_statistic = numeric(),
                        p_value = numeric(),
                        stringsAsFactors = FALSE)
  
  for (kmer in virus_values$Kmer) {
    virus_sizes <- unlist(virus_values$log10_segm_size[virus_values$Kmer == kmer])
    moustique_sizes <- unlist(moustique_values$log10_segm_size[moustique_values$Kmer == kmer])
    
    # Perform t-test only if sufficient data exists
    if (length(virus_sizes) > 1 && length(moustique_sizes) > 1) {
      test_result <- t.test(virus_sizes, moustique_sizes)
      
      # Extract means from the t-test results
      mean_virus <- test_result$estimate[1]  # Mean of virus sizes
      mean_moustique <- test_result$estimate[2]  # Mean of moustique sizes
      sd_virus <- sd(virus_sizes, na.rm = TRUE)  # Standard deviation of virus sizes
      sd_moustique <- sd(moustique_sizes, na.rm = TRUE)  # Standard deviation of moustique sizes
      
      # Calculate the difference using the means from the t-test results
      difference <- mean_virus - mean_moustique
      
    } else {
      mean_virus <- sd_virus <- mean_moustique <- sd_moustique <- difference <- NA
      test_result <- list(statistic = NA, p.value = NA)
    }
    
    # Append results to the data frame
    results <- rbind(results, data.frame(
      Kmer = kmer,
      mean_virus = mean_virus,
      sd_virus = sd_virus,
      mean_moustique = mean_moustique,
      sd_moustique = sd_moustique,
      difference = difference,
      t_statistic = test_result$statistic,
      p_value = test_result$p.value
    ))
  }
  
  return(results)
}


