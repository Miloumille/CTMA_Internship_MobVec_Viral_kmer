####### GENERATE REPORT  ########
library(KmerEnrich)

virus_genome_path <- "data_processed/virus/dengue/dengue_type_2/"
vector_genome_path <- c("data_processed/vector/aedes/aegypty/SRR14655772",
                        "data_processed/vector/aedes/albopictus/SRR8482204")
k <- c(5)
x <- 500

KmerEnrichFullReport(virus_genome_path, vector_genome_path, k, x,"dengue_type_2")

