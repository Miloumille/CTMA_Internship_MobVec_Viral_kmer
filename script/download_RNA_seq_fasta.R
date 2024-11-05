library(ShortRead)



accession_id <- "SRR23079314"

output_dir <- "../data/testing_data"

fastq_dump_command <- c(
  "run",
  "--rm",
  "-v", paste0(normalizePath(output_dir), ":/output"),
  "ncbi/sra-tools",
  "fastq-dump",
  accession_id,
  "--outdir", "/output",
  "--split-files",
  "-X 50000"
)

# Run the Docker command
system2("docker", args = fastq_dump_command, stdout = TRUE, stderr = TRUE)



create_kmers <- function(X) {
  # Create all combinations of 'A', 'C', 'G', 'T' of length X
  kmers <- do.call(paste0, expand.grid(rep(list(c("A", "C", "G", "T")), X)))
  return(kmers)
}


#### Check in Raw 

count_specific_kmer_in_fastq <- function(fastq_file, specific_kmer) {
  fq <- readFastq(fastq_file)
  sequences <- as.character(sread(fq))
  kmer_count <- 0
  total_size <- 0
  for (seq in sequences) {
    total_size <- total_size + nchar(seq)  # Add the length of the current sequence
    kmer_count <- kmer_count + (gregexpr(specific_kmer, seq, fixed = TRUE)[[1]] != -1) %>% sum()
  }
  ratio <- if (total_size > 0) {
    kmer_count / total_size
  } else {
    0
  }
  cat("Total size of all sequences:", total_size, "\n")
  cat("Count of specific k-mer ('", specific_kmer, "'):", kmer_count, "\n")
  cat("Ratio of count to size:", ratio, "\n")
}


# "SRR29420341" # D9 DENV bodies A
# "SRR29420335" # D9 Mock bodies A
# "SRR29420321" # D9 DENV midgut A
# "SRR29420355" # D9 Mock midgut A

# "SRR29420340" # D9 DENV bodies B
# "SRR29420334" # D9 Mock bodies B
# "SRR29420320" # D9 DENV midgut B
# "SRR29420354" # D9 Mock midgut B

# Other Experiment

# "SRR23079314" # D7 Naive All 


count_specific_kmer_in_fastq("../data/testing_data/SRR29420341_1.fastq","AGAG")

count_specific_kmer_in_fastq("../data/testing_data/SRR29420335_1.fastq","AGAG")

count_specific_kmer_in_fastq("../data/testing_data/SRR29420321_1.fastq","AGTCTA")

count_specific_kmer_in_fastq("../data/testing_data/SRR29420355_1.fastq","AGTCTA")



