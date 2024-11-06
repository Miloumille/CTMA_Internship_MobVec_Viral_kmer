# R/spe_kmers_list.R

#' Title of the Function
#'
#' Description of what the function does.
#' @param x Description of parameter x
#' @return Description of what the function returns
#' @export

library(biomartr)

download_virus <- function(path_csv_file,path_output,virus_name = "all"){
  Virus_csv <- read.csv(path_csv_file)
  if (virus_name == "all"){
    for (name in unique(Virus_csv$name)){
      print(paste("downloading",name))
      ids <- Virus_csv$accession_code[Virus_csv$name == name]
      combined_sequences <- entrez_fetch(db = "nuccore", id = ids, rettype = "fasta", retmode = "text")
      write(combined_sequences, file = paste(path_output,name,".fasta",sep = ""))
    }
  }
  else{
    ids <- Virus_csv$accession_code[Virus_csv$name == virus_name]
    combined_sequences <- entrez_fetch(db = "nuccore", id = ids, rettype = "fasta", retmode = "text")
    write(combined_sequences, file = paste(path_output,virus_name,".fasta",sep = ""))
  }
}


download_data <- function(accession_id, output_dir, n_reads) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  fastq_dump_command <- c(
    "run",
    "--rm",
    "-v", paste0(normalizePath(output_dir), ":/output"),
    "ncbi/sra-tools",
    "fastq-dump",
    accession_id,
    "--outdir", "/output",
    "--split-files",
    paste0("-X ", n_reads)
  )
  
  # Run the Docker command
  system2("docker", args = fastq_dump_command, stdout = TRUE, stderr = TRUE)
}