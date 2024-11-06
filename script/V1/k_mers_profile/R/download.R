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

download_vector <- function(path_csv_file,path_output){
  Vector <- read.csv(path_csv_file)
  for(assembly_code in Vector$NCBI.RefSeq.assembly){
    getGenome(db = "refseq",organism = assembly_code,path = path_output,reference = FALSE)
  }
}

