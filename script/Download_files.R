library(kmer.dist)


accession_id <- "SRR23079314"
output_dir <- "data_processed/vector/RNA_seq"
download_data(accession_id, output_dir, 10000)


#### Check in Raw 

# "SRR29420341" # D9 DENV bodies A
# "SRR29420335" # D9 Mock bodies A
# "SRR29420321" # D9 DENV midgut A
# "SRR29420355" # D9 Mock midgut A

# "SRR29420340" # D9 DENV bodies B
# "SRR29420334" # D9 Mock bodies B
# "SRR29420320" # D9 DENV midgut B
# "SRR29420354" # D9 Mock midgut B


download_virus("/data/virus/Virus_list.csv","/data/virus/reference/","all")




