This repository contains the different codes and pipelines used to obtain the results available in my thesis.

The repository is organized into several parts:

* **data\_processed** contains the various FASTQ and FASTA files corresponding to the reference viral genomes of dengue type 2 and the RNAseq files used for *Aedes albopictus*. This folder also contains a subfolder, **ref\_vector**, which holds the reference genome of *Aedes albopictus* used for mapping the reads.

* **data\_sequencing** contains the files generated and collected from the different sequencing runs. Due to their large size, these files are not available here. If needed, you can contact me via email at [milouovancau@gmail.com](mailto:milouovancau@gmail.com).

Only the folder for run 1 out of the 5 runs is available as a representation of the folder organization.

* **reports** contains the results obtained after running the different pipelines. It includes a section for the ICC results and another folder dedicated to reports generated with the KMerEnrich package, available on GitHub at https://github.com/Miloumille/KmerEnrich.

* **scripts** contains the various R scripts used to obtain the results presented in the thesis:

  * ICC analysis
  * Specific and generic pipelines to use the KMerEnrich package for pool or virus-specific purposes
  * The download input script, which allows downloading files directly from NCBI using rentrez
  * The script to generate the complete report from KMerEnrich
  * The script for using Kraken2 and Minimap2, which both associate reads with a reference genome.

More documentation about KMerEnrich is available on its GitHub page.
