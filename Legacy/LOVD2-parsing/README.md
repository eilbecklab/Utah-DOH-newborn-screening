## LOVD2 Variant Parsing

This program has been designed to web scrape the variant information from the Cincinnati Children's Medical Center (CCHMC) LOVD2 database for variants in specific genes relating to a disease.

Currently, the Cincinnati Children's Hospital Medical Center (CCHMC) database is the only one in LOVD2 format so no other databases are included in this parser. The LOVD2 format requires visiting the page of each variant to obtain all information about the variant.  

For information about preparing the Chromosome_name_to_accession.csv and Transcript_Info_For_Dictionaries.csv files, please see the [Obtaining-transcript-information](https://github.com/eilbecklab/Utah-DOH-newborn-screening/tree/master/Obtaining-transcript-information) folder.

Example usage: ./LOVD2_Variant_Parser.py --disease_gene_lists SCID_ny_panel.txt Metabolic_diseases_genes.txt --disease_names SCID Metabolic_Diseases --output_directory Total_Outputs --chr_to_accession Chromosome_name_to_accession.csv --transcript_info Transcript_Info_For_Dictionaries.csv

Optional Arguments:

--disease_gene_lists DISEASE_GENE_LISTS [DISEASE_GENE_LISTS ...]

                      This is a list of text files containing genes
                      associated with the disease of interest. The text
                      files should each contain a list of gene symbols for
                      one given disease, one gene per line. The name of the
                      disease will be specified by the argument
                      --disease_names.

--disease_names DISEASE_NAMES [DISEASE_NAMES ...]

                      This is a list of the disease names to accompany the
                      text files specified by the --disease_gene_lists
                      option. If you do not use this option, the file names
                      of the files specified in --disease_gene_lists
                      (without the extensions) will be used as the disease
                      names.

--output_directory OUTPUT_DIRECTORY

                      This is the name of the directory in which you want
                      your output results to be stored.

--chr_to_accession CHR_TO_ACCESSION

                      This is a csv file containing information about the
                      NCBI chromosome accession for each chromosome name. A
                      default will be generated if no file is provided. For
                      formatting, see Chromosome_name_to_accession.csv as an
                      example.

--transcript_info TRANSCRIPT_INFO

                      This is a csv file containing information about the
                      NCBI chromosome accession for each transcript and each
                      gene. This file will be used for variants where a gene
                      name or transcript accession is provided, but no
                      chromosome is provided. For formatting, see
                      Transcript_Info_For_Dictionaries.csv as an example.
