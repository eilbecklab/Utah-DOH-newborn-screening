## Generating the transcript tables for use with the LOVD parsers

This is the program that is used to parse transcript information from GFF files obtained from the NCBI websites. GFF files were chosen specifically because they contained all of the information necessary, while individual GTF files from various sources contained only some of the necessary information. To browse available GFF files from the NCBI, you can visit the [NCBI ftp directories](ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/).

For the programs included in this repository, GRCh37 was chosen as the reference because the majority of variant annotations contained GRCh37 coordinates, but many did not contain GRCh38 coordinates. To obtain the file used here, you can use `wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz` in a Linux kernel or use `curl -o GRCh37_latest_genomic.gff.gz ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz` in a Mac environment. Unzip the file with `gunzip GRCh37_latest_genomic.gff.gz` to prepare it for input to the python program.

The output files that will be produced from the GFF file downloaded on August 14, 2019 have been included with this repository. 

Example usage: ./LOVD3_Variant_Parser.py --gff_file GRCh37_latest_genomic.gff --output_transcripts_file Transcript_Info_For_Dictionaries.csv --output_chromosome_file Chromosome_name_to_accession.csv

Optional Arguments:

  --gff_file GFF_FILE   

                        This is a GFF file containing the information about
                        the transcripts of the genome. Default is
                        GRCh37_latest_genomic.gff.

  --output_transcripts_file OUTPUT_TRANSCRIPTS_FILE

                        This is the name of the directory in which you want
                        your output results to be stored. Default is
                        ranscript_Info_For_Dictionaries.csv

  --output_chromosome_file OUTPUT_CHROMOSOME_FILE

                        This is the name of the directory in which you want
                        your output results to be stored. Default is
                        Chromosome_name_to_accession.csv
