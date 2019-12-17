## Merge Transcripts

This program has been designed to merge the files that have been created for a
single gene due to multiple transcripts being present. This is specifically
designed to work with the output from LOVD3_Variant_Parser.py, which generates
an output file for each transcript present for specified genes. For input, this
program requires the config file that was used for LOVD3_Variant_Parser.py (the
json file), the list of disease names, the list of genes associated with the
given diseases, and the path to the directory containing all of the output files
from LOVD3_Variant_Parser.py. The individual transcript files will be moved to a
new folder 'Individual_Transcripts' and a new file with 'Combined_Transcripts'
in the name in place of the transcript names will be created. 

Example usage: ./merge_transcripts.py --config_file LOVD3_Databases.json
--disease_gene_lists SCID_ny_panel.txt Metabolic_diseases_genes.txt
--disease_names SCID Metabolic_Diseases --path LOVD3

optional arguments:

  -h, --help            

                        show this help message and exit

  -c CONFIG_FILE, --config_file CONFIG_FILE

                        This is the config file in JSON format that was used
                        for the input to LOVD3_Variant_Parser.py. For an
                        example of formatting, see the file
                        LOVD3_Databases.json. If this option is not used, the
                        program will attempt to find a json format file in
                        your present working directory to use.

  -g DISEASE_GENE_LISTS [DISEASE_GENE_LISTS ...], --disease_gene_lists DISEASE_GENE_LISTS [DISEASE_GENE_LISTS ...]

                        This is a list of text files containing genes
                        associated with the disease of interest. The text
                        files should each contain a list of gene symbols for
                        one given disease, one gene per line. The name of the
                        disease will be specified by the argument
                        --disease_names.

  -d DISEASE_NAMES [DISEASE_NAMES ...], --disease_names DISEASE_NAMES [DISEASE_NAMES ...]

                        This is a list of the disease names to accompany the
                        text files specified by the --disease_gene_lists
                        option. If you do not use this option, the file names
                        of the files specified in --disease_gene_lists
                        (without the extensions) will be used as the disease
                        names.

  -p PATH, --path PATH  

                        This is the path to the directory used for
                        LOVD3_Variant_Parser.py. If no path is specified, the
                        program will try to run from the current working
                        directory.
