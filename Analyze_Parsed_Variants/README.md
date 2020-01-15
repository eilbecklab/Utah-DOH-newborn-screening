## Analyze Parsed Variants

This program has been designed to gather information about the variants that
have been parsed with the parsers LOVD2_Variant_Parser.py,
LOVD3_Variant_Parser.py and ClinVar_Parser.py.

This program will count the number of variants of each type per each gene for
each database that passed validation with HGVS. Additionally, a file will be
created with all of the variants from all parsers and a file will be created
with all variants that passed validation but changed when they were validated
with HGVS.

There are many optional arguments that can be used with this parser, but most
of them are only used if you have not run all three of the parsers
(LOVD2_Variant_Parser.py, LOVD3_Variant_Parser.py, ClinVar_Parser.py) or if
you have used different disease gene lists for some of those parsers. Due to
multiple options starting with the same letter, many of the options for this
program do not have a short option.

Example usage: python Analyze_Parsed_Variants.py -c Output/ClinVar -2
Output/CCHMC -3 LOVD3 --config_file LOVD3_Databases.json --disease_gene_lists
SCID_ny_panel.txt Metabolic_diseases_genes.txt --disease_names SCID
Metabolic_Diseases --include_no_stars -p Parsed -d All_Parsed

Optional Arguments:

  -h, --help            

                        show this help message and exit

  -c CLINVAR, --clinvar CLINVAR

                        The path to the directory containing the output of
                        ClinVar_Parser.py. Note, if you specified 'Output'
                        using the --output_directory option then your
                        directory will be 'Output/ClinVar'.

  -2 LOVD2, --lovd2 LOVD2

                        The path to the directory containing the output of
                        LOVD2_Variant_Parser.py. Note, if you specified
                        'Output' using the --output_directory option then your
                        directory will be 'Output/CCHMC'.

  -3 LOVD3, --lovd3 LOVD3

                        The path to the directory containing the output of
                        LOVD3_Variant_Parser.py. Unlike the options -c and -2,
                        the directory used for this option is the one
                        specified in LOVD3_Variant_Parser.py. This is because
                        multiple databases are present for LOVD3. If you
                        specified 'Output' using the --output_directory option
                        then your directory will be 'Output'. Individual
                        databases can be determined by list of subdirectories
                        in the directory provided or they can be directly
                        stated using either the --config_file or -databases
                        options.

  --disease_gene_lists DISEASE_GENE_LISTS [DISEASE_GENE_LISTS ...]

                        Only use this option if you used the same gene lists
                        for all parsers. If you have used different lists for
                        the individual parsers, some of the genes may be
                        overlooked in the counting of variants of each type.
                        This is a list of text files containing genes
                        associated with the disease of interest. The text
                        files should each contain a list of gene symbols for
                        one given disease, one gene per line. The name of the
                        disease must be specified by the argument
                        --disease_names.

  --disease_names DISEASE_NAMES [DISEASE_NAMES ...]

                        This is a list of the disease names used for the
                        parsers. If neither this option nor the options for
                        specifying diseases used for individual parsers, then
                        the subdirectory names from the directories specified
                        with the -c, -2, and --config_file options will be
                        used. If more than one parser was used and different
                        disease names were used for the separate parsers, you
                        can specify the disease names for each parser with the
                        options --disease_names_clinvar,
                        --disease_names_lovd2, and --disease_names_lovd3.
                        Example: --disease_names SCID Metabolic_Diseases

  -p OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX

                        This is the prefix that you would like to use for the
                        output files. Default is Parsed_Variants

  -d OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY

                        This is the directory where you would like to store
                        the output files. Default is the current working
                        directory.

  --save_dictionary     

                        If this option is specified, then the program will
                        save the dictionary that contains the information
                        about the parsers, databases, and diseases specified
                        in the input along with the path to each of the files
                        associated with the input information. It will be
                        saved in a JSON file. By default, the program will not
                        save the dictionary.

The following two options are related to specifying which databases were used
for variant parsing with LOVD3_Variant_Parser.py.

  --databases DATABASES [DATABASES ...]

                        This is a list of the database names used in the
                        config file for LOVD3_Variant_Parser.py. This option
                        should not be used in conjunction with the
                        --config_file option. If neither option is used, all
                        subdirectories in the directory specified with the
                        --lovd3 option will be used. An example input for this
                        option would be -d BIPmed_SNPhg19 BIPmed_WES
                        Global_Variome Human_Variome MSeqDR-LSDB

  --config_file CONFIG_FILE

                        This is the config file in JSON format that was used
                        for the input to LOVD3_Variant_Parser.py. For an
                        example of formatting, see the file
                        LOVD3_Databases.json. This option should not be used
                        in conjunction with the --databases option. If neither
                        option is used, all subdirectories in the directory
                        specified with the --lovd3 option will be used.


The following options are for use when you have run differing gene lists for the
separate parsers but would like to analyze the variants together.

  --disease_names_clinvar DISEASE_NAMES_CLINVAR [DISEASE_NAMES_CLINVAR ...]

                        This is a list of the disease names used for
                        ClinVar_Parser.py. Use this option only if you have
                        used a different list of disease names for
                        ClinVar_Parser.py than were used for the other
                        parsers. If neither this option nor the
                        --disease_names option is used, the program will use
                        the subdirectories in the directory specified by the
                        --clinvar option as disease names. Example:
                        --disease_names_clinvar SCID Metabolic_Diseases

  --disease_names_lovd2 DISEASE_NAMES_LOVD2 [DISEASE_NAMES_LOVD2 ...]

                        This is a list of the disease names used for
                        LOVD2_Variant_Parser.py. Use this option only if you
                        have used a different list of disease names for
                        LOVD2_Variant_Parser.py than were used for the other
                        parsers. If neither this option nor the
                        --disease_names option is used, the program will use
                        the subdirectories in the directory specified by the
                        --lovd2 option as disease names. Example:
                        --disease_names_lovd2 SCID Metabolic_Diseases

  --disease_names_lovd3 DISEASE_NAMES_LOVD3 [DISEASE_NAMES_LOVD3 ...]

                        This is a list of the disease names used for
                        LOVD3_Variant_Parser.py. Use this option only if you
                        have used a different list of disease names for
                        LOVD3_Variant_Parser.py than were used for the other
                        parsers. If neither this option nor the
                        --disease_names option is used, the program will use
                        the subdirectories in the directory specified by the
                        --lovd2 option as disease names. Example:
                        --disease_names_lovd3 SCID Metabolic_Diseases


  --use_file_names      

                        Specify this option if you have used the same gene
                        lists for all parsers and would like to use the names
                        of the text files (without the final extension) for
                        the disease names. This is recommended only if you
                        allowed the parsers to determine the disease names based
                        upon the lists provided.

The following options are for telling the program whether to include the
variants that received no stars by ClinVar's review status. By default, only
variants that had at least one star by ClinVar's review status will be included.

  --include_no_stars    

                        Specify this option if you would like to include the
                        no-star variants from ClinVar_Parser.py. By default,
                        this program will ignore them.

  --exclude_stars       

                        Specify this option if you would like to include only
                        the no-star variants from ClinVar_Parser.py, but
                        exclude the variants that had 1 or more stars. By
                        default, this program will include only variants that
                        have 1 or more stars according to ClinVar.
