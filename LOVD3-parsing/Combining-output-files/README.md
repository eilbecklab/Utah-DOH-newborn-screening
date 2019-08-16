## Combining output files

These programs were written to facilitate combining the output files from the LOVD3 parser. The program combine_invalid_variants.py will combine the output files from all of the directories labeled "Invalid_Annotations" and the program combine_passed_variants.py will combine all of the output files containing variants that passed [biocommons HGVS normalization]((https://github.com/biocommons/hgvs)). These programs are useful for looking through all variants to determine the reasons for variants passing or failing normalization. Both programs work the same way, with the exception of which input files they take and the name of the output file. The output file from combine_invalid_variants.py is LOVD3_All_Invalid_HGVS_Annotations.csv and the output from combine_passed_variants.py is LOVD3_All_Valid_HGVS_Annotations.csv. Below is the help message displayed from combine_passed_variants.py.


Example usage: ./combine_invalid_variants.py --config_file
LOVD3_Databases.json --LOVD3_directory output_files --disease_names SCID
Metabolic_Diseases --output_directory output_files

Optional Arguments:


  --config_file CONFIG_FILE

                        This is the config file that was used for the LOVD3
                        parser. If this option is not provided, the program
                        will look for a json file in your present working
                        directory.

  --LOVD3_directory LOVD3_DIRECTORY

                        This is the name of the directory where the variants
                        from LOVD3 parser are stored. The final output file
                        will be stored in this directory if no output
                        directory is stored with the option
                        --output_directory. Default is output_files.

  --disease_names DISEASE_NAMES [DISEASE_NAMES ...]

                        This is the same list of disease names supplied into
                        LOVD3_Variant_Parser.py. If this is not specified, the
                        program will try to gather the information from the
                        directory specified with --LOVD3_directory.

  --output_directory OUTPUT_DIRECTORY [OUTPUT_DIRECTORY ...]

                        This is where your output file will be stored. If this
                        option is not used, the directory specified in
                        --LOVD3_directory will be used.
