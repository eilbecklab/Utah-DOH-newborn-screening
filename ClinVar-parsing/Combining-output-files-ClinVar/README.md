## Combining output files

These programs were written to facilitate combining the output files from the ClinVar parser. The program combine_variant_files_ClinVar.py will combine the output files from each group of variants (invalid HGVS, no star variants, variants with stars). This programs is useful for looking through all variants to determine the reasons for variants passing or failing normalization. This will result in three files being created: ClinVar_All_Invalid_HGVS_Annotations.csv, ClinVar_All_Valid_HGVS_Annotations_No_Stars.csv, and ClinVar_All_Valid_HGVS_Annotations_With_Stars.csv. Below is the help message displayed from combine_variant_files_ClinVar.py.


Example usage: ./combine_variant_files_ClinVar.py --ClinVar_directory output_files --disease_names SCID
Metabolic_Diseases --output_directory output_files

Optional Arguments:


  --ClinVar_directory CLINVAR_DIRECTORY

                        This is the name of the directory where the variants
                        from LOVD3 parser are stored. The final output file
                        will be stored in this directory if no output
                        directory is stored with the option
                        --output_directory. Default is output_files.

  --disease_names DISEASE_NAMES [DISEASE_NAMES ...]

                        This is the same list of disease names supplied into
                        ClinVar_Parser.py. If this is not specified, the
                        program will try to gather the information from the
                        directory specified with --ClinVar_directory.

  --output_directory OUTPUT_DIRECTORY [OUTPUT_DIRECTORY ...]

                        This is where your output file will be stored. If this
                        option is not used, the directory specified in
                        --ClinVar_directory will be used.
