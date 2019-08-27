## ClinVar XML File Parsing

This program has been designed to parse the XML file containing the entire ClinVar database to obtain variant information for specific genes relating to a disease.

The first step to running this program is obtaining the xml file containing all of the variant information from ClinVar. The ClinVar data files are stored on the NCBI server ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/. The latest release will always be labeled as "ClinVarFullRelease_00-latest.xml.gz".

This file can be obtained using the following commands:

On a Linux operating system: `wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz`

Or on MacOS: `curl -o ClinVarFullRelease_00-latest.xml.gz ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz`.

Unzip the file using `gunzip ClinVarFullRelease_00-latest.xml.gz` and the file will be ready for parsing.

Example usage: ./ClinVar_Parser.py --xml_file ClinVarFullRelease_00-latest.xml
--disease_gene_lists SCID_ny_panel.txt Metabolic_diseases_genes.txt
--disease_names SCID Metabolic_Diseases --output_directory Total_Outputs

Optional Arguements:

--xml_file XML_FILE   

                        This is the ClinVar xml file that can be obtained form
                        the NCBI database
                        (ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/). Default
                        file is ClinVarFullRelease_00-latest.xml

  --disease_gene_lists DISEASE_GENE_LISTS [DISEASE_GENE_LISTS ...]

                        This is a list of text files containing genes
                        associated with the disease of interest. The text
                        files should each contain a list of gene symbols for
                        one given disease, one gene per line. The name of the
                        disease will be specified by the arguement
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
