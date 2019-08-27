## ClinVar XML File Parsing

This program will soon be updated to read in arguements from the command line using the [argparse](https://docs.python.org/3/library/argparse.html) package in python.

The first step to running this program is obtaining the xml file containing all of the variant information from ClinVar. The page where we can see the different xml files available from ClinVar is [here](ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/). The latest release will always be labeled as "ClinVarFullRelease_00-latest.xml.gz". This can be obtained using wget with a Linux operating system `wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz` or curl on a Mac operating system `curl -o ClinVarFullRelease_00-latest.xml.gz ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz`. Unzip the file using `gunzip ClinVarFullRelease_00-latest.xml.gz` and the file will be ready for parsing.
