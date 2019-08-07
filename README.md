## Variant Parsing from ClinVar and LOVD websites

For this project we are obtaining the clinical variant information from [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) and the [Leiden Open Variation Databases (LOVD)](https://www.lovd.nl/). We are particularly interested in obtaining information about pathogenicity of variants that have previously been reported.

The ClinVar data is obtained from the XML file which can be downloaded directly from the [ClinVar site](ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/).

Multiple LOVD databases can be found, but we are specifically obtaining information from 6 separate LOVD databases.
1. [Global Variome LOVD3](https://www.lovd.nl/)
2. [Human Variome LOVD3](http://proteomics.bio21.unimelb.edu.au/lovd/status)
3. [Mitochondrial Disease Locus Specific Database (MSeqDR-LSDB) LOVD3](https://mseqdr.org/MITO/status)
4. [Brazilian Initiative on Precision Medicine (BIPMed) SNP Array LOVD3](http://bipmed.iqm.unicamp.br/snparray_hg19/docs/)
5. [Brazilian Initiative on Precision Medicine (BIPMed) Whole Exome Sequencing LOVD3](http://bipmed.iqm.unicamp.br/wes_hg19/docs/)
6. [Cincinnati Children's Hospital Medical Center (CCHMC) LOVD2](https://research.cchmc.org/LOVD2/home.php?action=switch_db)

After obtaining the variant information, all of the variants are normalized and validated through [Biocommons HGVS package](https://github.com/biocommons/hgvs). This normalization allows for comparison of variants that had been annotated in differing formats, such as VCF versus HGVS formatting.

After normalization, variants can be compared between databases to see which ones have been reported in multiple databases and which ones are unique to a given database.

For more information on the ClinVar parser, the LOVD3 parser or the LOVD2 parser, please see the directories dedicated to those programs and the accompanying documentation. 

### Dependencies

Each of the programs runs in Python and requires version 3.6 or higher (preferably [Python 3.7](https://www.python.org/downloads/)). Additionally, there are a number of dependencies. A requirements.txt file has been created for easy installation: `pip install -r requirememts.txt`
Below is a list of the dependencies along with the links to the documentation that can be used to give further information about installation if the pip installation cannot be used.
1. [pandas](https://pandas.pydata.org/)
2. [requests](https://pypi.org/project/requests/)
3. [json5](https://pypi.org/project/json5/)
4. [numpy](https://numpy.org/)
5. [HGVS from Biocommons](https://github.com/biocommons/hgvs)
