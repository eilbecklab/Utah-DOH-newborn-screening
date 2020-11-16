from argparse import ArgumentParser, FileType
import json5
import os
import glob
import pandas as pd
import copy


args = ArgumentParser('''./combine_variant_files_ClinVar.py', description='This program is a simple one to combine the output
	files of each group of variants (invalid HGVS, no star variants, variants with stars) into a single csv file. Example usage:
	./combine_variant_files_ClinVar.py --ClinVar_directory output_files --disease_names SCID Metabolic_Diseases --output_directory output_files''')

args.add_argument(
	'-c',
	'--ClinVar_directory',
	help="""This is the name of the directory where the variants from the ClinVar parser are stored. The final output file will be stored in this
	directory if no output directory is stored with the option --output_directory. Default is output_files.""",
	default = "output_files"
)

args.add_argument(
	'-d',
	'--disease_names',
	nargs='+',
	help="""This is the same list of disease names supplied into ClinVar_Parser.py. If this is not specified, the program will try to gather
	the information from the directory specified with --ClinVar_directory. """,
	default = None
)

args.add_argument(
	'-o',
	'--output_directory',
	help="This is where your output file will be stored. If this option is not used, the directory specified in --ClinVar_directory will be used.",
	default = None
)


args = args.parse_args()
directory = args.ClinVar_directory
output_directory = args.output_directory
if not output_directory :
	# If no output directory is specified, save the output in the same directory as the input
	output_directory = directory

disease_names = args.disease_names
if not disease_names:
	# If the disease names are not specified at the command line, you can find them from the subdirectories in the output directory
	disease_names = []
	subdirectories = glob.glob(directory+'/ClinVar/*/') # list all subdirectories in the first database directory
	for subdir in subdirectories:
		disease_names.append(subdir.rsplit('/')[-2])

columns_valid = ['Genome Assembly', 'Chr', 'Position Start', 'Position Stop', 'Ref', 'Alt',
			'Genomic Annotation', 'HGVS Normalized Genomic Annotation', 'Variant Type',
			'Variant Length', 'Pathogenicity', 'Disease', 'Genetic Origin', 'Inheritance Pattern',
			'Affected Genes', 'Gene Symbol', 'dbSNP ID', 'Compound Het Status', 'Transcript',
			'Transcript Notation', 'HGVS Transcript Notation', 'Protein Accession',
			'Protein Notation', 'HGVS Protein Annotation', 'Chr Accession', 'VCF Pos', 'VCF Ref',
			'VCF Alt', 'Database', 'Database Accession', 'Review Status', 'Star Level',
			'Submitter', 'Edited Date', 'Overall_MAF', 'Control_MAF', 'African_MAF',
			'NonFinish_Euro_MAF', 'Finnish_MAF', 'Ashkenazi_MAF','Latino_MAF', 'East_Asian_MAF',
			'South_Asian_MAF', 'Other_Ancestry_MAF', 'Transcript Normalization Failure Message',
			'Genomic Normalization Failure Message']
columns_Invalid = copy.deepcopy(columns_valid)
columns_Invalid.append('HGVS Normalization Failure Reason')

merged_invalid = pd.DataFrame(columns = columns_Invalid)

for disease_name in disease_names:
	base_directory = directory+"/ClinVar/"+disease_name
	files = [f for f in glob.glob(base_directory+"/Invalid_Annotations/*.csv")]
	for file in files:
		tmp_df = pd.read_csv(file, index_col = 0)
		merged_invalid = pd.concat([merged_invalid, tmp_df]).reset_index(drop=True)
merged_invalid.to_csv(output_directory+"/ClinVar_All_Invalid_HGVS_Annotations.csv")

merged_stars = pd.DataFrame(columns = columns_valid)

for disease_name in disease_names:
	base_directory = directory+"/ClinVar/"+disease_name
	files = [f for f in glob.glob(base_directory+"/*ClinVar_Results.csv")]
	for file in files:
		tmp_df = pd.read_csv(file, index_col = 0)
		merged_stars = pd.concat([merged_stars, tmp_df]).reset_index(drop=True)
merged_stars.to_csv(output_directory+"/ClinVar_All_Valid_HGVS_Annotations_With_Stars.csv")

merged_No_stars = pd.DataFrame(columns = columns_valid)

for disease_name in disease_names:
	base_directory = directory+"/ClinVar/"+disease_name
	files = [f for f in glob.glob(base_directory+"/*No_Star_Results.csv")]
	for file in files:
		tmp_df = pd.read_csv(file, index_col = 0)
		merged_No_stars = pd.concat([merged_No_stars, tmp_df]).reset_index(drop=True)
merged_No_stars.to_csv(output_directory+"/ClinVar_All_Valid_HGVS_Annotations_No_Stars.csv")
