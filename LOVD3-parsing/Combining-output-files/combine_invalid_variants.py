from argparse import ArgumentParser, FileType
import json5
import os
import glob
import pandas as pd


args = ArgumentParser('./combine_invalid_variants.py', description='This program is a simple one to combine the output files containing the invalid variants (anything that failed HGVS normalization) into one csv file. Example usage: ./combine_invalid_variants.py --config_file LOVD3_Databases.json --LOVD3_directory output_files --disease_names SCID Metabolic_Diseases --output_directory output_files')

args.add_argument(
	'--config_file',
	#type=FileType('r'),
	help="This is the config file that was used for the LOVD3 parser. If this option is not provided, the program will look for a json file in your present working directory.",
	default = None,
)

args.add_argument(
	'--LOVD3_directory',
	help="This is the name of the directory where the variants from LOVD3 parser are stored. The final output file will be stored in this directory if no output directory is stored with the option --output_directory. Default is output_files.",
	default = "output_files"
)

args.add_argument(
	'--disease_names',
	nargs='+',
	help="This is the same list of disease names supplied into LOVD3_Variant_Parser.py. If this is not specified, the program will try to gather the information from the directory specified with --LOVD3_directory. ",
	default = None
)

args.add_argument(
	'--output_directory',
	nargs='+',
	help="This is where your output file will be stored. If this option is not used, the directory specified in --LOVD3_directory will be used.",
	default = None
)


args = args.parse_args()
disease_names = args.disease_names
json_file = args.config_file
if not json_file:
	# If no input config file is present, look for one in the present working directory
	jsons_list = glob.glob("*.json")
	if len(jsons_list) == 1:
		json_file = jsons_list[0]
	try:
		if len(jsons_list) == 0:
			raise ValueError('No JSON file present to use as config file')
		if len(jsons_list) > 1:
			raise ValueError('More than one JSON file present and no file specified as config file. Please specify your config file with --config_file.')
	except ValueError as error:
		print(repr(error))
		raise SystemExit()

with open(json_file) as file:
	databases_list = json5.load(file)

directory = args.LOVD3_directory
output_directory = args.output_directory
if not output_directory :
	# If no output directory is specified, save the output in the same directory as the input
	output_directory = directory

if not disease_names:
	# If the disease names are not specified at the command line, you can find them from the subdirectories in the output directory
	disease_names = []
	first_database = databases_list[0]["name"]
	subdirectories = glob.glob(directory+"/"+first_database+'/*/') # list all subdirectories in the first database directory
	for subdir in subdirectories:
		disease_names.append(subdir.rsplit('/')[-2])

columns = ['Genome Assembly', 'Chr', 'Position Start', 'Position Stop', 'Ref', 'Alt',
			'Genomic Annotation', 'HGVS Normalized Genomic Annotation', 'Variant Type',
			'Variant Length', 'Pathogenicity', 'Disease', 'Genetic Origin', 'Inheritance Pattern',
			'Affected Genes', 'Gene Symbol', 'Compound Het Status', 'Transcript',
			'Transcript Notation', 'HGVS Transcript Notation', 'Protein Accession',
			'Protein Notation', 'HGVS Protein Annotation', 'Chr Accession', 'VCF Pos', 'VCF Ref',
			'VCF Alt', 'Database', 'ClinVar Accession', 'Review Status', 'Star Level',
			'Submitter', 'Edited Date', 'DNA change (genomic) (hg19)', 'Effect',
			'Exon','Reported', 'DB-ID', 'dbSNP ID', 'Published as', 'Variant remarks',
			'Reference', 'Frequency', 'Transcript Normalization Failure Message',
			'Genomic Normalization Failure Message', 'HGVS Normalization Failure Reason']

merged_df = pd.DataFrame(columns = columns)
for database in databases_list:
	database_name = database["name"]
	for disease_name in disease_names:
		base_directory = directory+"/"+database_name+"/"+disease_name
		files = [f for f in glob.glob(base_directory+"/Invalid_Annotations/*.csv")]
		for file in files:
			tmp_df = pd.read_csv(file, index_col = 0)
			merged_df = pd.concat([merged_df, tmp_df]).reset_index(drop=True)
merged_df.to_csv(output_directory+"/LOVD3_All_Invalid_HGVS_Annotations.csv")
