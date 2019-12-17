#!/usr/bin/env python
# coding: utf-8

import os
from argparse import ArgumentParser, FileType
import json5
import glob
from os.path import splitext
import sys
import pandas as pd


args = ArgumentParser('./merge_transcripts.py', description="""This program has been designed to merge the
files that have been created for a single gene due to multiple transcripts being present. For input, this
program requires the config file that was used for LOVD3_Variant_Parser.py, the list of disease names
and the list of genes associated with the given diseases.
Example usage: ./merge_transcripts.py --config_file LOVD3_Databases.json
--disease_gene_lists SCID_ny_panel.txt Metabolic_diseases_genes.txt --disease_names
SCID Metabolic_Diseases --path LOVD3""")

args.add_argument(
	'-c',
	'--config_file',
	help="""This is the config file in JSON format that was used for the input to
	LOVD3_Variant_Parser.py. For an example of formatting, see the file LOVD3_Databases.json.
	If this option is not used, the program will attempt to find a json format file in your
	present working directory to use.""",
	default=None,
)

args.add_argument(
	'-g',
	'--disease_gene_lists',
	nargs='+', # This tells the program that if they specify this flag, they have to give it at least one input. If they don't specify it, then the default will go in.
	help="""\
This is a list of text files containing genes associated with the disease of interest. The text files
should each contain a list of gene symbols for one given disease, one gene per line. The name of the disease will be specified
by the arguement --disease_names. """,
	default = None
)

args.add_argument(
	'-d',
	'--disease_names',
	nargs='+',
	help="""This is a list of the disease names to accompany the text files specified by the --disease_gene_lists
	option. If you do not use this option, the file names of the files specified in --disease_gene_lists
	(without the extensions) will be used as the disease names.""",
	default = None
)

args.add_argument(
	'-p',
	'--path',
	help="""This is the path to the directory used for LOVD3_Variant_Parser.py. If no path is
	specified, the program will try to run from the current working directory. """,
	default = "."
)

def print_help_message():
	print()
	print("""\tThis program has been designed to merge the files that have been created for a single
	gene due to multiple transcripts being present. For input, this program requires the
	config file that was used for LOVD3_Variant_Parser.py, the list of disease names
	and the list of genes associated with the given diseases. """)
	print()
	print("""\tExample usage: python merge_transcripts.py --config_file LOVD3_Databases.json
	--disease_gene_lists SCID_ny_panel.txt Metabolic_diseases_genes.txt
	--disease_names SCID Metabolic_Diseases --path LOVD3""")
	print()


args = args.parse_args()
config_file = args.config_file
if config_file == None: # If they did not specify a config file, then look for one in the directory
	json_files = glob.glob('*.json')
	if len(json_files) == 0:
		print_help_message()
		print("\tYou have no json file present in your current working directory.")
		print("\tPlease specify a config file with the --config_file option.")
		print()
		sys.exit(1) # Exit with a status of 1. They are probably trying to see what the program does.
	elif len(json_files) > 1:
		print_help_message()
		print("\tYou have multiple json files present in your current working directory. \n\tPlease specify a config file with the --config_file option.")
		print()
		print("\tPossible config files in your directory include:")
		print("\t"+', '.join(json_files))
		print()
		sys.exit(1) # Exit with a status of 1. They are probably trying to see what the program does.
	else:
		config_file = json_files[0]

with open(config_file, 'r') as file:
	databases_info = json5.load(file)

if args.disease_gene_lists == None:
	# If no disease gene files were given, exit the program and list the txt files in the directory that are likely to be the disease gene files
	print_help_message()

	print("""\tNo gene list files specified. Please create a text file for each disease of
	interest with a list of gene symbols, one on each line.""")
	print()
	txt_files_list = glob.glob("*.txt")
	if len(txt_files_list) != 0:
		print("\tPossible gene list files in your directory include: ")
		print("\t"+', '.join(txt_files_list))
		print()
	sys.exit(1) # Exit with a status of 1. They are probably trying to see what the program does.

input_gene_lists = []
for text_file in args.disease_gene_lists:
	input_gene_lists.append([line.rstrip('\n') for line in open(text_file)])

if args.disease_names == None:
	args.disease_names = [splitext(gene_list)[0] for gene_list in args.disease_gene_lists]
disease_names = args.disease_names

path = args.path

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
			'Genomic Normalization Failure Message']

for database in databases_info:
	database_name = database["name"]
	for disease, gene_list in zip(disease_names, input_gene_lists):
		for gene in gene_list:
			list_of_files = glob.glob(path+'/'+database_name+"/"+disease+"/"+gene+"_*.csv")
			if len(list_of_files) > 1:
				combined_df = pd.DataFrame(columns=columns)
				for file in list_of_files:
					combined_variants = combined_df["HGVS Normalized Genomic Annotation"].values
					current_df = pd.read_csv(file, index_col=0)
					current_variants = current_df["HGVS Normalized Genomic Annotation"].values
					unique_list = []
					for variant in current_variants:
						if variant not in combined_variants:
							unique_list.append(variant)
					unique_df = current_df[current_df["HGVS Normalized Genomic Annotation"].isin(unique_list)]
					combined_df = pd.concat([combined_df, unique_df]).reset_index(drop=True)
				combined_df.to_csv(path+'/'+database_name+'/'+disease+'/'+gene+'_Combined_Transcripts_'+database_name+'_results.csv')
				# Now that you have a new file with the variants from all transcripts, move the old files into a new directory
				os.makedirs(path+'/'+database_name+"/"+disease+"/Individual_Transcripts", exist_ok = True)
				# Make a directory to put the original files into instead of deleting them
				for file in list_of_files:
					file_name = file.split("/")[-1]
					os.rename(file, path+'/'+database_name+"/"+disease+"/Individual_Transcripts/"+file_name)

# I do not yet have a way to do this with the variants that did not pass validation because I use the normalized variant information for comparison.
