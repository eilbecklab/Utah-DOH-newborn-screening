#!/usr/bin/env python
# coding: utf-8

import os
from argparse import ArgumentParser
import json5
import glob
from os.path import splitext
import sys
import pandas as pd
import warnings
import json
import numpy as np
import matplotlib.pyplot as plt
import re

plt.rcParams.update({'figure.max_open_warning': 0})

args = ArgumentParser('./Analyze_Parsed_Variants.py', description="""This program has been
designed to gather information about the variants that have been parsed with the
parsers LOVD2_Variant_Parser.py, LOVD3_Variant_Parser.py and ClinVar_Parser.py.
This program will count the number of variants of each type per each gene for each
database that passed validation with HGVS. Additionally, a file will be created with
all of the variants from all parsers and a file will be created with all variants
that passed validation but changed when they were validated with HGVS.

Example usage: python Analyze_Parsed_Variants.py -c Output/ClinVar -2 Output/CCHMC
-3 LOVD3 --config_file LOVD3_Databases.json --disease_gene_lists SCID_ny_panel.txt
Metabolic_diseases_genes.txt --disease_names SCID Metabolic_Diseases
--include_no_stars -p Parsed -d All_Parsed""")

args.add_argument(
	'-c',
	'--clinvar',
	help="""The path to the directory containing the output of ClinVar_Parser.py.
	Note, if you specified 'Output' using the --output_directory option then your
	directory will be 'Output/ClinVar'. """,
	default=None,
)

args.add_argument(
	'-2',
	'--lovd2',
	help="""The path to the directory containing the output of LOVD2_Variant_Parser.py.
	Note, if you specified 'Output' using the --output_directory option then your
	directory will be 'Output/CCHMC'. """,
	default=None,
)

args.add_argument(
	'-3',
	'--lovd3',
	help="""The path to the directory containing the output of LOVD3_Variant_Parser.py.
	Unlike the options -c and -2, the directory used for this option is the one
	specified in LOVD3_Variant_Parser.py. This is because multiple databases are
	present for LOVD3. If you specified 'Output' using the --output_directory
	option then your directory will be 'Output'. Individual databases can be
	determined by list of subdirectories in the directory provided or they can
	be directly stated using either the --config_file or -databases options. """,
	default=None,
)

args.add_argument(
	'--disease_gene_lists',
	nargs='+', # This tells the program that if they specify this flag, they have to give it at least one input. If they don't specify it, then the default will go in.
	help="""Only use this option if you used the same gene lists for all parsers.
	If you have used different lists for the individual parsers, some of the
	genes may be overlooked in the counting of variants of each type. This is a
	list of text files containing genes associated with the disease of interest.
	The text files should each contain a list of gene symbols for one given
	disease, one gene per line. The name of the disease must be specified by the
	argument --disease_names. """,
	default = None
)

args.add_argument(
	'--disease_names',
	nargs='+',
	help="""This is a list of the disease names used for the parsers. If neither
	this option nor the options for specifying diseases used for individual
	parsers, then the subdirectory names from the directories specified with the
	-c, -2, and --config_file options will be used. If more than one parser was
	used and different disease names were used for the separate parsers, you can
	specify the disease names for each parser with the options
	--disease_names_clinvar, --disease_names_lovd2, and --disease_names_lovd3.
	Example: --disease_names SCID Metabolic_Diseases""",
	default = None
)

args.add_argument(
	'-p',
	'--output_prefix',
	help="""This is the prefix that you would like to use for the output files.
	Default is Parsed_Variants""",
	default = "Parsed_Variants"
)

args.add_argument(
	'-d',
	'--output_directory',
	help="""This is the directory where you would like to store the output files.
	Default is the current working directory.""",
	default = "./"
)

args.add_argument(
	'--databases',
	nargs='+',
	help="""This is a list of the database names used in the config file for
	LOVD3_Variant_Parser.py. This option should not be used in conjunction with
	the --config_file option. If neither option is used, all subdirectories in
	the directory specified with the --lovd3 option will be used. An example
	input for this option would be -d BIPmed_SNPhg19 BIPmed_WES Global_Variome
	Human_Variome MSeqDR-LSDB """,
	default = None
)

args.add_argument(
	'--config_file',
	help="""This is the config file in JSON format that was used for the input
	to LOVD3_Variant_Parser.py. For an example of formatting, see the file
	LOVD3_Databases.json. This option should not be used in conjunction with the
	--databases option. If neither option is used, all subdirectories in the
	directory specified with the --lovd3 option will be used.""",
	default=None,
)

args.add_argument(
	'--disease_names_clinvar',
	nargs='+',
	help="""This is a list of the disease names used for ClinVar_Parser.py. Use
	this option only if you have used a different list of disease names for
	ClinVar_Parser.py than were used for the other parsers. If neither this
	option nor the --disease_names option is used, the program will use the
	subdirectories in the directory specified by the --clinvar option as disease
	names. Example: --disease_names_clinvar	SCID Metabolic_Diseases""",
	default = None
)

args.add_argument(
	'--disease_names_lovd2',
	nargs='+',
	help="""This is a list of the disease names used for LOVD2_Variant_Parser.py.
	Use this option only if you have used a different list of disease names for
	LOVD2_Variant_Parser.py than were used for the other parsers. If neither this
	option nor the --disease_names option is used, the program will use the
	subdirectories in the directory specified by the --lovd2 option as disease
	names. Example: --disease_names_lovd2 SCID Metabolic_Diseases""",
	default = None
)

args.add_argument(
	'--disease_names_lovd3',
	nargs='+',
	help="""This is a list of the disease names used for LOVD3_Variant_Parser.py.
	Use this option only if you have used a different list of disease names for
	LOVD3_Variant_Parser.py than were used for the other parsers. If neither this
	option nor the --disease_names option is used, the program will use the
	subdirectories in the directory specified by the --lovd2 option as disease
	names. Example: --disease_names_lovd3 SCID Metabolic_Diseases""",
	default = None
)

args.add_argument(
	'--use_file_names',
	help="""Specify this option if you have used the same gene lists for all
	parsers and would like to use the names of the text files (without the final
	extension) for the disease names. This is recommended only if you allowed the
	parsers to determine the disease names based upon the lists provided. """,
	action= 'store_true'
)

args.add_argument(
	'--include_no_stars',
	help="""Specify this option if you would like to include the no-star variants
	from ClinVar_Parser.py. By default, this program will ignore them.""",
	action= 'store_true'
)

args.add_argument(
	'--exclude_stars',
	help="""Specify this option if you would like to include only the no-star
	variants from ClinVar_Parser.py, but exclude the variants that had 1 or more
	stars. By default, this program will include only variants that have 1 or
	more stars according to ClinVar.""",
	action= 'store_true'
)

args.add_argument(
	'--save_dictionary',
	help="""If this option is specified, then the program will save the dictionary
	that contains the information about the parsers, databases, and diseases
	specified in the input along with the path to each of the files associated
	with the input information. It will be saved in a JSON file. By default, the
	program will not save the dictionary.""",
	action= 'store_true'
)

args.add_argument(
	'--no_bar_charts',
	help="""By default, the program will save bar charts representing the number of
	variants of each type for each gene, database and parser. If this option is
	specified, then the program will not save these figures.""",
	action= 'store_true'
)

args.add_argument(
	'--ignore_invalid_hgvs',
	help="""By default, the program will go through the files containing variants
	that did not pass HGVS normalization through biocommons. If this option is
	specified, then the program will not save any files containing information
	about these variants.""",
	action= 'store_true'
)

def print_help_message():
	print()
	print("""\tWelcome to Analyze_Parsed_Variants.py. This program has been designed to count the number
	of variants of each type (SNV, deletion, indel, etc.) for each gene in each database that
	you have parsed. It will also create an output file that contains all HGVS validated
	variants and a file that contains variants that changed during validation. If you did not
	run one of the parsers, just do not use the options that point to the path for the parser. """)
	print()
	print("""\tExample usage: python Analyze_Parsed_Variants.py -c Output/ClinVar -2 Output/CCHMC
	-3 LOVD3 --config_file LOVD3_Databases.json --disease_gene_lists SCID_ny_panel.txt
	Metabolic_diseases_genes.txt --disease_names SCID Metabolic_Diseases
	--include_no_stars -p Parsed -d All_Parsed""")
	print()


args = args.parse_args()
clinvar_directory = args.clinvar
lovd2_directory = args.lovd2
lovd3_directory = args.lovd3
if all(not variable for variable in [clinvar_directory, lovd2_directory, lovd3_directory]):
	print_help_message()
	print("""\tYou have not entered a path for any variant files. Please use the -c, -2, or -3 option
	to specify the path to the outputs from ClinVar_Parser.py, LOVD2_Variant_Parser.py or
	LOVD3_Variant_Parser.py""")
	print()
	sys.exit(1) # Exit with a status of 1. They are probably trying to see what the program does.

lovd3_database_input = args.databases
lovd3_config = args.config_file
# If the lovd3_config is None, leave it alone
# Otherwise, load the json file
if lovd3_config:
	with open(lovd3_config, 'r') as file:
		lovd3_config = json5.load(file)

global_disease_names = args.disease_names
use_file_names = args.use_file_names
provided_gene_lists = args.disease_gene_lists
input_gene_lists = None
if provided_gene_lists:
	input_gene_lists = []
	for text_file in provided_gene_lists:
		invidual_disease_gene_list = [line.rstrip('\n') for line in open(text_file)]
		if '' in invidual_disease_gene_list:
			invidual_disease_gene_list.remove('') # This happens when they have extra lines at the end of the text file
		input_gene_lists.append(invidual_disease_gene_list)
	# If they did not specify global disease names but they said to use the file names, then obtain those
	if not global_disease_names:
		if use_file_names:
			global_disease_names = [splitext(gene_list)[0] for gene_list in provided_gene_lists]
		else:
			print_help_message()
			print("""\tYou have entered gene lists, but you have not entered the disease names and
	you have not specified to use the names of the input files (minus the extension)
	in place of the disease names. Please run using the --disease_names option or
	run without the --disease_gene_lists option and allow the program to use only
	the genes that are present in the directories.  """)
			sys.exit(1)
	if len(global_disease_names) != len(input_gene_lists):
		print_help_message()
		print("""\tYou have entered gene lists, but the length of your disease names does not
	match the length of your list of files containing gene lists. """)
		sys.exit(1)

clinvar_disease_names = args.disease_names_clinvar
lovd2_disease_names = args.disease_names_lovd2
lovd3_disease_names = args.disease_names_lovd3

include_no_stars = args.include_no_stars
exclude_stars = args.exclude_stars

output_prefix = args.output_prefix
output_dir = args.output_directory

save_dictionary = args.save_dictionary
no_bar_charts = args.no_bar_charts
ignore_invalid = args.ignore_invalid_hgvs

os.makedirs(output_dir, exist_ok = True)
# Now start building the dictionary with all of the path information
info_dict = {}
if not ignore_invalid:
	invalid_info_dict = {}
# Start with ClinVar
if clinvar_directory:
	# First add a slash to the end of this so that it will be consistent with the ones from LOVD3
	if not clinvar_directory.endswith('/'):
		clinvar_directory = clinvar_directory+'/'
	info_dict['ClinVar'] = {}
	# There is only one database for the ClinVar Parser, which is ClinVar
	info_dict['ClinVar']['ClinVar'] = {}
	if not ignore_invalid:
		invalid_info_dict['ClinVar'] = {}
		invalid_info_dict['ClinVar']['ClinVar'] = {}
	# Now add the disease names and the files associated with those diseases
	if not clinvar_disease_names:
		if global_disease_names:
			clinvar_disease_names = global_disease_names
		else:
			clinvar_disease_paths = glob.glob(clinvar_directory+'*/')
			clinvar_disease_names = []
			for path in clinvar_disease_paths:
				clinvar_disease_names.append(path.split('/')[-2])
	for disease in clinvar_disease_names:
		if exclude_stars:
			info_dict['ClinVar']['ClinVar'][disease] = glob.glob(clinvar_directory+disease+'/*No_Star_Results.csv')
		elif include_no_stars:
			info_dict['ClinVar']['ClinVar'][disease] = glob.glob(clinvar_directory+disease+'/*Results.csv')
		else:
			info_dict['ClinVar']['ClinVar'][disease] = glob.glob(clinvar_directory+disease+'/*ClinVar_Results.csv')
		if not ignore_invalid:
			invalid_info_dict['ClinVar']['ClinVar'][disease] = glob.glob(clinvar_directory+disease+'/Invalid_Annotations/*Results.csv')

# Now add the LOVD2 database
if lovd2_directory:
	# First add a slash to the end of this so that it will be consistent with the ones from LOVD3
	if not lovd2_directory.endswith('/'):
		lovd2_directory = lovd2_directory+'/'
	info_dict['LOVD2'] = {}
	if not ignore_invalid:
		invalid_info_dict['LOVD2'] = {}
	# This is set up to use only one database for LOVD2, the CCHMC database
	# Even though this should be consistent, I am going to use the final directory
	# name from the path specified in case they have changed it (perhaps spelling
	# out the full name)
	database_name = lovd2_directory.split('/')[-2]
	info_dict['LOVD2'][database_name] = {}
	if not ignore_invalid:
		invalid_info_dict['LOVD2'][database_name] = {}
	# Now add the disease names and the files associated with those diseases
	if not lovd2_disease_names:
		if global_disease_names:
			lovd2_disease_names = global_disease_names
		else:
			lovd2_disease_paths = glob.glob(lovd2_directory+'*/')
			lovd2_disease_names = []
			for path in clinvar_disease_paths:
				lovd2_disease_names.append(path.split('/')[-2])
	for disease in lovd2_disease_names:
		info_dict['LOVD2'][database_name][disease] = glob.glob(lovd2_directory+disease+'/*Results.csv')
		if not ignore_invalid:
			invalid_info_dict['LOVD2'][database_name][disease] = glob.glob(lovd2_directory+disease+'/Invalid_Annotations/*Results.csv')

# Now the tricky one, add LOVD3 databases

if lovd3_directory:
	if not lovd3_directory.endswith('/'):
		lovd3_directory = lovd3_directory+'/'
	info_dict['LOVD3'] = {}
	if not ignore_invalid:
		invalid_info_dict['LOVD3'] = {}
	# Now find all of the database names for the LOVD3
	# If the config file is present, use that as it will be consistent with
	# what was used for LOVD3_Variant_Parser.py
	lovd3_database_names = []
	if lovd3_config:
		if lovd3_database_input:
			warnings.warn("You have used both the --databases and --config_file options. The database names will be taken from the config_file provided.")
		for database in lovd3_config:
			lovd3_database_names.append(database['name'])
	elif lovd3_database_input:
		lovd3_database_names = lovd3_database_input
	else:
		# This means neither option was provided
		lovd3_database_paths = glob.glob(lovd3_directory+'*/')
		for path in lovd3_database_paths:
			lovd3_database_names.append(path.split('/')[-2])
	# Now that the database names are set, add them all to the dictionary
	for database in lovd3_database_names:
		info_dict['LOVD3'][database] = {}
		if not ignore_invalid:
			invalid_info_dict['LOVD3'][database] = {}
	# Now add the disease names.
	# If they have specified the disease names, use those for every directory
	if not lovd3_disease_names:
		if global_disease_names:
			lovd3_disease_names = global_disease_names
	if lovd3_disease_names:
		for database in info_dict['LOVD3']:
			for disease in lovd3_disease_names:
				info_dict['LOVD3'][database][disease] = glob.glob(lovd3_directory+database+'/'+disease+'/*results.csv')
				if not ignore_invalid:
					invalid_info_dict['LOVD3'][database][disease] = glob.glob(lovd3_directory+database+'/'+disease+'/Invalid_Annotations/*Results.csv')
	else:
		# This means no input was used to give the disease names
		# They might be different for each database, so use glob to find the names
		for database in info_dict['LOVD3']:
			disease_names = []
			disease_paths = glob.glob(lovd3_directory+database+'/*/')
			for path in disease_paths:
				disease_names.append(path.split('/')[-2])
			# Now we have the disease names for that specific database
			for disease in disease_names:
				info_dict['LOVD3'][database][disease] = glob.glob(lovd3_directory+database+'/'+disease+'/*results.csv')
				if not ignore_invalid:
					invalid_info_dict['LOVD3'][database][disease] = glob.glob(lovd3_directory+database+'/'+disease+'/Invalid_Annotations/*Results.csv')

# If they specified that they want to save the dictionary, save it to a json file
if not ignore_invalid:
	os.makedirs(output_dir+'/Invalid_HGVS', exist_ok = True)

if save_dictionary:
	dictionary_output = output_dir+'/'+output_prefix+'_info_dictionary.json'
	with open (dictionary_output, 'w') as file:
		json.dump(info_dict, file, indent='\t')
	if not ignore_invalid:
		invalid_dict_output = output_dir+'/Invalid_HGVS/'+output_prefix+'_Invalid_info_dictionary.json'
		with open (invalid_dict_output, 'w') as file:
			json.dump(invalid_info_dict, file, indent='\t')


# Now read the variants from all files into one dataframe.
ClinVar_Columns = ['Genome Assembly', 'Chr', 'Position Start', 'Position Stop', 'Ref',
					'Alt', 'Genomic Annotation', 'HGVS Normalized Genomic Annotation',
					'Variant Type', 'Variant Length', 'Pathogenicity', 'Disease',
					'Genetic Origin', 'Inheritance Pattern', 'Affected Genes',
					'Gene Symbol', 'dbSNP ID', 'Compound Het Status', 'Transcript',
					'Transcript Notation', 'HGVS Transcript Notation', 'Protein Accession',
					'Protein Notation', 'HGVS Protein Annotation', 'Chr Accession',
					'VCF Pos', 'VCF Ref', 'VCF Alt', 'Database', 'ClinVar Accession',
					'Review Status', 'Star Level', 'Submitter', 'Edited Date',
					'Transcript Normalization Failure Message', 'Genomic Normalization Failure Message']



combined_variants_df = pd.DataFrame(columns = ClinVar_Columns)

for parser in info_dict:
	for database in info_dict[parser]:
		for disease in info_dict[parser][database]:
			for file in info_dict[parser][database][disease]:
				tmp_df = pd.read_csv(file, index_col = 0)
				combined_variants_df = pd.concat([combined_variants_df, tmp_df[ClinVar_Columns]]).reset_index(drop=True)

if not ignore_invalid:
	invalid_columns = ClinVar_Columns+['HGVS Normalization Failure Reason']
	invalid_variants_df = pd.DataFrame(columns = invalid_columns)
	for parser in invalid_info_dict:
		for database in invalid_info_dict[parser]:
			for disease in invalid_info_dict[parser][database]:
				for file in invalid_info_dict[parser][database][disease]:
					tmp_df = pd.read_csv(file, index_col = 0)
					invalid_variants_df = pd.concat([invalid_variants_df, tmp_df[invalid_columns]], sort=False).reset_index(drop=True)

# For some reason Biocommons HGVS normalization uses "identity" in some places and "Identity" in others, so case sensitvity can be a problem
# The most common is lower case, so I am going to change all "Identity" to "identity"
def fix_case_identity(df):
	var_type = df['Variant Type']
	if var_type == 'Identity':
		var_type = "identity"
	df['Variant Type'] = var_type
	return df

combined_variants_df = combined_variants_df.apply(fix_case_identity, axis = 1)
# I don't need to fix the case for the invalid variants, none of them were entered as Identity
# Sometimes one gene is in multiple diseases, which would result in multiple entries here. Drop the duplicates to clean this up.
combined_variants_df = combined_variants_df.drop_duplicates()
combined_variants_output = output_dir+'/'+output_prefix+'_combined_variants.csv'
combined_variants_df.to_csv(combined_variants_output)
if not ignore_invalid:
	invalid_variants_df = invalid_variants_df.drop_duplicates()
	invalid_variants_output = output_dir+'/Invalid_HGVS/'+output_prefix+'_combined_invalid_variants.csv'
	invalid_variants_df.to_csv(invalid_variants_output)

# Further analysis (particularly counting the number of variants of each type)
# will exclude the ones from ClinVar that have multiple genes for the same variant
# Those ones will be saved in a separate file for easy access
# This is only applicable if they specified a ClinVar directory
# It isn't particularly helpful if only one disease was run, but is still makes it easier to find
if clinvar_directory:
	def find_multi_gene(df):
		symbol = df['Gene Symbol']
		multi_gene = False
		if ', ' in symbol:
			multi_gene = True
		df['Multiple Genes'] = multi_gene
		return df

	combined_variants_multigene = combined_variants_df.apply(find_multi_gene, axis = 1)
	multi_gene_df = combined_variants_multigene[combined_variants_multigene['Multiple Genes'] == True]
	multi_gene_df = multi_gene_df[multi_gene_df.columns[:-1]]
	if len(multi_gene_df) > 0:
		multi_gene_output = output_dir+'/'+output_prefix+'_ClinVar_MultiGene_Variants.csv'
		multi_gene_df.to_csv(multi_gene_output)


# Now for obtaining the counts of each variant type
rmdup_variants_df = combined_variants_df[['HGVS Normalized Genomic Annotation', 'Variant Type', 'Database', 'Gene Symbol']]
rmdup_variants_df = rmdup_variants_df.drop_duplicates()
# To deal with variants that affect multiple genes, I am going to add the variants with each individual gene.
multi_gene = rmdup_variants_df[rmdup_variants_df["Gene Symbol"].str.contains(", ")]
HGVS_list = []
Type_list = []
DB_List = []
Gene_List = []
for index, row in multi_gene.iterrows():
	HGVS = row['HGVS Normalized Genomic Annotation']
	var_type = row['Variant Type']
	gene_variable = row['Gene Symbol']
	database = row['Database']
	gene_list = gene_variable.split(', ')
	for gene in gene_list:
		HGVS_list.append(HGVS)
		Type_list.append(var_type)
		DB_List.append(database)
		Gene_List.append(gene)

split_multi_gene_df = pd.DataFrame({'HGVS Normalized Genomic Annotation': HGVS_list,
									'Variant Type': Type_list,
									'Database': DB_List,
									'Gene Symbol': Gene_List})

rmdup_variants_df = rmdup_variants_df.append(split_multi_gene_df, ignore_index=True)
rmdup_variants_df = rmdup_variants_df[~rmdup_variants_df["Gene Symbol"].str.contains(", ")]
# Now the variants that have multiple genes will be represented for each gene
count_object = rmdup_variants_df.groupby(['Database', 'Gene Symbol', 'Variant Type']).size()
# This is a pandas series object, but does not contain 0 for anything that is missing
# To obtain those 0 values, I need to make a dictionary that
# To account for the ones that have 0, first build a dictionary containing all
# genes and all variant types. Fill them with 0 and use the count_object to replace
# the 0 values with the counts when they are not 0 (present in the object)
# If they have specified gene lists and the global disease names, use those to build the dictionary
# If they have not provided it, use the info from the count_object

# If they have not provided gene lists, then try to get the gene names from the info_dict.
# We need a list of genes for each disease.
gene_lists_dict = {}
for parser in info_dict:
	for database in info_dict[parser]:
		for disease in info_dict[parser][database]:
			disease_gene_list = []
			for file in info_dict[parser][database][disease]:
				gene_name = file.split('/')[-1].split('_')[0]
				# ClinVar will have some genes twice if the no_stars files are inlcuded
				if gene_name not in disease_gene_list and "Multiple" not in gene_name:
					disease_gene_list.append(gene_name)
			# Database names should be unique, so we don't need the parser here
			if database not in gene_lists_dict:
				gene_lists_dict[database] = {}
			gene_lists_dict[database][disease] = disease_gene_list

counts_dictionary = {}
variant_types = ['single nucleotide variant', 'Deletion', 'Duplication', 'Insertion', 'Indel', 'identity', 'Inversion']
for parser in info_dict:
	counts_dictionary[parser] = {}
	for database in info_dict[parser]:
		counts_dictionary[parser][database] = {}
		for disease in info_dict[parser][database]:
			counts_dictionary[parser][database][disease] = {}
		# It has all been the same structure as the info_dict up until here,
		# but now we need to add the gene symbols and the variant types
		# If they have specified lists of genes, they must accompany global_disease_names


if input_gene_lists:
	for parser in counts_dictionary:
		for database in counts_dictionary[parser]:
			for gene_list, disease in zip(input_gene_lists, global_disease_names): # These two need to be paired together
				for gene in gene_list:
					counts_dictionary[parser][database][disease][gene] = {}
					for var_type in variant_types:
						counts_dictionary[parser][database][disease][gene][var_type] = 0
# If they did not put in gene lists, use the gene lists from the gene_lists_dict
else:
	for parser in counts_dictionary:
		for database in counts_dictionary[parser]:
			for disease in counts_dictionary[parser][database]:
				for gene in gene_lists_dict[database][disease]:
					counts_dictionary[parser][database][disease][gene] = {}
					for var_type in variant_types:
						counts_dictionary[parser][database][disease][gene][var_type] = 0

# The dictionary has been built and has all of the genes and variant types preset
# All counts are currently 0, so they need to be filled in from the count_object created above
# I have a lot of info from the count_object, but I do not have the parser name or the disease name
# I am going to have to grab those from the info_dict and the counts_dictionary
# There is only one parser per database, so this one is simple
db_parser_dict = {}
for parser in info_dict:
	for database in info_dict[parser]:
		db_parser_dict[database] = parser

# Now I need key, value pairs to go from the gene to the disease.
# This way if the disease is present in more than one disease, it will put the counts in both places
# I am going to use the counts_dictionary because it will contain the info from the input
# gene lists files if present or the genes from the files if input lists were not used
gene_disease_dict = {}
for parser in counts_dictionary:
	for database in counts_dictionary[parser]:
		for disease in counts_dictionary[parser][database]:
			for gene in counts_dictionary[parser][database][disease]:
				if gene not in gene_disease_dict:
					gene_disease_dict[gene] = []
				if disease not in gene_disease_dict[gene]:
					gene_disease_dict[gene].append(disease)

# Now populate the counts_dictionary
for index_tuple, var_count in zip(count_object.index, count_object):
	gene = index_tuple[1]
	database = index_tuple[0]
	var_type = index_tuple[2]
	parser = db_parser_dict[database]
	diseases_list = gene_disease_dict[gene]
	for disease in diseases_list:
		counts_dictionary[parser][database][disease][gene][var_type] = var_count


# Now run through each parser, database, disease, and gene and record the counts of
# each variant type into a csv file
header_line = ["Parser", "Database", "Disease", "Gene", "SNV", "Deletion", "Duplication", "Insertion", "Indel", "Identity", "Inversion", "Total Variants"]
counts_output = output_dir+'/'+output_prefix+'_variant_type_counts.csv'
counts_df = pd.DataFrame(columns = header_line)
for parser in counts_dictionary:
	for database in counts_dictionary[parser]:
		for disease in counts_dictionary[parser][database]:
			for gene in counts_dictionary[parser][database][disease]:
				snv = counts_dictionary[parser][database][disease][gene]['single nucleotide variant']
				deletion = counts_dictionary[parser][database][disease][gene]['Deletion']
				duplication = counts_dictionary[parser][database][disease][gene]['Duplication']
				insertion = counts_dictionary[parser][database][disease][gene]['Insertion']
				indel = counts_dictionary[parser][database][disease][gene]['Indel']
				identity = counts_dictionary[parser][database][disease][gene]['identity']
				inversion = counts_dictionary[parser][database][disease][gene]['Inversion']
				total_var = snv + deletion + duplication + insertion + indel + identity + inversion
				all_variables = [parser, database, disease, gene, snv, deletion, duplication, insertion, indel, identity, inversion, total_var]
				tmp_dict = {}
				for variable, value in zip(header_line, all_variables):
					tmp_dict[variable] = value
				counts_df = counts_df.append(tmp_dict, ignore_index=True)

counts_df.to_csv(counts_output)


# Unless they specified not to save the bar charts, make bar charts with the number of variants of each type for each disease, database, parser
def subset_df(df, specified_parser, specified_database, specified_disease):
	parser = df['Parser']
	db = df['Database']
	disease = df['Disease']
	output = False
	if parser == specified_parser and db == specified_database and disease == specified_disease:
		output = True
	df['Present_Disease'] = output
	return df

if not no_bar_charts:
	os.makedirs(output_dir+'/Variant_Count_Bar_Charts', exist_ok = True)
	# First print out the legend so that it doesn't cover up the bars in the charts
	# All will have the same legend because they all come from the same original DF
	# Give it only 30 lines so that it runs faster
	if len(counts_df) > 30:
		legend_df = counts_df[counts_df.columns[3:-2]].head(30)
	else:
		legend_df = counts_df[counts_df.columns[3:-2]]
	f, (ax1) = plt.subplots(1,1, sharex=True, figsize = (12,8)) # This is just to creat the object for extracting the legend
	legend_df.plot(kind="bar", stacked=True, ax = f.gca())
	figsize = (1.5, 1.5)
	fig_leg = plt.figure(figsize=figsize)
	ax_leg = fig_leg.add_subplot(111)
	# add the legend from the previous axes
	ax_leg.legend(*ax1.get_legend_handles_labels(), loc='center')
	# hide the axes frame and the x/y labels
	ax_leg.axis('off')
	fig_leg.savefig(output_dir+'/Variant_Count_Bar_Charts/Legend.png')
	# Now that the legend is saved, loop through each parser, database, and disease

	for parser in counts_dictionary:
		for database in counts_dictionary[parser]:
			for disease in counts_dictionary[parser][database]:
				tmp_df = counts_df.apply(subset_df, args=(parser, database, disease), axis = 1)
				tmp_df = tmp_df[tmp_df['Present_Disease'] == True]
				# If it has no counts, skip the plot
				if tmp_df['Total Variants'].max() == 0:
					print()
					print('\tThere are no variant counts for ')
					print('\t\tParser: '+parser)
					print('\t\tDatabase: '+ database)
					print('\t\tDisease: '+disease)
					print('\tThis plot will be skipped.')
					continue
				subset = tmp_df[tmp_df.columns[3:-3]].set_index('Gene')

				f = plt.figure(figsize = (12,9))
				plt.title(parser+' '+database+' '+disease)
				plt.ylabel('Variant Count')
				subset.plot(kind="bar", stacked=True, ax = f.gca())
				plt.legend().remove()
				plt.savefig(output_dir+'/Variant_Count_Bar_Charts/'+database+'_'+disease+'_Variant_Counts_Bar_Chart.png')
				plt.close('all')
				# It is common that there is one gene that has far more variants than all of the other ones.
				# Check to see if it is possible to make one with a split Y-axis.
				# The bottom part will always be SNV, so you need to split it in the middle of that one.
				highest_snv = subset['SNV'].max()
				second_highest_total = np.sort(tmp_df['Total Variants'].values)[-2]
				if second_highest_total == 0:
					continue
				if highest_snv > 2*second_highest_total:
					# First find appropriate cutoffs for the axes
					highest_value = tmp_df['Total Variants'].max()
					upper_split = 0.97*highest_snv
					lower_split = 1.05*second_highest_total
					upper_limit = 1.02*highest_value
					# Now plot it over two plots
					fig, (ax1, ax2) = plt.subplots(2,1, sharex=True, figsize = (12,9))
					ax1.spines['bottom'].set_visible(False)
					ax1.tick_params(axis='x',which='both',bottom=False)
					ax2.spines['top'].set_visible(False)
					# Set the split using the numbers calculated above
					ax2.set_ylim(0,lower_split)
					ax1.set_ylim(upper_split,upper_limit)
					subset.plot(ax=ax1,kind='bar', stacked = True)
					subset.plot(ax=ax2,kind='bar', stacked = True)
					plt.ylabel('Variant Count')
					# Add slashes on the axes to make it look clearer
					d = .005
					kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
					ax1.plot((-d, +d), (-d, +d), **kwargs)
					ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)
					kwargs.update(transform=ax2.transAxes)
					ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)
					ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)
					ax1.legend().remove()
					ax2.legend().remove()
					plt.savefig(output_dir+'/Variant_Count_Bar_Charts/'+database+'_'+disease+'_Variant_Counts_Split_Y_Bar_Chart.png')


####################
####################

# Now to look for the overlap between LOVD and ClinVar
# I have already made a dataframe with the duplicates filtered out, so now I just need to get the counts of overlap
def relabel_identity(df):
	type = df["Variant Type"]
	if type == "identity":
		type = "single nucleotide variant"
	df["Variant Type"] = type
	return df


if clinvar_directory and (lovd2_directory or lovd3_directory):

	clinvar_unique = rmdup_variants_df[rmdup_variants_df["Database"] == "ClinVar"]
	lovd_unique = rmdup_variants_df[rmdup_variants_df["Database"] != "ClinVar"]
	# Check the input_gene_lists and global_disease_names to make sure they are present. If not, create them.
	if not input_gene_lists:
		input_gene_lists = []
		if not global_disease_names:
			clinvar_disease_names = []
			lovd_disease_names = []
			global_disease_names = []
			for disease in gene_lists_dict["ClinVar"]:
				clinvar_disease_names.append(disease)
			for parser in gene_lists_dict:
				if parser != "ClinVar":
					for disease in gene_lists_dict[parser]:
						if disease not in lovd_disease_names:
							lovd_disease_names.append(disease)
			for disease in clinvar_disease_names:
				if disease in lovd_disease_names:
					global_disease_names.append(disease)
		# The global_disease_names list has now been built if it was not present to begin with
		for disease in global_disease_names:
			genes_list = gene_lists_dict["ClinVar"][disease]
			for parser in gene_lists_dict:
				if parser != "ClinVar":
					tmp_gene_list = gene_lists_dict[parser][disease]
					for gene in tmp_gene_list:
						if gene not in genes_list:
							genes_list.append(gene) # This will add any genes that were not in ClinVar, but present for the disease in the other databases
			input_gene_lists.append(genes_list)
	for gene_list, disease in zip(input_gene_lists, global_disease_names):
		os.makedirs(output_dir+'/Comparison_LOVD_Clinvar/'+disease, exist_ok = True)
		clinvar_only_counts = []
		lovd_only_counts = []
		overlap_counts = []
		clinvar_only_percentages = []
		lovd_only_percentages = []
		overlap_percentages = []
		disease_overlap_df = pd.DataFrame(columns = combined_variants_df.columns)
		disease_clinvar_df = pd.DataFrame(columns = combined_variants_df.columns)
		disease_lovd_df = pd.DataFrame(columns = combined_variants_df.columns)
		# don't need to recreate gene_list because it will go in order
		for gene in gene_list:
			if gene == '':
				continue # If they have left a blank line in the gene list file, this can happen
			clinvar_only = []
			lovd_only = []
			overlap = []
			clinvar_variants = clinvar_unique[clinvar_unique["Gene Symbol"] == gene]["HGVS Normalized Genomic Annotation"].values
			lovd_variants = lovd_unique[lovd_unique["Gene Symbol"] == gene][['HGVS Normalized Genomic Annotation']].drop_duplicates().values
			for variant in clinvar_variants:
				if variant in lovd_variants:
					overlap.append(variant)
				else:
					clinvar_only.append(variant)
			for variant in lovd_variants:
				if variant not in overlap:
					lovd_only.append(variant)
			clinvar_len = len(clinvar_only)
			lovd_len = len(lovd_only)
			overlap_len = len(overlap)
			total_counts = clinvar_len+lovd_len+overlap_len
			clinvar_only_counts.append(clinvar_len)
			lovd_only_counts.append(lovd_len)
			overlap_counts.append(overlap_len)

			if total_counts > 0: # can't calculate percentages if 0 variants
				clinvar_percent = clinvar_len/total_counts
				lovd_percent = lovd_len/total_counts
				overlap_percent = overlap_len/total_counts
				clinvar_only_percentages.append(clinvar_percent)
				lovd_only_percentages.append(lovd_percent)
				overlap_percentages.append(overlap_percent)
			else:
				clinvar_only_percentages.append(0)
				lovd_only_percentages.append(0)
				overlap_percentages.append(0)

			if len(overlap) > 0:
				gene_overlap_df = combined_variants_df[combined_variants_df['HGVS Normalized Genomic Annotation'].isin(overlap)].sort_values(by='HGVS Normalized Genomic Annotation')
				gene_overlap_df.to_csv(output_dir+'/Comparison_LOVD_Clinvar/'+disease+'/'+gene+'_LOVD_ClinVar_Overlap_Variants.csv')
				disease_overlap_df = pd.concat([disease_overlap_df, gene_overlap_df]).reset_index(drop=True)
			if len(clinvar_only) > 0:
				gene_clinvar_only_df = combined_variants_df[combined_variants_df['HGVS Normalized Genomic Annotation'].isin(clinvar_only)].sort_values(by='HGVS Normalized Genomic Annotation')
				gene_clinvar_only_df.to_csv(output_dir+'/Comparison_LOVD_Clinvar/'+disease+'/'+gene+'_ClinVar_Only_Variants.csv')
				disease_clinvar_df = pd.concat([disease_clinvar_df, gene_clinvar_only_df]).reset_index(drop=True)
			if len(lovd_only) > 0:
				gene_lovd_df = combined_variants_df[combined_variants_df['HGVS Normalized Genomic Annotation'].isin(lovd_only)].sort_values(by='HGVS Normalized Genomic Annotation')
				gene_lovd_df.to_csv(output_dir+'/Comparison_LOVD_Clinvar/'+disease+'/'+gene+'_LOVD_Only_Variants.csv')
				disease_lovd_df = pd.concat([disease_lovd_df, gene_lovd_df]).reset_index(drop=True)
		# Now save the combined dataframes
		disease_overlap_df.to_csv(output_dir+'/Comparison_LOVD_Clinvar/'+disease+'/'+disease+'_all_genes_LOVD_ClinVar_Overlap_Variants.csv')
		disease_clinvar_df.to_csv(output_dir+'/Comparison_LOVD_Clinvar/'+disease+'/'+disease+'_all_genes_ClinVar_Only_Variants.csv')
		disease_lovd_df.to_csv(output_dir+'/Comparison_LOVD_Clinvar/'+disease+'/'+disease+'_all_genes_LOVD_Only_Variants.csv')
		# Now make a dataframe out of the counts and save the file as a csv
		comparison_counts_df = pd.DataFrame({'Gene':gene_list,
											'ClinVar Only': clinvar_only_counts,
											'Overlap': overlap_counts,
											'LOVD Only': lovd_only_counts}).set_index('Gene')
		comparison_counts_df.to_csv(output_dir+'/Comparison_LOVD_Clinvar/'+disease+'_LOVD_ClinVar_Overlap_Counts.csv')
		comparison_percent_df = pd.DataFrame({'Gene':gene_list,
											'ClinVar Only': clinvar_only_percentages,
											'Overlap': overlap_percentages,
											'LOVD Only': lovd_only_percentages}).set_index('Gene')
		comparison_percent_df.to_csv(output_dir+'/Comparison_LOVD_Clinvar/'+disease+'_LOVD_ClinVar_Overlap_Percentages.csv')
		# now save the dataframe with the variants in each group
		# and append the counts (lenghts) to lists
		# Print the legend
		if "Legend.png" not in os.listdir(output_dir+'/Comparison_LOVD_Clinvar/'): # Only generate the legend once
			f, (ax1) = plt.subplots(1,1, sharex=True) # This is just to creat the object for extracting the legend
			comparison_counts_df.plot(kind="bar", stacked=True, ax = f.gca())
			figsize = (1.5, 1.5)
			fig_leg = plt.figure(figsize=figsize)
			ax_leg = fig_leg.add_subplot(111)
			# add the legend from the previous axes
			ax_leg.legend(*ax1.get_legend_handles_labels(), loc='center')
			# hide the axes frame and the x/y labels
			ax_leg.axis('off')
			fig_leg.savefig(output_dir+'/Comparison_LOVD_Clinvar/Legend.png')

		f = plt.figure(figsize = (12,9))
		plt.title(disease+ ' Variant Overlap between ClinVar and LOVD' )
		plt.ylabel('Variant Count')
		comparison_counts_df.plot(kind="bar", stacked=True, ax = f.gca())
		plt.legend().remove()
		plt.savefig(output_dir+'/Comparison_LOVD_Clinvar/'+disease+'_Comparison_Bar_Chart.png')

		f = plt.figure(figsize = (12,9))
		plt.title(disease+' Variant Overlap between ClinVar and LOVD')
		plt.ylabel('Proportion of Variants')
		comparison_percent_df.plot(kind="bar", stacked=True, ax = f.gca())
		plt.legend().remove()
		plt.savefig(output_dir+'/Comparison_LOVD_Clinvar/'+disease+'_Comparison_Proportion_Bar_Chart.png')

		## Now get the variant type counts per disease and per gene.
		var_types_dict = {"SNV": "single nucleotide variant",
						"Deletion": "Deletion",
						"Duplication": "Duplication",
						"Insertion": "Insertion",
						"Indel": "Indel",
						"Inversion": "Inversion"}
		lovd_within_disease = lovd_unique[lovd_unique["Gene Symbol"].isin(gene_list)]
		lovd_within_disease = lovd_within_disease.apply(relabel_identity, axis = 1)
		lovd_within_tmp = lovd_within_disease[["HGVS Normalized Genomic Annotation", "Variant Type"]]
		lovd_within_tmp = lovd_within_tmp.drop_duplicates() # If a variant is present in multiple LOVD databases, this will get rid of the duplicates
		lovd_within_count_object = lovd_within_tmp.groupby(["Variant Type"]).size()
		clinvar_within_disease = clinvar_unique[clinvar_unique["Gene Symbol"].isin(gene_list)]
		clinvar_within_disease = clinvar_within_disease.apply(relabel_identity, axis = 1)
		clinvar_within_tmp = clinvar_within_disease[["HGVS Normalized Genomic Annotation", "Variant Type"]]
		clinvar_within_tmp = clinvar_within_tmp.drop_duplicates() # If it has multiple genes covered by the variant, it will be in here multiple times
		clinvar_within_count_object = clinvar_within_tmp.groupby(["Variant Type"]).size()
		lovd_total_counts = 0
		clinvar_total_counts = 0
		disease_counts_dict = {}
		disease_counts_dict["Database"] = ["ClinVar", "LOVD Databases"]
		for key, value in zip(var_types_dict.keys(), var_types_dict.values()):
			lovd_count = 0
			clinvar_count = 0
			if value in clinvar_within_count_object:
				clinvar_count = clinvar_within_count_object[value]
				clinvar_total_counts += clinvar_count
			if value in lovd_within_count_object:
				lovd_count = lovd_within_count_object[value]
				lovd_total_counts += lovd_count
			disease_counts_dict[key] = [clinvar_count, lovd_count]
		disease_counts_dict["Total Variants"] = [clinvar_total_counts, lovd_total_counts]
		disease_var_counts_df = pd.DataFrame.from_dict(disease_counts_dict)
		disease_var_counts_df.to_csv(output_dir+'/Comparison_LOVD_Clinvar/'+disease+'_Variant_Type_Counts_by_disease.csv')

		plot_df = disease_var_counts_df[disease_var_counts_df.columns[:-1]].set_index('Database')
		if "Variant_Counts_Legend.png" not in os.listdir(output_dir+'/Comparison_LOVD_Clinvar/'): # Only generate the legend once
			f, (ax1) = plt.subplots(1,1, sharex=True) # This is just to creat the object for extracting the legend
			plot_df.plot(kind="bar", stacked=True, ax = f.gca())
			figsize = (2, 2.5)
			fig_leg = plt.figure(figsize=figsize)
			ax_leg = fig_leg.add_subplot(111)
			# add the legend from the previous axes
			ax_leg.legend(*ax1.get_legend_handles_labels(), loc='center')
			# hide the axes frame and the x/y labels
			ax_leg.axis('off')
			fig_leg.savefig(output_dir+'/Comparison_LOVD_Clinvar/Variant_Counts_Legend.png')

		f = plt.figure(figsize = (6,10))
		plt.title('Counts of Variants of Each Type for '+ disease)
		plt.ylabel('Variant Count')
		plot_df.plot(kind="bar", stacked=True, ax = f.gca())
		plt.legend().remove()
		plt.savefig(output_dir+'/Comparison_LOVD_Clinvar/'+disease+'_Variant_Counts_by_Condition_Bar_Chart.png')

if save_dictionary:
	gene_lists_output = output_dir+'/'+output_prefix+"_gene_lists_dict.json"
	with open (gene_lists_output, 'w') as file:
		json.dump(gene_lists_dict, file, indent='\t')


####################
####################

# Now make a dataframe with the variants that changed after the HGVS normalization
def label_hgvs_changed(df):
	original = df['Genomic Annotation']
	validated = df['HGVS Normalized Genomic Annotation']
	identical = True
	if original != validated:
		identical = False
	df['HGVS Identical'] = identical
	return df


hgvs_combined_variants_df = combined_variants_df.apply(label_hgvs_changed, axis = 1)
discordant = hgvs_combined_variants_df[hgvs_combined_variants_df['HGVS Identical'] == False]
discordant = discordant[discordant.columns[:-1]]
discordant.to_csv(output_dir+'/'+output_prefix+'_Nonidentical_HGVS.csv')
# Next thing to do is figure out why these ones changed and count the number that changed for each type


if not ignore_invalid:
	# The variants that did not pass HGVS normalization were read into a dataframe above.
	# Now they need to be parsed to count how many of these variants were present for each disease or database

	def simplify_failure_reason(df):
		hgvs = df['HGVS Normalization Failure Reason']
		if hgvs in ['Compound variant', 'Unknown breakpoint', 'Inserted unknown sequence, unknown breakpoint']:
			hgvs = "Complex HGVS Annotation"
		elif "No definitive failure reason" in hgvs:
			hgvs = "Other"
		elif hgvs == "Transcript Accession not in HGVS":
			hgvs = "Transcript Accession not in Biocommons"
		df['Failure Reason'] = hgvs
		return df

	failure_reasons = ['No variant information provided', 'Complex HGVS Annotation', 'Microsatellite',
						'Transcript Accession not in Biocommons','Incorrect reference base', 'Inserted unknown sequence', 'Other']

	invalid_variants_df = invalid_variants_df.apply(simplify_failure_reason, axis = 1)
	invalid_variants_df = invalid_variants_df.reset_index() # It has no unique identifier, but one will be needed later
	invalid_subset = invalid_variants_df[['index', 'Gene Symbol', 'Database', 'Failure Reason']]

	# First to count the failures of each type for all databases
	databases_list = []
	for parser in info_dict:
		for database in info_dict[parser]:
			databases_list.append(database)
	db_invalid_count_dict = {}
	for database in databases_list:
		db_invalid_count_dict[database] = {}
		for reason in failure_reasons:
			db_invalid_count_dict[database][reason] = 0
	invalid_counts = invalid_subset.groupby(['Database', 'Failure Reason']).size()
	for index_tuple, reason_count in zip(invalid_counts.index, invalid_counts):
		database = index_tuple[0]
		reason = index_tuple[1]
		db_invalid_count_dict[database][reason] = reason_count
	header_line = ["Database"] + failure_reasons
	db_invalid_count_df = pd.DataFrame.from_dict(db_invalid_count_dict, orient='index')
	db_invalid_count_df.to_csv(output_dir+'/Invalid_HGVS/'+output_prefix+'_Failure_reason_counts_per_database.csv')
	# Now to add graphs for the invalid counts
	if "Legend_Invalid.png" not in os.listdir(output_dir+'/Invalid_HGVS/'):
		f, (ax1) = plt.subplots(1,1, sharex=True)
		db_invalid_count_df.plot(kind="bar", stacked=True, ax = f.gca())
		figsize = (5, 2)
		fig_leg = plt.figure(figsize=figsize)
		ax_leg = fig_leg.add_subplot(111)
		# add the legend from the previous axes
		ax_leg.legend(*ax1.get_legend_handles_labels(), loc='center')
		# hide the axes frame and the x/y labels
		ax_leg.axis('off')
		fig_leg.savefig(output_dir+'/Invalid_HGVS/Legend_Invalid.png')

	f = plt.figure(figsize = (6,9))
	plt.title('Counts of Variants that Failed HGVS Normalization' )
	plt.ylabel('Variant Count')
	db_invalid_count_df.plot(kind="bar", stacked=True, ax = f.gca())
	plt.legend().remove()
	plt.savefig(output_dir+'/Invalid_HGVS/All_Databases_Invalid_HGVS_Bar_Chart.png')

	# There is no point in combining the LOVD databases if only one is used
	db_list = list(db_invalid_count_df.index)
	db_list.remove("ClinVar")
	if len(db_list) > 1:
		lovd_no_info = 0
		lovd_complex = 0
		lovd_micro = 0
		lovd_no_biocommons = 0
		lovd_incorrect = 0
		lovd_insert_unknown = 0
		lovd_other = 0
		for index, row in db_invalid_count_df.iterrows():
			no_info = row['No variant information provided']
			complex_hgvs = row['Complex HGVS Annotation']
			micro = row['Microsatellite']
			not_in_biocommons = row['Transcript Accession not in Biocommons']
			incorrect = row['Incorrect reference base']
			insert_unknown = row['Inserted unknown sequence']
			other = row['Other']
			if index != "ClinVar":
				lovd_no_info += no_info
				lovd_complex += complex_hgvs
				lovd_micro += micro
				lovd_no_biocommons += not_in_biocommons
				lovd_incorrect += incorrect
				lovd_insert_unknown += insert_unknown
				lovd_other += other
		lovd_tmp = pd.DataFrame({"Database": ["LOVD Databases"], "No variant information provided": lovd_no_info,
										'Complex HGVS Annotation': lovd_complex, 'Microsatellite': lovd_micro,
										'Transcript Accession not in Biocommons': lovd_no_biocommons,
										'Incorrect reference base': lovd_incorrect, 'Inserted unknown sequence': lovd_insert_unknown,
										'Other': lovd_other})
		combined_lovd_df = db_invalid_count_df[db_invalid_count_df.index == "ClinVar"].reset_index().rename(columns={"index":"Database"})
		combined_lovd_df = combined_lovd_df.append(lovd_tmp, ignore_index=True).set_index("Database")

		## This one uses the same legend as the previous one
		f = plt.figure(figsize = (6,9))
		plt.title('Counts of Variants that Failed HGVS Normalization' )
		plt.ylabel('Variant Count')
		combined_lovd_df.plot(kind="bar", stacked=True, ax = f.gca())
		plt.legend().remove()
		plt.savefig(output_dir+'/Invalid_HGVS/Combined_LOVD_Invalid_HGVS_Bar_Chart.png')

	## Now for counting the number of invalid variants of each type per gene
	## This one will be a little bit more complicated because many of these variants cover several genes
	multi_gene_invalid = invalid_subset[invalid_subset["Gene Symbol"].str.contains(", ")]
	invalid_index_list = []
	invalid_Gene_List = []
	invalid_database_List = []
	invalid_failure_List = []
	for index, row in multi_gene_invalid.iterrows():
		index = row['index']
		database = row['Database']
		failure_reason = row['Failure Reason']
		gene_list = row['Gene Symbol'].split(', ')
		for gene in gene_list:
			invalid_index_list.append(index)
			invalid_Gene_List.append(gene)
			invalid_database_List.append(database)
			invalid_failure_List.append(failure_reason)
	invalid_split_multi_gene_df = pd.DataFrame({'index': invalid_index_list,
												'Gene Symbol': invalid_Gene_List,
												'Database': invalid_database_List,
												'Failure Reason': invalid_failure_List})
	invalid_subset = invalid_subset.append(invalid_split_multi_gene_df, ignore_index=True)
	invalid_subset = invalid_subset[~invalid_subset["Gene Symbol"].str.contains(", ")]
	## Now each entry with multiple genes is listed with each gene, so it can loop through and count things per gene
	## I am going to set this up the same way
	invalid_counts_per_gene_dict = {} # Build a dictionary with 0 for each count
	# If they have entered the gene lists, that will be the easiest way to build this.
	if global_disease_names and input_gene_lists:
		for parser in info_dict:
			for database in info_dict[parser]:
				invalid_counts_per_gene_dict[database] = {}
				for disease, gene_list in zip(global_disease_names, input_gene_lists):
					invalid_counts_per_gene_dict[database][disease] = {}
					for gene in gene_list:
						invalid_counts_per_gene_dict[database][disease][gene] = {}
						for reason in failure_reasons:
							invalid_counts_per_gene_dict[database][disease][gene][reason] = 0
	# If they have not entered the gene lists, just use the gene_lists_dict that was built using
	else:
		for database in gene_lists_dict:
			invalid_counts_per_gene_dict[database] = {}
			for disease in gene_lists_dict[database]:
				invalid_counts_per_gene_dict[database][disease] = {}
				gene_list = gene_lists_dict[database][disease]
				for gene in gene_list:
					invalid_counts_per_gene_dict[database][disease][gene] = {}
					for reason in failure_reasons:
						invalid_counts_per_gene_dict[database][disease][gene][reason] = 0

	invalid_counts_per_gene_object = invalid_subset.groupby(['Gene Symbol', 'Database', 'Failure Reason']).size()
	for index_tuple, reason_count in zip(invalid_counts_per_gene_object.index, invalid_counts_per_gene_object):
		gene = index_tuple[0]
		database = index_tuple[1]
		reason = index_tuple[2]
		for disease in invalid_counts_per_gene_dict[database]:
			if gene in invalid_counts_per_gene_dict[database][disease]:
				invalid_counts_per_gene_dict[database][disease][gene][reason] = reason_count
	# Now change the format to a dataframe so that it can be saved as a csv file
	header_line = ["Database", "Disease", "Gene Symbol"] + failure_reasons + ["Total Count"]
	invalid_counts_per_gene_df = pd.DataFrame(columns = header_line)
	for database in invalid_counts_per_gene_dict:
		for disease in invalid_counts_per_gene_dict[database]:
			for gene in invalid_counts_per_gene_dict[database][disease]:
				no_variant = invalid_counts_per_gene_dict[database][disease][gene]["No variant information provided"]
				complex = invalid_counts_per_gene_dict[database][disease][gene]["Complex HGVS Annotation"]
				microsatellite = invalid_counts_per_gene_dict[database][disease][gene]["Microsatellite"]
				no_transcript = invalid_counts_per_gene_dict[database][disease][gene]["Transcript Accession not in Biocommons"]
				incorrect = invalid_counts_per_gene_dict[database][disease][gene]["Incorrect reference base"]
				unknown = invalid_counts_per_gene_dict[database][disease][gene]["Inserted unknown sequence"]
				other = invalid_counts_per_gene_dict[database][disease][gene]["Other"]
				total = no_variant + complex + microsatellite + no_transcript + incorrect + unknown + other
				all_variables = [database, disease, gene, no_variant, complex, microsatellite,
								no_transcript, incorrect, unknown, other, total]
				tmp_dict = {}
				for variable, value in zip(header_line, all_variables):
					tmp_dict[variable] = value
				invalid_counts_per_gene_df = invalid_counts_per_gene_df.append(tmp_dict, ignore_index=True)

	invalid_per_gene_output = output_dir+'/Invalid_HGVS/'+output_prefix+'_Failure_reason_counts_per_gene.csv'
	invalid_counts_per_gene_df.to_csv(invalid_per_gene_output)
	## Now I want the counts of invalid variants per disease/condition
	## Unfortunately, I have made it possible for different confitions per each parser

	databases = []
	for database in invalid_counts_per_gene_dict:
		if database not in databases:
			databases.append(database)
	disease_genes_dict = {}
	if global_disease_names and input_gene_lists:
		for condition, gene_list in zip(global_disease_names, input_gene_lists):
			disease_genes_dict[condition] = gene_list
	else:
		conditions = []
		for database in invalid_counts_per_gene_dict:
			for condition in invalid_counts_per_gene_dict[database]:
				if condition not in conditions:
					conditions.append(condition)
		for condition in conditions:
			condition_gene_list = []
			for database in invalid_counts_per_gene_dict:
				if condition in invalid_counts_per_gene_dict[database]:
					for gene in invalid_counts_per_gene_dict[database][condition]:
						if gene not in condition_gene_list:
							condition_gene_list.append(gene) # If the gene is paired with the condition for any of the databases, it will be added to this list
			disease_genes_dict[condition] = condition_gene_list

	# Now get the counts
	for condition in disease_genes_dict:
		os.makedirs(output_dir+'/Invalid_HGVS/'+condition, exist_ok = True)
		condition_subset = invalid_subset[invalid_subset["Gene Symbol"].isin(disease_genes_dict[condition])]
		invalid_counts_per_gene_dict = {}
		for database in databases+['LOVD Databases Combined']:
			invalid_counts_per_gene_dict[database] = {}
			for gene in disease_genes_dict[condition]:
				invalid_counts_per_gene_dict[database][gene] = {}
				for reason in failure_reasons+['Total']:
					invalid_counts_per_gene_dict[database][gene][reason] = 0
		per_gene_count_object = condition_subset.groupby(["Database", "Gene Symbol", "Failure Reason"]).size()
		for index_tuple, var_count in zip(per_gene_count_object.index, per_gene_count_object):
			gene = index_tuple[1]
			database = index_tuple[0]
			reason = index_tuple[2]
			invalid_counts_per_gene_dict[database][gene][reason] = var_count
		## Now to fill in the total counts and the combined LOVD
		for database in invalid_counts_per_gene_dict:
			for gene in invalid_counts_per_gene_dict[database]:
				total = 0
				for reason in failure_reasons:
					total += invalid_counts_per_gene_dict[database][gene][reason]
				invalid_counts_per_gene_dict[database][gene]['Total'] = total
		LOVD_databases = databases.copy()
		LOVD_databases.remove('ClinVar')
		for gene in disease_genes_dict[condition]:
			for database in LOVD_databases:
				for reason in failure_reasons+['Total']:
					count = invalid_counts_per_gene_dict[database][gene][reason]
					invalid_counts_per_gene_dict['LOVD Databases Combined'][gene][reason] += count
		# Now to put it to a dataframe
		header_line = ["Database", "Gene", "No variant information provided", "Complex HGVS Annotation",
						"Microsatellite", "Transcript Accession not in Biocommons", "Incorrect reference base",
						"Inserted unknown sequence", "Other", "Total"]
		invalid_counts_per_gene_df = pd.DataFrame(columns = header_line)
		for database in invalid_counts_per_gene_dict:
			for gene in invalid_counts_per_gene_dict[database]:
				no_info = invalid_counts_per_gene_dict[database][gene]["No variant information provided"]
				complex_hgvs = invalid_counts_per_gene_dict[database][gene]["Complex HGVS Annotation"]
				microsatellite = invalid_counts_per_gene_dict[database][gene]["Microsatellite"]
				not_in_biocommons = invalid_counts_per_gene_dict[database][gene]["Transcript Accession not in Biocommons"]
				incorrect_ref = invalid_counts_per_gene_dict[database][gene]["Incorrect reference base"]
				insert_unknown = invalid_counts_per_gene_dict[database][gene]["Inserted unknown sequence"]
				other = invalid_counts_per_gene_dict[database][gene]["Other"]
				total = invalid_counts_per_gene_dict[database][gene]["Total"]
				all_variables = [database, gene, no_info, complex_hgvs, microsatellite, not_in_biocommons,
								incorrect_ref, insert_unknown, other, total]
				tmp_dict = {}
				for variable, value in zip(header_line, all_variables):
					tmp_dict[variable] = value
				invalid_counts_per_gene_df = invalid_counts_per_gene_df.append(tmp_dict, ignore_index=True)
		gene_invalid_counts_output = output_dir+'/Invalid_HGVS/'+condition+'/'+condition+'_Invalid_HGVS_counts_per_gene.csv'
		invalid_counts_per_gene_df.to_csv(gene_invalid_counts_output)
		## Now for the bar charts
		## It has been pointed out that it is difficult to see the bars for LOVD variants because there are so many more for ClinVar,
		## so I am going to split them into separate graphs
		clinvar_invalid = invalid_counts_per_gene_df[invalid_counts_per_gene_df["Database"] == "ClinVar"]
		keep_columns = ['Gene', 'No variant information provided', 'Complex HGVS Annotation', 'Microsatellite',
						'Transcript Accession not in Biocommons', 'Incorrect reference base',
						'Inserted unknown sequence', 'Other']
		clinvar_invalid = clinvar_invalid[keep_columns].set_index('Gene')
		if "Legend_Invalid.png" not in os.listdir(output_dir+'/Invalid_HGVS/'+condition+'/'):
			f, (ax1) = plt.subplots(1,1, sharex=True)
			clinvar_invalid.plot(kind="bar", stacked=True, ax = f.gca())
			figsize = (5, 2)
			fig_leg = plt.figure(figsize=figsize)
			ax_leg = fig_leg.add_subplot(111)
			# add the legend from the previous axes
			ax_leg.legend(*ax1.get_legend_handles_labels(), loc='center')
			# hide the axes frame and the x/y labels
			ax_leg.axis('off')
			fig_leg.savefig(output_dir+'/Invalid_HGVS/'+condition+'/Legend_Invalid.png')

		f = plt.figure(figsize = (12,9))
		plt.title('Counts of Variants in ClinVar that Failed HGVS Normalization Per Gene' )
		plt.ylabel('Variant Count')
		clinvar_invalid.plot(kind="bar", stacked=True, ax = f.gca())
		plt.legend().remove()
		plt.savefig(output_dir+'/Invalid_HGVS/'+condition+'/'+condition+'_ClinVar_Invalid_HGVS_per_gene_Bar_Chart.png')

		lovd_invalid = invalid_counts_per_gene_df[invalid_counts_per_gene_df["Database"] == "LOVD Databases Combined"]
		keep_columns = ['Gene', 'No variant information provided', 'Complex HGVS Annotation', 'Microsatellite',
						'Transcript Accession not in Biocommons', 'Incorrect reference base',
						'Inserted unknown sequence', 'Other']
		lovd_invalid = lovd_invalid[keep_columns].set_index('Gene')

		f = plt.figure(figsize = (12,9))
		plt.title('Counts of Variants in LOVD Databases that Failed HGVS Normalization Per Gene' )
		plt.ylabel('Variant Count')
		lovd_invalid.plot(kind="bar", stacked=True, ax = f.gca())
		plt.legend().remove()
		plt.savefig(output_dir+'/Invalid_HGVS/'+condition+'/'+condition+'_LOVD_Invalid_HGVS_per_gene_Bar_Chart.png')


		invalid_counts_per_condition_dict = {}
		for database in databases:
			invalid_counts_per_condition_dict[database] = {}
			for reason in failure_reasons:
				invalid_counts_per_condition_dict[database][reason] = 0
		# Now I don't need the gene anymore
		condition_subset = condition_subset[['index', 'Database', 'Failure Reason']].drop_duplicates()
		count_object = condition_subset.groupby(['Database', 'Failure Reason']).size()
		for index_tuple, var_count in zip(count_object.index, count_object):
			reason = index_tuple[1]
			database = index_tuple[0]
			invalid_counts_per_condition_dict[database][reason] = var_count
		invalid_counts_per_condition_df = pd.DataFrame.from_dict(invalid_counts_per_condition_dict, orient='index')

		# The legend should match the one previously made
		f = plt.figure(figsize = (6,9))
		plt.title('Counts of Variants that Failed HGVS Normalization' )
		plt.ylabel('Variant Count')
		invalid_counts_per_condition_df.plot(kind="bar", stacked=True, ax = f.gca())
		plt.legend().remove()
		plt.savefig(output_dir+'/Invalid_HGVS/'+condition+'/'+condition+'_All_Databases_Invalid_HGVS_Bar_Chart.png')
		# Now to combine the counts from all of the LOVD databases
		lovd_no_info = 0
		lovd_complex = 0
		lovd_micro = 0
		lovd_no_biocommons = 0
		lovd_incorrect = 0
		lovd_insert_unknown = 0
		lovd_other = 0
		for index, row in invalid_counts_per_condition_df.iterrows():
			no_info = row['No variant information provided']
			complex_hgvs = row['Complex HGVS Annotation']
			micro = row['Microsatellite']
			not_in_biocommons = row['Transcript Accession not in Biocommons']
			incorrect = row['Incorrect reference base']
			insert_unknown = row['Inserted unknown sequence']
			other = row['Other']
			if index != "ClinVar":
				lovd_no_info += no_info
				lovd_complex += complex_hgvs
				lovd_micro += micro
				lovd_no_biocommons += not_in_biocommons
				lovd_incorrect += incorrect
				lovd_insert_unknown += insert_unknown
				lovd_other += other
		lovd_tmp = pd.DataFrame({"Database": ["LOVD Databases"], "No variant information provided": lovd_no_info,
								'Complex HGVS Annotation': lovd_complex, 'Microsatellite': lovd_micro,
								'Transcript Accession not in Biocommons': lovd_no_biocommons,
								'Incorrect reference base': lovd_incorrect, 'Inserted unknown sequence': lovd_insert_unknown,
								'Other': lovd_other})

		invalid_counts_per_condition_df = invalid_counts_per_condition_df.reset_index().rename(columns={"index":"Database"})
		invalid_counts_per_condition_df = invalid_counts_per_condition_df.append(lovd_tmp, ignore_index=True).set_index("Database")
		invalid_counts_per_condition_df.to_csv(output_dir+'/Invalid_HGVS/'+condition+'/'+condition+'_Invalid_Variant_Counts.csv')
		invalid_condition_tmp_df = invalid_counts_per_condition_df[invalid_counts_per_condition_df.index.isin(["ClinVar", "LOVD Databases"])]
		f = plt.figure(figsize = (6,9))
		plt.title('Counts of Variants that Failed HGVS Normalization' )
		plt.ylabel('Variant Count')
		invalid_condition_tmp_df.plot(kind="bar", stacked=True, ax = f.gca())
		plt.legend().remove()
		plt.savefig(output_dir+'/Invalid_HGVS/'+condition+'/'+condition+'_Combined_LOVD_Invalid_HGVS_Bar_Chart.png')



## All of the file paths are added to the dictionary.
####### The next steps are to read all of the files into one dataframe
####### save the dataframe
####### Count the number of variants of each type for each gene, disease, database, parser

####### Make a file of only the ones that have normalized annotations that aren't identical to the input genomic annotations
## Count the number of non-identical ones for each gene for each disease for each database
## Count the number of each variant type for those too
###### Then find the overlap between ClinVar and the rest of them

###### Count the number of ones that failed validation per each gene and why
## Then find if the ones that do overlap have the same pathogenicity rating

## Relabel "identity" as SNV for individual gene bar charts
## Simplify loop from 1044 to 1067 where everything is a 0 and then adds from there for LOVD variants
## Simplify failure reasons in all of the original 3 parsers and remove the simplify failure reason loop 
