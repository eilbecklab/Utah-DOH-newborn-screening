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
		input_gene_lists.append([line.rstrip('\n') for line in open(text_file)])
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

os.makedirs(output_dir, exist_ok = True)
# Now start building the dictionary with all of the path information
info_dict = {}
# Start with ClinVar
if clinvar_directory:
	# First add a slash to the end of this so that it will be consistent with the ones from LOVD3
	if not clinvar_directory.endswith('/'):
		clinvar_directory = clinvar_directory+'/'
	info_dict['ClinVar'] = {}
	# There is only one database for the ClinVar Parser, which is ClinVar
	info_dict['ClinVar']['ClinVar'] = {}
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

# Now add the LOVD2 database
if lovd2_directory:
	# First add a slash to the end of this so that it will be consistent with the ones from LOVD3
	if not lovd2_directory.endswith('/'):
		lovd2_directory = lovd2_directory+'/'
	info_dict['LOVD2'] = {}
	# This is set up to use only one database for LOVD2, the CCHMC database
	# Even though this should be consistent, I am going to use the final directory
	# name from the path specified in case they have changed it (perhaps spelling
	# out the full name)
	database_name = lovd2_directory.split('/')[-2]
	info_dict['LOVD2'][database_name] = {}
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

# Now the tricky one, add LOVD3 databases

if lovd3_directory:
	if not lovd3_directory.endswith('/'):
		lovd3_directory = lovd3_directory+'/'
	info_dict['LOVD3'] = {}
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
	# Now add the disease names.
	# If they have specified the disease names, use those for every directory
	if not lovd3_disease_names:
		if global_disease_names:
			lovd3_disease_names = global_disease_names
	if lovd3_disease_names:
		for database in info_dict['LOVD3']:
			for disease in lovd3_disease_names:
				info_dict['LOVD3'][database][disease] = glob.glob(lovd3_directory+database+'/'+disease+'/*results.csv')
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

# If they specified that they want to save the dictionary, save it to a json file
if save_dictionary:
	dictionary_output = output_dir+'/'+output_prefix+'_info_dictionary.json'
	with open (dictionary_output, 'w') as file:
		json.dump(info_dict, file, indent='\t')

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

# Sometimes one gene is in multiple diseases, which would result in multiple entries here. Drop the duplicates to clean this up.
combined_variants_df = combined_variants_df.drop_duplicates()
combined_variants_output = output_dir+'/'+output_prefix+'_combined_variants.csv'
combined_variants_df.to_csv(combined_variants_output)
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

count_object = combined_variants_df.groupby(['Database', 'Gene Symbol', 'Variant Type']).size()
################### How do I drop the ones that have multiple things per the same variant, like in MSeqDR
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
				disease_gene_list.append(file.split('/')[-1].split('_')[0])
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
	if ', ' in gene:
		continue # The dictionary will not populate with commas in the gene name
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
	for databse in counts_dictionary[parser]:
		for disease in counts_dictionary[parser][databse]:
			for gene in counts_dictionary[parser][databse][disease]:
				snv = counts_dictionary[parser][databse][disease][gene]['single nucleotide variant']
				deletion = counts_dictionary[parser][databse][disease][gene]['Deletion']
				duplication = counts_dictionary[parser][databse][disease][gene]['Duplication']
				insertion = counts_dictionary[parser][databse][disease][gene]['Insertion']
				indel = counts_dictionary[parser][databse][disease][gene]['Indel']
				identity = counts_dictionary[parser][databse][disease][gene]['identity']
				inversion = counts_dictionary[parser][databse][disease][gene]['Inversion']
				total_var = snv + deletion + duplication + insertion + indel + identity + inversion
				all_variables = [parser, databse, disease, gene, snv, deletion, duplication, insertion, indel, identity, inversion, total_var]
				tmp_dict = {}
				for variable, value in zip(header_line, all_variables):
					tmp_dict[variable] = value
				counts_df = counts_df.append(tmp_dict, ignore_index=True)

counts_df.to_csv(counts_output)


#####json_output = output_dir+'/'+output_prefix+'_testing_counts_dictionary.json'
#####with open (json_output, 'w') as file:
	#####json.dump(counts_dictionary, file, indent='\t')

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
					print('\tThere are no variant counts for ')
					print('\t\tParser: '+parser)
					print('\t\tDatabase: '+ database)
					print('\t\tDisease: '+disease)
					print('\tThis plot will be skipped.')
					continue
				subset = tmp_df[tmp_df.columns[3:-3]].set_index('Gene')

				f = plt.figure(figsize = (12,8))
				plt.title(parser+' '+database+' '+disease)
				plt.ylabel('Count')
				subset.plot(kind="bar", stacked=True, ax = f.gca())
				plt.legend().remove()
				plt.savefig(output_dir+'/Variant_Count_Bar_Charts/'+database+'_'+disease+'_Variant_Counts_Bar_Chart.png')
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
					fig, (ax1, ax2) = plt.subplots(2,1, sharex=True, figsize = (12,8))
					ax1.spines['bottom'].set_visible(False)
					ax1.tick_params(axis='x',which='both',bottom=False)
					ax2.spines['top'].set_visible(False)
					# Set the split using the numbers calculated above
					ax2.set_ylim(0,lower_split)
					ax1.set_ylim(upper_split,upper_limit)
					subset.plot(ax=ax1,kind='bar', stacked = True)
					subset.plot(ax=ax2,kind='bar', stacked = True)
					plt.ylabel('Count')
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






## All of the file paths are added to the dictionary.
####### The next steps are to read all of the files into one dataframe
####### save the dataframe
####### Count the number of variants of each type for each gene, disease, database, parser

####### Make a file of only the ones that have normalized annotations that aren't identical to the input genomic annotations
## Count the number of non-identical ones for each gene for each disease for each databse
## Count the number of each variant type for those too
## Then find the overlap between ClinVar and the rest of them

## Count the number of ones that failed validation per each gene and why
## Then find if the ones that do overlap have the same pathogenicity rating
