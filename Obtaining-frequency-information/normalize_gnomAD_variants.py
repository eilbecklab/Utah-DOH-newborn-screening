#!/usr/bin/env python
# coding: utf-8

from argparse import ArgumentParser, FileType

args = ArgumentParser('./normalize_gnomAD_variants.py', description='''This program is designed to take a vcf file with a subset of
variants from the gnomAD vcf file and normalize the variants to HGVS format rather than VCF format. The data will be stored as a JSON
file that can be read in as a dictionary for the parsers.
Example usage: ./normalize_gnomAD_variants.py -v gene_variants.vcf -j frequency_info.json -c Chromosome_name_to_accession.csv -o frequency_info.csv''')

args.add_argument(
	'-v',
	'--vcf_file',
	help='''This is a variant vall format (VCF) file containing the frequency information about a subset of variants
from the gnomAD VCF file. No default is specified.''',
	default=None,
)

args.add_argument(
	'-j',
	'--output_json_file',
	help="""This is the name of the output json file that contains normalized vairants with the frequency information for
those variants. Default gnomad_frequencies.json""",
	default = "gnomad_frequencies.json"
)

args.add_argument(
	'-c',
	'--chr_to_accession',
	help="""This is a csv file containing information about the NCBI chromosome accession for each chromosome name. A default
will be generated if no file is provided. For formatting, see Chromosome_name_to_accession.csv as an example.""",
	default = None
)

args.add_argument(
	'-o',
	'--output_csv',
	help="""This is the name of the output csv file that contains normalized vairants with the frequency information for
those variants. No csv file will be saved by default.""",
	default = None
)


args = args.parse_args()

import csv
import pandas as pd
import re
import glob
import sys
import json

vcf_file = args.vcf_file
def print_error_message():
	print("""Welcome to normalize_gnomAD_variants.py. This program has been designed to normalize variants to HGVS format from a vcf
	file that contains a subset of variants from the gnomAD vcf file. The normalized variants will be stored to a json file along with
	frequency information for several populations.
	Example usage: ./normalize_gnomAD_variants.py -v gene_variants.vcf -j frequency_info.json""")


if not vcf_file:
	vcf_files_list = glob.glob('*.vcf')
	print()
	print_error_message()
	print()
	print("You have not specified a vcf file to be used with this program. Please specify a vcf file using the --vcf_file or -v option.")
	print()
	if len(vcf_files_list) == 0:
		print("No vcf files were found in your current working directory.")
	else:
		print("The vcf file(s) that were found in your current working directory are:")
		print(vcf_files_list)
	sys.exit(1)


chr_to_accession_file = args.chr_to_accession
if chr_to_accession_file != None:
	chr_to_accession_df = pd.read_csv(chr_to_accession_file, index_col = 0)
	chr_accession_dict = dict(zip(list(chr_to_accession_df["Chromosome"].values), list(chr_to_accession_df["Chromosome_Accession"].values)))
else:
	# If no file was put into the command line, then use this hard coded dictionary
	# It is more difficult to update than one that will come in from the command line, but is just a fallback option
	chr_accession_dict = {"1":"NC_000001.10", "2":"NC_000002.11", "3":"NC_000003.11", "4":"NC_000004.11", "5":"NC_000005.9",
							"6":"NC_000006.11", "7":"NC_000007.13", "8":"NC_000008.10", "9":"NC_000009.11", "10":"NC_000010.10",
							"11":"NC_000011.9", "12":"NC_000012.11", "13":"NC_000013.10", "14":"NC_000014.8", "15":"NC_000015.9",
							"16":"NC_000016.9", "17":"NC_000017.10", "18":"NC_000018.9", "19":"NC_000019.9", "20":"NC_000020.10",
							"21":"NC_000021.8", "22":"NC_000022.10", "X":"NC_000023.10", "Y":"NC_000024.9", "M":"NC_012920.1", "MT":"NC_012920.1"}


output_json = args.output_json_file
if ".json" not in output_json:
	output_json = output_json+".json"

output_csv = args.output_csv

skip_count = 0
with open(vcf_file, 'r') as file:
	variant_file = csv.reader(file, delimiter='\t')
	for line in variant_file:
		if line[0].startswith('#'):
			skip_count += 1
		else:
			break

vcf_df = pd.read_csv(vcf_file, index_col = None, header = None, skiprows = skip_count, sep='\t', dtype={0:'str', 1:'int', 2:'str', 3:'str', 4:'str', 5:'float', 6:'str', 7:'str'})
vcf_df.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

def obtain_frequency(df):
	info = df['INFO']
	info_dict = {}
	for piece in info.split(';'):
		second_split = piece.split('=')
		if len(second_split) == 2:
			info_dict[second_split[0]] = second_split[1]
	overall_MAF = 0
	control_MAF = 0
	African_MAF = 0
	NonFinish_Euro_MAF = 0
	Finnish_MAF = 0
	Ashkenazi_MAF = 0
	Latino_MAF = 0
	East_Asian_MAF = 0
	South_Asian_MAF = 0
	Other_Ancestry_MAF = 0
	allele_type = '-'

	if 'AF' in info_dict:
		overall_MAF = info_dict['AF']
	if 'controls_AF_raw' in info_dict:
		control_MAF = info_dict['controls_AF_raw']
	if 'AF_afr' in info_dict:
		African_MAF = info_dict['AF_afr']
	if 'AF_nfe' in info_dict:
		NonFinish_Euro_MAF = info_dict['AF_nfe']
	if 'AF_fin' in info_dict:
		Finnish_MAF = info_dict['AF_fin']
	if 'AF_asj' in info_dict:
		Ashkenazi_MAF = info_dict['AF_asj']
	if 'AF_amr' in info_dict:
		Latino_MAF = info_dict['AF_amr']
	if 'AF_eas' in info_dict:
		East_Asian_MAF = info_dict['AF_eas']
	if 'AF_sas' in info_dict:
		South_Asian_MAF = info_dict['AF_sas']
	if 'AF_oth' in info_dict:
		Other_Ancestry_MAF = info_dict['AF_oth']
	if 'allele_type' in info_dict:
		allele_type = info_dict['allele_type']

	df['Overall_MAF'] = overall_MAF
	df['Control_MAF'] = control_MAF
	df['African_MAF'] = African_MAF
	df['NonFinish_Euro_MAF'] = NonFinish_Euro_MAF
	df['Finnish_MAF'] = Finnish_MAF
	df['Ashkenazi_MAF'] = Ashkenazi_MAF
	df['Latino_MAF'] = Latino_MAF
	df['East_Asian_MAF'] = East_Asian_MAF
	df['South_Asian_MAF'] = South_Asian_MAF
	df['Other_Ancestry_MAF'] = Other_Ancestry_MAF
	df['allele_type'] = allele_type

	if len(info) > 32000:
		info = info[:32000]
	df['INFO'] = info

	return df


# Now build an HGVS style annotation from the vcf info and normalize it using the biocommons tool

import hgvs.normalizer
import hgvs.variantmapper
import hgvs.dataproviders.uta
import hgvs.assemblymapper
import hgvs.parser
import hgvs.validator
import hgvs.exceptions

hdp = hgvs.dataproviders.uta.connect()
hp = hgvs.parser.Parser()
hn = hgvs.normalizer.Normalizer(hdp)
vr = hgvs.validator.Validator(hdp=hdp)

def build_hgvs(df):
	chr_short = df['CHROM']
	if chr_short in chr_accession_dict:
		chr_accession = chr_accession_dict[chr_short]
	else:
		chr_accession = chr_short
	pos = df['POS']
	ref = df['REF']
	alt = df['ALT']
	allele_type = df['allele_type']
	hgvs_not_norm = '-'
	# The only three types that were present in the filtered vcf for newborn genes were snv, indel, del
	# There are some in the vcf listed as 'mixed', but none were present so I will need to dig further to find those
	# Duplicates are listed as insertions
	# Inversions and indels will be listed as 'mixed'
	if allele_type == 'snv':
		hgvs_not_norm = chr_accession+':g.'+str(pos)+ref+'>'+alt
	elif allele_type == 'ins':
		# All ins in gnomAD have one base in the 'ref' that is included as the first base of the 'alt' column
		pos_stop = pos+1
		inserted_sequence = alt[1:]
		hgvs_not_norm = chr_accession+':g.'+str(pos)+'_'+str(pos_stop)+'ins'+inserted_sequence
	elif allele_type == 'del':
		# The vcf format includes an extra base in the reference, so the start is one later than the listed start position
		# The stop position will not be altered as long as the length is the specified ref length (without removing the one base)
		pos_start = pos+1
		ref_length = len(ref)
		pos_stop = pos+ref_length
		hgvs_not_norm = chr_accession+':g.'+str(pos_start)+'_'+str(pos_stop)+'del'
	# the hgvs has been built, now normalize it using biocommons
	normalization_error = '-'
	hgvs_start = '-'
	hgvs_stop = '-'
	var_length = '-'
	hgvs_ref = '-'
	hgvs_alt = '-'
	normalized_hgvs = '-'
	hgvs_allele_type = '-'

	if allele_type == 'snv':
		#converted_object = None
		try:
			validated = vr.validate(hp.parse_hgvs_variant(hgvs_not_norm))
		except hgvs.exceptions.HGVSError as e:
			normalization_error = str(e)
			validated = False
		if validated:
			normalized_hgvs = hgvs_not_norm
			hgvs_allele_type = 'snv'
			var_length = 1
			hgvs_start = pos
			hgvs_stop = pos
			hgvs_ref = ref
			hgvs_alt = alt
	else:
		try:
			converted_object = hn.normalize(hp.parse_hgvs_variant(hgvs_not_norm))
		except hgvs.exceptions.HGVSError as e:
			normalization_error = str(e)
			converted_object = None
		if converted_object:
			normalized_hgvs = str(converted_object)
			try:
				hgvs_allele_type = converted_object.posedit.edit.type
			except AttributeError:
				hgvs_allele_type = '-' # this really does nothing
			if hgvs_allele_type == 'sub':
				hgvs_allele_type = 'snv' # just to be consistent with the gnomAD vcf
			# find the start and stop positions
			try:
				hgvs_start = converted_object.posedit.pos.start.base
			except AttributeError:
				hgvs_start = '-'
			try:
				hgvs_stop = converted_object.posedit.pos.end.base
			except AttributeError:
				hgvs_stop = '-'
			# Find the reference allele
			if hgvs_allele_type != 'ins':
				try:
					hgvs_ref = converted_object.posedit.edit.ref
				except AttributeError:
					hgvs_ref = '-'
			# find the alternate allele
			if hgvs_allele_type == 'dup':
				hgvs_alt = hgvs_ref+hgvs_ref
			elif hgvs_allele_type != 'inv' and hgvs_allele_type != 'del':
				try:
					hgvs_alt = converted_object.posedit.edit.alt
				except AttributeError:
					hgvs_alt = '-'
			# find the variant length
			if hgvs_allele_type == 'inv' or hgvs_allele_type == 'dup':
				var_length = str(int(hgvs_stop)-int(hgvs_start)+1)
			elif hgvs_allele_type == 'ins':
				var_length = str(len(hgvs_alt))
			else:
				try:
					var_length = converted_object.posedit.edit.ref_n
				except AttributeError:
					var_length = '-'
	df['Unnormalized HGVS'] = hgvs_not_norm
	df['Normalized HGVS Annotation'] = normalized_hgvs
	df['HGVS Start'] = hgvs_start
	df['HGVS Stop'] = hgvs_stop
	df['HGVS Ref'] = hgvs_ref
	df['HGVS Alt'] = hgvs_alt
	df['HGVS Allele Type'] = hgvs_allele_type
	df['HGVS Variant Length'] = var_length
	df['HGVS Normalization Error'] = normalization_error

	return df


vcf_df = vcf_df.apply(obtain_frequency, axis = 1)
vcf_df = vcf_df.apply(build_hgvs, axis = 1)

# Now make a dictionary, then save it as a json so that it can be used by the other programs for obtaining the frequency info
frequency_dict = {}
for index, row in vcf_df.iterrows():
	frequency_dict[row['Normalized HGVS Annotation']] = {}
	frequency_dict[row['Normalized HGVS Annotation']]['Overall_MAF'] = row['Overall_MAF']
	frequency_dict[row['Normalized HGVS Annotation']]['Control_MAF'] = row['Control_MAF']
	frequency_dict[row['Normalized HGVS Annotation']]['African_MAF'] = row['African_MAF']
	frequency_dict[row['Normalized HGVS Annotation']]['NonFinish_Euro_MAF'] = row['NonFinish_Euro_MAF']
	frequency_dict[row['Normalized HGVS Annotation']]['Finnish_MAF'] = row['Finnish_MAF']
	frequency_dict[row['Normalized HGVS Annotation']]['Ashkenazi_MAF'] = row['Ashkenazi_MAF']
	frequency_dict[row['Normalized HGVS Annotation']]['Latino_MAF'] = row['Latino_MAF']
	frequency_dict[row['Normalized HGVS Annotation']]['East_Asian_MAF'] = row['East_Asian_MAF']
	frequency_dict[row['Normalized HGVS Annotation']]['South_Asian_MAF'] = row['South_Asian_MAF']
	frequency_dict[row['Normalized HGVS Annotation']]['Other_Ancestry_MAF'] = row['Other_Ancestry_MAF']

with open(output_json, 'w')as file:
	json.dump(frequency_dict, file, indent='\t')

if output_csv:
	if ".csv" not in output_csv:
		output_csv = output_csv+'.csv'
	vcf_df.to_csv(output_csv)
