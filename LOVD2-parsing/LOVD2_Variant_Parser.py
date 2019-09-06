#!/usr/bin/env python
# coding: utf-8

import os
import re
from bs4 import BeautifulSoup
import warnings
import requests
import hgvs.normalizer
import hgvs.variantmapper
import hgvs.dataproviders.uta
import hgvs.assemblymapper
import hgvs.parser
import hgvs.validator
import hgvs.exceptions
import pandas as pd
import numpy as np
from itertools import count
#import json

from argparse import ArgumentParser, FileType
import json5
from os.path import splitext


args = ArgumentParser('./LOVD2_Variant_Parser.py', description='This program has been designed to web scrape the variant information from the Cincinnati Children\'s Medical Center (CCHMC) LOVD2 databases for variants in specific genes relating to a disase. For details on preparing input files, see README.md. Example usage: ./LOVD2_Variant_Parser.py --disease_gene_lists SCID_ny_panel.txt Metabolic_diseases_genes.txt --disease_names SCID Metabolic_Diseases --output_directory Total_Outputs --chr_to_accession Chromosome_name_to_accession.csv --transcript_info Transcript_Info_For_Dictionaries.csv')

args.add_argument(
	'--disease_gene_lists',
	nargs='+', # This tells the program that if they specify this flag, they have to give it at least one input. If they don't specify it, then the default will go in.
	help="""\
This is a list of text files containing genes associated with the disease of interest. The text files
should each contain a list of gene symbols for one given disease, one gene per line. The name of the disease will be specified
by the arguement --disease_names. """,
	default=["SCID_ny_panel.txt", "Metabolic_diseases_genes.txt"]
)

args.add_argument(
	'--disease_names',
	nargs='+',
	help="This is a list of the disease names to accompany the text files specified by the --disease_gene_lists option. If you do not use this option, the file names of the files specified in --disease_gene_lists (without the extensions) will be used as the disease names.",
	default = None
)

args.add_argument(
	'--output_directory',
	help="This is the name of the directory in which you want your output results to be stored.",
	default = "output_files"
)

args.add_argument(
	'--chr_to_accession',
	help="This is a csv file containing information about the NCBI chromosome accession for each chromosome name. A default will be generated if no file is provided. For formatting, see Chromosome_name_to_accession.csv as an example.",
	default = None
)

args.add_argument(
	'--transcript_info',
	help="This is a csv file containing information about the NCBI chromosome accession for each transcript and each gene. This file will be used for variants where a gene name or transcript accession is provided, but no chromosome is provided. For formatting, see Transcript_Info_For_Dictionaries.csv as an example.",
	default = None
)

args = args.parse_args()
# If no input was added for the disease_names, then use the disease_gene_lists to fill it in
if args.disease_names == None:
	args.disease_names = [splitext(gene_list)[0] for gene_list in args.disease_gene_lists]

disease_names = args.disease_names
input_gene_lists = []
for text_file in args.disease_gene_lists:
	input_gene_lists.append([line.rstrip('\n') for line in open(text_file)])

output_directory = args.output_directory
for disease in disease_names:
	os.makedirs(output_directory+"/CCHMC/"+disease+"/Invalid_Annotations/", exist_ok=True)

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

transcript_info_file = args.transcript_info
if transcript_info_file == None:
	# If they do not provide an input file, set these to blank dictionaries so that they are passed by in the web_scrape function
	transcript_dictionary = {}
	gene_to_transcript_dict = {}
	transcript_to_chr_accession_dict = {}
	gene_to_chr_accession_dict = {}
else:
	transcripts_df = pd.read_csv(transcript_info_file, index_col = 0)
	transcript_dictionary = dict(zip(list(transcripts_df["Transcript"].values), list(transcripts_df["Transcript_ID_Full"].values)))
	gene_to_transcript_dict = dict(zip(list(transcripts_df["Gene"].values), list(transcripts_df["Transcript_ID_Full"].values)))

	correct_chromosome_accessions = transcripts_df[transcripts_df["Chromosome_Accession"].str.contains("NC_")]
	unique_df = correct_chromosome_accessions[["Gene", "Chromosome_Accession"]].drop_duplicates()
	gene_to_chr_accession_dict = dict(zip(list(unique_df["Gene"].values), list(unique_df["Chromosome_Accession"].values)))

# These take quite some time to connect, so they have been moved below the parser arguements.
hdp = hgvs.dataproviders.uta.connect()
hp = hgvs.parser.Parser()
hn = hgvs.normalizer.Normalizer(hdp)
am = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name='GRCh37', alt_aln_method='splign', replace_reference=True)
vr = hgvs.validator.Validator(hdp=hdp)

CCHMC_home_link = "https://research.cchmc.org/LOVD2/home.php?select_db={gene}"
CCHMC_data_link = "https://research.cchmc.org/LOVD2/variants.php?select_db={gene}&action=search_all&order=Variant%2FDNA%2CASC&hide_col=&show_col=&limit=1000&search_pathogenic_=&search_Variant%2FDNA=&search_Variant%2FProtein=&search_Variant%2FRemarks=&search_Variant%2FdbSNP=&search_Patient%2FPatient%2FRace=&search_Patient%2FReference=&search_Patient%2FOrigin%2FEthnic="

def fail_request_message(status):
	if status == 301:
		response_status = "Response status is 301. This means the server is redirecting you to a different endpoint. Please double check that the website has not been moved."
	elif status == 401:
		response_status = "Response status is 401. This means the server is not authenticated, which indicates that the website has been updated to need a password."
	elif status == 400:
		response_status = "Response status is 400. This means the request was not a proper request. Please check source code."
	elif status == 403:
		response_status = "Response status is 403. This means you need special permissions to access the website. This should not occur for LOVD websites. Please double check the gene name and the website."
	elif status == 404:
		response_status = "Response status is 404. This means the request did not find the associated website. It is possible that the page for this gene simply has not been created."
	else:
		response_status = "Response status not coded. We suggest checking the error code online."
	return response_status


def web_scrape(disease_names, input_gene_lists):
	# Create an empty dictionary that can be filled with the variant information
	combined_diseases_dictionary = {}

	for disease, input_gene_list in zip(disease_names, input_gene_lists):
		combined_diseases_dictionary[disease] = {} # For each disease, add it to the dictionary
		for gene in input_gene_list:
			home_link = CCHMC_home_link.format(gene=gene)
			response = requests.get(home_link)
			if response.status_code != 200:
				response_status = fail_request_message(response.status_code)
				warnings.warn("Failed to read in web page for gene ", gene, " for the CCHMC LOVD2 Database.", sep = "")
				warnings.warn(home_link)
				warnings.warn(response_status)
				warnings.warn("Response code:", response.status_code)
				continue

			# If the loop continues to run, then it means the response status was 200 and the page downloaded correctly
			response_html = response.content.decode("utf-8")
			bs_format = BeautifulSoup(response_html, "html.parser")
			info_variable = bs_format.select('th')
			potential_error_message = info_variable[0].get_text() # If the gene is present, this should say "General information"
			if potential_error_message == 'Error: Init':
				print("The CCHMC LOVD2 database does not contain any variants for the gene "+gene+".", sep = "")
				# I could change this out for a warning, but the CCHMC database doesn't contain any info for the majority of the genes
				# so I thought this would be more appropriate to standard output instead of a warning to standard error
				continue
			# If it has made it this far, add the gene to the dictionary, with valid and invalid as options for the individual variant
			combined_diseases_dictionary[disease][gene] = {'Valid': [], 'Invalid': []}

			updated_date = "-"
			chromosome = "-"
			transcript = "-"
			transcript_new = "-"
			Chr_Accession = "-"
			for term in info_variable:
				term_text = term.get_text()
				if term_text == "Last update" :
					updated_date = term.parent.select('td')[0].get_text()
				elif term_text == "Chromosome Location":
					chromosome_search_result = re.search(r"([0-9,X,Y,M,T]+)[p|q]\d+", term.parent.select('td')[0].get_text())
					if chromosome_search_result != None:
						chromosome = chromosome_search_result.group(1)
				elif "Transcript" in term_text:
					transcript = term.parent.select('td')[0].get_text()
			# Some of the transcripts are using old versions. We need to find the current version. The only changes will be in the UTR, which may cause them to not find the location but it will not lead to erroneous normalization
			transcript_seach = re.search(r"NM_\d+", transcript)
			if transcript_seach != None:
				transcript_parent = transcript_seach.group(0)
				transcript_new = transcript_dictionary[transcript_parent]
			##
			## The ones that had no transcript listed on the home page only have one transcript
			## so I am going to use the dictionary to look up the transcript for those ones
			if transcript_new == "-":
				transcript_new = gene_to_transcript_dict[gene]

			if chromosome in chr_accession_dict:
				Chr_Accession = chr_accession_dict[chromosome]
			elif gene in gene_to_chr_accession_dict:
				Chr_Accession = gene_to_chr_accession_dict[gene]

			database_link = CCHMC_data_link.format(gene=gene)
			response = requests.get(database_link)
			if response.status_code != 200:
				response_status = fail_request_message(response.status_code)
				warnings.warn("Failed to read in variant page for gene ", gene, " for the CCHMC LOVD2 Database, despite there being a home page for the gene.", sep = "")
				warnings.warn(database_link)
				warnings.warn(response_status)
				warnings.warn(home_link)
				warnings.warn("Response code:", response.status_code)
				continue

			response_html = response.content.decode("utf-8")
			bs_format = BeautifulSoup(response_html, "html.parser")
			potential_error_message = bs_format.select('img[src="./gfx/header_variant_listings.png"]')[0].parent.parent.select('td')[0].get_text().strip()
			if "currently no public variants" in potential_error_message:
				print("The CCHMC LOVD2 database does not contain any variants for the gene "+gene+".", sep = "")
				continue
				# This should not happen because it should have happened with the first link if this were the case.
			# Now get all of the links to the individual variants
			links = []
			anchor_data = bs_format.select('a[class="data"]')
			for anchor in anchor_data:
				links.append("https://research.cchmc.org"+anchor.attrs['href'])
			# Now we have the links to all of the individual pages, one per variant
			for link in links:
				Individual_Variant_List = [] # Clear this out just to make sure no mistakes are made carrying over data from previous variant
				response = requests.get(link)
				if response.status_code != 200:
					response_status = fail_request_message(response.status_code)
					warnings.warn("Individual variant link for gene ", gene, " for the CCHMC LOVD2 Database failed to load properly.", sep = "")
					warnings.warn(link)
					warnings.warn(response_status)
					warnings.warn("Response code:", response.status_code)
					continue

				response_html = response.content.decode("utf-8")
				bs_format = BeautifulSoup(response_html, "html.parser")
				ths = bs_format.select('th[valign="top"]')
				individual_results = []
				individual_labels = []
				for th in ths:
					# This can get all of the column values for the things we care about in a list format
					# Unfortunately, some variables are missing for some variants
					individual_results.append(th.parent.td.get_text().replace(u'\xa0', u' ').replace('(View in UCSC Genome Browser, Ensembl)', ''))
				for label in bs_format.select('th[valign="top"]'):
					# This will be the labels that correspond with the values obtained above.
					individual_labels.append(label.get_text().replace(u'\xa0', u' '))
				# I have the values and their labels in separate lists at this point.
				# Now I am going to put them into a dictionary so that I can have key, value pairs
				# and I can put in a dash for anything that isn't present
				individual_dictionary = dict(zip(individual_labels, individual_results))
				terms_list = ["Race", "Reference", "Ethnic origin", "Submitter", "Allele", "Reported pathogenicity",
								"Concluded pathogenicity", "DNA change", "Protein", "Variant remarks", "dbSNP ID"]
				for term in terms_list:
					if term not in individual_dictionary:
						individual_dictionary[term] = "-"
				# Strip off the extra white space from the DNA change
				individual_dictionary["DNA change"] = individual_dictionary["DNA change"].strip()
				# Now I am going to convert the pathogenicity to the ClinVar terms
				clinvar_path_dict = {"No known pathogenicity": "Benign", "Probably no pathogenicity": "Likely benign",
									"Unknown": "Uncertain significance", "Probably pathogenic": "Likely pathogenic",
									"Pathogenic": "Pathogenic"}
				pathogenicity_clinvar = "-"
				if individual_dictionary["Reported pathogenicity"] in clinvar_path_dict:
					pathogenicity_clinvar = clinvar_path_dict[individual_dictionary["Reported pathogenicity"]]

				# The next step is to normalize the variant and take it from transcript to genomic coordinates in GRCh37
				#
				hgvs_transcript = "-"
				transcript_error = "-"
				Genomic_annotation = "-"
				Var_Type = "-"
				Position_g_start = "-"
				Position_g_stop = "-"
				Ref_allele = "-"
				Alt_allele = "-"
				Var_Length = "-"
				variant = individual_dictionary["DNA change"]
				if transcript_new != "-" and variant != "-":
					hgvs_transcript = transcript_new+":"+variant
				elif transcript != "-" and variant != "-":
					# if there is no new transcript, give it a test with the old one
					hgvs_transcript = transcript+":"+variant
					# Leave it as a dash if the DNA change is a dash, which is common
				converted_object = None
				if hgvs_transcript != "-":
					# The ClinVar parser used some regex lines to make sure that duplications, insertions and indels didn't end in a number
					# but I haven't seen any of those in the CCHMC. It can be copied from the ClinVar_Parser.py if needed
					# see lines 413-423 in ClinVar_Parser.py
					try:
						converted_object = am.c_to_g(hp.parse_hgvs_variant(hgvs_transcript))
					except hgvs.exceptions.HGVSError as e:
						transcript_error = str(e)
				if converted_object != None:
					Genomic_annotation = str(converted_object)
					# If a converted object is saved, then it can be used to obtain lots of information about the variant.
					try:
						Var_variable = converted_object.posedit.edit.type
					except AttributeError:
						Var_variable = "-"
						# Change these to be consistent with ClinVar
					if Var_variable == "ins":
						Var_Type = "Insertion"
					elif Var_variable == "del":
						Var_Type = "Deletion"
					elif Var_variable == "sub":
						Var_Type = "single nucleotide variant"
					elif Var_variable == "inv":
						Var_Type = "Inversion"
					elif Var_variable == "delins":
						Var_Type = "Indel"
					elif Var_variable == "dup":
						Var_Type = "Duplication"
					elif Var_variable == "identity":
						Var_Type = "Identity"
					# If the result was anything else, including nothing or a dash, then leave it as it was
					try:
						Position_g_start = converted_object.posedit.pos.start.base
					except AttributeError:
						Position_g_start = "-"
					try:
						Position_g_stop = converted_object.posedit.pos.end.base
					except AttributeError:
						Position_g_stop = "-"
					# The other variables depend on the variant type
					if Var_Type == "Insertion":
						Ref_allele = "-"
					else:
						try:
							Ref_allele = converted_object.posedit.edit.ref
						except AttributeError:
							Ref_allele = "-"

					if Var_Type == "Duplication" and Ref_allele != "-":
						Alt_allele = Ref_allele+Ref_allele
					elif Var_Type != "Inversion" and Var_Type != "Deletion":
						try:
							Alt_allele = converted_object.posedit.edit.alt
						except AttributeError:
							Alt_allele = "-"

					if Var_Type == "Insertion" or Var_Type == "Inversion" or Var_Type == "Duplication":
						if Position_g_stop != "-" and Position_g_start != "-":
							Var_Length = str(int(Position_g_stop)-int(Position_g_start)+1)
					elif Var_Type == "single nucleotide variant":
						Var_Length = 1
					else:
						try:
							Var_Length = converted_object.posedit.edit.ref_n
						except AttributeError:
							Var_Length = "-"
					# It turns out some of the reference and alternate alleles from the large variants (>30,000 bp) return valid answers, but they are too big to be opened in Excel
					# Here I am just going to change those to something Excel can handle
					if len(Ref_allele) > 1000:
						Ref_allele = "Reference allele too long to include in table"
					if len(Alt_allele) > 1000:
						Alt_allele = "Alternate allele too long to include in table"
				# Some of the variants that didn't pass validation still have some data that can be taken from the transcript info
				if Var_Type == "-":
					if ">" in hgvs_transcript:
						Var_Type = "single nucleotide variant"
					elif "delins" in hgvs_transcript:
						Var_Type = "Indel"
					elif "del" in hgvs_transcript and "ins" in hgvs_transcript:
						Var_Type = "Indel" # this will pick up ones like 123delCATinsATTG
					elif "del" in hgvs_transcript:
						Var_Type = "Deletion"
					elif "ins" in hgvs_transcript:
						Var_Type = "Insertion"
					elif "dup" in hgvs_transcript:
						Var_Type = "Duplication"
					elif "inv" in hgvs_transcript:
						Var_Type = "Inversion"

				if Var_Type == "single nucleotide variant":
					Var_Length = 1
					if Ref_allele == "-" or Alt_allele == "-": # Only run the regex loop if something is missing
						ref_allele_search = re.search(r"\d+([A,C,T,G])>([A,C,T,G])", hgvs_transcript)
						if ref_allele_search != None:
							if Ref_allele == "-":
								Ref_allele = ref_allele_search.group(1)
							if Alt_allele == "-":
								Alt_allele = ref_allele_search.group(2)
				elif Var_Type == "Indel" or Var_Type == "Insertion":
					# Can't obtain reference allele without looking it up, and it already failed so ref will stay as a dash
					if Alt_allele == "-": # Only run the regex loop if Alt_allele is missing
						alt_allele_search = re.search(r"ins\(?([A,C,T,G]+)", hgvs_transcript) # Some of these have the inserted variant in parentheses
						if alt_allele_search != None:
							Alt_allele = alt_allele_search.group(1)
					if Var_Type == "Insertion" and Var_Length == "-" and Alt_allele != "-": # Var_Length for indels is inconsistent in ClinVar, so I will leave these as a dash
						Var_Length = len(Alt_allele)
				elif Var_Type == "Duplication" or Var_Type == "Deletion" or Var_Type == "Inversion":
					# Can't get reference or alternate for any of these
					if Var_Length == "-": # Only run the regex loop if something is missing
						multiple_location_search = re.search(r"\d+_\d+", input_for_filling_blanks)
						if multiple_location_search == None:
							Var_Length = 1
						elif Position_g_stop != "-" and Position_g_start != "-":
							Var_Length = int(Position_g_stop) - int(Position_g_start) + 1

				# If the genomic normalization did not work, find out the reason from the failure from the error messages
				failure_reason = "-"
				if Genomic_annotation == "-":
					# Now find the failure reason if it didn't normalize properly, use the transcript_error and genomic_error messages
					# Use the genomic error message first, because many have transcript messages that will say "out of region"
					microsatellite_search = re.search(r"\[\d+\]", transcript_error)
					inserted_unknown_search = re.search(r"ins\(?.?(\d+|\?)", transcript_error)
					if "-:" in transcript_error:
						failure_reason = "No variant information provided"
					elif microsatellite_search != None:
						failure_reason = "Microsatellite"
					elif inserted_unknown_search != None:
						# May also contain an unknown breakpoint
						if "(?_" in transcript_error or "_?)" in transcript_error:
							failure_reason = "Inserted unknown sequence, unknown breakpoint"
						else:
							failure_reason = "Inserted unknown sequence"
					elif "?" in transcript_error:
						failure_reason = "Unknown breakpoint"
					elif ")_(" in transcript_error:
						failure_reason = "Compound variant"
					else:
						failure_reason = "No definitive failure reason detected, likely compound variant with nontraditional formatting"
				#
				# Now put the variant information obtained into a list, and append it to the dictionary in the correct location
				if Genomic_annotation == "-":
					variant_category = "Invalid"
				else:
					variant_category = "Valid"
				Individual_Variant_List = ['GRCh37', chromosome, Position_g_start, Position_g_stop, Ref_allele,
											Alt_allele, Genomic_annotation, Genomic_annotation, Var_Type, Var_Length,
											pathogenicity_clinvar, "-", individual_dictionary["Allele"], "-", gene, gene,
											"-", transcript_new,  individual_dictionary["DNA change"], hgvs_transcript, "-",
											individual_dictionary["Protein"], "-", Chr_Accession, "-", "-", "-", "CCHMC", "-",
											"-", "-", individual_dictionary["Submitter"], updated_date,
											individual_dictionary["Reported pathogenicity"], individual_dictionary["Concluded pathogenicity"],
											individual_dictionary["Race"], individual_dictionary["Ethnic origin"],
											individual_dictionary["Reference"], individual_dictionary["Variant remarks"],
											transcript, individual_dictionary["dbSNP ID"], transcript_error, "-"]
				if variant_category == "Invalid":
					Individual_Variant_List.append(failure_reason)
				# Now that the individual_variant_list has been created, append it to the dictionary where it belongs
				combined_diseases_dictionary[disease][gene][variant_category].append(Individual_Variant_List)
	# Now that everything has been parsed and put into the super dictionary, make the labels, put each list into a dataframe and save them out to csv files
	column_labels = ["Genome Assembly", "Chr", "Position Start", "Position Stop", "Ref", "Alt",
						"Genomic Annotation", "HGVS Normalized Genomic Annotation", "Variant Type",
						"Variant Length", "Pathogenicity", "Disease", "Genetic Origin", "Inheritance Pattern",
						"Affected Genes" , "Gene Symbol", "Compound Het Status", "Transcript",
						"Transcript Notation", "HGVS Transcript Notation","Protein Accession",
						"Protein Notation", "HGVS Protein Annotation", "Chr Accession", "VCF Pos",
						"VCF Ref", "VCF Alt", "Database", "ClinVar Accession", "Review Status",
						"Star Level", "Submitter", "Edited Date", "CCHMC Reported Pathogenicity",
						"CCHMC Concluded Pathogenicity", "Race", "Ethnic Origin", "Reference Paper",
						"Variant Remarks", "CCHMC Reported Transcript", "dbSNP ID", "Transcript Normalization Failure Message",
						"Genomic Normalization Failure Message"]
	column_labels_invalid = ["Genome Assembly", "Chr", "Position Start", "Position Stop", "Ref", "Alt",
						"Genomic Annotation", "HGVS Normalized Genomic Annotation", "Variant Type",
						"Variant Length", "Pathogenicity", "Disease", "Genetic Origin", "Inheritance Pattern",
						"Affected Genes" , "Gene Symbol", "Compound Het Status", "Transcript",
						"Transcript Notation", "HGVS Transcript Notation","Protein Accession",
						"Protein Notation", "HGVS Protein Annotation", "Chr Accession", "VCF Pos",
						"VCF Ref", "VCF Alt", "Database", "ClinVar Accession", "Review Status",
						"Star Level", "Submitter", "Edited Date", "CCHMC Reported Pathogenicity",
						"CCHMC Concluded Pathogenicity", "Race", "Ethnic Origin", "Reference Paper",
						"Variant Remarks", "CCHMC Reported Transcript", "dbSNP ID", "Transcript Normalization Failure Message",
						"Genomic Normalization Failure Message", "HGVS Normalization Failure Reason"]
	for disease_key in combined_diseases_dictionary.keys():
		for gene_key in combined_diseases_dictionary[disease_key].keys():
			for output_file_key in combined_diseases_dictionary[disease_key][gene_key].keys():
				individual_results_list = combined_diseases_dictionary[disease_key][gene_key][output_file_key]
				if len(individual_results_list) > 0:
					individual_df = pd.DataFrame(individual_results_list)
					if output_file_key == "Invalid":
						individual_df.columns = column_labels_invalid
						individual_df.to_csv(output_directory+"/CCHMC/"+disease_key+"/Invalid_Annotations/"+gene_key+"_CCHMC_InvalidResults.csv")
					elif output_file_key == "Valid":
						individual_df.columns = column_labels
						individual_df.to_csv(output_directory+"/CCHMC/"+disease_key+"/"+gene_key+"_CCHMC_Results.csv")


# Now run the function
web_scrape(disease_names, input_gene_lists)
