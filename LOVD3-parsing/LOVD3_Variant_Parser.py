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


from argparse import ArgumentParser, FileType
import json5
from os.path import splitext


args = ArgumentParser('./LOVD3_Variant_Parser.py', description='This program has been designed to web scrape the LOVD3 databases to obtain the information listed for variants in specific genes relating to a disase. For details on preparing input files, see LOVD3_README.md. Example usage: ./LOVD3_Variant_Parser.py --config_file LOVD3_Databases.json --disease_gene_lists SCID_ny_panel.txt Metabolic_diseases_genes.txt --disease_names SCID Metabolic_Diseases --output_directory Total_Outputs --chr_to_accession Chromosome_name_to_accession.csv --transcript_info Transcript_Info_For_Dictionaries.csv')

args.add_argument(
	'--config_file',
	type=FileType('r'),
	help="This is a config file in JSON format with the information about the databases to be read in. For an example of formatting, see the file LOVD3_Databases.json",
	default="LOVD3_Databases.json",
)

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
if args.disease_names == None:
	args.disease_names = [splitext(gene_list)[0] for gene_list in args.disease_gene_lists]


databases_list = json5.load(args.config_file)


disease_names = args.disease_names
input_gene_lists = []
for text_file in args.disease_gene_lists:
	input_gene_lists.append([line.rstrip('\n') for line in open(text_file)])


output_directory = args.output_directory
for variable in databases_list:
	database_name = variable["name"]
	for disease in disease_names:
		os.makedirs(output_directory+"/"+database_name+"/"+disease+"/Invalid_Annotations/", exist_ok=True)

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
	transcript_to_chr_accession_dict = {}
	gene_to_chr_accession_dict = {}
else:
	transcripts_df = pd.read_csv(transcript_info_file, index_col = 0)
	correct_chromosome_accessions = transcripts_df[transcripts_df["Chromosome_Accession"].str.contains("NC_")]
	transcript_to_chr_accession_dict = dict(zip(list(correct_chromosome_accessions["Transcript_ID_Full"].values), list(correct_chromosome_accessions["Chromosome_Accession"].values)))
	unique_df = correct_chromosome_accessions[["Gene", "Chromosome_Accession"]].drop_duplicates()
	gene_to_chr_accession_dict = dict(zip(list(unique_df["Gene"].values), list(unique_df["Chromosome_Accession"].values)))

# These take quite some time to connect, so they have bene moved below the parser arguements.
hdp = hgvs.dataproviders.uta.connect()
hp = hgvs.parser.Parser()
hn = hgvs.normalizer.Normalizer(hdp)
am = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name='GRCh37', alt_aln_method='splign', replace_reference=True)
vr = hgvs.validator.Validator(hdp=hdp)

def get_data_from_tables(database_data_link, gene, transcript_accession, LOVD_id, database_name):
	# need database_data_link, LOVD ID, gene, page_num, transcript_accessions[i]
	output_results = []
	for page_number in count(1):
		data_link = database_data_link.format(gene=gene, LOVD_id=LOVD_id, page_num = str(page_number))
		response = requests.get(data_link)
		if response.status_code != 200:
			warnings.warn("The link to the data for transcript "+transcript_accession+" in the database "+database_name+" did not work. This might indicate a change in the web address. Please double check the data link.")
			warnings.warn(data_link)
			# re-scrape the first page to get the headers and return whatever data we have up to this point
			# If the count is 1, then there is no data to really return
			if count == 1:
				header_list = []
			else:
				data_link = database_data_link.format(gene=gene, LOVD_id=LOVD_id, page_num = "1")
				response = requests.get(data_link)
				response_html = response.content.decode("utf-8")
				bs_format = BeautifulSoup(response_html, "html.parser")
				header_tmp_list = [tag.get_text().strip() for tag in bs_format.select('div')]
				header_list = [i.replace(u'\xa0', u' ') for i in header_tmp_list[1:]]

			return output_results, header_list # return what you have so far if one of the pages fails to open properly
		response_html = response.content.decode("utf-8")
		bs_format = BeautifulSoup(response_html, "html.parser")
		# From this get the variant information from each line in the database (split by newline character)
		tmp_variable = bs_format.select('tr.data')
		results_list = [tag.get_text().replace("\r","").strip().split("\n") for tag in tmp_variable]
		output_results = output_results+results_list
		if len(results_list) == 0:
			if count == 1:
				# this means the first page pulled up but had no results.
				header_list = []
				return output_results, header_list
			# This will only happen if there are a multiple of 1000 variants, meaning there wasn't really a next page to check
			# Re-scrape the first page to get the headers. Not ideal to scrape again, but current bs_format will not have any headers.
			data_link = database_data_link.format(gene=gene, LOVD_id=LOVD_id, page_num = "1")
			response = requests.get(data_link)
			response_html = response.content.decode("utf-8")
			bs_format = BeautifulSoup(response_html, "html.parser")
			header_tmp_list = [tag.get_text().strip() for tag in bs_format.select('div')]
			header_list = [i.replace(u'\xa0', u' ') for i in header_tmp_list[1:]]
			return output_results, header_list
		elif len(results_list) < 1000:
			header_tmp_list = [tag.get_text().strip() for tag in bs_format.select('div')]
			header_list = [i.replace(u'\xa0', u' ') for i in header_tmp_list[1:]]
			return output_results, header_list


def normalize_variant(row):
	transcript_notation = row["HGVS Transcript Notation"]
	Genomic_annotation = row["Genomic Annotation"]
	transcript_error = "-"
	genomic_error = "-"
	Var_Type = "-"
	Position_g_start = "-"
	Position_g_stop = "-"
	Var_Length = "-"
	Ref_allele = "-"
	Alt_allele = "-"
	Genomic_Normalized = "-"

	try:
		# this should work for the majority of the variants
		converted_object = am.c_to_g(hp.parse_hgvs_variant(transcript_notation))
	except hgvs.exceptions.HGVSError as e:
		transcript_error = str(e)
		# If this doesn't run using the transcript_notation, then we need to run either the normalization
		# or the validation on the genomic annotation

		if ">" in Genomic_annotation:
			# These are SNVs, they just need to run validation
			converted_object = None # Needs to be present later in loop
			try:
				validated = vr.validate(hp.parse_hgvs_variant(Genomic_annotation))
			except hgvs.exceptions.HGVSError as e:
				genomic_error = str(e)
				validated = False
			if validated == True:
				Genomic_Normalized = Genomic_annotation
				Var_Type = "single nucleotide variant"
				Var_Length = 1
				pos_search = re.search(r"g.(\d+)", Genomic_Normalized)
				if pos_search != None:
					Position_g_start = pos_search.group(1)
				Position_g_stop = Position_g_start # May be a dash
				ref_and_alt_allele_search = re.search(r"g.\d+([A,C,T,G])>([A,C,T,G])", Genomic_Normalized)
				if ref_and_alt_allele_search != None:
					Ref_allele = ref_and_alt_allele_search.group(1)
					Alt_allele = ref_and_alt_allele_search.group(2)
			# No else statement needed, all will stay as dashes if failed validation (incorrect base)


		else:
			# for everything that isn't a SNV, it needs to be normalized
			try:
				converted_object = hn.normalize(hp.parse_hgvs_variant(Genomic_annotation))
			except hgvs.exceptions.HGVSError as e:
				genomic_error = str(e)
				Genomic_Normalized = "-"
				converted_object = None
	# Now use the converted object to fill everything in that you possibly can

	if converted_object != None:
		# If a converted object is saved, then it can be used to obtain lots of information about the variant.
		Genomic_Normalized = str(converted_object)
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
		elif Var_variable == None:
			Var_Type = "-"
		else:
			Var_Type = Var_variable # Anything that was not in this list will just end up putting the variable so it can be seen in the output file

		try:
			Position_g_start = converted_object.posedit.pos.start.base
		except AttributeError:
			Position_g_start = "-"
		try:
			Position_g_stop = converted_object.posedit.pos.end.base
		except AttributeError:
			Position_g_stop = "-"

		if Var_Type == "Insertion":
			Ref_allele = "-"
		else:
			try:
				Ref_allele = converted_object.posedit.edit.ref
			except AttributeError:
				Ref_allele = "-"

		if Var_Type == "Duplication":
			Alt_allele = Ref_allele+Ref_allele
		elif Var_Type == "Inversion" or Var_Type == "Deletion":
			Alt_allele = "-"
		else:
			try:
				Alt_allele = converted_object.posedit.edit.alt
			except AttributeError:
				Alt_allele = "-"

		if Var_Type == "Insertion" or Var_Type == "Inversion" or Var_Type == "Duplication":
			Var_Length = str(int(Position_g_stop)-int(Position_g_start)+1)
		else:
			Var_Length = converted_object.posedit.edit.ref_n
	# At this point if something is missing (predominantly just the ones that didn't pass validation)
	# Then we can try to fill in the blanks from what is provided
	input_for_filling_blanks = "-"
	using_transcript = False # I can't get the variant position if I use the transcript for info
	if Genomic_Normalized != "-":
		input_for_filling_blanks = Genomic_Normalized
	elif "g.?" not in Genomic_annotation:
		input_for_filling_blanks = Genomic_annotation
	elif "c.?" not in transcript_notation:
		input_for_filling_blanks = transcript_notation
		using_transcript = True
	if input_for_filling_blanks != "-":
		if Var_Type == "-":
			if ">" in input_for_filling_blanks:
				Var_Type = "single nucleotide variant"
			elif "delins" in input_for_filling_blanks:
				Var_Type = "Indel"
			elif "del" in input_for_filling_blanks and "ins" in input_for_filling_blanks:
				Var_Type = "Indel" # this will pick up ones like 123delCATinsATTG
			elif "del" in input_for_filling_blanks:
				Var_Type = "Deletion"
			elif "ins" in input_for_filling_blanks:
				Var_Type = "Insertion"
			elif "dup" in input_for_filling_blanks:
				Var_Type = "Duplication"
			elif "inv" in input_for_filling_blanks:
				Var_Type = "Inversion"

		if Position_g_start == "-" and not using_transcript: # Only do the regex loop if start position is missing
			pos_start_search = re.search(r"g.(\d+)", input_for_filling_blanks)
			if pos_start_search != None:
				Position_g_start = pos_start_search.group(1)

		## Filling in the rest of the information about the variant depends on what the variant type is

		if Var_Type == "single nucleotide variant":
			if Position_g_stop == "-":
				Position_g_stop = Position_g_start
			Var_Length = 1
			if Ref_allele == "-" or Alt_allele == "-": # Only run the regex loop if something is missing
				ref_allele_search = re.search(r"g.\d+([A,C,T,G])>([A,C,T,G])", input_for_filling_blanks)
				if ref_allele_search != None:
					if Ref_allele == "-":
						Ref_allele = ref_allele_search.group(1)
					if Alt_allele == "-":
						Alt_allele = ref_allele_search.group(2)
		elif Var_Type == "Indel":
			# Var_Length for indels is inconsistent in ClinVar, so I will leave these as a dash
			# Can't obtain reference allele without looking it up, which would take a request and be slow so ref will stay as a dash
			if Alt_allele == "-": # Only run the regex loop if Alt_allele is missing
				alt_allele_search = re.search(r"ins\(?([A,C,T,G]+)", input_for_filling_blanks) # Some of these have the inserted variant in parentheses
				if alt_allele_search != None:
					Alt_allele = alt_allele_search.group(1)
			if Position_g_stop == "-" and not using_transcript: # Only run the regex loop if Position_g_stop is missing
				multiple_location_search = re.search(r"\d+_(\d+)del", input_for_filling_blanks)
				if multiple_location_search == None:
					Position_g_stop = Position_g_start
				else:
					Position_g_stop = multiple_location_search.group(1)
		elif Var_Type == "Insertion":
			# Reference will always be a dash
			# I think the stop should always be start plus one, but in case it is funky I am going to use the regex anyway
			if Position_g_stop == "-" or Alt_allele == "-": # Only run the regex loop if something is missing
				info_search = re.search(r"_(\d+)ins\(?([A,C,T,G]+)", input_for_filling_blanks) # Some of these have the inserted variant in parenthese
				if info_search != None:
					if Position_g_stop == "-" and not using_transcript:
						Position_g_stop = info_search.group(1)
					if Alt_allele == "-":
						Alt_allele = info_search.group(2)
				# If the stop position is still a dash because the info search didn't return a result
				if Position_g_stop == "-" and not using_transcript: # This happens when they end in an insertion of a number or a question mark
					info_search = re.search(r"_(\d+)ins", input_for_filling_blanks)
					if info_search != None:
						Position_g_stop = info_search.group(1)
			if Alt_allele != "-" and Var_Length == "-":
				Var_Length = len(Alt_allele)
		elif Var_Type == "Duplication" or Var_Type == "Deletion" or Var_Type == "Inversion":
			# Can't get reference or alternate for any of these
			if Position_g_stop == "-": # Only run the regex loop if something is missing
				multiple_location_search = re.search(r"\d+_\d+", input_for_filling_blanks)
				if multiple_location_search == None:
					Position_g_stop = Position_g_start # If only one of these was a dash, the other one will be replaced too but should be replaced with the same thing
					Var_Length = 1
				elif not using_transcript:
					stop_search = re.search(r"_(\d+)[a-z]+", input_for_filling_blanks) # This only yields valid results if not using the transcript info
					if stop_search != None:
						Position_g_stop = stop_search.group(1)
			if Var_Length == "-" and Position_g_stop != "-" and Position_g_start != "-":
				Var_Length = int(Position_g_stop) - int(Position_g_start) + 1

	# It turns out some of the reference and alternate alleles from the large variants (>30,000 bp) return valid answers, but they are too big to be opened in Excel
	# Here I am just going to change those to something Excel can handle
	if len(Ref_allele) > 1000:
		Ref_allele = "Reference allele too long to include in table"
	if len(Alt_allele) > 1000:
		Alt_allele = "Alternate allele too long to include in table"

	row["HGVS Normalized Genomic Annotation"] = Genomic_Normalized
	row["Transcript Normalization Failure Message"] = transcript_error
	row["Genomic Normalization Failure Message"] = genomic_error
	row['Variant Type'] = Var_Type
	row['Position Start'] = Position_g_start
	row['Position Stop'] = Position_g_stop
	row['Ref'] = Ref_allele
	row['Alt'] = Alt_allele
	row['Variant Length'] = Var_Length

	return row


def find_failure_reason(row):
	transcript_err = row["Transcript Normalization Failure Message"]
	genome_err = row["Genomic Normalization Failure Message"]
	reason = "-" # None should return this, but just in case
	if ":g.?:" not in genome_err:
		microsatellite_search = re.search(r"\[\d+\]", genome_err)
		inserted_unknown_search = re.search(r"ins\(?.?(\d+|\?)", genome_err)
		if "does not agree with reference sequence" in genome_err:
			reason = "Incorrect reference base"
		elif microsatellite_search != None:
			reason = "Microsatellite"
		elif inserted_unknown_search != None:
			# May also contain an unknown breakpoint
			if "(?_" in genome_err or "_?)" in genome_err:
				reason = "Inserted unknown sequence, unknown breakpoint"
			else:
				reason = "Inserted unknown sequence"
		elif "?" in genome_err:
			reason = "Unknown breakpoint"
		elif ")_(" in genome_err:
			reason = "Compound variant"
		else:
			reason = "No definitive failure reason detected, likely compound variant with nontraditional formatting"
	else:
		# These ones have no information for genomic, but may have something for the transcript annotation
		microsatellite_search = re.search(r"\[\d+\]", transcript_err)
		inserted_unknown_search = re.search(r"ins\(?.?(\d+|\?)", transcript_err)
		if ":c.?:" in transcript_err:
			reason = "No variant information provided"
		elif microsatellite_search != None:
			reason = "Microsatellite"
		elif inserted_unknown_search != None:
			# May also contain an unknown breakpoint
			if "(?_" in transcript_err or "_?)" in transcript_err:
				reason = "Inserted unknown sequence, unknown breakpoint"
			else:
				reason = "Inserted unknown sequence"
		elif "?" in transcript_err:
			reason = "Unknown breakpoint"
		elif ")_(" in transcript_err:
			reason = "Compound variant"
		else:
			reason = "No definitive failure reason detected, likely compound variant with nontraditional formatting"
	row["HGVS Normalization Failure Reason"] = reason
	return row

def web_scrape(disease_names, input_gene_lists, databases_list):
	for database in databases_list:
		database_name = database["name"]
		database_home_link = database["home link"] # will still have variables in it
		database_data_link = database["data link"] # will still have variables in it
		database_contains_unique = database["contains unique variants"] # will be boolean true or false
		for num, input_gene_list in enumerate(input_gene_lists):
			disease_name = disease_names[num]
			for gene in input_gene_list:
				home_link = database_home_link.format(gene=gene)
				response = requests.get(home_link)
				if response.status_code != 200:
					if response.status_code == 301:
						response_status = "Response status is 301. This means the server is redirecting you to a different endpoint. Please double check that the website has not been moved."
					elif response.status_code == 401:
						response_status = "Response status is 401. This means the server is not authenticated, which indicates that the website has been updated to need a password."
					elif response.status_code == 400:
						response_status = "Response status is 400. This means the request was not a proper request. Please check source code."
					elif response.status_code == 403:
						response_status = "Response status is 403. This means you need special permissions to access the website. This should not occur for LOVD websites. Please double check the gene name and the website."
					elif response.status_code == 404:
						response_status = "Response status is 404. This means the request did not find the associated website. It is possible that the page for this gene simply has not been created."
					else:
						response_status = "Response status not coded. We suggest checking the error code online."

					warnings.warn("Failed to read in web page for gene ", gene, " for the Human Variome Database.", sep = "")
					warnings.warn(link)
					warnings.warn(response_status)
					warnings.warn("Response code:", response.status_code)
					continue

				# If the loop continues to run, then it means the response status was 200 and the page downloaded correctly
				response_html = response.content.decode("utf-8")
				bs_format = BeautifulSoup(response_html, "html.parser")

				no_such_gene_variable = bs_format.select('td[valign="middle"]')
				if len(no_such_gene_variable) > 0:
					warning_message = no_such_gene_variable[0].get_text()
					if warning_message == "No such ID!":
						warnings.warn("Gene "+gene+" is not present in the Database "+database_name+".")
						continue
					else:
						warnings.warn("Gene "+gene+" has a message that is not standard. This might mean that the gene is no present in the database "+database_name+".")
						warnings.warn("Warning message: "+warning_message)
						continue
						# I am going to use continue here because this page is formatted differently from the ones that have a message
				# If the loop continues to run for the gene, it means that the page is present and formatted properly
				updated_date = "-"
				for term in bs_format.select('th'):
					if "updated" in term.get_text() :
						updated_date = term.parent.select('td')[0].get_text()
				if updated_date == "-": # If no updated date is present, use the created date
					for term in bs_format.select('th'):
						if "created" in term.get_text() :
							updated_date = term.parent.select('td')[0].get_text()


				chromosome = bs_format.select('tr[class="data"]')[0].select('td')[1].get_text()
				header_variable = bs_format.select('th[class="order"]')
				header_list = []
				for child in header_variable:
					header_list.append(child.get_text().strip().replace('\xa0',' '))
				LOVD_id_index = header_list.index("ID")
				transcript_accession_index = header_list.index("NCBI ID")
				protein_accession_index = header_list.index("NCBI Protein ID")

				LOVD_ids = []
				transcript_accessions = []
				protein_accessions = []

				tr_data = bs_format.select('tr[class="data"]')
				## If there is more than one transcript, there will be more than one child in the tr data
				for child in tr_data:
					LOVD_ids.append(child.select('td')[LOVD_id_index].get_text())
					transcript_accessions.append(child.select('td')[transcript_accession_index].get_text())
					protein_accessions.append(child.select('td')[protein_accession_index].get_text())
				# Look up the chromosome accession by gene. If that isn't there, try the transcript and then try working from the chromsome web scraped
				# If none of these work, use a dash
				if gene in gene_to_chr_accession_dict:
					Chr_accession = gene_to_chr_accession_dict[gene]
				elif transcript_accessions[0] in transcript_to_chr_accession_dict:
					Chr_accession = transcript_to_chr_accession_dict[transcript_accessions[0]]
				elif chromosome in chr_accession_dict:
					Chr_accession = chr_accession_dict[chromosome]
				else:
					Chr_accession = "-"
				###
				###	Now that we have the important information from the home page,
				###	we can move on to getting the info from the pages from the individual transcripts
				###
				for i, LOVD_id in enumerate(LOVD_ids):
					transcript_accession = transcript_accessions[i]

					results_list, label_list = get_data_from_tables(database_data_link, gene, transcript_accession, LOVD_id, database_name)
					if len(results_list) == 0:
						warnings.warn("Transcript "+transcript_accession+" in the database "+database_name+" has a page, but does not contain any variants.")
						continue

					# Some of the results have a \n (newline) and a \r (carriage return from Excel) in the middle of a column, which will make it
					# This will shift a couple of the variables at the end (columns after the comments) but those are more comment fields that don't affect the variant information
					necessary_length = len(label_list)
					for i, sub_list in enumerate(results_list):
						if len(sub_list) > necessary_length:
							tmp_list = sub_list[:necessary_length]
							results_list[i] = tmp_list



					gene_df = pd.DataFrame(results_list)
					gene_df.columns = label_list
					gene_df['Database'] = database_name
					gene_df["Transcript"] = transcript_accession
					gene_df["Chr"] = chromosome
					gene_df["Edited Date"] = updated_date
					gene_df["Chr Accession"] = Chr_accession
					column_names = gene_df.columns
					if "DNA change (cDNA)" in column_names:
						gene_df = gene_df.rename(index = str, columns = {"DNA change (cDNA)": "Transcript Notation"})
					else:
						gene_df["Transcript Notation"] = "-"
					if "DNA change (genomic)" in column_names:
						gene_df = gene_df.rename(index = str, columns = {"DNA change (genomic)": "DNA change (genomic) (hg19)"})


					if database_contains_unique == False: # This is one of the things coded in the JSON config file
						counts_df = pd.DataFrame(gene_df.groupby("Transcript Notation").size())
						counts_df.columns = ["Reported"]
						gene_df = pd.merge(gene_df, counts_df, on = "Transcript Notation").drop_duplicates(["Transcript Notation", "DNA change (genomic) (hg19)"], keep='first')

					### Not all of the LOVD databases give the same output, so I am going to fill some of those in with a dash if they are not present

					# The column_names list needs to be remade to match the new column names
					column_names = gene_df.columns
					if "DNA change (genomic) (hg19)" not in column_names:
						gene_df["DNA change (genomic) (hg19)"] = "-"

					if "Effect" not in column_names:
						gene_df["Effect"] = "-"
					if "Exon" not in column_names:
						gene_df["Exon"] = "-"
					if "Reported" not in column_names:
						gene_df["Reported"] = "-"
					if "DB-ID" not in column_names:
						gene_df["DB-ID"] = "-"
					if "dbSNP ID" not in column_names:
						gene_df["dbSNP ID"] = "-"
					if "Protein" not in column_names:
						gene_df["Protein"] = "-"
					if "Published as" not in column_names:
						gene_df["Published as"] = "-"
					if "Variant remarks" not in column_names:
						gene_df["Variant remarks"] = "-"
					if "Reference" not in column_names:
						gene_df["Reference"] = "-"
					if "Frequency" not in column_names:
						gene_df["Frequency"] = "-"


					if "ClassClinical" in column_names:
						gene_df = gene_df.rename(index = str, columns = {"ClassClinical": "Pathogenicity"})
					else:
						gene_df["Pathogenicity"] = "-"
					if "ClinVar ID" in column_names:
						gene_df = gene_df.rename(index = str, columns = {"ClinVar ID": "ClinVar Accession"})
					else:
						gene_df["ClinVar Accession"] = "-"
					if "Origin" in column_names:
						gene_df = gene_df.rename(index = str, columns = {"Origin": "Genetic Origin"})
					else:
						gene_df["Genetic Origin"] = "-"
					if "Protein" in column_names:
						gene_df = gene_df.rename(index = str, columns = {"Protein": "Protein Notation"})
					else:
						gene_df["Protein Notation"] = "-"
					if "Owner" in column_names:
						gene_df = gene_df.rename(index = str, columns = {"Owner": "Submitter"})
					else:
						gene_df["Submitter"] = "-"

					###
					###	Now everything that I can get out of the database is present in a pandas format
					###	The next step is to normalize or validate the variant
					###
					gene_df["Genomic Annotation"] = Chr_accession+":"+gene_df["DNA change (genomic) (hg19)"]
					gene_df["HGVS Transcript Notation"] = transcript_accession+":"+gene_df["Transcript Notation"]
					gene_df = gene_df.apply(normalize_variant, axis = 1)

					## Now I am going to try to format this to match the output from ClinVar as close as possible
					## but most of the things I am adding will either be a dash or the same word for every column
					gene_df["Genome Assembly"] = "GRCh37"
					gene_df["Disease"] = "-"
					gene_df["Inheritance Pattern"] = "-"
					gene_df["Affected Genes"] = gene
					gene_df["Gene Symbol"] = gene
					gene_df["Protein Accession"] = "-"
					gene_df["HGVS Protein Annotation"] = "-"
					gene_df['Compound Het Status'] = "-"
					gene_df["VCF Pos"] = "-"
					gene_df["VCF Ref"] = "-"
					gene_df["VCF Alt"] = "-"
					gene_df["Review Status"] = "-"
					gene_df["Star Level"] = "-"

					gene_df = gene_df[['Genome Assembly', 'Chr', 'Position Start', 'Position Stop', 'Ref', 'Alt',
										'Genomic Annotation', 'HGVS Normalized Genomic Annotation', 'Variant Type',
										'Variant Length', 'Pathogenicity', 'Disease', 'Genetic Origin', 'Inheritance Pattern',
										'Affected Genes', 'Gene Symbol', 'Compound Het Status', 'Transcript',
										 'Transcript Notation', 'HGVS Transcript Notation', 'Protein Accession',
										'Protein Notation', 'HGVS Protein Annotation', 'Chr Accession', 'VCF Pos', 'VCF Ref',
										'VCF Alt', 'Database', 'ClinVar Accession', 'Review Status', 'Star Level',
										'Submitter', 'Edited Date', 'DNA change (genomic) (hg19)', 'Effect',
										 'Exon','Reported', 'DB-ID', 'dbSNP ID', 'Published as', 'Variant remarks',
										'Reference', 'Frequency', 'Transcript Normalization Failure Message', 'Genomic Normalization Failure Message']]
					failed_HGVS_df = gene_df[gene_df['HGVS Normalized Genomic Annotation'] == "-"]
					successful_HGVS_df = gene_df[gene_df['HGVS Normalized Genomic Annotation'] != "-"]
					if len(successful_HGVS_df) > 0:
						# There have been a few of these that only have a couple of annotations and all of them fail normalization
						successful_HGVS_df.to_csv(output_directory+"/"+database_name+"/"+disease_name+"/"+gene+"_"+transcript_accession+"_"+database_name+"_results.csv")
					if len(failed_HGVS_df) > 0:
						## There are so few that are failing currently that I think I will just put this out as one csv
						failed_HGVS_df = failed_HGVS_df.apply(find_failure_reason, axis = 1)
						failed_HGVS_df.to_csv(output_directory+"/"+database_name+"/"+disease_name+"/Invalid_Annotations/"+gene+"_"+transcript_accession+"_"+database_name+"_InvalidResults.csv")
					print("Finished "+transcript_accession+" for gene "+gene+" in database "+database_name+".", sep = "")

# Now run the function
web_scrape(disease_names, input_gene_lists, databases_list)
