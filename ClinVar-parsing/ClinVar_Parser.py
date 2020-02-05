#!/usr/bin/env python
# coding: utf-8

import os
from argparse import ArgumentParser, FileType
from os.path import splitext
import pandas as pd
import re
from xml.etree import ElementTree as ET
import hgvs.normalizer
import hgvs.variantmapper
import hgvs.dataproviders.uta
import hgvs.assemblymapper
import hgvs.parser
import hgvs.validator
import hgvs.exceptions
import pandas as pd
import numpy as np
import glob
from itertools import chain


args = ArgumentParser('./ClinVar_Parser.py', description="""This program has been designed to parse the XML file
containing the entire ClinVar database to obtain variant information for specific genes relating to a disease.
For details on obtaining the ClinVar XML file, see README.md.
Example usage: ./ClinVar_Parser.py --xml_file ClinVarFullRelease_00-latest.xml
--disease_gene_lists SCID_ny_panel.txt Metabolic_diseases_genes.txt --disease_names SCID Metabolic_Diseases
--output_directory Total_Outputs""")

args.add_argument(
	'--xml_file',
	help="This is the ClinVar xml file that can be obtained form the NCBI database (ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/). Default file is ClinVarFullRelease_00-latest.xml",
	default="ClinVarFullRelease_00-latest.xml",
)

args.add_argument(
	'--disease_gene_lists',
	nargs='+', # This tells the program that if they specify this flag, they have to give it at least one input. If they don't specify it, then the default will go in.
	help="""\
This is a list of text files containing genes associated with the disease of interest. The text files
should each contain a list of gene symbols for one given disease, one gene per line. The name of the disease will be specified
by the arguement --disease_names. """,
	default = None
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


args = args.parse_args()
ClinVar_File = args.xml_file

if args.disease_gene_lists == None:
	# If no disease gene files were given, exit the program and list the txt files in the directory that are likely to be the disease gene files
	txt_files_list = glob.glob("*.txt")
	print("No gene list files specified. Please create a text file with a list of gene symbols, one on each line. \nExample usage: ./ClinVar_Parser.py --xml_file ClinVarFullRelease_00-latest.xml --disease_gene_lists SCID_ny_panel.txt Metabolic_diseases_genes.txt --disease_names SCID Metabolic_Diseases --output_directory Total_Outputs \nPossible gene list files in your directory include: "+', '.join(txt_files_list))
	exit(-1)
input_gene_lists = []
for text_file in args.disease_gene_lists:
	input_gene_lists.append([line.rstrip('\n') for line in open(text_file)])

if args.disease_names == None:
	args.disease_names = [splitext(gene_list)[0] for gene_list in args.disease_gene_lists]
disease_names = args.disease_names

output_directory = args.output_directory
for disease in disease_names:
	os.makedirs(output_directory+"/ClinVar/"+disease+"/Invalid_Annotations/", exist_ok=True)
	# This will make the output_directory, ClinVar, disease name, and Invalid_Annotations directories if they are not present

# Moved these below all argparse statements because they take a long time to execute
hdp = hgvs.dataproviders.uta.connect()
hp = hgvs.parser.Parser()
hn = hgvs.normalizer.Normalizer(hdp)
am = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name='GRCh37', alt_aln_method='splign', replace_reference=True)
vr = hgvs.validator.Validator(hdp=hdp)


def parse_clinvar_xml(disease_names, input_gene_lists, ClinVar_File, output_directory):
	# Create an object that has all of the genes in all of the lists. If the variant is in any of the lists, we want to parse out the information
	all_genes_list = list(chain(*input_gene_lists)) # The star function unpacks the list within input_gene_lists and chain combines them

	# Create an empty dictionary that can be filled with the variant information
	combined_diseases_dictionary = {}
	for disease, gene_list in zip(disease_names, input_gene_lists):
		combined_diseases_dictionary[disease] = {}
		for gene in gene_list:
			combined_diseases_dictionary[disease][gene] = {'Stars': [], 'No_Stars': [], 'Invalid': []}
		# Add another element to the dictionary to account for any variants that affect multiple genes within the list
		combined_diseases_dictionary[disease]['Multiple_Genes'] = {'Stars': [], 'No_Stars': [], 'Invalid': []}

	# It is unlikely that anything will say it has a gene for the RCV but then have variants with no gene listed in the annotation,
	# But if they do arise (like a compound het where one variant is a loss of an entire chromosome arm), the variant with no gene
	# needs a place to be stored. It will be saved out without any gene name associated with it.
	not_annotated_variants = []
	no_gene_symbol = []

	# Parse the variants
	for event, elem in ET.iterparse(ClinVar_File, events=['end']):
		if elem.tag != 'ClinVarSet':
			continue
		# No else statement needed here because it will go back to the start (continue) if it isn't a ClinVarSet
		genes_in_rcv = []
		for child in (elem.findall('ReferenceClinVarAssertion//MeasureSet/Measure/MeasureRelationship/Symbol/ElementValue')):
			genes_in_rcv.append(child.text)
			# It is possible for this to end up being a blank list. That is just fine.

		if not any(x in genes_in_rcv for x in all_genes_list):
			elem.clear() # Clear the memory once we realize that the RCV does not contain variants in genes within the input list
			continue

		# Set the variables to a dash. If we can find information for the variable, the dashes will end up being replaced
		Pathogenicity = "-"
		Disease = "-"
		Genetic_Origin = "-"
		Compound_Het = "-"
		Inheritance_Pattern = "-"
		RCV_num = "-"
		Edited_Date = "-"
		Review_Status = "-"
		Star_Level = "-"
		Submitter = "-"

		genotype_set_elem = elem.find('ReferenceClinVarAssertion/GenotypeSet')
		if genotype_set_elem != None:
			Compound_Het = genotype_set_elem.get("Type", "-")

		submitter_elem = elem.findall('ClinVarAssertion/ClinVarSubmissionID')
		submitter_len = len(submitter_elem)
		if submitter_len == 1:
			Submitter = submitter_elem[0].get("submitter", "-")
		elif submitter_len > 1:
			submitter_list = []
			for child in submitter_elem:
				submitter_list.append(child.get("submitter", "-"))
			Submitter = ("Multiple Submitters: "+', '.join(submitter_list))

		review_status_elem = elem.find('ReferenceClinVarAssertion/ClinicalSignificance/ReviewStatus')
		if review_status_elem != None:
			Review_Status = review_status_elem.text

		if Review_Status in ['criteria provided, single submitter', 'criteria provided, conflicting interpretations']:
			Star_Level = 1
		elif Review_Status == 'criteria provided, multiple submitters, no conflicts':
			Star_Level = 2
		elif Review_Status == 'reviewed by expert panel':
			Star_Level = 3
		elif Review_Status == 'practice guideline':
			Star_Level = 4
		else:
			Star_Level = 0 # Ones I see here are 'no assertion criteria provided' and 'no assertion provided'

		pathogenicity_elem = elem.find('ReferenceClinVarAssertion/ClinicalSignificance/Description')
		if pathogenicity_elem != None:
			Pathogenicity = pathogenicity_elem.text

		condition_elem = elem.findall('ReferenceClinVarAssertion/TraitSet/Trait/Name/ElementValue[@Type = "Preferred"]')
		condition_length = len (condition_elem)
		if condition_length == 1:
			Disease = condition_elem[0].text
		elif condition_length > 1:
			condition_list = []
			for child in condition_elem:
				condition_list.append(child.text)
			Disease = ("Multiple conditions listed as primary: "+', '.join(condition_list))

		genetic_origin_elem = elem.find('ReferenceClinVarAssertion/ObservedIn/Sample/Origin')
		if genetic_origin_elem != None:
			Genetic_Origin = genetic_origin_elem.text

		attribute_elem = elem.find('ReferenceClinVarAssertion/AttributeSet/Attribute')
		if attribute_elem != None:
			Inheritance_Pattern = attribute_elem.text

		RCV_elem = elem.find('ReferenceClinVarAssertion/ClinVarAccession')
		if RCV_elem != None:
			RCV_num = RCV_elem.get('Acc', "-")

		Date_elem = elem.find('ReferenceClinVarAssertion/ClinicalSignificance')
		if Date_elem != None:  # ReferenceClinVarAssertion
			Edited_Date = Date_elem.get('DateLastEvaluated', "-")

		measure_set_elem = elem.findall('ReferenceClinVarAssertion//MeasureSet/Measure')
		if measure_set_elem == None:
			# If nothing is here, that means there is no information about actual variant.
			# I hope this never happens, but if it does then we can't assign it to a gene, so it must go to its own output

			# Most of these things can't be obtained without variant information, so they will be marked with dashes, These include:
			# Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession, Protein_Notation, Protein_HGVS, Assembly,
			# Chromosome, Position_g_start, Position_g_stop, Ref_allele, Alt_allele, Genomic_annotation, Genomic_Normalized,
			# Var_Type, Var_Length, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt, transcript_error, genomic_error
			Individual_Variant_List = []
			Gene_Symbol = "-"
			Affected_Genes = ', '.join(genes_in_rcv)
			found_disease = False
			# If there are no genes in RCV, then just skip this variant
			if len(genes_in_rcv) == 0:
				elem.clear()
				continue
			for disease, gene_list in zip(disease_names, input_gene_lists):
				# If the RCV has something, find out which gene(s) is(are) within this disease
				genes_from_diseases_list = []
				for gene in genes_in_rcv:
					if gene in gene_list:
						genes_from_diseases_list.append(gene)

				num_disease_genes = len(genes_from_diseases_list)
				if num_disease_genes >= 1:
					# If at least one gene from the variant overlaps with the disease genes, then change variable to say it has found a disease associated with a gene from this variant
					found_disease = True
				if num_disease_genes == 1:
					Gene_Symbol = genes_from_diseases_list[0]
				elif num_disease_genes > 1:
					Gene_Symbol = ', '.join(genes_from_diseases_list)

				Individual_Variant_List = ["-", "-", "-", "-", "-", "-", "-", "-", "-", "-", Pathogenicity, Disease, Genetic_Origin,
											Inheritance_Pattern, Affected_Genes, Gene_Symbol, Compound_Het,"-", "-", "-", "-", "-",
											"-", "-", "-", "-", "-", "ClinVar", RCV_num, Review_Status, Star_Level, Submitter,
											Edited_Date, "-", "-", "-", "No Variant Information Present"]
				if num_disease_genes == 1:
					# Append to dictionary for individual gene
					combined_diseases_dictionary[disease][Gene_Symbol]['Invalid'].append(Individual_Variant_List)
				elif num_disease_genes > 1:
					# Append to dictionary for multiple genes
					combined_diseases_dictionary[disease]['Multiple_Genes']['Invalid'].append(Individual_Variant_List)

			if not found_disease:
				Individual_Variant_List = ["-", "-", "-", "-", "-", "-", "-", "-", "-", "-", Pathogenicity, Disease, Genetic_Origin,
											Inheritance_Pattern, Affected_Genes, Gene_Symbol, "-", Compound_Het,"-", "-", "-", "-", "-",
											"-", "-", "-", "-", "-", "ClinVar", RCV_num, Review_Status, Star_Level, Submitter,
											Edited_Date, "-", "No Variant Information Present"]
				not_annotated_variants.append(Individual_Variant_List)
				print(RCV_num+" contains at least one variant that did not have any genes listed that were associated with any of the diseases in the input. The information associated with this variant will be stored to "+output_directory+"/Not_Annotated_Variants.csv")
				print("Please check this variant on the ClinVar website: https://www.ncbi.nlm.nih.gov/clinvar/ ")

			elem.clear()
			continue

		for child in measure_set_elem:
			# Reset the variables so that they can't carry over from a previous child term
			Individual_Variant_List = []
			Assembly = "-"
			Transcript = "-"
			Gene_Symbol = "-"
			Transcript_notation = "-"
			Transcript_HGVS = "-"
			Protein_Notation = "-"
			Protein_HGVS = "-"
			Protein_Accession = "-"
			Chromosome = "-"
			Position_g_start = "-"
			Position_g_stop = "-"
			Ref_allele = "-"
			Alt_allele = "-"
			Chr_Accession = "-"
			Pos_VCF = "-"
			VCF_Ref = "-"
			VCF_Alt = "-"
			Var_Type = ""
			Var_Length = "-"
			Genomic_annotation = "-"
			Genomic_Normalized = "-"
			Affected_genes_list = []
			Affected_Genes = "-"
			transcript_error = "-"
			genomic_error = "-"
			failure_reason = "-"
			dbSNP = "-"

			Var_Type = child.get('Type', '-')

			for grandchild in child.findall('MeasureRelationship/Symbol/ElementValue'):
				Affected_genes_list.append(grandchild.text)
			Affected_Genes = ', '.join(Affected_genes_list) # This says to separate them by a comma and a space

			for grandchild in child.findall('XRef'):
				if grandchild.get('DB', '-') == "dbSNP":
					dbSNP_id = grandchild.get("ID", "-")
					if dbSNP_id != "-":
						dbSNP = "rs"+dbSNP_id

			for grandchild in child.findall('Name/ElementValue[@Type="Preferred"]'): # Again, always comes back with something, just sometimes a blank list
				# There should only be one "Preferred" term per variant listed, so I shouldn't have to reset this.
				if "NM_" in grandchild.text:
					match = re.search(r"NM_\d+\.\d", grandchild.text)
					if match != None:
						Transcript = match.group(0)

					match = re.search(r"\(([A-Z]+.*)\):", grandchild.text)
					if match != None:
						Gene_Symbol = match.group(1)
					match = re.search(r"c\.[^\s]+", grandchild.text)
					if match != None:
						Transcript_notation = match.group(0).replace("[", "").replace("]", "").replace("=", "").replace("/", "")
					match = re.search(r"p\.[^\s\)]+", grandchild.text)
					if match != None:
						Protein_Notation = match.group(0).replace("[", "").replace("]", "").replace("=", "").replace("/", "")

			## Get the official HGVS genomic annotations if there are any present
			for grandchild in child.findall('AttributeSet/Attribute'):
				if grandchild.attrib['Type'] == 'HGVS, genomic, top level, previous':
					Genomic_annotation = grandchild.text
				elif grandchild.attrib['Type'] == 'HGVS, coding, RefSeq':
					Transcript_HGVS = grandchild.text
				elif grandchild.attrib['Type'] == 'HGVS, protein, RefSeq':
					Protein_HGVS = grandchild.text
			### Now get the protein accession if I can from the protein HGVS
			if Protein_HGVS != "-":
				potential_protein = re.search(r"NP_\d+\.\d+", Protein_HGVS)
				if potential_protein != None:
					Protein_Accession = potential_protein.group(0)

			for grandchild in (child.findall('SequenceLocation[@Assembly="GRCh37"]')):
				# Again, there should be only one term for each of these, I simply need this for loop to be able to get the attributes
				Assembly = grandchild.get('Assembly', '-') # These ones are instead of grandchild.attrib['Assembly'] with a try except
				Chromosome = grandchild.get('Chr', '-')
				Position_g_start = grandchild.get('start', '-')
				Position_g_stop = grandchild.get('stop', '-')
				Ref_allele = grandchild.get('referenceAllele', '-')
				Alt_allele = grandchild.get('alternateAllele', '-')
				Var_Length = grandchild.get('variantLength', '-')
				Chr_Accession = grandchild.get('Accession', '-')
				Pos_VCF = grandchild.get('positionVCF', '-')
				VCF_Ref = grandchild.get('referenceAlleleVCF', '-')
				VCF_Alt = grandchild.get('alternateAlleleVCF', '-')

			# Some of these don't have a genomic annotation yet, but all of the information needed to build the genomic annotation is present
			# This next little bit is just to build the genomic annotation for any of those ones, because it could be used for HGVS normalization
			if Genomic_annotation == "-":
				if Var_Type == "single nucleotide variant":
					if Chr_Accession != "-" and Position_g_start != "-" and Ref_allele != "-" and Alt_allele != "-":
						Genomic_annotation = Chr_Accession+":g."+str(Position_g_start)+Ref_allele+">"+Alt_allele

				elif Var_Type == "Deletion":
					if Var_Length == "-":
						pass
					elif Var_Length == "1": # These are strings, and some are dashes so I can't convert them to ints
						if Chr_Accession != "-" and Position_g_start != "-":
							Genomic_annotation = Chr_Accession+":g."+str(Position_g_start)+"del"
					else:
						if Chr_Accession != "-" and Position_g_start != "-" and Position_g_stop != "-":
							Genomic_annotation = Chr_Accession+":g."+str(Position_g_start)+"_"+str(Position_g_stop)+"del"

				elif Var_Type == "Duplication":
					if Var_Length == "-":
						pass
					elif Var_Length == "1":
						if Chr_Accession != "-" and Position_g_start != "-":
							Genomic_annotation = Chr_Accession+":g."+str(Position_g_start)+"dup"
					else:
						if Chr_Accession != "-" and Position_g_start != "-" and Position_g_stop != "-":
							Genomic_annotation = Chr_Accession+":g."+str(Position_g_start)+"_"+str(Position_g_stop)+"dup"

				elif Var_Type == "Insertion":
					if Chr_Accession != "-" and Position_g_start != "-" and Alt_allele != "-" and Position_g_stop != "-":
						Genomic_annotation = Chr_Accession+":g."+str(Position_g_start)+"_"+str(Position_g_stop)+"ins"+Alt_allele

				elif Var_Type == "Indel":
					if Ref_allele != "-":
						ref_len = len(Ref_allele)
						if ref_len == 1:
							if Chr_Accession != "-" and Position_g_start != "-" and Alt_allele != "-":
								Genomic_annotation = Chr_Accession+":g."+str(Position_g_start)+"delins"+Alt_allele
						if ref_len > 1:
							if Chr_Accession != "-" and Position_g_start != "-" and Alt_allele != "-" and Position_g_stop != "-":
								Genomic_annotation = Chr_Accession+":g."+str(Position_g_start)+"_"+str(Position_g_stop)+"delins"+Alt_allele

			## The next step is to normalize the genomic annotation
			## This runs faster from the transcript than it does from the genomic annotation, so I am going to start with that
			## But if the transcript fails, then I am going to use the genomic annotation
			biocommons_format = Transcript_HGVS
			if "[=/" in Transcript_HGVS:
				biocommons_format = Transcript_HGVS.replace("[=/", "").replace("]", "") # Multiple RCVs have these put it in square brackets and have a =/ before the c.
			try:
				converted_object = am.c_to_g(hp.parse_hgvs_variant(biocommons_format))
			except hgvs.exceptions.HGVSError as e:
				transcript_error = str(e)
				# If the transcript notation didn't work, now we have to use the genomic annotation
				if ">" in Genomic_annotation:
					converted_object = None # Needed later in loop
					# the SNVs don't need to be normalized, instead they need to be validated
					try:
						validated = vr.validate(hp.parse_hgvs_variant(Genomic_annotation))
					except hgvs.exceptions.HGVSError as e:
						genomic_error = str(e)
						validated = False
					if validated == True:
						Genomic_Normalized = Genomic_annotation
					# Now fill in anything that has a dash with the information in the genomic annotation
					# If something was present, I want to leave it, in case something doesn't match
					if Var_Type == "-":
						Var_Type = "single nucleotide variant"
					if Var_Length == "-":
						Var_Length = 1
					if Position_g_start == "-":
						pos_search = re.search(r"g.(\d+)", Genomic_Normalized)
						if pos_search != None:
							Position_g_start = pos_search.group(1)
					if Position_g_stop == "-":
						Position_g_stop = Position_g_start # May be a dash
					if Ref_allele == "-" or Alt_allele == "-":
						ref_and_alt_allele_search = re.search(r"g.\d+([A,C,T,G])>([A,C,T,G])", Genomic_Normalized)
						if ref_and_alt_allele_search != None:
							Ref_allele = ref_and_alt_allele_search.group(1)
							Alt_allele = ref_and_alt_allele_search.group(2)
				else:
					# This is all non-SNV variants. They need to be normalized rather than validated
					# There are a number of variants that contain the information, but it is in a format that needs to be updated before HGVS normalization
					use_for_biocommons = Genomic_annotation
					match_ins = re.search(r'(.+ins)\d+', Genomic_annotation)
					if match_ins != None and Alt_allele != "-":
						if Ref_allele == Alt_allele[0]:
							# Many of the insertions keep one base for the Ref_allele and then include it in the alternate
							use_for_biocommons = match_ins.group(1)+Alt_allele[1:]
						else:
							# This will catch ones with Ref_allele == "-" and the indels that have an actual sequence for the ref allele
							use_for_biocommons = match_ins.group(1)+Alt_allele
					match_dup = re.search(r'(NC_\d+\.\d+:g\.\d+_?\d*dup)\d+', Genomic_annotation)
					if match_dup != None:
						use_for_biocommons = match_dup.group(1)
					try:
						converted_object = hn.normalize(hp.parse_hgvs_variant(use_for_biocommons))
					except hgvs.exceptions.HGVSError as e:
						genomic_error = str(e)
						converted_object = None

			if converted_object != None:
				Genomic_Normalized = str(converted_object)
				# If a converted object is saved, then it can be used to obtain lots of information about the variant.
				# The variant type sometimes changes after normalization, like an insertion may be changed to a duplicate
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
				if Position_g_start == "-":
					try:
						Position_g_start = converted_object.posedit.pos.start.base
					except AttributeError:
						Position_g_start = "-"
				if Position_g_stop == "-":
					try:
						Position_g_stop = converted_object.posedit.pos.end.base
					except AttributeError:
						Position_g_stop = "-"
				if Var_Type != "Insertion":
					if Ref_allele == "-":
						try:
							Ref_allele = converted_object.posedit.edit.ref
						except AttributeError:
							Ref_allele = "-"
				if Alt_allele == "-":
					if Var_Type == "Duplication" and Ref_allele != "-":
						Alt_allele = Ref_allele+Ref_allele
					elif Var_Type != "Inversion" and Var_Type != "Deletion":
						try:
							Alt_allele = converted_object.posedit.edit.alt
						except AttributeError:
							Alt_allele = "-"
				if Var_Length == "-":
					if Var_Type == "Insertion" or Var_Type == "Inversion" or Var_Type == "Duplication":
						if Position_g_stop != "-" and Position_g_start != "-":
							Var_Length = str(int(Position_g_stop)-int(Position_g_start)+1)
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

			## If there was no converted object and some information was missing, there is still a possibility to get the information from other parts of the variant information
			using_transcript = False
			input_for_filling_blanks = "-"
			if Genomic_Normalized != "-":
				input_for_filling_blanks = Genomic_Normalized
			elif Genomic_annotation != "-":
				input_for_filling_blanks = Genomic_annotation
			elif Transcript_notation != "-":
				input_for_filling_blanks = Transcript_HGVS
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

				if Position_g_start == "-" and not using_transcript: # Only do the regex loop if start position is missing, only valid result when not using transcript annotation
					pos_start_search = re.search(r"g.(\d+)", input_for_filling_blanks)
					if pos_start_search != None:
						Position_g_start = pos_start_search.group(1)
				## Filling in the rest of the information about the variant depends on what the variant type is
				if Var_Type == "single nucleotide variant" and not using_transcript:
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
					# Can't obtain reference allele without looking it up, and it already failed so ref will stay as a dash
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
					if Position_g_stop == "-" or Var_Length == "-": # Only run the regex loop if something is missing
						multiple_location_search = re.search(r"\d+_\d+", input_for_filling_blanks)
						if multiple_location_search == None:
							if Position_g_stop == "-" and not using_transcript:
								Position_g_stop = Position_g_start # If only one of these was a dash, the other one will be replaced too but should be replaced with the same thing
							Var_Length = 1
						else:
							stop_search = re.search(r"_(\d+)[a-z]+", input_for_filling_blanks)
							if stop_search != None:
								if Position_g_stop == "-" and not using_transcript:
									Position_g_stop = stop_search.group(1)
							if Var_Length == "-" and Position_g_stop != "-" and Position_g_start != "-":
								Var_Length = int(Position_g_stop) - int(Position_g_start) + 1

			# There are still some things that can be filled in
			if Assembly == "-":
				if Genomic_annotation != "-":
					Assembly = "GRCh37" # If the variant is still being processed by these methods, then the information obtained is using GRCh37
			if Chr_Accession == "-":
				potential_genomic = re.search(r"NC_\d+\.\d+", Genomic_annotation)
				if potential_genomic != None:
					Chr_Accession = potential_genomic.group(0)
			if Chromosome == "-":
				if "NC_" in Chr_Accession:
					potential_genomic = re.search('NC_0+(\d+)\.\d+', Chr_Accession)
					if potential_genomic != None:
						Chromosome = potential_genomic.group(1)
			# Some of them from this method will give us numbers when they should be X, Y, or MT
			if Chromosome == "23": # These are from if I got the chromosome from the chromosome accession
				Chromosome = "X"
			elif Chromosome == "24":
				Chromosome = "Y"
			elif Chromosome == "12920":
				Chromosome = "MT"

			# If the genomic normalization did not work, find out the reason from the failure from the error messages
			if Genomic_Normalized == "-":
				# Now find the failure reason if it didn't normalize properly, use the transcript_error and genomic_error messages
				# Use the genomic error message first, because many have transcript messages that will say "out of region"
				indicate_no_variant_info = [":g.?:", "-:"]
				if not any(l in genomic_error for l in indicate_no_variant_info):
					microsatellite_search = re.search(r"\[\d+\]", genomic_error)
					inserted_unknown_search = re.search(r"ins\(?.?(\d+|\?)", genomic_error)
					if "does not agree with reference sequence" in genomic_error:
						failure_reason = "Incorrect reference base"
					elif microsatellite_search != None:
						failure_reason = "Microsatellite"
					elif inserted_unknown_search != None:
						# May also contain an unknown breakpoint
						if "(?_" in genomic_error or "_?)" in genomic_error:
							failure_reason = "Inserted unknown sequence, unknown breakpoint"
						else:
							failure_reason = "Inserted unknown sequence"
					elif "?" in genomic_error:
						failure_reason = "Unknown breakpoint"
					elif ")_(" in genomic_error:
						failure_reason = "Compound variant"
					else:
						failure_reason = "No definitive failure reason detected, likely compound variant with nontraditional formatting"
				else:
					# These ones have no information for genomic, but may have something for the transcript annotation
					microsatellite_search = re.search(r"\[\d+\]", transcript_error)
					inserted_unknown_search = re.search(r"ins\(?.?(\d+|\?)", transcript_error)
					if any(l in transcript_error for l in indicate_no_variant_info):
						failure_reason = "No variant information provided"
					elif microsatellite_search != None:
						failure_reason = "Microsatellite"
					elif inserted_unknown_search != None:
						# May also contain an unknown breakpoint
						if "(?_" in transcript_err or "_?)" in transcript_error:
							failure_reason = "Inserted unknown sequence, unknown breakpoint"
						else:
							failure_reason = "Inserted unknown sequence"
					elif "?" in transcript_error:
						failure_reason = "Unknown breakpoint"
					elif ")_(" in transcript_error:
						failure_reason = "Compound variant"
					else:
						failure_reason = "No definitive failure reason detected, likely compound variant with nontraditional formatting"

			## Each variant may need to be appended to the dictionary in multiple places,
			## because a gene might be in multiple disease panels, or multiple genes might be affected
			saved_to_dictionary = False
			for disease, gene_list in zip(disease_names, input_gene_lists):
				multi_gene = False
				# If no gene symbol has been obtained yet, try to see if anything can be taken from the list of genes in the RCV
				if Gene_Symbol == "-":
					if len(genes_in_rcv) == 1:
						Gene_Symbol= genes_in_rcv[0]
					else:
						genes_from_diseases_list = []
						for gene in genes_in_rcv:
							if gene in gene_list and gene not in genes_from_diseases_list:
								genes_from_diseases_list.append(gene)
						num_disease_genes = len(genes_from_diseases_list)
						if num_disease_genes == 1:
							Gene_Symbol = genes_from_diseases_list[0]
						elif num_disease_genes > 1:
							Gene_Symbol = ', '.join(genes_from_diseases_list)
							multi_gene = True

				if Genomic_Normalized == "-":
					variant_category = "Invalid"
				else:# Determine if it has stars or not for the ones that passed validation
					if Star_Level > 0:
						variant_category = "Stars"
					else:
						variant_category = "No_Stars"
				Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
											Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
											Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol, dbSNP,
											Compound_Het, Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
											Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
											"ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date,
											transcript_error, genomic_error]
				if Genomic_Normalized == "-":
					Individual_Variant_List.append(failure_reason)

				if multi_gene == True:
					saved_to_dictionary = True
					combined_diseases_dictionary[disease]['Multiple_Genes'][variant_category].append(Individual_Variant_List)
				elif Gene_Symbol in gene_list:
					saved_to_dictionary = True
					combined_diseases_dictionary[disease][Gene_Symbol][variant_category].append(Individual_Variant_List)

			if saved_to_dictionary == False:
				# This means the variant wasn't saved anywhere
				# These are probably large copy number variants where the "preferred" variant gene was not one in the input disease lists
				# Need to check the gene list and pick one that is within a disease
				for disease, gene_list in zip(disease_names, input_gene_lists):
					multi_gene = False
					genes_from_diseases_list = []
					for gene in genes_in_rcv:
						if gene in gene_list and gene not in genes_from_diseases_list:
							genes_from_diseases_list.append(gene)
					num_disease_genes = len(genes_from_diseases_list)
					if num_disease_genes == 1:
						Gene_Symbol = genes_from_diseases_list[0]
					elif num_disease_genes > 1:
						Gene_Symbol = ', '.join(genes_from_diseases_list)
						multi_gene = True
					# If num_disease_genes == 0, then just leave it alone

					Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
												Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
												Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol, dbSNP,
												Compound_Het, Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
												Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
												"ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date,
												transcript_error, genomic_error]

					if multi_gene == True:
						saved_to_dictionary = True
						combined_diseases_dictionary[disease]['Multiple_Genes'][variant_category].append(Individual_Variant_List)
					elif Gene_Symbol in gene_list:
						saved_to_dictionary = True
						combined_diseases_dictionary[disease][Gene_Symbol][variant_category].append(Individual_Variant_List)

			if saved_to_dictionary == False:
				# If after that previous loop, there is still no gene associated with a disease, then put these into an output with no associated disease
				if variant_category != "Invalid":
					Individual_Variant_List.append("-")
				no_gene_symbol.append(Individual_Variant_List)
				print(RCV_num +" contains at least one variant that did not have a gene symbol associated with any of the diseases in the input. The information associated with this variant will be stored to "+output_directory+"/ClinVar/Not_Associated_With_Specified_Disease.csv")
				print("Please check this variant on the ClinVar website: https://www.ncbi.nlm.nih.gov/clinvar/ ")

	# Now save out the collected lists
	column_labels = ["Genome Assembly", "Chr", "Position Start", "Position Stop", "Ref", "Alt",
						"Genomic Annotation", "HGVS Normalized Genomic Annotation", "Variant Type",
						"Variant Length", "Pathogenicity", "Disease", "Genetic Origin", "Inheritance Pattern",
						"Affected Genes" , "Gene Symbol", "dbSNP ID", "Compound Het Status", "Transcript",
						"Transcript Notation", "HGVS Transcript Notation","Protein Accession",
						"Protein Notation", "HGVS Protein Annotation", "Chr Accession", "VCF Pos",
						"VCF Ref", "VCF Alt", "Database", "ClinVar Accession", "Review Status",
						"Star Level", "Submitter", "Edited Date", "Transcript Normalization Failure Message",
						"Genomic Normalization Failure Message"]
	column_labels_invalid = ["Genome Assembly", "Chr", "Position Start", "Position Stop", "Ref", "Alt",
						"Genomic Annotation", "HGVS Normalized Genomic Annotation", "Variant Type",
						"Variant Length", "Pathogenicity", "Disease", "Genetic Origin", "Inheritance Pattern",
						"Affected Genes" , "Gene Symbol", "dbSNP ID", "Compound Het Status", "Transcript",
						"Transcript Notation", "HGVS Transcript Notation","Protein Accession",
						"Protein Notation", "HGVS Protein Annotation", "Chr Accession", "VCF Pos",
						"VCF Ref", "VCF Alt", "Database", "ClinVar Accession", "Review Status",
						"Star Level", "Submitter", "Edited Date", "Transcript Normalization Failure Message",
						"Genomic Normalization Failure Message", "HGVS Normalization Failure Reason"]
	for disease_key in combined_diseases_dictionary.keys():
		for gene_key in combined_diseases_dictionary[disease_key].keys():
			for output_file_key in combined_diseases_dictionary[disease_key][gene_key].keys():
				individual_results_list = combined_diseases_dictionary[disease_key][gene_key][output_file_key]
				if len(individual_results_list) > 0:
					individual_df = pd.DataFrame(individual_results_list)
					if output_file_key == "Invalid":
						individual_df.columns = column_labels_invalid
						individual_df.to_csv(output_directory+"/ClinVar/"+disease_key+"/Invalid_Annotations/"+gene_key+"_ClinVar_InvalidResults.csv")
					elif output_file_key == "No_Stars":
						individual_df.columns = column_labels
						individual_df.to_csv(output_directory+"/ClinVar/"+disease_key+"/"+gene_key+"_ClinVar_No_Star_Results.csv")
					elif output_file_key == "Stars":
						individual_df.columns = column_labels
						individual_df.to_csv(output_directory+"/ClinVar/"+disease_key+"/"+gene_key+"_ClinVar_Results.csv")
	if len(not_annotated_variants) > 0:
		not_annotated_df = pd.DataFrame(not_annotated_variants)
		not_annotated_df.columns = column_labels_invalid
		not_annotated_df.to_csv(output_directory+"/ClinVar/Not_Annotated_Variants.csv")
	if len(no_gene_symbol) > 0:
		no_gene_df = pd.DataFrame(no_gene_symbol)
		no_gene_df.columns = column_labels_invalid
		no_gene_df.to_csv(output_directory+"/ClinVar/Not_Associated_With_Specified_Disease.csv")


parse_clinvar_xml(disease_names, input_gene_lists, ClinVar_File, output_directory)
