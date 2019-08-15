#!/usr/bin/env python
# coding: utf-8

# The first thing needed here is to obtain the command line arguements

from argparse import ArgumentParser, FileType

args = ArgumentParser('./LOVD3_Variant_Parser.py', description='This program is for obtaining information about transcripts and chromosomes from a GFF file. Example usage: ./LOVD3_Variant_Parser.py --gff_file GRCh37_latest_genomic.gff --output_transcripts_file Transcript_Info_For_Dictionaries.csv --output_chromosome_file Chromosome_name_to_accession.csv')

args.add_argument(
	'--gff_file',
	type=FileType('r'),
	help="This is a GFF file containing the information about the transcripts of the genome. Default is GRCh37_latest_genomic.gff.",
	default="GRCh37_latest_genomic.gff",
)

args.add_argument(
	'--output_transcripts_file',
	help="This is the name of the directory in which you want your output results to be stored. Default is ranscript_Info_For_Dictionaries.csv",
	default = "Transcript_Info_For_Dictionaries.csv"
)

args.add_argument(
	'--output_chromosome_file',
	help="This is the name of the directory in which you want your output results to be stored. Default is Chromosome_name_to_accession.csv",
	default = "Chromosome_name_to_accession.csv"
)

import csv
import pandas as pd
import re

args = args.parse_args()

output_transcripts = args.output_transcripts_file
if ".csv" not in output_transcripts:
	output_transcripts = output_transcripts+".csv"

output_chromosome = args.output_chromosome_file
if ".csv" not in output_chromosome:
	output_chromosome = output_chromosome+".csv"

file = args.gff_file
gff_list = []
gff_file = csv.reader(file, delimiter="\t")
for line in gff_file:
	gff_list.append(line)

# Now grab the rows that contain info about transcripts, don't use the ones from the mitochondria becaue they don't have transcript accessions
transcript_types = ['ncRNA', 'transcript', 'primary_transcript', 'lnc_RNA', 'mRNA', 'snoRNA', 'antisense_RNA', 'snRNA', 'rRNA', 'telomerase_RNA', 'vault_RNA', 'Y_RNA', 'RNase_MRP_RNA', 'RNase_P_RNA', 'SRP_RNA']
transcript_data_list = []
for line in gff_list:
	if len(line) == 9 and line[2] in transcript_types and line[0] != "NC_012920.1":
		transcript_data_list.append(line)
transcripts_df = pd.DataFrame(transcript_data_list)
transcripts_df.columns = ["Chromosome_Accession", "Source", "Type", "Start", "Stop", "-", "Sense", "-", "Info"]

# Now use regex functions to get the data out of the info column
def get_transcript_id(row):
	data_frame_column = row['Info']
	transcript_result = re.search(r"transcript_id=((N(M|R)_\d+)\.\d)", data_frame_column)
	if transcript_result != None:
		transcript_full = transcript_result.group(1)
		transcript_partial = transcript_result.group(2)
	else:
		transcript_full = "-"
		transcript_partial = "-"
	row["Transcript_ID_Full"] = transcript_full
	row["Transcript"] = transcript_partial
	return row

transcripts_df = transcripts_df.apply(get_transcript_id, axis = 1)

def get_gene_name(data_frame_column):
	gene_result = re.search(r"gene=([A-Za-z0-9\-\._]+);", data_frame_column)
	if gene_result == None:
		return "-"
	else:
		return (gene_result.group(1))

transcripts_df["Gene"] = transcripts_df["Info"].apply(get_gene_name)

# Now the info about the transcripts has been obtained, but we are missing anything from the mitochondria
# The mitochondrial genes do not have accessions, so instead this program will use the protein IDs
# They will not be useful for the dictionaries, but they will be helpful in determining which on chromosome the gene is located

mito_data_list = []
for line in gff_list:
	if len(line) == 9 and line[2] == "CDS" and line[0] == "NC_012920.1":
		mito_data_list.append(line)

mito_df = pd.DataFrame(mito_data_list)
mito_df.columns = ["Chromosome_Accession", "Source", "Type", "Start", "Stop", "-", "Sense", "-", "Info"]

# Now use regex functions to get the information
def get_protein_id(row):
	data_frame_column = row["Info"]
	protein_result = re.search(r"protein_id=((YP_\d+)\.\d)", data_frame_column)
	if protein_result != None:
		protein_full = protein_result.group(1)
		protein_partial = protein_result.group(2)
	else:
		protein_full = "-"
		protein_partial = "-"
	row["Transcript_ID_Full"] = protein_full # Must still be transcript ID because that is the column it must merge with
	row["Transcript"] = protein_partial
	return row

mito_df["Gene"] = mito_df["Info"].apply(get_gene_name)
mito_df = mito_df.apply(get_protein_id, axis = 1)

# The LOVD databases use gene names with "MT-" prepending the gene names, but the GFF file does not
# I could just change them, but I want to be able to usae this file for either LOVD or ClinVar, both will be included
mito_df["Gene2"] = mito_df["Gene"].apply(lambda x: "MT-"+x)

df1 = transcripts_df[['Chromosome_Accession', 'Source', 'Type', 'Start', 'Stop', 'Sense', 'Info', 'Transcript', 'Transcript_ID_Full', 'Gene']]
df2 = mito_df[['Chromosome_Accession', 'Source', 'Type', 'Start', 'Stop', 'Sense', 'Info', 'Transcript', 'Transcript_ID_Full', 'Gene']]
df3 = mito_df[['Chromosome_Accession', 'Source', 'Type', 'Start', 'Stop', 'Sense', 'Info', 'Transcript', 'Transcript_ID_Full', 'Gene2']].rename(index = str, columns = {"Gene2": "Gene"})

merged_transcript_info = pd.concat([df1, df2, df3]).reset_index()

# Now remove any of the ones that returned a dash (hopefully only a couple, in example input there are only 2 lncRNAs that had dashes)
merged_transcript_info = merged_transcript_info[(merged_transcript_info["Transcript"] != "-") & (merged_transcript_info["Transcript_ID_Full"] != "-") & (merged_transcript_info["Gene"] != "-")]
# Save out the transcript information
info_for_dictionaries = merged_transcript_info[["Gene", "Chromosome_Accession", "Transcript", "Transcript_ID_Full"]]
info_for_dictionaries.to_csv(output_transcripts)

# Now for the chromosome name and accession
# This is easy to fill in manually, but this will be helpful for obtaining the
# last digit after the decimal point if things have been updated_dated
# Chr 1-22 all are NC_0000 followed by the chromosome number (2 digits, add a 0 in front of 0-9)
# Chr X is NC_000023, Chr Y is NC_000024, Chr M (or MT) is NC_012920

all_chr_accessions = list(info_for_dictionaries[["Chromosome_Accession"]].values.flatten())
chromosome_accessions_list = []
for term in all_chr_accessions:
	if "NC" in term and term not in chromosome_accessions_list:
		chromosome_accessions_list.append(term)

# Read through the list of chromosome accessions and append a csv file with the chromosome name and accession, including both mitochondrial chromosome symbols
with open(output_chromosome, "w") as f:
	writer = csv.writer(f)
	writer.writerow(["","Chromosome","Chromosome_Accession"])
	i = 0
	for term in chromosome_accessions_list:
		nuclear_chromosome_search = re.search(r"NC_0000(\d+)\.\d+", term)
		mt_search = re.search(r"NC_012920\.\d+", term)
		if nuclear_chromosome_search != None:
			if int(nuclear_chromosome_search.group(1)) <= 22:
				writer.writerow([i,nuclear_chromosome_search.group(1).lstrip("0"),term])
			elif nuclear_chromosome_search.group(1) == "23":
				writer.writerow([i,"X",term])
			elif nuclear_chromosome_search.group(1) == "24":
				writer.writerow([i,"Y",term])
			else:
				writer.writerow([i,"0",term]) # This should not happen, but will pick up mistaken terms
			i += 1
		if mt_search != None:
			writer.writerow([i,"M",term])
			i += 1
			writer.writerow([i,"MT",term]) # I have to add this one twice to account for the differences between the American and European annotation of mitochondrial chromosome
			i += 1
