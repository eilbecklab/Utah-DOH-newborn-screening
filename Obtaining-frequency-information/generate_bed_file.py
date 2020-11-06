from argparse import ArgumentParser, FileType

args = ArgumentParser('./generate_bed_file.py', description='''This program is designed to generate a bed file from a GFF file.
Example usage: ./generate_bed_file.py -g GRCh37_latest_genomic.gff -o filtered_genes.bed -l SCID_ny_panel.txt
Metabolic_diseases_genes.txt MLD_genes.txt -a 1000''')

args.add_argument(
	'-g',
	'--gff_file',
	type=FileType('r'),
	help="This is a GFF file containing the information about the transcripts of the genome. Default is GRCh37_latest_genomic.gff.",
	default="GRCh37_latest_genomic.gff",
)

args.add_argument(
	'-o',
	'--output_bed_file',
	help="This is the name of the output bed file that can be produced to be used for filtering the gnomAD vcf file. The default name is 'genes.bed'.",
	default = "genes.bed"
)

args.add_argument(
	'-l',
	'--disease_gene_lists',
	nargs='+', # This tells the program that if they specify this flag, they have to give it at least one input. If they don't specify it, then the default will go in.
	help="""\
This is a list of text files containing genes associated with the disease of interest. The text files
should each contain a list of gene symbols for one given disease, one gene per line. If nothing is
specified, the bed file will contain info from all genes in the gff file.""",
	default = None
)

args.add_argument(
	'-a',
	'--added_distance',
	help="This is the distance which will be subtracted from the start position and added to the stop position for each gene. By default, no distance will be added.",
	default = 0
)

import csv
import pandas as pd
import re
from itertools import chain

args = args.parse_args()

output_bed = args.output_bed_file
if ".bed" not in output_bed:
	output_bed = output_bed+".bed"

add_distance = int(args.added_distance)
gene_lists = args.disease_gene_lists

file = args.gff_file
gff_list = []
gff_file = csv.reader(file, delimiter="\t")
for line in gff_file:
	gff_list.append(line)

gene_lines = []
for line in gff_list:
	if len(line) == 9 and line[2] == 'gene':
		gene_lines.append(line)
gene_df = pd.DataFrame(gene_lines)
gene_df.columns = ["Chromosome_Accession", "Source", "Type", "Start", "Stop", "Score", "Sense", "-", "Info"]

def get_gene_name(row):
	info = row['Info']

	gene_result = re.search(r"gene=([A-Za-z0-9\-\._]+);", info)
	gene_name = "-"
	if gene_result != None:
		gene_name = gene_result.group(1)

	row["Gene"] = gene_name
	return row


gene_df = gene_df.apply(get_gene_name, axis = 1)

#
chromosome_accessions_list = []
for term in gene_df['Chromosome_Accession'].values:
	if "NC" in term and term not in chromosome_accessions_list:
		chromosome_accessions_list.append(term)

accession_to_chr_dict = {}
for term in chromosome_accessions_list:
	nuclear_chromosome_search = re.search(r"NC_0000(\d+)\.\d+", term)
	mt_search = re.search(r"NC_012920\.\d+", term)
	if nuclear_chromosome_search:
		if int(nuclear_chromosome_search.group(1)) <= 22:
			accession_to_chr_dict[term] = nuclear_chromosome_search.group(1).lstrip("0")
		elif nuclear_chromosome_search.group(1) == "23":
			accession_to_chr_dict[term] = 'X'
		elif nuclear_chromosome_search.group(1) == "24":
			accession_to_chr_dict[term] = 'Y'
	elif mt_search:
		# MT is not in gnomad, but I am going to use MT for in case they update it later
		accession_to_chr_dict[term] = 'MT'

# If they have provided input gene lists, filter the dataset to only hvae the genes from these
if gene_lists:
	input_gene_lists = []
	for text_file in gene_lists:
		input_gene_lists.append([line.rstrip('\n') for line in open(text_file)])
	all_genes_list = list(chain(*input_gene_lists))
	gene_df = gene_df[gene_df['Gene'].isin(all_genes_list)]

# Now get the chromosome name that will be used by gnomAD
def format_for_bed(df):
	chr_accession = df['Chromosome_Accession']
	chr_name = 'contig'
	if chr_accession in accession_to_chr_dict:
		chr_name = accession_to_chr_dict[chr_accession]
	df['Chr'] = chr_name
	pos_start = df['Start']
	new_start = int(pos_start) - add_distance
	if new_start <1:
		new_start = 1
	pos_end = df['Stop']
	new_stop = int(pos_end) + add_distance
	df['New_Start'] = new_start
	df['New_Stop'] = new_stop
	return df


gene_df = gene_df.apply(format_for_bed, axis = 1)
# Now get the columns for the bed file
# Chr, start, stop, name, score (dot), strand (+ or -)
bed_format = gene_df[['Chr', 'New_Start', 'New_Stop', 'Gene','Score', 'Sense']]
# Some of the genes have another gene listed on one of the contigs.
# Remove those ones because gnomAD can't find anything not on the cannnical chromosomes
bed_format = bed_format[bed_format['Chr'] != 'contig']
# Save it to a tab delimited file
bed_format.to_csv(output_bed, sep='\t', index=False, header=False)
