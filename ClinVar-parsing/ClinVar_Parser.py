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


args = ArgumentParser('./ClinVar_Parser.py', description="""This program has been designed to parse the XML file
containing the entire ClinVar database to obtain variant information for specific genes relating to a disase.
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


def parse_clinvar_xml(disease, list_of_genes, ClinVar_File, output_file):

    ClinVar_Variants = []
    not_annotated_variants = []
    compound_variants = []
    copy_number_variants = []
    microsatellites = []
    unknown_breakpoint = []
    no_change_from_reference = []
    no_stars = []
    inserted_unknown_sequence = []
    invalid_hgvs = []
    incorrect_base = []


    for event, elem in ET.iterparse(ClinVar_File, events=['end']):
        if elem.tag != 'ClinVarSet':
            continue
# No else statement needed here because it will go back to the start (continue) if it isn't a ClinVarSet
        genes_in_rcv = []
        for child in (elem.findall('ReferenceClinVarAssertion//MeasureSet/Measure/MeasureRelationship/Symbol/ElementValue')):
            genes_in_rcv.append(child.text)
            # It is possible for this to end up being a blank list. That is just fine.


        if not any(x in genes_in_rcv for x in list_of_genes): # May change this out for something else looking for a disease instead of the genes
            elem.clear() # Clear the memory once we realize that the RCV has genes, just not genes we are interested in
            continue

        # This is where I need to run the loop that I made last time
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
            # If nothing is here, that means there is no information about actual variants.
            # I hope this never happens, but I am worried this might happen for all of the copy number variants
            # Set each thing to a dash so I can save it in the same format as the normal, annotated variants
            Individual_Variant_List = []
            Transcript = "-"
            Gene_Symbol = "-"
            Transcript_notation = "-"
            Transcript_HGVS = "-"
            Protein_Accession = "-"
            Protein_Notation = "-"
            Protein_HGVS = "-"
            Assembly = "-"
            Chromosome = "-"
            Position_g_start = "-"
            Position_g_stop = "-"
            Ref_allele = "-"
            Alt_allele = "-"
            Var_Length = "-"
            Chr_Accession = "-"
            Pos_VCF = "-"
            VCF_Ref = "-"
            VCF_Alt = "-"
            Var_Type = "-"
            Genomic_annotation = "-"
            Affected_Genes = "-"
            Genomic_Normalized = "-"
            ### Double check that everything is here
            Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                       Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                       Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                       Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession, Protein_Notation,
                                       Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt, "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]




            not_annotated_variants.append(Individual_Variant_List)

            elem.clear()
		#del elem # I hope there aren't very many of these
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






            Var_Type = child.get('Type', '-') # This one seems to pick up something, even if I tell it the path includes Monkey

            for grandchild in child.findall('MeasureRelationship/Symbol/ElementValue'):
                Affected_genes_list.append(grandchild.text)
            Affected_Genes = ', '.join(Affected_genes_list) # This says to separate them by a comma and a space

            for grandchild in child.findall('Name/ElementValue[@Type="Preferred"]'): # Again, always comes back with something, just sometimes a blank list
                # There should only be one "Preferred" term per variant listed, so I shouldn't have to reset this.

                if "NM_" in grandchild.text:
                    match = re.search(r"NM_\d+\.\d", grandchild.text)
                    if match != None:
                        Transcript = match.group(0)

                    match = re.search(r"\([A-Z]+.*\):", grandchild.text)
                    if match != None:
                        Gene_Symbol = match.group(0).replace("(", "").replace(")", "").replace(":", "")
                    match = re.search(r"c\.[^\s]+", grandchild.text)
                    if match != None:
                        Transcript_notation = match.group(0).replace("[", "").replace("]", "").replace("=", "").replace("/", "")
                    match = re.search(r"p\.[^\s\)]+", grandchild.text)
                    if match != None:
                        Protein_Notation = match.group(0).replace("[", "").replace("]", "").replace("=", "").replace("/", "")


##################### Dave updated section to not use biocommons so much for speed purposes
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

####################
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

######### this section is for updating things that didn't get the right result from above
            # Need to get genomic annotation from transcript if that is There
            if Genomic_annotation == "-":
                if Transcript_HGVS != "-":
                    if "[=/" in Transcript_HGVS:
                        try:
                            biocommons_format = Transcript_HGVS.replace("[=/", "").replace("]", "") #Multiple RCVs have these put it in square brackets and have a =/ before the c.
                            biocommons_result = am.c_to_g(hp.parse_hgvs_variant(biocommons_format))
                            Genomic_annotation = str(biocommons_result)
                        except:
                            Genomic_annotation = "-"
                    else:
                        try:
                            biocommons_result = am.c_to_g(hp.parse_hgvs_variant(Transcript_HGVS))
                            Genomic_annotation = str(biocommons_result)
                        except:
                            Genomic_annotation = "-"  # In the first 10% of variants this would fix 6 of them

            if Chromosome == "-" or Chromosome == "":
                if Chr_Accession == "-" or Chr_Accession == "":
                    potential_genomic = re.search('NC_0+(\d+)\.\d+', Genomic_annotation)
                    if potential_genomic != None:
                        Chromosome = potential_genomic.group(1)
                elif "NC_" in Chr_Accession:
                    potential_genomic = re.search('NC_0+(\d+)\.\d+', Chr_Accession)
                    if potential_genomic != None:
                         Chromosome = potential_genomic.group(1)


            if Chr_Accession == "-" or Chr_Accession == "":
                potential_genomic = re.search(r"NC_\d+\.\d+", Genomic_annotation)
                if potential_genomic != None:
                    Chr_Accession = potential_genomic.group(0)

            ### Fill in anything blank with information elsewhere if possible
            if Var_Type == "single nucleotide variant":
                if Position_g_start == "-":
                    if Genomic_annotation != "-":
                        search_result = re.search(r"g\.\d+", Genomic_annotation)
                        if search_result != None:
                            Position_g_start = search_result.group(0).replace("g.", "")
                if Position_g_stop == "-":
                    Position_g_stop = Position_g_start ## this might be a blank or a - but that is ok
                if Ref_allele == "-":
                    if Genomic_annotation != "-":
                        search_result = re.search(r"g\.\d+([A,C,G,T])", Genomic_annotation)
                        if search_result != None:
                            Ref_allele = search_result.group(1)
                if Alt_allele == "-":
                    if Genomic_annotation != "" and Genomic_annotation != "-":
                        search_result = re.search(r">([A,C,G,T])", Genomic_annotation)
                        if search_result != None:
                            Alt_allele = search_result.group(1)
                if Var_Length == "-":
                    Var_Length = 1
                if Pos_VCF == "-":
                    Pos_VCF = Position_g_start #May be a dash or blank, but that is fine.
                if VCF_Ref == "-":
                    VCF_Ref = Ref_allele #May be a dash or blank, but that is fine.
                if VCF_Alt == "-":
                    VCF_Alt = Alt_allele #May be a dash or blank, but that is fine.

            elif Var_Type == "Deletion":
                # This part depends on variant length
                if Var_Length == "-":
                    if Ref_allele != "-":
                        Var_Length = len(Ref_allele)
                    elif Position_g_start != "-" and Position_g_stop != "-":
                        try: # this is the only way I can think of to make sure these are actual numbers
                            Var_Length = (int(Position_g_stop)-int(Position_g_start) + 1)
                        except:
                            Var_Length = "-"
                if Var_Length == "-":
                    pass
                if Var_Length == "1":
                    if Position_g_start == "-":
                        if Genomic_annotation != "-":
                            search_result = re.search(r"g\.(\d+)", Genomic_annotation)
                            if search_result != None:
                                Position_g_start = search_result.group(1)
                    if Position_g_stop == "-":
                        Position_g_stop = Position_g_start ## may be blank
                    ### Can't fill in ref allele this way,
                    ### don't want to slow it down by using biocommons to get the ref allele
                    ### VCF info would need normalization
                else:
                    ### These are deletions of more than one base
                    if Position_g_start == "-":
                        if Genomic_annotation != "-":
                            search_result = re.search(r"g\.(\d+)", Genomic_annotation)
                            if search_result != None:
                                Position_g_start = search_result.group(1)
                    if Position_g_stop == "-":
                        if Genomic_annotation != "" and Genomic_annotation != "-":
                            search_result = re.search(r"(\d+)del", Genomic_annotation)
                            if search_result != None:
                                Position_g_stop = search_result.group(1)

                    ### I can't get normalized VCF info here
            elif Var_Type == "Duplication":
                # This part depends on variant length
                if Var_Length == "-":
                    if Ref_allele != "-":
                        Var_Length = len(Ref_allele)
                    # I have found that the genomic positions from ClinVar weren't completely reliable for duplications (didn't match HGVS perfectly)
                    # so I am not going to include the next part.

                if Var_Length == "-":
                    pass
                elif Var_Length == "1":
                    if Position_g_start == "-":
                        if Genomic_annotation != "" and Genomic_annotation != "-":
                            search_result = re.search(r"g\.(\d+)", Genomic_annotation)
                            if search_result != None:
                                Position_g_start = search_result.group(1)
                    if Position_g_stop == "-":
                        Position_g_stop = Position_g_start ## may be blank
                        ### Can't fill in ref allele this way,
                        ### don't want to slow it down by using biocommons to get the ref allele
                    if Alt_allele == "-" and Ref_allele != "-":
                        Alt_allele = Ref_allele+Ref_allele
                    ### VCF info would need normalization
                else:
                    ### These are deletions of more than one base
                    if Position_g_start == "-" or Position_g_start == "":
                        if Genomic_annotation != "" and Genomic_annotation != "-":
                            search_result = re.search(r"g\.(\d+)", Genomic_annotation)
                            if search_result != None:
                                Position_g_start = search_result.group(1)
                    if Position_g_stop == "-" or Position_g_stop == "":
                        if Genomic_annotation != "-":
                            search_result = re.search(r"(\d+)dup", Genomic_annotation)
                            if search_result != None:
                                Position_g_stop = search_result.group(1)
                    ### I can't get normalized VCF info here

            elif Var_Type == "Insertion":
                if Position_g_start == "-":
                    if Genomic_annotation != "-":
                        search_result = re.search(r"g\.(\d+)", Genomic_annotation)
                        if search_result != None:
                            Position_g_start = search_result.group(1)
                if Position_g_stop == "-" :
                    if Genomic_annotation != "-":
                        search_result = re.search(r"(\d+)ins", Genomic_annotation)
                        if search_result != None:
                            Position_g_stop = search_result.group(1)
                if Alt_allele == "-":
                    search_result = re.search(r"ins([A,C,T,G]+)", Genomic_annotation)
                    if search_result != None:
                        Alt_allele = search_result.group(1)
                if Var_Length == "-":
                    if Alt_allele != "-":
                        Var_Length = len(Alt_allele)

            elif Var_Type == "Indel":
                if Position_g_start == "-" or Position_g_start == "":
                    if Genomic_annotation != "-":
                        search_result = re.search(r"g\.(\d+)", Genomic_annotation)
                        if search_result != None:
                            Position_g_start = search_result.group(1)
                if Position_g_stop == "-":
                    if Genomic_annotation != "" and Genomic_annotation != "-":
                        search_result = re.search(r"(\d+)delins", Genomic_annotation)
                        if search_result != None:
                            Position_g_stop = search_result.group(1)
                if Alt_allele == "-":
                    search_result = re.search(r"ins([A,C,T,G]+)", Genomic_annotation) # because it ends in ins, we can use just this
                    if search_result != None:
                        Alt_allele = search_result.group(1)
                # I don't want to use biocommons to get the reference allele or the vcf info
                # I don't know how they determine the variant length
            # I don't need an else statement here because it would just be pass


            ###### One last shot to try to manually build the genomic annotation if it hasn't been found yet
            if Genomic_annotation == "-":
                if Var_Type == "single nucleotide variant":
                    if Chr_Accession != "-" and Position_g_start != "-" and Ref_allele != "-" and Alt_allele != "-":
                        Genomic_annotation = Chr_Accession+":g."+str(Position_g_start)+Ref_allele+">"+Alt_allele
                elif Var_Type == "Deletion":
                    if Var_Length == "-":
                        pass
                    elif Var_Length == "1":
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
            ###
            ### Fixing up just a couple of things that might have useful info hidden in another location
            if Gene_Symbol == "-":
                Gene_Symbol = Affected_Genes
            if Assembly == "-":
                if Genomic_annotation != "-":
                    Assembly = "GRCh37" # If something was picked up by these methods, then it has to be GRCh37
            if Chromosome == "23": # These are from if I got the chromosome from the chromosome accession
                Chromosome = "X"
            if Chromosome == "24":
                Chromosome = "Y"

            ###
            ### Now to normalize or validate the genomic annotation if I have one
            ###
            if Var_Type == "Complex":
                # No need to normalize, these ones have funky or missing genomic annotations
                Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                           Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                           Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                           Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                           Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                           "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                compound_variants.append(Individual_Variant_List)
            elif "copy number" in Var_Type:
                Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                           Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                           Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                           Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                           Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                           "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                copy_number_variants.append(Individual_Variant_List)
            elif Var_Type == "Microsatellite":
                Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                           Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                           Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                           Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                           Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                           "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                microsatellites.append(Individual_Variant_List)
            elif Genomic_annotation == "-":
                Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                           Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                           Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                           Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                           Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                           "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                not_annotated_variants.append(Individual_Variant_List)
            elif "?" in Genomic_annotation:
                Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                           Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                           Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                           Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                           Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                           "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                unknown_breakpoint.append(Individual_Variant_List)
            elif ")_(" in Genomic_annotation:
                Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                           Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                           Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                           Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                           Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                           "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                compound_variants.append(Individual_Variant_List)
            elif "=" in Genomic_annotation:
                Genomic_Normalized = Genomic_annotation
                Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                           Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                           Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                           Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                           Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                           "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                no_change_from_reference.append(Individual_Variant_List)
            elif "[" in Genomic_annotation or "(" in Genomic_annotation:
                # At this point anything with a square bracket or a parentheses should be a Microsatellite
                Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                           Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                           Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                           Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                           Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                           "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                microsatellites.append(Individual_Variant_List)
            elif Star_Level == 0:
                ## There will be a lot of these, especially those submitted by OMIM
                ## but I want to validate/normalize them anyway so that they can be used if we decide to include them
                if Var_Type == "single nucleotide variant":
                    try:
                        validated = vr.validate(hp.parse_hgvs_variant(Genomic_annotation))
                    except hgvs.exceptions.HGVSError:
                        validated = False
                    if validated == True:
                        Genomic_Normalized = Genomic_annotation
                        Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                   Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                   Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                   Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                   Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                   "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]

                        no_stars.append(Individual_Variant_List)
                    else:
                        ## This will leave Genomic_Normalized as a dash and save it to a different output
                        Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                   Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                   Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                   Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                   Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                   "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]

                        incorrect_base.append(Individual_Variant_List) # These ones have no stars, but will still be included in the incorrect base output
                else:
                    ## All variant types that are not SNV need to be normalized
                    match_ins = re.search(r'ins\d+', Genomic_annotation)
                    match_dup = re.search(r'(NC_\d+\.\d+:g\.\d+_?\d*dup)\d+', Genomic_annotation)
                    if match_ins != None:
                        # These are ones that end in ins followed by a number not a sequence
                        # These aren't the most common, but there are still enough of them to make this step necessary
                        search_alt_allele = re.search(r"([ACTG])([ACTG]*)", Alt_allele)
                        if search_alt_allele != None:
                            if search_alt_allele.group(1) == Ref_allele:
                                first_part_of_genomic = re.search(r"(.*ins)\d+", Genomic_annotation)
                                if first_part_of_genomic != None:
                                    new_genomic = first_part_of_genomic.group(1) + search_alt_allele.group(2)
                                    try:
                                        Genomic_annotation = new_genomic
                                        Genomic_Normalized = str(hn.normalize(hp.parse_hgvs_variant(Genomic_annotation)))
                                        Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                                   Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                                   Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                                   Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                                   Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                                   "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                                        no_stars.append(Individual_Variant_List)
                                    except hgvs.exceptions.HGVSError:
                                        # failed normalization
                                        Genomic_annotation = new_genomic
                                        Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                                   Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                                   Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                                   Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                                   Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                                   "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                                        invalid_hgvs.append(Individual_Variant_List)
                                else:
                                    # No info from first part of genomic. Don't know what else to do with it.
                                    Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                               Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                               Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                               Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                               Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                               "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                                    inserted_unknown_sequence.append(Individual_Variant_List)
                            else:
                                # These ones have a different first base and ref allele, probably a dash for ref allele.
                                first_part_of_genomic = re.search(r"(.*ins)\d+", Genomic_annotation)
                                if first_part_of_genomic != None:
                                    new_genomic = first_part_of_genomic.group(1) + Alt_allele
                                    try:
                                        Genomic_annotation = new_genomic
                                        Genomic_Normalized = str(hn.normalize(hp.parse_hgvs_variant(Genomic_annotation)))
                                        Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                                   Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                                   Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                                   Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                                   Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                                   "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                                        no_stars.append(Individual_Variant_List)
                                    except hgvs.exceptions.HGVSError:
                                        # failed normalization
                                        Genomic_annotation = new_genomic
                                        Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                                   Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                                   Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                                   Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                                   Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                                   "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                                        invalid_hgvs.append(Individual_Variant_List)
                                else:
                                    # No returned search result from the first_part_of_genomic regex
                                    Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                               Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                               Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                               Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                               Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                               "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                                    inserted_unknown_sequence.append(Individual_Variant_List)

                        else:
                            # These are the ones that have no alternate allele
                            Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                       Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                       Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                       Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                       Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                       "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                            inserted_unknown_sequence.append(Individual_Variant_List)

                    elif match_dup != None:
                        # These are the duplicates that end in a number, which biocmmons can't handle
                        try:
                            Genomic_Normalized = str(hn.normalize(hp.parse_hgvs_variant(match_dup.group(1))))
                            Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                       Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                       Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                       Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                       Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                       "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                            no_stars.append(Individual_Variant_List)
                        except hgvs.exceptions.HGVSError:
                            Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                       Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                       Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                       Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                       Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                       "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                            invalid_hgvs.append(Individual_Variant_List)
                    else:
                        # This is everything else, should be the majority of non-SNV variants with no stars
                        try:
                            Genomic_Normalized = str(hn.normalize(hp.parse_hgvs_variant(Genomic_annotation)))
                            Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                       Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                       Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                       Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                       Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                       "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                            no_stars.append(Individual_Variant_List)
                        except hgvs.exceptions.HGVSError:
                            # Many of these are the insertions that list the start and stop positions as the same thing
                            # Biocommons can't handle that, but it can be fixed
                            if Var_Type == "Insertion":
                                match_failed_ins = re.search(r'(NC_\d+\.\d+:g\.)(\d+)_(\d+)(ins[A,C,T,G]+)', Genomic_annotation)
                                if match_failed_ins != None:
                                    if match_failed_ins.group(2) == match_failed_ins.group(3):
                                        new_stop = int(match_failed_ins.group(2)) + 1
                                        biocommons_input = match_failed_ins.group(1)+str(match_failed_ins.group(2))+"_"+str(new_stop)+match_failed_ins.group(4)
                                        try:
                                            Genomic_Normalized = str(hn.normalize(hp.parse_hgvs_variant(biocommons_input)))
                                            Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                                       Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                                       Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                                       Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                                       Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                                       "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                                            no_stars.append(Individual_Variant_List)
                                        except hgvs.exceptions.HGVSError:
                                            Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                                       Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                                       Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                                       Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                                       Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                                       "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                                            invalid_hgvs.append(Individual_Variant_List)
                                    else:
                                        # If the numbers were not the same, it is an invalid_hgvs
                                        Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                                   Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                                   Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                                   Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                                   Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                                   "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                                        invalid_hgvs.append(Individual_Variant_List)
                                else:
                                    # If there was no regex match, it is an invalid_hgvs
                                    Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                               Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                               Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                               Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                               Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                               "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                                    invalid_hgvs.append(Individual_Variant_List)
                            else:
                                # If the variant type is not insertion, then I don't know how to save it at this point
                                Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                           Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                           Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                           Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                           Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                           "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                                invalid_hgvs.append(Individual_Variant_List)




                    #  check for funky duplicates with a number at the ends
                    # insertions where start and stop position are the same
                    # check for insertions or indels that end in a number
            else:
                ## At this point all of the ones that have something odd in the genomic annotation should be weeded out
                ## And the ones with no stars have been weeded out
                ## SNVs don't need to normalize, just validate
                if Var_Type == "single nucleotide variant":
                    try:
                        validated = vr.validate(hp.parse_hgvs_variant(Genomic_annotation))
                    except hgvs.exceptions.HGVSError:
                        validated = False
                    if validated == True:
                        Genomic_Normalized = Genomic_annotation
                        Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                   Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                   Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                   Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                   Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                   "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]

                        ClinVar_Variants.append(Individual_Variant_List)
                    else:
                        ## This will leave Genomic_Normalized as a dash and save it to a different output
                        Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                   Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                   Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                   Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                   Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                   "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]

                        invalid_hgvs.append(Individual_Variant_List) # These ones have no stars, but will still be included in the incorrect base output
                else:
                    ## All variant types that are not SNV need to be normalized
                    match_ins = re.search(r'ins\d+', Genomic_annotation)
                    match_dup = re.search(r'(NC_\d+\.\d+:g\.\d+_?\d*dup)\d+', Genomic_annotation)
                    if match_ins != None:
                        # These are ones that end in ins followed by a number not a sequence
                        # These aren't the most common, but there are still enough of them to make this step necessary
                        search_alt_allele = re.search(r"([ACTG])([ACTG]*)", Alt_allele)
                        if search_alt_allele != None:
                            if search_alt_allele.group(1) == Ref_allele:
                                first_part_of_genomic = re.search(r"(.*ins)\d+", Genomic_annotation)
                                if first_part_of_genomic != None:
                                    new_genomic = first_part_of_genomic.group(1) + search_alt_allele.group(2)
                                    try:
                                        Genomic_annotation = new_genomic
                                        Genomic_Normalized = str(hn.normalize(hp.parse_hgvs_variant(Genomic_annotation)))
                                        Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                                   Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                                   Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                                   Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                                   Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                                   "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                                        ClinVar_Variants.append(Individual_Variant_List)
                                    except hgvs.exceptions.HGVSError:
                                        # failed normalization
                                        Genomic_annotation = new_genomic
                                        Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                                   Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                                   Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                                   Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                                   Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                                   "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                                        invalid_hgvs.append(Individual_Variant_List)
                                else:
                                    # No info from first part of genomic. Don't know what else to do with it.
                                    Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                               Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                               Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                               Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                               Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                               "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                                    inserted_unknown_sequence.append(Individual_Variant_List)
                            else:
                                # These ones have a different first base and ref allele, probably a dash for ref allele.
                                first_part_of_genomic = re.search(r"(.*ins)\d+", Genomic_annotation)
                                if first_part_of_genomic != None:
                                    new_genomic = first_part_of_genomic.group(1) + Alt_allele
                                    try:
                                        Genomic_annotation = new_genomic
                                        Genomic_Normalized = str(hn.normalize(hp.parse_hgvs_variant(Genomic_annotation)))
                                        Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                                   Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                                   Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                                   Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                                   Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                                   "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                                        ClinVar_Variants.append(Individual_Variant_List)
                                    except hgvs.exceptions.HGVSError:
                                        # failed normalization
                                        Genomic_annotation = new_genomic
                                        Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                                   Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                                   Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                                   Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                                   Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                                   "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                                        invalid_hgvs.append(Individual_Variant_List)
                                else:
                                    # No returned search result from the first_part_of_genomic regex
                                    Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                               Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                               Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                               Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                               Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                               "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                                    inserted_unknown_sequence.append(Individual_Variant_List)

                        else:
                            # These are the ones that have no alternate allele
                            Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                       Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                       Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                       Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                       Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                       "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                            inserted_unknown_sequence.append(Individual_Variant_List)

                    elif match_dup != None:
                        # These are the duplicates that end in a number, which biocmmons can't handle
                        # It must still be valid HGVS because there are a lot of these
                        try:
                            Genomic_Normalized = str(hn.normalize(hp.parse_hgvs_variant(match_dup.group(1))))
                            Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                       Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                       Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                       Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                       Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                       "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                            ClinVar_Variants.append(Individual_Variant_List)
                        except hgvs.exceptions.HGVSError:
                            Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                       Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                       Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                       Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                       Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                       "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                            invalid_hgvs.append(Individual_Variant_List)
                    else:
                        # This is everything else, should be the majority of non-SNV variants with no stars
                        try:
                            Genomic_Normalized = str(hn.normalize(hp.parse_hgvs_variant(Genomic_annotation)))
                            Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                       Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                       Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                       Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                       Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                       "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                            ClinVar_Variants.append(Individual_Variant_List)
                        except hgvs.exceptions.HGVSError:
                            # Many of these are the insertions that list the start and stop positions as the same thing
                            # Biocommons can't handle that, but it can be fixed
                            if Var_Type == "Insertion":
                                match_failed_ins = re.search(r'(NC_\d+\.\d+:g\.)(\d+)_(\d+)(ins[A,C,T,G]+)', Genomic_annotation)
                                if match_failed_ins != None:
                                    if match_failed_ins.group(2) == match_failed_ins.group(3):
                                        new_stop = int(match_failed_ins.group(2)) + 1
                                        biocommons_input = match_failed_ins.group(1)+str(match_failed_ins.group(2))+"_"+str(new_stop)+match_failed_ins.group(4)
                                        try:
                                            Genomic_Normalized = str(hn.normalize(hp.parse_hgvs_variant(biocommons_input)))
                                            Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                                       Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                                       Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                                       Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                                       Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                                       "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                                            ClinVar_Variants.append(Individual_Variant_List)
                                        except hgvs.exceptions.HGVSError:
                                            Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                                       Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                                       Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                                       Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                                       Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                                       "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                                            invalid_hgvs.append(Individual_Variant_List)
                                    else:
                                        # If the numbers were not the same, it is an invalid_hgvs
                                        Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                                   Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                                   Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                                   Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                                   Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                                   "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                                        invalid_hgvs.append(Individual_Variant_List)
                                else:
                                    # If there was no regex match, it is an invalid_hgvs
                                    Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                               Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                               Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                               Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                               Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                               "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                                    invalid_hgvs.append(Individual_Variant_List)
                            else:
                                # If the variant type is not insertion, then I don't know how to save it at this point
                                Individual_Variant_List = [Assembly, Chromosome, Position_g_start, Position_g_stop, Ref_allele,
                                                           Alt_allele, Genomic_annotation, Genomic_Normalized, Var_Type, Var_Length, Pathogenicity,
                                                           Disease, Genetic_Origin, Inheritance_Pattern, Affected_Genes, Gene_Symbol,
                                                           Compound_Het,Transcript, Transcript_notation, Transcript_HGVS, Protein_Accession,
                                                           Protein_Notation, Protein_HGVS, Chr_Accession, Pos_VCF, VCF_Ref, VCF_Alt,
                                                           "ClinVar", RCV_num, Review_Status, Star_Level, Submitter, Edited_Date]
                                invalid_hgvs.append(Individual_Variant_List)


        elem.clear()
#del elem # Clear the memory once an RCV of an important gene is processed


    column_labels = ["Genome Assembly", "Chr", "Position Start", "Position Stop", "Ref", "Alt",
                    "Genomic Annotation", "HGVS Normalized Genomic Annotation", "Variant Type",
                    "Variant Length", "Pathogenicity", "Disease", "Genetic Origin", "Inheritance Pattern",
                    "Affected Genes" , "Gene Symbol", "Compound Het Status", "Transcript",
                    "Transcript Notation", "HGVS Transcript Notation","Protein Accession",
                    "Protein Notation", "HGVS Protein Annotation", "Chr Accession", "VCF Pos",
                    "VCF Ref", "VCF Alt", "Database", "ClinVar Accession", "Review Status",
                    "Star Level", "Submitter", "Edited Date"]
    if len(not_annotated_variants) > 0 :
        not_annotated_variants_DF = pd.DataFrame(not_annotated_variants)
        not_annotated_variants_DF.columns = column_labels
        not_annotated_variants_DF.to_csv(output_directory+"/ClinVar/"+disease+"/Invalid_Annotations/"+disease+"_not_annotated_variants_ClinVar.csv")

    if len(compound_variants) > 0 :
        compound_variants_DF = pd.DataFrame(compound_variants)
        compound_variants_DF.columns = column_labels
        compound_variants_DF.to_csv(output_directory+"/ClinVar/"+disease+"/Invalid_Annotations/"+disease+"_compound_variants_ClinVar.csv")

    if len(copy_number_variants) > 0 :
        copy_number_variants_DF = pd.DataFrame(copy_number_variants)
        copy_number_variants_DF.columns = column_labels
        copy_number_variants_DF.to_csv(output_directory+"/ClinVar/"+disease+"/Invalid_Annotations/"+disease+"_copy_number_variants_ClinVar.csv")

    if len(microsatellites) > 0 :
        microsatellites_DF = pd.DataFrame(microsatellites)
        microsatellites_DF.columns = column_labels
        microsatellites_DF.to_csv(output_directory+"/ClinVar/"+disease+"/Invalid_Annotations/"+disease+"_microsatellites_ClinVar.csv")

    if len(unknown_breakpoint) > 0 :
        unknown_breakpoint_DF = pd.DataFrame(unknown_breakpoint)
        unknown_breakpoint_DF.columns = column_labels
        unknown_breakpoint_DF.to_csv(output_directory+"/ClinVar/"+disease+"/Invalid_Annotations/"+disease+"_unknown_breakpoint_variants_ClinVar.csv")

    if len(no_change_from_reference) > 0 :
        no_change_from_reference_DF = pd.DataFrame(no_change_from_reference)
        no_change_from_reference_DF.columns = column_labels
        no_change_from_reference_DF.to_csv(output_directory+"/ClinVar/"+disease+"/Invalid_Annotations/"+disease+"_no_change_from_reference_ClinVar.csv")

    if len(inserted_unknown_sequence) > 0 :
        inserted_unknown_sequence_DF = pd.DataFrame(inserted_unknown_sequence)
        inserted_unknown_sequence_DF.columns = column_labels
        inserted_unknown_sequence_DF.to_csv(output_directory+"/ClinVar/"+disease+"/Invalid_Annotations/"+disease+"_inserted_unknown_sequence_ClinVar.csv")

    if len(invalid_hgvs) > 0 :
        invalid_hgvs_DF = pd.DataFrame(invalid_hgvs)
        invalid_hgvs_DF.columns = column_labels
        invalid_hgvs_DF.to_csv(output_directory+"/ClinVar/"+disease+"/Invalid_Annotations/"+disease+"_invalid_hgvs_variants_ClinVar.csv")

    if len(no_stars) > 0 :
        no_stars_DF = pd.DataFrame(no_stars)
        no_stars_DF.columns = column_labels
        no_stars_DF.to_csv(output_directory+"/ClinVar/"+disease+"/Invalid_Annotations/"+disease+"_no_star_variants_ClinVar.csv")


    # ClinVar_Variants should always have something, so I am not checking for length here
    ClinVar_DF = pd.DataFrame(ClinVar_Variants)
    ClinVar_DF.columns= column_labels
    ClinVar_DF.to_csv(output_directory+"/ClinVar/"+disease+"/"+disease+"_Variants_ClinVar.csv")

for num, disease in enumerate(disease_names):
    parse_clinvar_xml(disease, input_gene_lists[num], ClinVar_File, output_directory)
