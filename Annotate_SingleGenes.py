#Created on 11/5/15
author = "Juan A. Ugalde"

from tools.data_input import AnnotationData
import argparse
from collections import defaultdict
from tools.AnnotationInfo import cog_definitions
import os

program_description = "This script annotates a list of genes, where the first column is the genome name" \
                      "and the second is the gene ID. The annotation is taken from the folder generated by" \
                      "the import scripts"

parser = argparse.ArgumentParser(description=program_description)

#Arguments
parser.add_argument("-d", "--data_folders", type=str,
                    help="Folder with the processed genomes", required=True)
parser.add_argument("-c", "--gene_file", type=str,
                    help="Gene file", required=True)
parser.add_argument("-o", "--output_file", type=str,
                    help="Output file", required=True)
parser.add_argument("-g", "--genome_list", type=str,
                    help="Genome list", required=True)

args = parser.parse_args()

#Get the genome information
species_name_dict = defaultdict()
for line in open(args.genome_list, 'r'):
    line = line.rstrip()
    img_id, genome_name, short_name = line.split("\t")
    species_name_dict[short_name] = genome_name

#Read the gene list

genome_gene_info = defaultdict(list)

for line in open(args.gene_file):
    if line.strip():
        line = line.rstrip()
        genome_gene_info[line.split("\t")[0]].append(line.split("\t")[1])


#Get the coords for all the proteins and store in a dictionary
gene_coords = defaultdict(tuple)
coord_folder = args.data_folders + "/coords/"

for coord_file in os.listdir(coord_folder):
    for line in open(coord_folder + coord_file, 'r'):
        line = line.rstrip()
        contig_id, protein_id, start, stop = line.split("\t")
        gene_coords[protein_id] = (contig_id, start, stop)

#Get the annotation information
annotation_folder = args.data_folders + "/annotation"

protein_annotation, function_definitions = \
    AnnotationData.parse_annotation_folder(genome_gene_info.keys(), annotation_folder)

#Print output table
output_file = open(args.output_file, 'w')

#Get the COG definitions
cog_one_letter, desc_cog_letter, desc_cog_number = cog_definitions()

for genome in genome_gene_info:
    for protein in genome_gene_info[genome]:

        try:
            product = protein_annotation[protein]["Product"]
        except KeyError:
            product = None

        try:
            COG_number = protein_annotation[protein]["COG"]
        except KeyError:
            COG_number = None

        COG_description = None

        if not COG_number is None:
            try:
                COG_description = function_definitions[COG_number]
            except KeyError:
                COG_description = None

        try:
            PFAM_number = protein_annotation[protein]["PFAM"]
        except KeyError:
            PFAM_number = None

        PFAM_description = None
        if not PFAM_number is None:
            PFAM_description = function_definitions[PFAM_number]

        #Get COG one letter code
        try:
            letter_cog = cog_one_letter[COG_number]
        except KeyError:
            letter_cog = None

        #Description
        try:
            letter_description = desc_cog_letter[letter_cog[0]]
        except IndexError:
            letter_description = None

        #Get the locus tag
        try:
            locus_tag = protein_annotation[protein]["locus_tag"]
        except KeyError:
            locus_tag = None

        #Get the protein coords
        try:
            contig_id, start, stop = gene_coords[protein]
        except KeyError:
            contig_id, start, stop = None, None, None

        output_line = [species_name_dict[genome], genome, contig_id, protein, start, stop, locus_tag, product,
                       COG_number, COG_description, letter_cog, letter_description, PFAM_number, PFAM_description]

        output_file.write("\t".join(str(x) for x in output_line) + "\n")

#print protein_annotation
#print function_definitionsw