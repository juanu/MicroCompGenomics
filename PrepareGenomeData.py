#Created on 1/27/2015

__author__ = 'Juan A. Ugalde'

#TODO Implement RAST and NCBI importers

import os
import argparse

from tools.data_input import GenomeData

program_description = "This script takes a list of genomes, provided on a genome list file, and a folder with genome" \
                      "information derived from different sources. The output are four folders:\n" \
                      "protein, nucleotide, genome and annotation\n" \
                      "\n" \
                      "The format of the genome list is a three column file with: genome ID, genome name, " \
                      "genome prefix\n" \
                      "Annotation sources are:\n" \
                      "jgi: Data dump obtained from the JGI website (directly or using Globus)\n" \
                      "img_single: Data obtained by downloading each genome directly from IMG, using the generate\n" \
                      "genbank option\n" \
                      "rast: Genome annotated with RAST (not implemented yet)\n" \
                      "ncbi: Genome obtained from ncbi (not implemented yet)\n"

parser = argparse.ArgumentParser(description=program_description)

parser.add_argument("-g", "--genome_list", type=str,
                    help="Tabular file with the list of genomes to include and their prefix\n", required=True)
parser.add_argument("-i", "--input_folder", type=str,
                    help="input folder with the genomes. Within this folder each genome needs to have its own folder,"
                         "with the name being the genome ID used on the genome list", required=True)
parser.add_argument("-s", "--genome_source", help="Source of the genome annotation (jgi, img_single, rast, ncbi)",
                    required=True, type=str)
parser.add_argument("-o", "--output_directory", type=str, help="Output folder", required=True)

args = parser.parse_args()

#Create output directory

if not os.path.exists(args.output_directory):
    os.makedirs(args.output_directory)
    os.makedirs(args.output_directory + "/nucleotide")
    os.makedirs(args.output_directory + "/protein")
    os.makedirs(args.output_directory + "/genome")
    os.makedirs(args.output_directory + "/annotation")
    os.makedirs(args.output_directory + "/coords")

#Read the genome list and create a dictionary with the information, and the total genome count
genome_dictionary, total_genome_count = GenomeData.read_genome_list(args.genome_list)


#Read the source option, and process the genome data

processed_genomes = 0

if args.genome_source == "jgi":
    processed_genomes = GenomeData.parse_jgi_dump(genome_dictionary, args.input_folder, args.output_directory)

if args.genome_source == "img_single":
    processed_genomes = GenomeData.parse_single_img(genome_dictionary, args.input_folder, args.output_directory)

print "Done processing %d genomes" % total_genome_count
print "A total of %d genomes were found in the input folder" % processed_genomes

