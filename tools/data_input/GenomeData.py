#Created on 1/27/2015
__author__ = 'Juan A. Ugalde'

import sys

def read_genome_list(input_file):
    """
    This function reads the genome list, and process that information. The format is:

    First column, ID of the genome depending the annotation source (IMG, JGI, Genbank)

    Second column: Full species name

    Third column: Prefix to use in the analysis

    :param input_file: Genome list
    :return: Dictionary with genome information, and total count
    """

    genome_count = 0
    genome_info = {}  # Dictionary with the genome information
    genome_name_list = []  # List of genome names
    genome_prefix_list = []  # List of genome prefix

    for line in open(input_file, 'r'):
        if line.strip():
            line = line.rstrip()
            element = line.split("\t")

            #Check that there are 3 elements in the file. It should only have three columns
            if not len(element) == 3:
                error_mesage = "There are more than 3 columns in one line of the file. Double check the file\n " \
                               "The error line is" + line
                sys.exit(error_mesage)

            #Check for duplicates
            if element[0] in genome_info.keys():  # Duplicate genome ID
                print "Duplicate genome_id found: " + line
                sys.exit("Check for duplicates")

            elif element[1] in genome_name_list:  # Duplicate genome name
                print "Duplicate genome name found: " + line
                sys.exit("Check for duplicates")

            elif element[2] in genome_prefix_list:  # Duplicate prefix name
                print "Duplicate prefix found: " + line
                sys.exit("Check for duplicates")

            else:
                genome_info[element[2]] = element[0]
                genome_count += 1
                genome_name_list.append(element[1])
                genome_prefix_list.append(element[2])

    return genome_info, genome_count

def modify_fasta(input_fasta_file, new_fasta_prefix, output_file):
    """
The input is a fasta file. This will rename the file with the new ID,
and modify each ID in each entry of the fasta file
"""

    from Bio import SeqIO  # Import from tools to read fasta, from Biopython

    #Rename the fasta
    for record in SeqIO.parse(input_fasta_file, "fasta"):
        output_file.write(">" + new_fasta_prefix + "|" + record.id + "\n")
        output_file.write(str(record.seq) + "\n")