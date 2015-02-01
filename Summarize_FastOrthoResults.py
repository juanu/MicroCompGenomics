#Created on 1/2/2015
__author__ = 'Juan A. Ugalde'

import os
import sys
import argparse
from tools.data_input import GenomeData
from collections import defaultdict
import re


def get_protein_info(gl, fd):
    """
    This function is used to parse the fasta files. It will read a folder with fasta files (extension .fasta)
    and it will get the information for each protein in the genomes, including protein ID and length.
    :param gl: List of genomes to use.
    :param fd: Fasta directory with the protein files (.fasta)
    :return:
    """

    from Bio import SeqIO
    from collections import defaultdict

    proteins_in_genomes = defaultdict(list)  # output list of the proteins that are present in each genome
    protein_length = defaultdict(int)  # length of each protein. Assumes a unique ID for each protein

    file_counter = 0  # Count the number of files being read

    fasta_files = [fd + "/" + input_fasta + ".fasta" for input_fasta in gl]

    for input_fasta in fasta_files:
        for record in SeqIO.parse(input_fasta, "fasta"):
            proteins_in_genomes[record.id.split('|')[0]].append(record.id)
            protein_length[record.id] = int(len(record.seq))

        file_counter += 1

    return proteins_in_genomes, protein_length, file_counter


def parse_fastortho(cf, gl):
    """
    Function used to parse the FastOrtho results
    :param cf: cluster file (usually with extension .end)
    :param gl: List of genomes to use
    :return:
    """

    fastortho_results = open(cf, 'r')

    cluster_dictionary = defaultdict(list)  # Dictionary with the clusters. Stores only the clusters that will be used
    unique_proteins_genome_count = defaultdict(int)  # Proteins that are not in any cluster
    proteins_in_cluster = set()  # Proteins in cluster. Needed to search for proteins that are absent
    total_cluster_count = 0
    clusters_removed = 0

    for line in fastortho_results:
        total_cluster_count += 1
        line = line.strip('\n')
        ids_proteins = re.split(":\s+", line)
        proteins = ids_proteins[1].split(" ")

        clean_protein_list = []  # This is used to remove proteins from genomes not in the list

        for genome in gl:  # Adding the proteins that we need
            [clean_protein_list.append(protein) for protein in [x for x in proteins if x.startswith(genome)]]

        #Now I need to evaluate those clusters that now are unique and zero
        if len(clean_protein_list) == 0:
            clusters_removed += 1
            continue

        if len(clean_protein_list) == 1:
            unique_proteins_genome_count[clean_protein_list[0].split("|")[0]] += 1
            clusters_removed += 1
            continue

        for protein in clean_protein_list:
            cluster_dictionary[ids_proteins[0]].append(protein)
            proteins_in_cluster.add(protein)

    return cluster_dictionary, proteins_in_cluster, unique_proteins_genome_count, total_cluster_count, clusters_removed