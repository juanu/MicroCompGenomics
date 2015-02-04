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
            cluster_id = ids_proteins[0].split(" (")[0]

            cluster_dictionary[cluster_id].append(protein)
            proteins_in_cluster.add(protein)

    return cluster_dictionary, proteins_in_cluster, unique_proteins_genome_count, total_cluster_count, clusters_removed


def read_group_files(gf):
    """
    Reads a file with the group list and returns a dictionary containing the name of the group
    and a list of the genomes that are present in that group.
    The input file requires that the first column is the prefix of the genome and the second the name
    of the group
    :param gf: Group file
    :return: Genome groups
    """

    from collections import defaultdict

    genome_groups = defaultdict(list)

    for line in open(gf, 'r'):
        line = line.rstrip()
        element = line.split("\t")

        genome_groups[element[0]].append(element[1])

    return genome_groups


def get_unique_seqs_genome(protein_genomes, protein_clusters, protein_length, min_length):
    """
    This module takes a dictionary with their genome and proteins, and a set of proteins
    and look for proteins that are not in the total set
    """

    from collections import defaultdict

    large_unique_sequences = defaultdict(list)
    short_unique_sequences = defaultdict(list)
    processed_sequences = 0

    for genome in protein_genomes:
        for seq_name in protein_genomes[genome]:
            processed_sequences += 1

            if seq_name in protein_clusters:
                continue
            else:
                if protein_length[seq_name] > min_length:
                    large_unique_sequences[genome].append(seq_name)
                else:
                    short_unique_sequences[genome].append(seq_name)

    return large_unique_sequences, short_unique_sequences, processed_sequences


def seqs_shared_clusters(cd, gd):
    """
    :param cd: Cluster dictionary
    :param gd: Genome dictionary
    :return: Genome groups
    """

    unique_clusters = defaultdict(list)  # Dictionary with the unique clusters
    shared_single_clusters = []  # List with clusters that are shared and single copy
    shared_multiple_clusters = []  # List with clusters that are shared and in multiple copies

    genomes_in_matrix = sorted(gd.keys())
    header = ["Cluster_ID"]
    header.extend(sorted(gd.keys()))
    all_clusters_matrix = [header]

    for cluster in cd:

        genome_list = [protein.split("|")[0] for protein in cd[cluster]]  # Create a list with the genomes

        count = {x: genome_list.count(x) for x in genome_list}  # Count the occurences

        #Create the matrix
        cluster_matrix = [cluster]

        for genome in genomes_in_matrix:

            if genome in genome_list:
                cluster_matrix.append(count[genome])
            else:
                cluster_matrix.append(0)

        #print cluster_matrix
        all_clusters_matrix.append(cluster_matrix)

        if len(count) == 1:
            unique_clusters[genome_list[0]].append(cluster)

        elif len(count) == len(gd.keys()):
            if sum(count.itervalues()) == len(gd.keys()):
                shared_single_clusters.append(cluster)
            else:
                shared_multiple_clusters.append(cluster)

    return unique_clusters, shared_single_clusters, shared_multiple_clusters, all_clusters_matrix


def clusters_in_groups(clusters, groups):
    """

    """
    from collections import defaultdict

    import itertools
    unique_group_clusters = defaultdict(list)  # Count the unique clusters in each group

    combination_clusters = defaultdict(list)

    #Create inverted dictionary
    genome_group_info = defaultdict()

    for group in groups:
        for genome in groups[group]:
            genome_group_info[genome] = group

    for cluster in clusters:

        group_count = defaultdict(lambda: defaultdict(int))

        for protein in clusters[cluster]:
            genome_id = protein.split("|")[0]
            group_for_genome = genome_group_info[genome_id]

            group_count[group_for_genome][genome_id] += 1

        ##Unique clusters for each group

        if len(group_count) == 1:
            for group in group_count:
                if len(group_count[group]) == len(groups[group]):
                    unique_group_clusters[group].append(cluster)

                else:  # I could add something here to count the number of proteins not unique
                    pass

        #Shared, all possible combinations
        else:
            for combination_count in range(2, len(group_count) + 1):
                for group_combinations in itertools.combinations(groups.keys(), combination_count):

                    # Check that all the genomes in the group are represented
                    group_check = 0

                    for i in range(0, len(group_combinations)):
                        if len(group_count[group_combinations[i]]) == len(groups[group_combinations[i]]):
                            continue
                        else:
                            group_check = 1

                    if group_check == 0:
                        combination_clusters[group_combinations].append(cluster)

    return unique_group_clusters, combination_clusters

