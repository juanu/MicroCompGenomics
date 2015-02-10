#Created on 1/2/2015
__author__ = 'Juan A. Ugalde'

from collections import defaultdict


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

    gg = defaultdict(list)

    for line in open(gf, 'r'):
        line = line.rstrip()
        element = line.split("\t")

        gg[element[0]].append(element[1])

    return gg


def get_unique_seqs_genome(protein_genomes, protein_clusters, protein_length, min_length):
    """
    This module takes a dictionary with their genome and proteins, and a set of proteins
    and look for proteins that are not in the total set
    """

    from collections import defaultdict

    large_unique_sequences = defaultdict(list)
    short_unique_sequences = defaultdict(list)
    processed_sequences = 0

    for genome_entry in protein_genomes:
        for seq_name in protein_genomes[genome_entry]:
            processed_sequences += 1

            if seq_name in protein_clusters:
                continue
            else:
                if protein_length[seq_name] > min_length:
                    large_unique_sequences[genome_entry].append(seq_name)
                else:
                    short_unique_sequences[genome_entry].append(seq_name)

    return large_unique_sequences, short_unique_sequences, processed_sequences


def seqs_shared_clusters(cd, gd):
    """
    :param cd: Cluster dictionary
    :param gd: Genome dictionary
    :return: Genome groups
    """

    u_clusters = defaultdict(list)  # Dictionary with the unique clusters
    ss_clusters = []  # List with clusters that are shared and single copy
    sm_clusters = []  # List with clusters that are shared and in multiple copies

    genomes_in_matrix = sorted(gd.keys())
    header = ["Cluster_ID"]
    header.extend(sorted(gd.keys()))
    total_cluster_matrix = [header]

    for single_cluster in cd:

        genome_list = [prot.split("|")[0] for prot in cd[single_cluster]]  # Create a list with the genomes

        occurence_count = {x: genome_list.count(x) for x in genome_list}  # Count the occurences

        #Create the matrix
        cluster_matrix = [single_cluster]

        for genome_entry in genomes_in_matrix:

            if genome_entry in genome_list:
                cluster_matrix.append(occurence_count[genome_entry])
            else:
                cluster_matrix.append(0)

        #print cluster_matrix
        total_cluster_matrix.append(cluster_matrix)

        if len(occurence_count) == 1:
            u_clusters[genome_list[0]].append(single_cluster)

        elif len(occurence_count) == len(gd.keys()):
            if sum(occurence_count.itervalues()) == len(gd.keys()):
                ss_clusters.append(single_cluster)
            else:
                sm_clusters.append(single_cluster)

    return u_clusters, ss_clusters, sm_clusters, total_cluster_matrix


def clusters_in_groups(clusters, groups):
    """

    """
    from collections import defaultdict

    import itertools
    ug_clusters = defaultdict(list)  # Count the unique clusters in each group

    comb_clusters = defaultdict(list)

    #Create inverted dictionary
    genome_group_info = defaultdict()

    for entry in groups:
        for genome_entry in groups[entry]:
            genome_group_info[genome_entry] = entry

    for single_cluster in clusters:

        group_count = defaultdict(lambda: defaultdict(int))

        for protein_entry in clusters[single_cluster]:
            genome_id = protein_entry.split("|")[0]
            group_for_genome = genome_group_info[genome_id]

            group_count[group_for_genome][genome_id] += 1

        ##Unique clusters for each group

        if len(group_count) == 1:
            for grp in group_count:
                if len(group_count[grp]) == len(groups[grp]):
                    ug_clusters[grp].append(single_cluster)

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
                        comb_clusters[group_combinations].append(single_cluster)

    return ug_clusters, comb_clusters

if __name__ == '__main__':
    import os
    import sys
    import argparse
    from tools.data_input.GenomeData import read_genome_list
    from tools.data_input.ClusterInput import parse_fastortho, parse_tribemcl

    #Create the options and program description
    program_description = "This script summarize the results of orthoMCL, and create several summary files." \
                          " The inputs are:" \
                          "-List of clusters, generated by orthoMCL" \
                          "-A genome list" \
                          "-A folder with fasta files" \
                          "- An optional group file, to group genomes. For example, all genomes from the same species "

    parser = argparse.ArgumentParser(description=program_description)

    parser.add_argument("-l", "--genome_list_index", type=str,
                        help="File with the genome list. Format GenomeID, FullName, ShortName", required=True)
    parser.add_argument("-c", "--cluster_file", type=str, help="Ortholog file, generated by OrthoMCL", required=True)
    parser.add_argument("-f", "--fasta_aa_directory", type=str, help="Directory with the fasta files", required=True)
    parser.add_argument("-t", "--clustering_algorithm", type=str, help="Type of clustering used. Options are fastortho,"
                                                                       "tribemcl", required=True)
    parser.add_argument("-g", "--group_information", type=str, help="Group file")
    parser.add_argument("-o", "--output_directory", type=str, help="Output directory", required=True)

    args = parser.parse_args()

    #Create the output directory
    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)

    #Create a log file
    run_summary = open(args.output_directory + "/logfile.txt", 'w')

    #####Read the genome list
    genome_id_dictionary, genome_count = read_genome_list(args.genome_list_index)

    run_summary.write("Genomes in the genome list: %d" % genome_count + "\n")

    ######Read the cluster information, and check that everything is ok
    #cluster_information, set_of_proteins_in_clusters, unique_cluster_count, total_clusters, removed_clusters = \
    #    get_orthomcl_results(args.cluster_file, [i for i in genome_id_dictionary.itervalues()])

    #Create the output variables
    cluster_information = None
    set_of_proteins_in_clusters = None
    unique_cluster_count = None
    total_clusters = None
    removed_clusters = None

    if args.clustering_algorithm == "fastortho":
        cluster_information, set_of_proteins_in_clusters, unique_cluster_count, total_clusters, removed_clusters = \
            parse_fastortho(args.cluster_file, genome_id_dictionary.keys())

    elif args.clustering_algorithm == "tribemcl":
           cluster_information, set_of_proteins_in_clusters, unique_cluster_count, total_clusters, removed_clusters = \
               parse_tribemcl(args.cluster_file, genome_id_dictionary.keys())


    run_summary.write("Total number of clusters: %d" % len(cluster_information) + "\n")
    run_summary.write("Total number of protein in clusters: %d" % len(set_of_proteins_in_clusters) + "\n")
    run_summary.write("Total number of removed clusters (not present in the genome file): %d" % removed_clusters + "\n")

    #Check the counts, to see if everything is going ok
    if total_clusters - removed_clusters != len(cluster_information):
        sys.exit("The number of removed clusters clusters plus the retained clusters, "
                 "doesn't match the total of original clusters in the file")

    #####Read the fasta file
    #dic_protein_in_genomes, dic_protein_length, files_read_counter = \
    #    get_protein_info([i for i in genome_id_dictionary.itervalues()], args.fasta_aa_directory)

    dic_protein_in_genomes, dic_protein_length, files_read_counter = \
        get_protein_info(genome_id_dictionary.keys(), args.fasta_aa_directory)

    run_summary.write("Total fasta files: %d" % files_read_counter + "\n")
    run_summary.write("Total number of proteins in the fasta files: %d" % len(dic_protein_length) + "\n")

    #####Read the genome groups (if present)
    genome_groups = {}
    if args.group_information:
        genome_groups = read_group_files(args.group_information)
        run_summary.write("Total defined groups: %d" % len(genome_groups) + "\n")

    #####################################################
    #Look for unique protein in each genome
    selected_unique_proteins, removed_unique_proteins, total_number_proteins = \
        get_unique_seqs_genome(dic_protein_in_genomes, set_of_proteins_in_clusters, dic_protein_length, 50)

    #Check that everything looks ok
    check_number = 0
    for value in selected_unique_proteins.itervalues():
        check_number += len(value)
    for value in removed_unique_proteins.itervalues():
        check_number += len(value)

    if total_number_proteins - len(set_of_proteins_in_clusters) - check_number != 0:
        print "Total number of proteins:" + str(total_number_proteins)
        print "Total number of proteins in clusters:" + str(len(set_of_proteins_in_clusters))
        print "Removed proteins:" + str(check_number)
        sys.exit("Failed checkpoint. The number of unique proteins and proteins in "
                 "clusters does not match the total number of proteins")

    #Print the output files
    count_unique_proteins = open(args.output_directory + "/count_unique_sequences.txt", 'w')
    list_unique_proteins = open(args.output_directory + "/list_unique_sequences.txt", 'w')

    count_unique_proteins.write("Genome\tSelected\tTooShort\n")

    for genome in selected_unique_proteins:
        count_unique_proteins.write(genome + "\t" + str(len(selected_unique_proteins[genome])) + "\t" +
                                    str(len(removed_unique_proteins[genome])) + "\n")

        for protein in selected_unique_proteins[genome]:
            list_unique_proteins.write(genome + "\t" + protein.split("|")[1] + "\n")

    count_unique_proteins.close()
    list_unique_proteins.close()

    ############################
    ##Get the clusters shared between genomes and unique clusters to each genome
    matrix_output = open(args.output_directory + "/matrix_output.txt", 'w')
    list_unique_clusters = open(args.output_directory + "/list_unique_clusters.txt", 'w')
    count_unique_clusters = open(args.output_directory + "/count_unique_clusters.txt", 'w')
    list_shared_single_copy_clusters = open(args.output_directory + "/list_single_copy_clusters.txt", 'w')
    list_shared_multiple_copy_clusters = open(args.output_directory + "/list_shared_multiple_copy.txt", 'w')

    unique_clusters, shared_single_clusters, shared_multiple_clusters, all_clusters_matrix = \
        seqs_shared_clusters(cluster_information, genome_id_dictionary)

    #Print counters
    run_summary.write("Number of shared single copy clusters: %d" % len(shared_single_clusters) + "\n")
    run_summary.write("Number of shared multiple copy clusters: %d" % len(shared_multiple_clusters) + "\n")

    #Print the outputs
    matrix_output.write("\n".join(["\t".join(map(str, r)) for r in all_clusters_matrix]))  # Matrix output

    # Unique clusters per genome (duplicate or paralogs?)
    count_unique_clusters.write("Genome\tNumber of Clusters\n")
    for genome in unique_clusters:
        count_unique_clusters.write(genome + "\t" + str(len(unique_clusters[genome])) + "\n")

        for cluster in unique_clusters[genome]:
            list_unique_clusters.write(genome + "\t" + cluster + "\t"
                                       + ",".join(protein for protein in cluster_information[cluster]) + "\n")

    # Single copy shared clusters

    for cluster in shared_single_clusters:
        list_shared_single_copy_clusters.write(cluster + "\t" + ",".join(cluster_information[cluster]) + "\n")

    # Multiple copy shared clusters

    for cluster in shared_multiple_clusters:
        list_shared_multiple_copy_clusters.write(cluster + "\t" + ",".join(cluster_information[cluster]) + "\n")

    matrix_output.close()
    list_unique_clusters.close()
    count_unique_clusters.close()
    list_shared_single_copy_clusters.close()
    list_shared_multiple_copy_clusters.close()

    ###Save the cluster information
    list_all_clusters = open(args.output_directory + "/list_all_clusters.txt", 'w')
    for cluster in cluster_information:
        list_all_clusters.write(cluster + "\t" + ",".join(cluster_information[cluster]) + "\n")

    list_all_clusters.close()

    ###############
    ##Get clusters shared by groups
    if args.group_information:

        #print genome_groups
        unique_group_clusters, combination_clusters = clusters_in_groups(cluster_information, genome_groups)

        list_unique_clusters_group = open(args.output_directory + "/list_unique_clusters_group.txt", 'w')
        list_all_group_combinations = open(args.output_directory + "/list_all_group_combinations.txt", 'w')
        count_group_results = open(args.output_directory + "/count_groups.txt", 'w')
        protein_count_group_results = open(args.output_directory + "/count_proteins_group.txt", 'w')

        for group in unique_group_clusters:
            protein_count = sum(len(cluster_information[cluster]) for cluster in unique_group_clusters[group])

            count_group_results.write(group + "\t" +
                                      str(len(unique_group_clusters[group])) + "\t" + str(protein_count) + "\n")

            for cluster in unique_group_clusters[group]:
                list_unique_clusters_group.write(group + "\t" + cluster + "\t" + ",".join(cluster_information[cluster])
                                                 + "\n")

        count_group_results.write("\n")

        for combination in combination_clusters:
            combination_name = "-".join(combination)
            protein_count = sum(len(cluster_information[cluster]) for cluster in combination_clusters[combination])

            count_group_results.write(combination_name + "\t" +
                                      str(len(combination_clusters[combination])) + "\t" + str(protein_count) + "\n")

            for cluster in combination_clusters[combination]:
                list_all_group_combinations.write(combination_name + "\t"
                                                  + cluster + "\t" + ",".join(cluster_information[cluster]) + "\n")

            protein_count_group_results.write(combination_name)

            for group in combination:

                count = 0

                genomes = genome_groups[group]

                for cluster in combination_clusters[combination]:
                    for genome in genomes:
                        for protein in cluster_information[cluster]:
                            if protein.startswith(genome):
                                count += 1

                protein_count_group_results.write("\t" + group + ":" + str(count))

            protein_count_group_results.write("\n")

        count_group_results.close()
        list_unique_clusters_group.close()
        list_all_group_combinations.close()
        run_summary.close()