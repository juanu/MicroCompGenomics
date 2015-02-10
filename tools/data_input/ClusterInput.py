from collections import defaultdict
import re

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
        prots = ids_proteins[1].split(" ")

        prots = [x.split("(")[0] for x in prots]  # Remove the genome names on the protein IDS (Fastortho does that)

        clean_protein_list = []  # This is used to remove proteins from genomes not in the list

        for genome_entry in gl:  # Adding the proteins that we need
            [clean_protein_list.append(prot) for prot in [x for x in prots if x.startswith(genome_entry)]]

        #Now I need to evaluate those clusters that now are unique and zero
        if len(clean_protein_list) == 0:
            clusters_removed += 1
            continue

        if len(clean_protein_list) == 1:
            unique_proteins_genome_count[clean_protein_list[0].split("|")[0]] += 1
            clusters_removed += 1
            continue

        for prot in clean_protein_list:
            cluster_id = ids_proteins[0].split(" (")[0]

            cluster_dictionary[cluster_id].append(prot)
            proteins_in_cluster.add(prot)

    return cluster_dictionary, proteins_in_cluster, unique_proteins_genome_count, total_cluster_count, clusters_removed


def parse_tribemcl(cf, gl):
    """
    Function used to parse TribeMCL results
    :param cf:
    :param gl:
    :return:
    """

    tribemcl_results = open(cf, 'r')

    cluster_dictionary = defaultdict(list)  # Dictionary with the clusters. Stores only the clusters that will be used
    unique_proteins_genome_count = defaultdict(int)  # Proteins that are not in any cluster
    proteins_in_cluster = set()  # Proteins in cluster. Needed to search for proteins that are absent
    total_cluster_count = 0
    clusters_removed = 0
    cluster_number = 1

    for line in tribemcl_results:
        total_cluster_count += 1
        line = line.strip('\n')
        prots = line.split("\t")

        prots = [x.split("(")[0] for x in prots]  # Remove the genome names on the protein IDS (Fastortho does that)

        clean_protein_list = []  # This is used to remove proteins from genomes not in the list

        for genome_entry in gl:  # Adding the proteins that we need
            [clean_protein_list.append(prot) for prot in [x for x in prots if x.startswith(genome_entry)]]

        #Now I need to evaluate those clusters that now are unique and zero
        if len(clean_protein_list) == 0:
            clusters_removed += 1
            continue

        if len(clean_protein_list) == 1:
            unique_proteins_genome_count[clean_protein_list[0].split("|")[0]] += 1
            clusters_removed += 1
            continue

        cluster_id = "Cluster" + str(cluster_number)

        for prot in clean_protein_list:

            cluster_dictionary[cluster_id].append(prot)
            proteins_in_cluster.add(prot)

        cluster_number += 1

    return cluster_dictionary, proteins_in_cluster, unique_proteins_genome_count, total_cluster_count, clusters_removed