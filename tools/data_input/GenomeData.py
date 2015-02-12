#Created on 1/27/2015
__author__ = 'Juan A. Ugalde'

import sys
import os
from collections import defaultdict
from Bio import SeqIO


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
    This script modifies a fasta file to transform it into a standard format that can be used in the rest of the
    pipeline, and is compatible with OrthoMCL.

    :param input_fasta_file: Fasta file to modify
    :param new_fasta_prefix: Prefix for the fasta file, this will be used to replace the IDs on each entry
    :param output_file: A new fasta file, suitable to use for OrthoMCL
    :return:
    """

    #Rename the fasta
    for record in SeqIO.parse(input_fasta_file, "fasta"):
        output_file.write(">" + new_fasta_prefix + "|" + record.id + "\n")
        output_file.write(str(record.seq) + "\n")

    ##I should add something to check that the fasta is in the right format.##


def parse_jgi_dump(genome_dictionary, input_path, output_path):
    """
    :param genome_dictionary:
    :param input_path:
    :param output_path:
    :return:
    """

    genome_count = 0

    for prefix in genome_dictionary:
        img_id = genome_dictionary[prefix]

        organism_folder = input_path + "/" + img_id
        file_prefix = organism_folder + "/" + img_id

        if not os.path.exists(organism_folder):
            continue

        nucleotide_file = open(output_path + "/nucleotide/" + prefix + ".fna", 'w')
        aminoacid_file = open(output_path + "/protein/" + prefix + ".fasta", 'w')
        genome_file = open(output_path + "/genome/" + prefix + ".fna", 'w')
        annotation_file = open(output_path + "/annotation/" + prefix + ".txt", 'w')
        coords_file = open(output_path + "/coords/" + prefix + ".txt", 'w')

        ##Generate the output files with the sequence data

        #First the input files that are present in the JGI folder
        input_nucleotide = file_prefix + ".genes.fna"
        input_aa = file_prefix + ".genes.faa"
        input_genome = file_prefix + ".fna"

        #Modify the files and output to the right folder
        modify_fasta(input_nucleotide, prefix, nucleotide_file)
        modify_fasta(input_aa, prefix, aminoacid_file)
        modify_fasta(input_genome, prefix, genome_file)

        ##Now to generate the annotation outputs
        #Generate the output of the annotation
        input_cog = file_prefix + ".cog.tab.txt"
        input_pfam = file_prefix + ".pfam.tab.txt"
        input_ko = file_prefix + ".ko.tab.txt"
        input_gff = file_prefix + ".gff"

        #Store the COG, PFAM, KO, and product annotation
        #Format Gene_ID Description(COG, PFAM, KO) Annotation
        annotation_summary = defaultdict(lambda: defaultdict(str))
        annotation_defs = defaultdict()

        for line in open(input_cog, 'r'):
            if line.strip():
                if line.startswith("gene_oid"):
                    continue
                else:
                    line = line.rstrip()
                    elements = line.split("\t")
                    annotation_summary[elements[0]]["COG"] = elements[9]
                    annotation_defs[elements[9]] = elements[10]

        if os.path.exists(input_pfam):

            for line in open(input_pfam, 'r'):
                if line.strip():
                    if line.startswith("gene_oid"):
                        continue
                    else:
                        line = line.rstrip()
                        elements = line.split("\t")
                        annotation_summary[elements[0]]["PFAM"] = elements[8]
                        annotation_defs[elements[8]] = elements[9]

        for line in open(input_ko, 'r'):
            if line.strip():
                if line.startswith("gene_oid"):
                    continue
                else:
                    line = line.rstrip()
                    elements = line.split("\t")
                    annotation_summary[elements[0]]["KO"] = elements[9]

                    try:
                        annotation_defs[elements[9]] = elements[10]
                    except IndexError:
                        continue

        for line in open(input_gff, 'r'):
            if line.strip():
                if line.startswith("#"):
                    continue
                else:
                    line = line.rstrip()
                    elements = line.split("\t")

                    feature_type = elements[2]

                    #annotation_summary[protein_id]["feature_type"] = feature_type

                    if not feature_type == "CDS":  # only look at CDS
                        continue

                    try:
                        features = elements[8]
                    except IndexError:
                        continue

                    protein_id = None
                    product_desc = None
                    locus_tag = None

                    for feature in features.split(";"):
                        if feature.startswith("ID"):
                            name, protein_id = feature.split("=")
                        if feature.startswith("product"):
                            info = feature.split("=")
                            product_desc = info[1]
                        if feature.startswith("locus_tag"):
                            ltag, locus_tag = feature.split("=")

                    if protein_id is not None:
                        annotation_summary[protein_id]["Feature_type"] = elements[2]
                        annotation_summary[protein_id]["start"] = elements[3]
                        annotation_summary[protein_id]["end"] = elements[4]
                        annotation_summary[protein_id]["contig"] = elements[0]
                        annotation_summary[protein_id]["locus_tag"] = locus_tag

                    if product_desc is not None:
                        annotation_summary[protein_id]["Product"] = product_desc

        for protein in annotation_summary:
            for annotation_type in annotation_summary[protein]:

                function = annotation_summary[protein][annotation_type]

                if annotation_type == "Product" or annotation_type == "Feature_type":
                    output_line = protein + "\t" + annotation_type + "\t" + function
                    annotation_file.write(output_line + "\n")

                else:
                    if annotation_defs.has_key(function):
                        output_line = protein + "\t" + annotation_type + "\t" + function + "\t" + annotation_defs[function]
                    else:
                        output_line = protein + "\t" + annotation_type + "\t" + function
                    annotation_file.write(output_line + "\n")

            coords_file.write(annotation_summary[protein]["contig"] + "\t" + protein + "\t" +
                              annotation_summary[protein]["start"] + "\t" + annotation_summary[protein]["end"] + "\n")

        nucleotide_file.close()
        aminoacid_file.close()
        genome_file.close()
        annotation_file.close()
        coords_file.close()

        genome_count += 1

    return genome_count

def parse_single_img(genome_dictionary, input_path, output_path):
    """
    :param genome_dictionary:
    :param input_path: Input path with the genome information. Two files need to be present here, the genome annotation
    (usually with the extension .info.xls, and a GBK file, usually with the .gbf extension
    :param output_path: Output for the processed genomes
    :return: A folder with genomes ready to be analyzed
    """

    genome_count = 0

    for prefix in genome_dictionary:
        img_id = genome_dictionary[prefix]
        print img_id, prefix

        organism_folder = input_path + "/" + img_id
        file_prefix = organism_folder + "/" + img_id

        if not os.path.exists(organism_folder):
            continue

        nucleotide_file = open(output_path + "/nucleotide/" + prefix + ".fna", 'w')
        aminoacid_file = open(output_path + "/protein/" + prefix + ".fasta", 'w')
        genome_file = open(output_path + "/genome/" + prefix + ".fna", 'w')
        annotation_file = open(output_path + "/annotation/" + prefix + ".txt", 'w')
        coords_file = open(output_path + "/coords/" + prefix + ".txt", 'w')

        ##Generate the output files with the sequence data

        #First the input files that are present in the IMG folder
        input_gbk = file_prefix + ".gbf"
        input_annotation = file_prefix + ".info.xls"

        locustag_to_id = defaultdict()

        annotation_dictionary = defaultdict(lambda: defaultdict())

        #I have to process the annotation first, to go from the locus_tag to the IMG gene number
        for line in open(input_annotation, 'r'):
            if line.strip():
                if line.startswith("gene_oid"):
                    continue
                else:
                    line = line.rstrip()
                    elements = line.split("\t")
                    locustag_to_id[elements[1]] = elements[0]

                    annotation_dictionary[elements[0]][elements[2]] = elements[3:]

        #Process the annotation dictionary
        for gene_entry in annotation_dictionary:

            #Only process CDS
            if not annotation_dictionary[gene_entry]["Locus_type"][1] == "CDS":
                continue

            #Get all the PFAM entries

            pfam_entries = [k for k in annotation_dictionary[gene_entry].keys() if k.startswith("pfam")]

            for pfam in pfam_entries:
                pfam_description = annotation_dictionary[gene_entry][pfam][0]
                annotation_file.write("%s\tPFAM\t%s\t%s\n" % (gene_entry, pfam, pfam_description))


            #Get all the COG entries (maybe there is more than one)

            cog_entries = [k for k in annotation_dictionary[gene_entry].keys() if k.startswith("COG") and k is not "COG_category"]

            for cog in cog_entries:
                cog_description = annotation_dictionary[gene_entry][cog][0]
                annotation_file.write("%s\tCOG\t%s\t%s\n" % (gene_entry, cog, cog_description))

            #Get all the KO entries (maybe there is more than one)

            ko_entries = [k for k in annotation_dictionary[gene_entry].keys() if k.startswith("KO")]

            for ko in ko_entries:
                ko_description = annotation_dictionary[gene_entry][ko][0]
                annotation_file.write("%s\tKO\t%s\t%s\n" % (gene_entry, ko, ko_description))

            #Print locus type, for consistency with other data. Probably not necessary

            annotation_file.write("%s\tFeature_type\t%s\n" % (gene_entry, annotation_dictionary[gene_entry]["Locus_type"][1]))

            #Print product name
            try:
                product_name = annotation_dictionary[gene_entry]["Product_name"][1]

            except IndexError:
                product_name = None

            annotation_file.write("%s\tProduct\t%s\n" % (gene_entry, product_name))

            #Print the coordinate data

            scaf_id = annotation_dictionary[gene_entry]["Scaffold"][1]
            start, stop = annotation_dictionary[gene_entry]["Coordinates"][1].split("..")
            stop = stop[:-3]
            coords_file.write("%s\t%s\t%s\t%s\n" % (scaf_id, gene_entry, start, stop))

        #Process the GBK file to create the aminoacid, nucleotide and genome files

        #First, I need to check that the gbk file is on the right format. There is issue where the LOCUS line
        #on the gbk files generated by JGI is compressing the size information. For example:
        #LOCUS       JCM14758DRAFT_BALL01000001_1.1233475 bp   DNA linear
        #should be
        #LOCUS       JCM14758DRAFT_BALL01000001_1.1  233475 bp   DNA linear


        for record in SeqIO.parse(input_gbk, "genbank"):
            scaf_id = prefix + "|" + record.id
            #Write the genome sequence
            genome_file.write(">" + scaf_id + "\n" + str(record.seq) + "\n")

            locus_id = None

            for feature in record.features:

                if feature.type == "gene":
                    locus_id = feature.qualifiers["locus_tag"][0]

                elif feature.type == "CDS":

                    protein_id = locustag_to_id[locus_id]
                    cds_id = prefix + "|" + protein_id
                    nucleotide_sequence = feature.extract(record.seq)
                    aminoacid_sequence = feature.qualifiers["translation"][0]

                    nucleotide_file.write(">" + cds_id + "\n" + str(nucleotide_sequence) + "\n")
                    aminoacid_file.write(">" + cds_id + "\n" + str(aminoacid_sequence) + "\n")

        #input_gbk.close()

        #Close output files
        nucleotide_file.close()
        aminoacid_file.close()
        genome_file.close()
        annotation_file.close()
        coords_file.close()

        genome_count += 1

    return genome_count