#Created on 1/31/2015
__author__ = 'Juan A. Ugalde'

#TODO Do iteration over several inflation values

import os
import argparse


def create_option_file(mcl, od, inflation, blast, pd):
    """

    :param mcl: Path of MCL
    :param od: Output directory
    :param inflation: Inflation value to use
    :param blast: Input blast file
    :param pd: Protein directory
    :return:
    """
    #Open the parameter file, and print some of the data
    filename = od + "/options_I" + str(inflation) + ".txt"
    parameter_file = open(filename, 'w')

    working_directory = od + "/fastortho_I" + str(inflation)

    parameter_file.write("--mcl_path %s\n" % mcl)
    parameter_file.write("--working_directory %s\n" % working_directory)
    parameter_file.write("--project_name %s\n" % ("fastortho_I" + str(inflation)))
    parameter_file.write("--blast_file %s\n" % blast)
    parameter_file.write("--inflation %s\n" % inflation)

    #Create the working directory
    if not os.path.exists(working_directory):
        os.makedirs(working_directory)

    #Add the protein info to the parameter file
    for root, dirs, files in os.walk(pd):
        for name in files:
            parameter_file.write("--single_genome_fasta %s\n" % (pd + "/" + name))

    parameter_file.close()

    return filename

program_description = "This script will run the FastOrtho analysis. FastOrtho is a reimplementation of OrthoMCL" \
                      "that does not require the use of a MySQL database. More detailes on " \
                      "http://enews.patricbrc.org/fastortho/"

parser = argparse.ArgumentParser(description=program_description)

parser.add_argument("-f", "--fastortho", type=str, required=True, help="Location of the FastOrtho binary")
parser.add_argument("-w", "--output_directory", type=str, required=True)
parser.add_argument("-p", "--protein_directory", type=str, required=True, help="Path with the protein files,"
                                                                                      "one per genome")
parser.add_argument("-m", "--mcl", type=str, required=True, help="Location of the MCL binary")
parser.add_argument("-b", "--input_blast", type=str, required=True, help="Input BLAST file, tabular format")
parser.add_argument("-i", "--iteration", action="store_true", help="If this flag is used, the FastOrtho will be run with"
                                                                   "several inflation values for testing. "
                                                                   "The output will also be a table and a plot showing"
                                                                   "the number of clusters found at different inflation"
                                                                   "values")
args = parser.parse_args()

#Create the output directory
if not os.path.exists(args.output_directory):
    os.makedirs(args.output_directory)


#Create an option file. If the iteration flag is not set, we will use the default value of 1.5
#After creation, run FastOrtho

#Move to the output directory
#os.chdir(args.output_directory)


if not args.iteration:

    path_option_file = create_option_file(args.mcl, args.output_directory, 1.1, args.input_blast, args.protein_directory)

    #Run FastOrtho
    run_fast_ortho = args.fastortho + " --option_file " + path_option_file
    print run_fast_ortho
    os.system(run_fast_ortho)

    #Create an option file with iteration. The range will go from 1 to 10, in 0.5 intervals. This could
    #change, depending on some results.

else:
    i = 1.1
    cluster_inflation = []

    while i < 5.1:
        path_option_file = create_option_file(args.mcl, args.output_directory, i, args.input_blast, args.protein_directory)
        run_fast_ortho = args.fastortho + " --option_file " + path_option_file
        os.system(run_fast_ortho)

        #Count number of clusters
        cluster_file = args.output_directory + "/fastortho_I" + str(i) + "/fastortho_I" + str(i) + ".end"
        number_cluster = 0
        with open(cluster_file) as f:
            number_cluster = sum(1 for _ in f)

        cluster_inflation.append([i, number_cluster])

        i += 0.2

    cluster_data = open(args.output_directory + "/Inflation_ClusterSize.txt", 'w')
    for value in cluster_inflation:
        cluster_data.write(str(value[0]) + "\t" + str(value[1]) + "\n")

    cluster_data.close()
