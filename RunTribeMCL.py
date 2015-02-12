#Created on 1/31/2015
__author__ = 'Juan A. Ugalde'


import os
import argparse


program_description = "This script will run a TribeMCL analysis, using MCL and iterating over a series of different" \
                      "inflation values (if requested)"

parser = argparse.ArgumentParser(description=program_description)

parser.add_argument("-o", "--output_directory", type=str, required=True)
parser.add_argument("-b", "--input_blast", type=str, required=True, help="Input BLAST file, tabular format")
parser.add_argument("-i", "--iteration", action="store_true", help="If this flag is used, the FastOrtho will be run"
                                                                   "with several inflation values for testing. "
                                                                   "The output will also be a table and a plot showing"
                                                                   "the number of clusters found at different inflation"
                                                                   "values")
parser.add_argument("-t", "--threads", type=int, help="Number of threads to use for MCL")
args = parser.parse_args()

#Create the output directory
if not os.path.exists(args.output_directory):
    os.makedirs(args.output_directory)


#Create an option file. If the iteration flag is not set, we will use the default value of 1.5
#After creation, run FastOrtho

#Move to the output directory
#os.chdir(args.output_directory)

threads = 1

if args.threads:
    threads = args.threads

    #Run FastOrtho
    #Run TribeMCL
abc_file = args.output_directory + "/seq.abc"
tab_file = args.output_directory + "/seq.tab"
mci_file = args.output_directory + "/seq.mci"

make_abc_file = "cut -f 1,2,11 " + args.input_blast + " > " + abc_file
run_mcxload = "mcxload -abc " + abc_file + " --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o " \
              + mci_file + " -write-tab " + tab_file

os.system(make_abc_file)
os.system(run_mcxload)

if not args.iteration:

    output_file = args.output_directory + "/out.seq.mci.I.1.4"
    run_mcl = "mcl " + mci_file + " -I 1.4 -te " + str(threads) + " -o " + output_file + " -use-tab " + tab_file

    os.system(run_mcl)

    #Create an option file with iteration. The range will go from 1 to 10, in 0.5 intervals. This could
    #change, depending on some results.

else:
    i = 1.1
    cluster_inflation = []

    while i < 5.1:

        output_file = args.output_directory + "/out.seq.mci.I" + str(i)
        run_mcl = "mcl " + mci_file + " -I " + str(i) + " -te " + str(threads) + \
                  " -o " + output_file + " -use-tab " + tab_file

        os.system(run_mcl)

        #Count number of clusters
        number_cluster = 0
        with open(output_file) as f:
            number_cluster = sum(1 for _ in f)

        cluster_inflation.append([i, number_cluster])

        i += 0.1

    cluster_data = open(args.output_directory + "/Inflation_ClusterSize.txt", 'w')
    for value in cluster_inflation:
        cluster_data.write(str(value[0]) + "\t" + str(value[1]) + "\n")

    cluster_data.close()
