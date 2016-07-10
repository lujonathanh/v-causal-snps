__author__ = 'jlu96'

import sys
import geneTSmunging as gtm
import itertools
import pickle
import os


def partition_inputs(input, number):
    num_inputs = len(input)
    return [input[num_inputs * i/number:num_inputs * (i+1)/number] for i in range(number)]

def get_parser():
    # Parse arguments
    import argparse

    description = 'Apply a pre-specified causal test to an input dataset where each row is a gene' \
                  'and its time points'
    parser = argparse.ArgumentParser(description=description)


    parser.add_argument('-d', '--data_file', required=True)

    parser.add_argument('-a', '--args_file', required=True)

    parser.add_argument('-t', '--test', required=True)

    parser.add_argument('-plp', '--pairlist_fileprefix', default=None)

    parser.add_argument('-o', '--output_name', required=True)

    parser.add_argument('-n', '--job_num', type=int, default=3)

    return parser

def run(args):

    df = gtm.load_file_and_avg(args.data_file)

    genes = df['gene'].values

    n = len(genes)

    pairs = list(itertools.permutations(range(n), 2))

    pairlists = partition_inputs(pairs, args.job_num)

    if args.pairlist_fileprefix == None:
        args.pairlist_fileprefix = args.output_name + "-pairs-"


    script_filenames = []
    output_filenames = []

    for pairlist, i in zip(pairlists, range(len(pairlists))):
        pairlist_filename = args.pairlist_fileprefix + str(i) + ".p"

        pickle.dump(pairlist, open(pairlist_filename, 'w'))

        print "Pair lists written to", pairlist_filename


        script_filename = args.output_name + "-script-" + str(i) + ".sh"
        script_filenames.append(script_filename)

        output_filename = args.output_name + "-" + str(i) + ".p"
        output_filenames.append(output_filename)
        # prepare the job associated with this

        command_string = "python run_causal.py -d " + args.data_file.split('/')[-1] + " -a " + args.args_file.split('/')[-1] + " -t " + args.test + " -pl " + \
                         pairlist_filename + " -o " + output_filename

        with open(script_filename, 'w') as outputfile:
            outputfile.write("#!/bin/bash\n")
            outputfile.write("module load python/2.7\n")
            outputfile.write("module load python/2.7/scipy-mkl\n")
            outputfile.write("module load python/2.7/numpy-mkl\n")
            outputfile.write("module load anaconda\n")
            outputfile.write(command_string)
            outputfile.write("\n")
        os.chmod(script_filename, 0777)

        print "Script written to ", script_filename

    # submit the jobs soon


    with open("script_list.txt", 'w') as scriptfile:
        for script_filename in script_filenames:
            scriptfile.write(script_filename + "\n")
        print "Script list written to script_list.txt"

    with open("output_list.txt", 'w') as outputfile:
        for output_filename in output_filenames:
            outputfile.write(output_filename + "\n")
        print "Output list written to output_list.txt"

    with open("integrate_outputs.sh", 'w') as ifile:
        integrated_filename = args.output_name + ".p"
        ifile.write("python integrate_outputs.py -i output_list.txt -o " + integrated_filename + "\n")
        print "Integration script written to integrate_outputs.sh"
        os.chmod("integrate_outputs.sh", 0777)

def main():
    run(get_parser().parse_args(sys.argv[1:]))

if __name__ == '__main__':
    main()
