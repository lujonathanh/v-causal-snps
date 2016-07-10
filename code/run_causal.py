__author__ = 'jlu96'

import sys
import CausalTests as ct
import geneTSmunging as gtm
import pandas as pd
import numpy as np
import pickle

def load_kwargs_file(argsfile):

    df = pd.read_csv(argsfile, sep='\t')

    args = df['Argument'].values
    values = df['Value'].values

    newdict = {}
    for arg, value in zip(args, values):
        try:
          newdict[arg] = eval(value)
        except NameError:
            newdict[arg] = value

    return newdict


def get_parser():
    # Parse arguments
    import argparse

    description = 'Apply a pre-specified causal test to an input dataset where each row is a geene' \
                  'and its tim points, write output to a pickle file.'
    parser = argparse.ArgumentParser(description=description)


    parser.add_argument('-d', '--data_file', required=True)

    parser.add_argument('-a', '--args_file', required=True)

    parser.add_argument('-t', '--test', required=True)

    parser.add_argument('-pl', '--pairlist_file', default=None)

    parser.add_argument('-o', '--output_name', required=True)


    return parser



def run(args):

    # load the data

    df = gtm.load_file_and_avg(args.data_file)

    genes = df['gene'].values

    found_genes, geneTS = gtm.get_gene_TS(df, genes)

    args_dict = load_kwargs_file(argsfile=args.args_file)

    if args.pairlist_file == None:
        pairlist = None
    else:
        pairlist = np.load(open(args.pairlist_file))

    print args_dict

    if args.test == 'g':
        output = ct.pairwise_granger_causality_all(geneTS, pairlist, **args_dict)
        with open(args.output_name, 'w') as outfile:
            pickle.dump(output, outfile)

    print "HELLOOOOOOOO"
    print "Output written to ", args.output_name


def main():
    run(get_parser().parse_args(sys.argv[1:]))

if __name__ == '__main__':
    main()
