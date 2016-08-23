__author__ = 'jlu96'


__author__ = 'jlu96'

import sys
import csv
import pickle
import numpy as np
import pandas as pd



def get_parser():
    # Parse arguments
    import argparse

    description = 'Integrate output to a pickle file.'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-i', '--output_dfname', required=True)

    parser.add_argument('-t', '--type', required=True, help="Args.type must be m (matrix) or d (pandas dataframe)")
    # if "m", then pickle
    # if "d", then pandas

    parser.add_argument('-a', '--axis', type=int, default = -1)

    parser.add_argument('-o', '--int_name_dfname', required=True)
    return parser




def run(args):

    # load the data

    output_df = pd.read_csv(args.output_dfname, sep="\t")
    int_name_df = pd.read_csv(args.int_name_dfname, sep="\t")
    print int_name_df.head()

    assert set(output_df.columns.values) == set(int_name_df.columns.values)

    for x in output_df:
        filenames = output_df[x].values

        print
        print x
        print "Files to integrate are ", filenames
        integrated_filename = int_name_df[x].values[0]
        print "Integrated_file is ", integrated_filename


        if args.type == "m":

            if args.axis == -1:
                raise ValueError("You must set the axis of concatenation for the matrices.")

            all_matrs = []
            for filename in filenames:
                all_matrs.append(pickle.load(open(filename, 'rU')))

            integrated_matr = np.concatenate(tuple(all_matrs), axis=args.axis)

            print "Final matrix shape", integrated_matr.shape

            with open(integrated_filename, 'w') as outfile:
                pickle.dump(integrated_matr, outfile)

            print "Integrated matr written to ", integrated_filename


        elif args.type == "d":
            all_dfs = []
            for filename in filenames:
                all_dfs.append(pd.read_csv(filename, sep="\t"))

            integrated_df = pd.concat(all_dfs)

            print "Final df shape:", integrated_df.shape

            integrated_df.to_csv(integrated_filename, sep="\t", index=False)

            print "Integrated df written to ", integrated_filename

        else:
            raise ValueError("Args.type must be m (matrix) or d (pandas dataframe)")



def main():
    run(get_parser().parse_args(sys.argv[1:]))

if __name__ == '__main__':
    main()



# __author__ = 'jlu96'
#
#
# __author__ = 'jlu96'
#
# import sys
# import csv
# import pickle
# import numpy as np
#
#
#
# def get_parser():
#     # Parse arguments
#     import argparse
#
#     description = 'Integrate output to a pickle file.'
#     parser = argparse.ArgumentParser(description=description)
#
#     parser.add_argument('-i', '--output_listname', default="output_list.txt")
#
#     parser.add_argument('-o', '--integrated_filename', required=True)
#
#     parser.add_argument('-p_pos', type=int, default=0)
#
#     parser.add_argument('-p_add', type=bool, default=True)
#
#     return parser
#
#
#
#
# def run(args):
#
#     # load the data
#
#     filenames = []
#     with open(args.output_listname, 'rU') as csvfile:
#         reader = csv.reader(csvfile, delimiter='\t')
#         for row in reader:
#             filenames.append(row[0])
#
#     print "Files to integrate are ", filenames
#
#     output_tuples = []
#
#     for filename in filenames:
#         output_tuple = pickle.load(open(filename, 'rU'))
#         output_tuples.append(output_tuple)
#
#     output_lists = zip(*output_tuples)
#
#     integrated_tuple = tuple(np.sum(output_list, axis=0) for output_list in output_lists)
#
#     if args.p_add:
#         print "Adding ones to p-matr at location ", args.p_pos
#         p_matr = integrated_tuple[args.p_pos]
#         p_matr += np.diagflat(np.ones(p_matr.shape[0]))
#         integrated_tuple = integrated_tuple[0:args.p_pos] + (p_matr,) + integrated_tuple[args.p_pos + 1:]
#
#     with open(args.integrated_filename, 'w') as outfile:
#
#         pickle.dump(integrated_tuple, outfile)
#
#     print "Integrated written to ", args.integrated_filename
#
# def main():
#     run(get_parser().parse_args(sys.argv[1:]))
#
# if __name__ == '__main__':
#     main()
