__author__ = 'jlu96'

import prep_jobs as pj
import sys
import os
import geneTSmunging as gtm
import pandas as pd
import math

def get_parser():
    # Parse arguments
    import argparse

    description = 'Prepare cluster jobs by partitioning tests by rows.'
    parser = argparse.ArgumentParser(description=description)


    parser.add_argument('-d', '--data_file', required=True)

    parser.add_argument('-d2', '--rand_data_file', required=True, help="The effect genes")

    parser.add_argument('-a', '--args_file', required=True)

    parser.add_argument('-t', '--test', required=True)

    parser.add_argument('-o', '--output_name', required=True)

    parser.add_argument('-n', '--job_num', type=int, default=3)

    parser.add_argument('-p', '--parallel_num', type=int, default=0)

    parser.add_argument('-f', '--fdr', type=float, default=0.05)

    parser.add_argument('-c', '--coef_num', type=int, required=True)

    return parser

def run(args):


    data_file = args.data_file.split('/')[-1]
    rand_data_file = args.rand_data_file.split('/')[-1]


    df = gtm.load_file_and_avg(data_file)

    genes = df['gene'].values

    n = len(genes)

    script_filenames = []
    output_filenames = []
    output_rand_filenames = []

    if args.test == "e":
        all_res_filenames = []
        use_filenames = []
        all_res_rand_filenames = []
        use_rand_filenames = []
    else:
        all_res_filenames = None
        use_filenames = None
        all_res_rand_filenames = None
        use_rand_filenames = None

    partition_rows = pj.partition_inputs(range(n), args.job_num)


    for partition_row, i in zip(partition_rows, range(len(partition_rows))):

        script_filename = args.output_name + "-script-" + str(i) + ".sh"
        script_filenames.append(script_filename)


        output_filename = args.output_name + "-" + str(i) + ".p"
        output_filenames.append(output_filename)

        output_rand_filename = args.output_name + "-randomized-" + str(i) + ".p"
        output_rand_filenames.append(output_rand_filename)

        # prepare the job associated with this

        row_filename = args.output_name + "-row-" + str(i) + ".txt"

        command_string = "python run_causal_rand_row.py -d " + data_file +  " -rd " + rand_data_file + \
                         " -a " + args.args_file.split('/')[-1] + " -t " + args.test + " -rl " + \
                         str(row_filename) + " -o " + output_filename + " -or " + output_rand_filename

        if args.test == "e":
            all_res_filename = args.output_name + "-all-params-" + str(i) + ".txt"
            all_res_filenames.append(all_res_filename)

            use_filename = args.output_name + "-used-params-" + str(i) + ".txt"
            use_filenames.append(use_filename)

            all_res_rand_filename = args.output_name + "-all-params-randomized-" + str(i) + ".txt"
            all_res_rand_filenames.append(all_res_rand_filename)

            use_rand_filename = args.output_name + "-used-params-randomized-" + str(i) + ".txt"
            use_rand_filenames.append(use_rand_filename)

            command_string += " -oa " + all_res_filename + " -ou " + use_filename + " -ora " + all_res_rand_filename + " -oru " + use_rand_filename


        with open(row_filename, 'w') as rowfile:
            rowfile.write(str(partition_row) + "\n")

        print "Partition row written to ", row_filename


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

    integrated_name_dict = {}
    integrated_name_dict["Output"] = args.output_name + ".p"
    integrated_name_dict["Rand-Output"] = args.output_name + "-randomized.p"
    integrated_name_dict["All-Params"] = args.output_name + "-all-params.txt"
    integrated_name_dict["Use-Params"] = args.output_name + "-use-params.txt"
    integrated_name_dict["All-Rand-Params"] = args.output_name + "-all-params-randomized.txt"
    integrated_name_dict["Use-Rand-Params"] = args.output_name + "-use-params-randomized.txt"


    with open("script_list.txt", 'w') as scriptfile:
        for script_filename in script_filenames:
            scriptfile.write(script_filename + "\n")
        print "Script list written to script_list.txt"


    # list of matrices to integrate
    output_matr_dict = {"Output": output_filenames, "Rand-Output": output_rand_filenames}
    output_matr_df = pd.DataFrame(output_matr_dict)
    output_matr_df.to_csv("output_matr_list.txt", sep="\t", index=False)
    print "Output matrices written to output_matr_list.txt"

    int_matr_dict = dict([(x, integrated_name_dict[x]) for x in ["Output", "Rand-Output"]])
    int_matr_df = pd.DataFrame(int_matr_dict, index=[0])
    int_matr_df.to_csv("int_matr_list.txt", sep="\t", index=False)
    print "integrated matrices written to int_matr_list.txt"


    if args.test == "e":
        # lists of dataframes (param files) to integrate
        # These will only be integrated if
        output_df_dict = {}
        output_df_lists = [all_res_filenames, use_filenames, all_res_rand_filenames, use_rand_filenames]
        output_df_names = ["All-Params", "Use-Params", "All-Rand-Params", "Use-Rand-Params"]
        for out_list, out_name in zip(output_df_lists, output_df_names):
            if out_list != None:
                output_df_dict[out_name] = out_list

        output_df_df = pd.DataFrame(output_df_dict)
        output_df_df.to_csv("output_df_list.txt", sep="\t", index=False)
        print "output dfs written to output_df_list.txt"


        int_df_dict = dict([(x, integrated_name_dict[x]) for x in set(output_df_names).intersection(output_df_dict.keys())])
        int_df_df = pd.DataFrame(int_df_dict, index=[0])
        int_df_df.to_csv("int_df_list.txt", sep="\t", index=False)
        print "Integrated dfs written to int_df_list.txt"


    with open("integrate_outputs.sh", 'w') as ifile:

        if args.test == "e":
            # here , "a" means the axis to integrate by
            ifile.write("python integrate_outputs_rand_row.py -i output_matr_list.txt -t m -o int_matr_list.txt -a 1 && " + \
                        "python integrate_outputs_rand_row.py -i output_df_list.txt -t d -o int_df_list.txt\n")

        else:
            ifile.write("python integrate_outputs_rand_row.py -i output_matr_list.txt -t m -o int_matr_list.txt -a 1\n")

        print "Integration script written to integrate_outputs.sh"
        os.chmod("integrate_outputs.sh", 0777)

    with open("fdr_control.sh", 'w') as ffile:
        fdr_string = "python fdr_control.py -m " + integrated_name_dict["Output"] + " -rm " + integrated_name_dict["Rand-Output"] + \
                    " -d " + data_file + " -rd " + rand_data_file + " -n " + args.output_name + " -f \"" + str(args.fdr) + "\" " + \
                    " -c " + str(args.coef_num) + " -mn " + str(1) + " -pp " + args.output_name + "-all-beta-histogram "
        ffile.write(fdr_string + " -sb e && " + fdr_string + " -sb n\n")
        print "FDR CONTROL script written to fdr_control.sh"
        os.chmod("fdr_control.sh", 0777)


    if args.parallel_num > 0:
        print "Parallel Number (# processes per job): " + str(args.parallel_num)

        script_groups = pj.partition_inputs(script_filenames, number=int(math.ceil(len(script_filenames) * 1.0/args.parallel_num)))

        print "Number of script groups ", len(script_groups)


        parallel_scripts = []
        for i, script_group in zip(range(len(script_groups)), script_groups):
            appended_script_filenames = ["./" + script_filename for script_filename in script_group]
            parallel_script = " & ".join(appended_script_filenames)
            print "Parallel Script ", i, ":", parallel_script
            parallel_scripts.append(parallel_script)

        with open("parallel_script_list.txt", 'w') as scriptfile:
            for parallel_script in parallel_scripts:
                scriptfile.write(parallel_script + "\n")
            print "Parallel script list written to parallel_script_list.txt"


def main():
    run(get_parser().parse_args(sys.argv[1:]))

if __name__ == '__main__':
    main()
