__author__ = 'jlu96'

import sys
import numpy as np
import collections
import pandas as pd
import geneTSmunging as gtm
import pickle
import network_helpers as nh
import scipy.stats as stats
import matplotlib.pyplot as plt


def get_num_above(betas, threshold):
    return len(np.where(betas >= threshold)[0])

def FDR_above_threshold(orig, null, FDR):
    pos_values = np.sort(orig[np.where(orig > 0)])

    for pos_value in pos_values:
        origP = get_num_above(orig, pos_value)
        nullP = get_num_above(null, pos_value)
        fdr = nullP * 1.0 / (nullP + origP)
        if fdr < FDR:
            return pos_value

    return None


def get_num_below(betas, threshold):
    return len(np.where(betas <= threshold)[0])

def FDR_below_threshold(orig, null, FDR):
    neg_values = (orig[np.where(orig < 0)])
    neg_values.sort()
    neg_values = neg_values[::-1]

    for neg_value in neg_values:
        origP = get_num_below(orig, neg_value)
        nullP = get_num_below(null, neg_value)
        fdr = nullP * 1.0 / (nullP + origP)
        if fdr < FDR:
            return neg_value

    return None





def get_thresh(beta_matr, rand_beta_matr, fdr, stratify_by="effect"):
    """
    :param beta_matr: a cause x effect matrix
    :param rand_beta_matr: a cause x effect matrix where causes were randomized by time
    :param fdr: the false discovery rate, treating the causes as randomized by time
    :param stratify_by: col: control the FDR by stratifying by this
    :return:
    """
    print "Calculating thresholds"
    print "Stratifying by ", stratify_by
    thresh_matr = beta_matr.copy()
    beta_threshes = []

    if stratify_by not in {"effect", "none"}:
        raise ValueError("Need to stratify thresholding by effect or over none")

    if stratify_by == "effect":
        for j in range(beta_matr.shape[1]):
            beta_vec = beta_matr[:, j]
            rand_beta_vec = rand_beta_matr[:, j]

            beta_thresh = FDR_above_threshold(beta_vec, rand_beta_vec, fdr)

            beta_threshes.append(beta_thresh)

            if beta_thresh == None:
                thresh_matr[:, j] = np.zeros(thresh_matr.shape[0])
            else:
                c = thresh_matr[:, j]
                c[np.where(c < beta_thresh)] = 0
        return thresh_matr, beta_threshes

    elif stratify_by == "none":
        beta_thresh = FDR_above_threshold(beta_matr.flatten(), rand_beta_matr.flatten(), fdr)
        if beta_thresh == None:
            thresh_matr = np.zeros(thresh_matr.shape)
        else:
            thresh_matr[np.where(thresh_matr < beta_thresh)] = 0
        beta_threshes.append(beta_thresh)
        return thresh_matr, beta_threshes



def get_abs_thresh(beta_matr, rand_beta_matr, fdr, stratify_by="effect"):

    copy_beta_matr = beta_matr.copy()

    abs_beta_matr = np.absolute(beta_matr)
    abs_rand_beta_matr = np.absolute(rand_beta_matr)

    thresh_matr, beta_threshes = get_thresh(abs_beta_matr, abs_rand_beta_matr, fdr, stratify_by = stratify_by)

    copy_beta_matr[np.where(thresh_matr == 0)] = 0

    return copy_beta_matr, beta_threshes


def get_pos_thresh(beta_matr, rand_beta_matr, fdr, stratify_by="effect"):

    pos_beta_matr = beta_matr.copy()
    pos_beta_matr[np.where(pos_beta_matr < 0)] = 0
    pos_rand_beta_matr = rand_beta_matr.copy()
    pos_rand_beta_matr[np.where(pos_rand_beta_matr < 0)] = 0


    thresh_matr, beta_threshes = get_thresh(pos_beta_matr, pos_rand_beta_matr, fdr, stratify_by = stratify_by)

    pos_beta_matr[np.where(thresh_matr == 0)] = 0

    return pos_beta_matr, beta_threshes

def get_neg_thresh(beta_matr, rand_beta_matr, fdr, stratify_by="effect"):

    neg_beta_matr = beta_matr.copy()
    neg_beta_matr[np.where(neg_beta_matr > 0)] = 0
    neg_rand_beta_matr = rand_beta_matr.copy()
    neg_rand_beta_matr[np.where(neg_rand_beta_matr > 0)] = 0


    thresh_matr, beta_threshes = get_thresh(-1 * neg_beta_matr, -1 * neg_rand_beta_matr, fdr, stratify_by = stratify_by)

    neg_beta_matr[np.where(thresh_matr == 0)] = 0

    return neg_beta_matr, beta_threshes


def get_pos_neg_thresh(beta_matr, rand_beta_matr, fdr, stratify_by="effect"):
    pos_thresh_matr, pos_beta_threshes = get_pos_thresh(beta_matr, rand_beta_matr, fdr, stratify_by=stratify_by)
    neg_thresh_matr, neg_beta_threshes = get_neg_thresh(beta_matr, rand_beta_matr, fdr, stratify_by=stratify_by)

    thresh_matr = pos_thresh_matr + neg_thresh_matr

    beta_threshes = zip(pos_beta_threshes, neg_beta_threshes)


    return thresh_matr, beta_threshes

def cap_matr(matr, cap):
    print "Before cap: Num significant betas ", len(np.where(matr != 0)[0])

    matr[np.where(np.absolute(matr) > cap)] = 0
    print "After cap: Num significant betas ", len(np.where(matr != 0)[0])

    return matr


def write_readme(matr, name, fdr, readme_name, filename):
    rd = collections.OrderedDict()

    rd["Name"] = name
    rd["Filename"] = filename
    rd["Genes"] = len(matr)
    rd["Pairs"] = matr.shape[0] * matr.shape[1]
    rd["FDR"] = fdr
    rd["Significant"] = len(np.where(matr)[0])
    rd["% Significant"] = len(np.where(matr)[0]) * 100.0 / (matr.shape[0] * matr.shape[1])

    rd_df = pd.DataFrame(rd, index=["Value"])

    rd_df.transpose().to_csv(readme_name, sep="\t")

    print "Readme for", name, "written to", readme_name

def plot_betas(unshuffled, shuffled, filename=None, zoom_in_percentile=None, xlabel="Beta", ylabel="Count", title="Histogram of Causal Coefficients", nbins=30):


    both = np.concatenate((unshuffled, shuffled))
    bins = np.linspace(min(both), max(both), nbins)

    fig = plt.figure(figsize=(12,8))
    plt.hist(unshuffled, alpha=0.5, color='red', label="Unshuffled", bins=bins)
    plt.hist(shuffled, alpha=0.5, color='blue', label="Shuffled", bins=bins)
    if zoom_in_percentile != None:
        extreme_perc = 0.5 * (100 - zoom_in_percentile)
        top_percentile = stats.scoreatpercentile(both, 100 - extreme_perc)
        bottom_percentile = stats.scoreatpercentile(both, extreme_perc)
        plt.xlim(bottom_percentile, top_percentile)

    plt.legend()
    plt.xlabel(xlabel, fontsize=20)
    plt.ylabel(ylabel, fontsize=20)
    plt.title(title, fontsize=20)

    if filename != None:
        print "Overlaid betas saved to ", filename
        fig.savefig(filename)
    plt.show()
    plt.close()



def get_parser():
    # Parse arguments
    import argparse

    description = 'Apply a pre-specified causal test to an input dataset where each row is a geene' \
                  'and its tim points, specifying which rows to test as effect,'\
                    'Save the results (and parameters if needed), write output coefficients to a pickle file.'
    parser = argparse.ArgumentParser(description=description)


    parser.add_argument('-d', '--original_data', required=True)

    parser.add_argument('-rd', '--randomized_data', required=True)

    parser.add_argument('-m', '--original_matrix', required=True)

    parser.add_argument('-rm', '--randomized_matrix', required=True)

    parser.add_argument('-n', '--name', required=True)

    parser.add_argument('-c', '--coef_num', type=int, required=True)

    parser.add_argument('-f', '--fdr', type=float, required=True)

    parser.add_argument('-sb', '--stratify_by', type=str, required=True)

    parser.add_argument('-mn', '--make_network', type=int, required=True)

    parser.add_argument("-pp", '--plot_prefix', type=str, default=None)

    return parser



def run(args):
    data = gtm.load_file_and_avg(args.original_data)
    rand_data = gtm.load_file_and_avg(args.randomized_data)

    matr = pickle.load(open(args.original_matrix, 'rB'))[:, :, args.coef_num - 1]
    rand_matr = pickle.load(open(args.randomized_matrix, 'rB'))[:, :, args.coef_num - 1]

    if args.stratify_by not in {"e", "n"}:
        raise ValueError("Stratify_by must be either 'e' for effect or 'n' for none")
    else:
        if args.stratify_by == "e":
            stratify_by = "effect"
        elif args.stratify_by == "n":
            stratify_by = "none"

    genes = data["gene"]
    rand_genes = rand_data["gene"]

    if (genes != rand_genes).any():
        raise ValueError("Genes are not the same!")


    print "Original matrix for ", args.name, "saved to", args.name + "-unshuffled-matrix.txt"
    gtm.save_gene_matrix(matrix=matr, filename=args.name + "-unshuffled-matrix.txt", genes=genes)

    print "Randomized matrix for ", args.name, "saved to", args.name + "-shuffled-matrix.txt"
    gtm.save_gene_matrix(matrix=rand_matr, filename=args.name + "-shuffled-matrix.txt", genes=rand_genes)


    if args.plot_prefix != None:
        plot_betas(matr.flatten(), rand_matr.flatten(), filename=args.plot_prefix)
        plot_betas(matr.flatten(), rand_matr.flatten(), filename=args.plot_prefix + "_zoom-in-95", zoom_in_percentile=95)



    print "Using original"
    print "Trying to have an FDR of ", args.fdr
    print args.name


    functions = [get_abs_thresh, get_pos_thresh, get_neg_thresh, get_pos_neg_thresh]
    types = ["abs-thresh", "pos-thresh", "neg-thresh", "pos-neg-thresh"]
    # whether to take absolute value of given matrices
    absoluted = [True, False, False, True]

    for function, t, a in zip(functions, types, absoluted):
        out_prefix = args.name + "-unshuffled-" + t + "-FDR-" + str(args.fdr) + "-stratby-" + stratify_by


        thresh_matr, threshes = function(matr, rand_matr, args.fdr, stratify_by = stratify_by)


        matr_df = gtm.save_gene_matrix(out_prefix + "-matrix.txt", thresh_matr, genes)
        pickle.dump(threshes, open(out_prefix + "-threshes.p", 'w'))

        print "Matrix written to ", out_prefix + "-matrix.txt"
        print "Threshes written to ", out_prefix + "-threshes.p"

        write_readme(thresh_matr, out_prefix, args.fdr, out_prefix + '-README.txt', out_prefix + "-matrix")

        if args.make_network:
            net_df = nh.matr_to_net(matr_df, args.name, make_pair=False)

            net_df.to_csv(out_prefix + "-network.txt", sep="\t", index=False)

            print "Network written to ", out_prefix + "-network.txt"

        if absoluted:
            abs_matr = np.absolute(thresh_matr)

            abs_prefix = args.name + "-unshuffled-" + t + "-absoluted-FDR-" + str(args.fdr) + "-stratby-" + stratify_by

            abs_df = gtm.save_gene_matrix(abs_prefix + "-matrix", abs_matr, genes)

            write_readme(abs_matr, abs_prefix, args.fdr, abs_prefix + '-README.txt', abs_prefix + "-matrix")

            if args.make_network:
                abs_net_df = nh.matr_to_net(abs_df, args.name, make_pair=False)

                abs_net_df.to_csv(abs_prefix + "-network.txt", sep="\t", index=False)

                print "Network written to ", abs_prefix + "-network.txt"



def main():
    run(get_parser().parse_args(sys.argv[1:]))

if __name__ == '__main__':
    main()