__author__ = 'jlu96'


import matplotlib.pyplot as plt
import scipy.stats as stats
import pandas as pd
import numpy as np

def load_file_and_avg(filename, keys = ["t00", "t05", "t1_", "t2_", "t3_", "t4_", "t5_", "t6_", "t7_", "t8_", "t10_", "t12_"], fold_change=True,
                      diff=True):
    data = pd.read_csv(filename, sep="\t")
    num_per_key = []
    for key in keys:
        cols = [col for col in list(data.columns.values) if col[:len(key)] == key and col != key]
        num_per_key.append(len(cols))
        print cols
        data[key] = pd.Series(sum([data[col] for col in cols]) * 1.0 / len(cols), index=data.index)
        std = np.std([data[col] for col in cols], axis=0)
        data[key + 'std'] = std


    # normalized
    # look at differences

    if fold_change:
        fold_keys = []
        for i in range(len(keys) - 1):
            key1 = keys[i]
            key2 = keys[i + 1]
            fold_key = key1 + "-" + key2 + " fold"
            data[fold_key] = data[key2] * 1.0/ data[key1]
            fold_keys.append(fold_key)
    if diff:
        diff_keys = []
        for i in range(len(keys) - 1):
            key1 = keys[i]
            key2 = keys[i + 1]
            diff_key = key1 + "-" + key2 + " diff"
            data[diff_key] = data[key2] - data[key1]
            diff_keys.append(diff_key)



    if std_diff:
        std_diff_keys = []
        data["Mean_Diff"] = data[diff_keys].mean(axis=1)
        data["Std_Diff"] = data

        for i in range(len(keys) - 1):
            key1 = keys[i]
            key2 = keys[i + 1]
            std_diff_key = key1 + "-" + key2 + " std_diff"
            data[std_diff_key] = data[key2] - data[key1]
            std_diff_keys.append(std_diff_key)



        

    for key, key_num in zip(keys, num_per_key):
        print key, "has", key_num, "data points"

    return data

def get_gene_TS(data, genes, keys=["t00", "t05", "t1_", "t2_", "t3_", "t4_", "t5_", "t6_", "t7_", "t8_", "t10_", "t12_"]):
    """Return genes found
    and geneTS an n x t matrix of expression levels where n is number of genes and t is number of timepoints, i.e. keys
    """

    gene_set = set(genes)

    key_indices = np.where([value in keys for value in data.columns.values])[0]


    geneTS = data[data.gene.isin(gene_set)][key_indices].values
    found_genes = data[data.gene.isin(gene_set)]['gene'].values

    return found_genes, geneTS


def plot_genes(data, genes, keys = ["t00", "t05", "t1_", "t2_", "t3_", "t4_", "t5_", "t6_", "t7_", "t8_", "t10_", "t12_"], title="Gene Expression From a few species",
               num_per_keys  = np.array([4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4])):
    """Assumes dataframe has gene and time point averages` in first row."""
    plt.figure(figsize=(12,8))
    ax = plt.subplot(111)
    for gene in genes:
        # time points, average for the index, std for index
        gene_avg = data[data['gene'] == gene][keys].values.flatten()
        try:
            gene_std = data[data['gene'] == gene][[key + 'std' for key in keys]].values.flatten()

            plt.errorbar(range(len(keys)), gene_avg,
                     yerr=stats.t.ppf(0.95, num_per_keys - 1) * gene_std, label=gene)
        except KeyError:
            plt.errorbar(range(len(keys)), gene_avg, label=gene)
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.xlabel("Time points", fontsize=20)
    plt.ylabel("Expression level", fontsize=20)
    plt.title(title, fontsize=20)
    plt.show()

def plot_gene_pairs(data, gene_pairs, keys = ["t00", "t05", "t1_", "t2_", "t3_", "t4_", "t5_", "t6_", "t7_", "t8_", "t10_", "t12_"],
                    title="Paired Gene Expression From a few species",
               num_per_keys  = np.array([4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4]), no_error = False):
    """Assumes dataframe has gene and time point averages` in first row."""
    colors = ['red', 'blue', 'green', 'cyan', 'magenta', 'yellow']
    plt.figure(figsize=(12,8))
    ax = plt.subplot(111)
    for gene_pair, i in zip( gene_pairs, range(len(gene_pairs))):
        # time points, average for the index, std for index

        gene1, gene2 = gene_pair

        print gene1, gene2
        gene1_avg = data[data['gene'] == gene1][keys].values.flatten()
        gene2_avg = data[data['gene'] == gene2][keys].values.flatten()
        if no_error:
            plt.errorbar(range(len(keys)), gene1_avg, label=gene1, color=colors[i])
        else:
            try:
                gene1_std = data[data['gene'] == gene1][[key + 'std' for key in keys]].values.flatten()

                plt.errorbar(range(len(keys)), gene1_avg,
                         yerr=stats.t.ppf(0.95, num_per_keys - 1) * gene1_std, label=gene1, color=colors[i])
            except KeyError:
                plt.errorbar(range(len(keys)), gene1_avg, label=gene1, color=colors[i])

        if no_error:
            plt.errorbar(range(len(keys)), gene2_avg, label=gene2, color=colors[i], linestyle='dashed')
        else:
            try:
                gene2_std = data[data['gene'] == gene2][[key + 'std' for key in keys]].values.flatten()

                plt.errorbar(range(len(keys)), gene2_avg,
                         yerr=stats.t.ppf(0.95, num_per_keys - 1) * gene2_std, label=gene2, color=colors[i], linestyle='dashed')
            except KeyError:
                plt.errorbar(range(len(keys)), gene2_avg, label=gene2, color=colors[i], linestyle='dashed')


    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.xlabel("Time points", fontsize=20)
    plt.ylabel("Expression level", fontsize=20)
    plt.title(title, fontsize=20)
    plt.show()



def get_sig_gene_pairs(sig_matr, genes):
    """
    :param sig_matr: Matrix of indices where significant
    :param genes: List of genes
    :return: List of pairs of gene by indices
    """

    rows, cols = np.where(sig_matr)

    grows, gcols = [genes[row] for row in rows], [genes[col] for col in cols]

    return zip(grows, gcols)



def compare_sig_matr(sig_matr_list):
    """
    Returns:
    true-values across all matrices
    tuples of # in, same, sorted
    total # sig
    """

    num = len(sig_matr_list)

    sig_matr_sum = np.sum(np.array(sig_matr_list), axis=0)

    all_sig_matr = sig_matr_sum >= num

    all_sig_num = len(np.where(all_sig_matr.flatten())[0])

    not_sig_num = len(np.where(sig_matr_sum.flatten())[0]) - all_sig_num

    return all_sig_matr, all_sig_num, not_sig_num

def randomize_geneTS(geneTS, replace=False, num_genes=None, by_time=True):
    """
    Given an n x T (T is timepoints) matrix of gene expression levels over time, return a new matrix with randomization.

    Replace=True means don't sample with replacement the new values
    """

    if num_genes==None:
        num_genes = geneTS.shape[0]
        T = geneTS.shape[1]

    if by_time:

        newgeneTS = np.array([np.random.choice(geneTS[i], size=T, replace=replace) for i in range(num_genes)])

        return newgeneTS

    else:
        TbyG = geneTS.T


        newTbyG = np.array([np.random.choice(TbyG[i], size=num_genes, replace=replace) for i in range(TbyG.shape[0])])

        newgeneTS = newTbyG.T

        return newgeneTS