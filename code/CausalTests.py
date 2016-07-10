__author__ = 'jlu96'

#GRANGER CAUSALITY
import sys
import time
from sklearn import linear_model
import scipy.stats as stats
import itertools
import geneTSmunging as gtm
import parallel_compute_working as pac

from scipy.spatial import cKDTree
import numpy as np
import pickle

import CausalParams as args

def get_lagged_vector(X_list, index, tau, E):
    """Generate lagged vector of X: (X[index], X[index - tau], ..., X[index - tau*E - 1]).

    :param X_list: list of X observations at equally spaced time intervals.
    :param index: index to start from
    :param tau: time lag
    :param E: dimension of lagged vector

    :return: (X[index], X[index - tau], ..., X[index - tau*E - 1])
    """
    assert index >= tau * (E - 1)
    return [X_list[index - e * tau] for e in range(E)]


def pairwise_granger_causality(i, j, model_order):
    """
    Computes p-value of fit to determine if i granger causes j.
    i and j must have length at least model_order (to construct lagged vectors)
    Returns: pvalue
    """
    if len(i) != len(j):
        raise ValueError("i and j must have same length")

    m = len(i) - model_order
    if m <= 2 * model_order:
        raise ValueError("i and j must have length > 2*model_order + model_order, so that target has legnth 2model_order" +
                "with 2model_order predictors")

    T = len(i)
    #lagged_i goes from p up to T, the length of i. It is a T-p matrix
    lagged_i = np.array([get_lagged_vector(i, index, 1, model_order) for index in range(model_order-1, T-1)])

    lagged_j = np.array([get_lagged_vector(j, index, 1, model_order) for index in range(model_order-1, T-1)])

    target = np.array(j[model_order:T])


    # Restricted classification
    clf_j = linear_model.LinearRegression()
    clf_j.fit(lagged_j, target)
    pred_j = clf_j.predict(lagged_j)

    SSE_j = np.sum(np.power(target - pred_j, 2))

    # Unrestricted classification
    clf_ij = linear_model.LinearRegression()
    lagged_ij = np.concatenate((lagged_i, lagged_j), axis=1)
    clf_ij.fit(lagged_ij, target)
    pred_ij = clf_ij.predict(lagged_ij)

    SSE_ij = np.sum(np.power(target - pred_ij, 2))


    beta_i = clf_ij.coef_[0:lagged_i.shape[1]]

    # Formula is taken from: http://pages.uoregon.edu/aarong/teaching/G4075_Outline/node4.html
    # The variance explained by i itself compared to variance of just j
    # IF i explains more variance then we reject
    #DOF of first is p = 2p - p = # in ij fit - # in j fit (# of params added)
    # DOF of bottom is is n - p - 1 (where p is the number of parameters in the bottom model)
    F = ((SSE_j - SSE_ij) * 1.0 / model_order)/(SSE_ij/(len(target) - 2 * model_order - 1))

    p_value = 1.0 - stats.f.cdf(F, model_order, len(target) - 2 * model_order - 1)

    return p_value, beta_i


def pairwise_granger_causality_all(matr, pairs, model_order=2, use_processes=False, procnum=32,
                                   add_ones=False):
    """
    Assume matr is an n by T matrix, where n is number of genes, T is # timepoints.
    """
    n, T = np.array(matr).shape

    p_matr = np.zeros(shape=(n,n))

    beta_matr = np.zeros(shape=(n,n,model_order))

    if pairs == None:
       pairs = itertools.permutations(range(n), 2)

    if use_processes:
        pairs = list(pairs)

        print pairs[0:10], pairs[-10:]

        function = pairwise_granger_causality_process
        args = [matr, model_order, pairs]
        input_list = pairs
        input_index = 2
        partition_input_function = pac.partition_inputs
        join_functions_dict = {0: join_pvalue_matrices, 1: join_pvalue_matrices}

        p_matr, beta_matr = pac.parallel_compute_new(function, args, input_list, input_index, partition_input_function,
                                           join_functions_dict, number=procnum, multi_return_values=True)





    else:
        for pair in pairs:
            i, j = pair
            p_matr[i][j], beta_matr[i][j] = pairwise_granger_causality(matr[i], matr[j], model_order)

    if add_ones:
        p_matr += np.diagflat(np.ones(p_matr.shape[0]))

    return (p_matr, beta_matr)

def pairwise_granger_causality_process(matr, model_order, pairs):
    """Return a matrix of zeros with the causal p-values return for the pairs.
    """

    n, T = np.array(matr).shape

    p_matr = np.zeros(shape=(n,n))
    beta_matr = np.zeros(shape=(n,n,model_order))


    for pair in pairs:
        i, j = pair
        p_matr[i][j], beta_matr[i][j] = pairwise_granger_causality(matr[i], matr[j], model_order)

    return p_matr, beta_matr

def join_pvalue_matrices(p_matr_a, p_matr_b):
    return p_matr_a + p_matr_b

def get_CCM(X_list, Y_list, L, tau=1, E=3, test_indices=None, num_test=100, use_same=True):
    """
    Compute the correlation coeffiicent of using X's cross-mapped estimates to
    predict Y at the test_indices.

    :param X_list: list of X observations
    :param Y_list: list of Y observations
    :param L: Library length, i.e. number of X/Y observations to use to construct manifold for estiamtion
    :param  tau: Time lag of lagged vecotr
    :param E: dimesnion of lagged vector
    :param test_indices: indices of Y to estimate. All shou be > than L. Default is random.
    :param num_test: number of indeces to test if random
    :param use_same: If we find an X that has the exact value of the test X, then use only those Xs to generate estiamtes of Y

    :return: rho, the correlation coefficient of the estimates and true alues

    Workflow
    1 ) make all lagged vectors
    2) assign first L - first_lag to train, in kdtree
    3) predict by taking test_indices, convert
    4) find closets tin lag tree
    5) calculate weights
    6) get indices in old
    7) make estimate
    """

    length = len(X_list)
    first_lag = tau * (E - 1)

    train_indices = np.arange(first_lag, L)
    other_indices = np.arange(L, length)
    all_indices = np.arange(first_lag, length)
    if test_indices == None:
        test_indices = np.random.choice(other_indices, num_test, replace=False)

    # make all lagged vectors
    x_lag_list = np.array([get_lagged_vector(X_list, i, tau, E) for i in all_indices])
    y_lag_list = np.array([get_lagged_vector(X_list, i, tau, E) for i in all_indices])


    # put training X and Y (used for estimation) into kdtree for nearest neighbors seaerch
    x_lag_train = x_lag_list[np.ix_(train_indices - first_lag)]

    x_lagtree = cKDTree(x_lag_train, leafsize=100)




    Y_target = Y_list[test_indices]
    Y_ests = []

    # generate each estimate
    for k in range(len(test_indices)):
        test_index = test_indices[k]
        # for each t, find contemporaneous x[t]
        x_lag = x_lag_list[test_index - first_lag]



        # Find e+1 nearest neighbors, Calculate distances
        distances, indices = x_lagtree.query(x_lag, k=E+1)
        min_dist = min(distances)


        # Case 1: we find an X that has the exact same value as the test X
        # In this case, use only those Xs that have same value as test X, and take
        # the average of Ys
        if (use_same and min_dist == 0) or all([dist == 0 for dist in distances]):
            zipped = zip(*[(dist, i) for dist, i in zip(distances, indices) if dist == 0])
            distances, indices = zipped[0], zipped[1]
            weights = [1.0 / len(distances)] * len(distances)
            Y_est = sum([weight * Y_list[i + first_lag] for weight, i in zip(weights, indices)])

        # Case 2: all Xs are different
        # Use exponential weighting as in paper
        else:
            zipped = zip(*[(dist, i) for dist, i in zip(distances, indices) if dist != 0])
            distances, indices = zipped[0], zipped[1]
            min_dist = min(distances)

            # Calculate weights
            weights = np.array([np.exp(-dist * 1.0/ min_dist) for dist in distances])

            weights /= sum(weights)

            # generate estimate y^
            Y_est = sum([weight * Y_list[i + first_lag] for weight, i in zip(weights, indices)])

        Y_ests.append(Y_est)

    rho = np.corrcoef(Y_ests, Y_target)[0][1]


    return rho


def get_CCM_complete(X_list, Y_list, Ls, tau=1, E=3, test_indices=None, num_test=100, use_same=True):
    """
    :param X_list: list of X observations
    :param Y_list: list of Y observations
    :param Ls: Library lengths, i.e. number of X/Y observations to use to construct manifold for estiamtion
    :param  tau: Time lag of lagged vecotr
    :param E: dimesnion of lagged vector
    :param test_indices: indices of Y to estimate. All shou be > than L. Default is random.
    :param num_test: number of indeces to test if random
    :param use_same: If we find an X that has the exact value of the test X, then use only those Xs to generate estiamtes of Y

    :return: rhos, the correlation coefficient for each value of L
    """


    rhos = [get_CCM(X_list, Y_list, L, tau=tau, E=E, test_indices=test_indices, num_test=num_test, use_same=use_same)
           for L in Ls]

    for i in range(len(rhos)):
        if np.isnan(rhos[i]):
            rhos[i] = 0

    return rhos, Ls




def main():
    tstart = time.time()


    input_file = args.input_file
    out_file_prefix = args.out_file_prefix


    start_index = args.start_index
    end_index = args.end_index


    df = gtm.load_file_and_avg(input_file)

    genes = df['gene'][start_index:end_index].values

    found_genes, geneTS = gtm.get_gene_TS(df, genes)

    cause_type = args.cause_type

    if cause_type == 'g':
        model_orders = range(args.model_order_min, args.model_order_max + 1)

        threshold = args.p_threshold

        p_matr_list = []
        sig_matr_list = []

        for model_order in model_orders:
            t_gc = time.time()
            p_matr = pairwise_granger_causality_all(geneTS, model_order=model_order, use_processes=args.use_processes, procnum=args.procnum)
            print "Time for granger causality", time.time() - t_gc


            sig_matr = p_matr < threshold

            p_matr_list.append(p_matr)
            sig_matr_list.append(sig_matr)



        all_sig_matr, all_sig_num, not_sig_num = gtm.compare_sig_matr(sig_matr_list=sig_matr_list)

        print "Total number of significant pairs ", all_sig_num + not_sig_num
        print "Pairs significant across all matrices ", all_sig_num, all_sig_num * 1.0 / (all_sig_num + not_sig_num)


        out_file_name = out_file_prefix + "_GC.p"
        pickle.dump([model_orders, p_matr_list, sig_matr_list, (all_sig_matr, all_sig_num, not_sig_num)], open(out_file_name, "w"))

        print "Results written  to", out_file_name



    # compare the significant matrices

    # save the output p matrices

    print "Total time used ", time.time() - tstart

if __name__ == '__main__':
    main()

