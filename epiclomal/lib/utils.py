'''
Created on 2017-01-17

@author: Mirela Andronescu
'''
from __future__ import division

from scipy.special import logsumexp as log_sum_exp
from scipy.special import gammaln as log_gamma, psi

import numpy as np
import pandas as pd

def compute_e_log_dirichlet(x):
    # E[log(x)] = psi(alpha_i) - psi(alpha_0), where alpha_0 is sum_i=1^K alpha_i and x is Dirichlet - distributed
    # psi(z) is the digamma function, the derivative of the logarithm of the gamma function evaluated at z
    return psi(x) - psi(np.expand_dims(x.sum(axis=-1), axis=-1))


# THis function obsolete, moved inside the class
def get_indicator_matrix(states, X):
# for every position that is not missing in X, return 0 if it has state s or 1 otherwise
# all the missing positions are replaced with 0
    Y = np.zeros((len(states), X.shape[0], X.shape[1]))   # size SxNxM, now SxNxRxL

    for i, s in enumerate(states):
        Y[i, :, :] = (X == s).astype(int)

    return Y

def compute_e_log_p_dirichlet(posterior, prior):
# For p we use both the prior and the hidden variable, such as gamma_star and gamma_0
# This computes the E log of Dirichlet (x, alpha), where alpha=prior x=posterior
# Dirichlet(alpha) = gamma(sum(alpha))/prod(gamma(alpha)) prod x_i^{alpha_i-1}
# log(Dirichlet(alpha)) = log_gamma(sum(alpha)) - sum(log_gamma(alpha)) + sum(alpha_i-1 * log(x_i))
# Elog(Dirichlet(alpha)) = log_gamma(sum(alpha)) - sum(log_gamma(alpha)) + sum(alpha_i-1 * Elog(x_i))
    log_p = log_gamma(prior.sum()) - \
            log_gamma(prior).sum() + \
            safe_multiply(prior - 1, compute_e_log_dirichlet(posterior)).sum()

    return log_p

def compute_e_log_q_dirichlet(x):
    # for q, we don't use the priors of the Dirichlets (alpha_0, beta_0, gamma_0), just the hidden variables
    # x is the prior

    a_0 = x.sum()

    K = len(x)
    # I checked, and the following give exactly the same values
    # print ('Val1 ', log_gamma(a_0) - log_gamma(x).sum() + safe_multiply(x - 1, compute_e_log_dirichlet(x)).sum())
    # print ('Val2 ', log_gamma(a_0) - log_gamma(x).sum() + safe_multiply(x - 1, psi(x)).sum() - safe_multiply(a_0 - K, psi(a_0)))
    return log_gamma(a_0) - log_gamma(x).sum() + safe_multiply(x - 1, psi(x)).sum() - safe_multiply(a_0 - K, psi(a_0))


def compute_e_log_q_discrete(log_x):
    return np.sum(safe_multiply(np.exp(log_x), log_x))

def safe_multiply(x, y):
    # MA: if x or y are 0, then log doesn't apply. So I add 1e-10.
    return np.sign(x) * np.sign(y) * np.exp(np.log(np.abs(x+1e-20)) + np.log(np.abs(y+1e-20)))
    # return np.sign(x) * np.sign(y) * np.exp(np.log(np.abs(x)) + np.log(np.abs(y)))
    # without adding anything, I get [ 0.47323602  0.52676398]
    # adding 1e-20, I get [ 0.47323602  0.52676398]

def log_space_normalise(log_X, axis=0):
    return log_X - np.expand_dims(log_sum_exp(log_X, axis=axis), axis=axis)

def init_log_pi_star(K, N, initial_clusters_data):
    # This function assigns random or given clusters to the Z variable, and then takes the log
    log_pi_star = np.zeros((N, K))

    if (initial_clusters_data is not None):
        # Epiclomal clusters start from 0, making the initial clusters to start from 0 too
        # Checking to see what is the min
        mincl = initial_clusters_data.iloc[:,0].min()
        # print 'mincl ', mincl
        initial_clusters_data = initial_clusters_data - mincl

        # initialize to this and exit
        labels = np.reshape(initial_clusters_data.values, N)

    else:
        labels = np.random.random(size=(N, K)).argmax(axis=1)
        # np.random.random returns random float in the interval [0,1)
        # np.random.random(size=(N, K)) is a matrix of size NxK
        # labels will be the cluster that has the largest sample for every row 1..N
        # print 'Labels inside init_pi_star: {0}'.format(labels)


        # print 'Range of len of set of labels: {0}'.format(range(len(set(labels))))

        # MA: I think this is wrong. For example for K=5, N=10, if there wasn't any label 3, as in the example
        #  labels inside init_pi_star: [2 2 0 1 0 0 2 1 0 4], then Range of len of set of labels is [0, 1, 2, 3]
        #  then the last cell doesn't get any 1
        # for i, s in enumerate(range(len(set(labels)))):
        #     log_pi_star[:, i] = (labels == s).astype(int)
        # Replacing len(set(labels)) with K
    print ('Initial labels:')
    print (*labels)

    for i, s in enumerate(range(K)):
        log_pi_star[:, i] = (labels == s).astype(int)


    # now log_pi_star  of size NxK has, for each row, one 1 for the randomly chosen cluster, and 0 for the other clusters

    # print 'pi_star before log'
    # print(log_pi_star)
    log_pi_star = np.log(log_pi_star + 1e-10)       # adding 1e-10 to avoid undetermined log(0)
    # I tried adding 1e-20 instead of 1e-10 and it seems to make no difference
    # with 1e-20 0.526767468326
    # with 1e-10 0.526767468326
    # It makes a very slight difference (of 0.000005) if we have K=10 or K=2

    # Now make sure they sum up to 1 over columns
    return log_space_normalise(log_pi_star, axis=1)


def load_labels(cell_ids, file_name):
    if file_name is None:
        labels = None

    else:
        labels = pd.read_csv(file_name, compression='gzip', index_col='cell_id', sep='\t')

        labels = labels.loc[cell_ids, 'cluster']

    return labels
