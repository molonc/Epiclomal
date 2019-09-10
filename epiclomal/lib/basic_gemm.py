'''
Created on 2017-01-17

@author: Mirela Andronescu
'''

# Very nice blog about einsum and how to convert it into nested sums
# https://obilaniu6266h16.wordpress.com/2016/02/04/einstein-summation-in-numpy/

from __future__ import division

import numpy as np
import pandas as pd
from numba import njit, prange

from lib.utils import compute_e_log_dirichlet, compute_e_log_q_dirichlet, compute_e_log_p_dirichlet, \
                      compute_e_log_q_discrete, init_log_pi_star, log_space_normalise, safe_multiply

from scipy.stats import dirichlet
from sklearn.metrics import mean_absolute_error

from collections import defaultdict
from random import random, shuffle

from math import exp

class BasicGeMM(object):
    # @profile
    def __init__(self,
                 gamma_prior,
                 alpha_prior,
                 beta_prior,
                 X,
                 regions,
                 initial_clusters_data,
                 mu_has_k,
                 Bishop_model_selection = False,
                 slsbulk_data = None):

        self.K = len(alpha_prior)               # max number of clusters
        self.mu_has_k = mu_has_k
        self.Bishop_model_selection = Bishop_model_selection
        self.regions = regions
        self.slsbulk_data = slsbulk_data
        print ('Max number of clusters: ', self.K)
        print ('Including k in mu: ', self.mu_has_k)
        print ('Using Bishop model selection: ', self.Bishop_model_selection)

        # In python3 I am replacing X.keys() to list(X.keys()) (list(X) would work too)
        self.N = X[list(X.keys())[0]].shape[0]        # number of cells

        self.gamma_prior = gamma_prior

        self.alpha_prior = alpha_prior

        # this is only used in the case of Bishop model selection
        self.pi = np.ones(self.K)
        for k in range(self.K):
            self.pi[k] = 1/self.K
        # print(*self.pi)

        if (self.Bishop_model_selection is False):
            self._init_alpha_star()
        # MA: this just initializes with 1 for every clone. if 4 clones, it'll be [1 1 1 1]

        self.data_types = list(X.keys())

        # I am actually implementing it as a one region, that way we can extend to regions more easily
        self.M = {}
        self.R = {}
        self.maxL = {}    # maxL = max_r(l)  = the maximum region size

        self.S = {}

        self.T = {}

        self.IX = {}

        self.X = {}

        self.gamma_star = {}

        self.beta_prior = {}
        self.include_bulk = {}      # true if bulk is used as prior

        self.beta_star = {}

        self.log_mu_star = {}

        for data_type in self.data_types:
            if X[data_type].shape[0] != self.N:
                raise Exception('All data types must have the same number of rows (cells).')

            self.M[data_type] = X[data_type].shape[1]
            # M is the number of loci

            self._set_region_params(data_type, regions)

            self.S[data_type] = gamma_prior[data_type].shape[0]
            # S is 2 for 2 states

            self.T[data_type] = gamma_prior[data_type].shape[1]
            # T is 2 when gamma_prior is [99, 1; 1, 99]

            self.X[data_type] = self._region_data_matrix(data_type, X[data_type])

            self.IX[data_type] = self._get_indicator_X(data_type, self.X[data_type])
            # size TxNxM, now TxNxRxmaxL
            # self.IX is actually  the indicator matrix I(X_nm=t)
            # Now self.IX contains 2 matrices:
            # the first is I(Xnm=0); this is the opposite of the input matrix
            #   (0 instead of 1 and 1 instead of 0, and with 0 instead of NaNs)
            # the second matrix is I(Xnm=1) and in this case (2 states) it is the same as the input matrix, and with 0 instead of NaNs

            self._init_gamma_star(data_type)
            # This just initializes gamma_star with gamma_prior

            self._init_beta_prior(data_type, beta_prior)
            # This just initializes beta_prior, one vector of size S if bulk is not used,
            #   or if bulk is used as prior, a vector of size S for each locus, matrix RxmaxLxS

            self._init_beta_star(data_type)
            # This just initializes beta_star for each cluster k with beta_prior, matrix KxS, now KxRxS

        self.log_pi_star = init_log_pi_star(self.K, self.N, initial_clusters_data)
        # This function assigns random or given clusters to the Z variable.
        # For each cell (row), one of the k values will be 1 and the others 0.
        # Then, return the log of this matrix (log_pi_star) of size N rows x K columns

        self.pi_star = self.set_pi_star()

        self.lower_bound = [float('-inf')]
        self.log_likelihood = float('-inf')
        self.log_posterior = [float('-inf')]

        self._debug_lower_bound = [float('-inf')]

        self.converged = False


    def _set_region_params(self, data_type, regions):
        self.R[data_type] = 1       # 1 region for the basic model
        self.maxL[data_type] = self.M[data_type]

    def _region_data_matrix(self, data_type, X):
    # transform data from unregioned (index m) to regioned (index rl)
    # this will be done differently in the regions class
        return X.values.reshape(self.N, self.R[data_type], self.maxL[data_type])

    # NOT used
    def _unregion_data_matrix(self, data_type, matrix):
    # transform data back from regioned (index rl) to unregioned (index m)
    # this will be done differently in the regions class
        return matrix.reshape(self.N, self.M[data_type])

    def unregion_mu_star(self, data_type):
    # transform data back from regioned (index rl) to unregioned (index m)
        mu_star = self.get_mu_star(data_type)
        states = range(self.S[data_type])
        # this also takes the argmax of mu_star, we don't want that
        # return np.argmax(mu_star.reshape(self.K, self.M[data_type], len(states)), axis=2)
        return mu_star.reshape(len(states), self.K, self.M[data_type])

    def _get_indicator_X(self, data_type, X):
    # for every position that is not missing in X, return 0 if it has state s or 1 otherwise
    # all the missing positions are replaced with 0
        Y = np.zeros((self.S[data_type], self.N, self.R[data_type], int(self.maxL[data_type])))   # size SxNxM, now SxNxRxmaxL
        # Y = np.zeros((self.S[data_type], self.N, self.M[data_type]))
        total = 0
        for i, s in enumerate(range(self.S[data_type])):
            Y[i, :, :, :] = (X == s).astype(int)
            ti = np.sum(Y[i,:,:,:])
            print("Number of observed ", i, "'s: ", ti)
            total = total + ti
        entries = self.N*self.M[data_type]
        print ("Total number of loci in all cells: ", entries)
        print ("Missing proportion: ", (entries-total)/entries)
        return Y


    def _init_gamma_star(self, data_type):
        self.gamma_star[data_type] = self.gamma_prior[data_type].copy()


    def _init_beta_prior(self, data_type, beta_prior):
        if (len(beta_prior[data_type]) == self.S[data_type]):
            # not including bulk as prior
            self.include_bulk[data_type] = False
            self.beta_prior[data_type] = beta_prior[data_type]
        else:
            # yes including bulk as prior
            self.include_bulk[data_type] = True
            # beta_prior will be a matrix of size RxmaxLxS
            self.beta_prior[data_type] = np.zeros((self.R[data_type], int(self.maxL[data_type]), self.S[data_type]))

            # this should be done differently in the regions class
            # NOTE: I am adding 1 to all the bulk reads to avoid 0 values. This is equivalent to adding a very small extra prior
            self.beta_prior[data_type] = 1 + beta_prior[data_type].values.reshape(self.R[data_type], int(self.maxL[data_type]), self.S[data_type])

            #print self.beta_prior[data_type]


    def _init_beta_star(self, data_type):
        if (self.mu_has_k):
            # include k
            if (self.include_bulk[data_type]):
                # this includes a mu for each locus
                # MA: this has to be a matrix of size KxS, with bulk it is KxRxmaxLxS
                self.beta_star[data_type] = np.zeros((self.K, self.R[data_type], int(self.maxL[data_type]), self.S[data_type]))
                for k in range(self.K):
                     self.beta_star[data_type][k] = self.beta_prior[data_type].copy()
            else:
                # MA: this has to be a matrix of size KxS, without bulk it is KxRxS
                self.beta_star[data_type] = np.zeros((self.K, self.R[data_type], self.S[data_type]))
                for k in range(self.K):
                    for r in range(self.R[data_type]):
                        self.beta_star[data_type][k][r] = self.beta_prior[data_type].copy()
        else:
            # do not include k
            if (self.include_bulk[data_type]):
                # this includes a mu for each locus
                # MA: this has to be a matrix of size S, with bulk it is RxmaxLxS
                self.beta_star[data_type] = np.zeros((self.R[data_type], int(self.maxL[data_type]), self.S[data_type]))
                self.beta_star[data_type] = self.beta_prior[data_type].copy()
            else:
                # MA: this has to be a matrix of size S, without bulk it is RxS
                self.beta_star[data_type] = np.zeros((self.R[data_type], self.S[data_type]))
                for r in range(self.R[data_type]):
                    self.beta_star[data_type][r] = self.beta_prior[data_type].copy()



    def _init_alpha_star(self):
        self.alpha_star = np.ones(self.K)


    def get_e_log_epsilon(self, data_type):
        # computes expectation(log(epsilon)) using the digamma function
        return compute_e_log_dirichlet(self.gamma_star[data_type])

    def get_e_log_mu(self, data_type):
        # computes expectation(log(mu)) using the digamma function
        # I tested this, and it returns K rows, as if I had a for loop
        # NOTE: the values of self.beta_star[data_type] shouldn't be 0, otherwise it returns inf
        return compute_e_log_dirichlet(self.beta_star[data_type])

    def get_e_log_pi(self):
        # computes expectation(log(pi)) using the digamma function
        return compute_e_log_dirichlet(self.alpha_star)

    def get_mu_star(self, data_type):
        return np.exp(self.log_mu_star[data_type])

    @property
    def mu_star(self):
        mu_star = {}

        for data_type in self.data_types:
            mu_star[data_type] = self.get_mu_star(data_type)

        return mu_star

    def set_pi_star(self):
        return np.exp(self.log_pi_star)

    # @profile
    def fit(self, convergence_tolerance=1e-4, debug=False, num_iters=100):
        print ("Iter  ELBO difference")

        # print 'pi_star for cell 14 initially: ', self.pi_star[13,]

        for i in range(num_iters):

            # update E(I(Gmk=s))
            self._update_mu_star()

            if debug:
                print ('ELBO, diff after update_mu_star', self._diff_lower_bound())

            if (self.Bishop_model_selection is False):
                # update alpha_star
                self._update_alpha_star()
                if debug:
                    print ('ELBO, diff after update_alpha_star', self._diff_lower_bound())

            # update beta_star
            self._update_beta_star()

            if debug:
                print ('ELBO, diff after update_beta_star', self._diff_lower_bound())

            # update gamma_star
            self._update_gamma_star()

            if debug:
                print ('ELBO, diff after update_gamma_star', self._diff_lower_bound())

            # update pi_star
            self._update_pi_star()

            if debug:
                print ('ELBO, diff after update_pi_star', self._diff_lower_bound())

            # print 'pi_star for cell 14 after pi_star update: ', self.pi_star[13,]

            if (self.Bishop_model_selection is True):
                pi_star_sum = self.pi_star.sum(axis=0)

                for k in range(self.K):
                    self.pi[k] = pi_star_sum[k] / self.N


            # update rho_star, but here in the basic model nothing happens
            self._update_rho_star()
            if debug:
                print ('ELBO, diff after update_rho_star', self._diff_lower_bound())


            # From Blei 2016: computing the ELBO for the full data set may be too expensive, we can compute on a held-out set
            (ELBO, ElogP, ElogP_X, ElogQ) = self._compute_lower_bound()
            self.lower_bound.append(ELBO)

            diff = (self.lower_bound[-1] - self.lower_bound[-2]) / np.abs(self.lower_bound[-1])

            print (i, self.lower_bound[-1], diff)

            if abs(diff) < convergence_tolerance:
                print ('Converged')

                self.converged = True
                self.ElogP = ElogP
                self.ElogQ = ElogQ

                self._unregion_rho_star()

                print ('ELogP:', self.ElogP)
                print ('ELogQ:', self.ElogQ)
                print ('ELBO:', ELBO)

                self.mean_or_mode = "mean"
                loglik = self._compute_log_likelihood()
                print ('Log likelihood mean: ', loglik)
                self.log_likelihood_mean = loglik
                self.log_likelihood = loglik
                (loglik1, loglik2) = self._compute_log_likelihood_times_priors()
                print ('Log posterior mean unnormalized cluster Ks: ', loglik1)
                print ('Log posterior mean unnormalized all Ks: ', loglik2)
                self.log_posterior_mean_clusterK = loglik1
                self.log_posterior_mean_allK = loglik2
                self.DIC_pd = -2 * ElogP_X + 2 * loglik
                self.DIC_measure = -4 * ElogP_X + 2 * loglik
                print('DIC_pd: ', self.DIC_pd)
                print('ElogP_X: ', ElogP_X)
                print('loglik: ', loglik)
                print('DIC_measure: ', self.DIC_measure)

                self.nclusters = self._compute_num_clusters()   # this one take the clusters with highest probability; sometimes clusters 0 and 1 are the same clusters but with probability 0.51 and 0.49
                self.prevalences = self._compute_cluster_prevalence()
                print('Obtained ', self.nclusters, ' clusters with prevalences ')
                print(*self.prevalences)
                self.nnonzeropi = self._compute_num_nonzero_pi()
                print ('Obtained ', self.nnonzeropi, ' clusters with nonzero pi')

                break

            elif diff < 0:
                print ('Lower bound decreased')

                if not debug:
                    self.converged = False

                    break

    ####################

    def _update_rho_star(self):
        # Nothing happens here
        # Different in basic_miss_gemm
        return

    def _unregion_rho_star(self):
        # Nothing happens here
        # Different in basic_miss_gemm
        return

    ####################

    def _update_mu_star(self):
        for data_type in self.data_types:
            self._update_mu_star_d(data_type)

    def _update_mu_star_d(self, data_type):

        log_mu_star = self._update_mu_star_big_sum(data_type)

        # Now log_mu_star has 3 dimensions: s (2 groups), k (total number of clusters), m (num loci)
        # With 1 region, log_mu_star has 4 dimensions: s, k, r, maxl

        # KxS
        # Here I am just transposing the matrix
        if (self.mu_has_k):
            if (self.include_bulk[data_type]):
                e_log_mu = np.einsum('krls -> skrl',self.get_e_log_mu(data_type))
                # Here we just need a sum, not einsum
                log_mu_star = e_log_mu[:, :, :, :] + log_mu_star
            else:
                e_log_mu = np.einsum('krs -> skr',self.get_e_log_mu(data_type))
                log_mu_star = e_log_mu[:, :, :, np.newaxis] + log_mu_star
        else:
            if (self.include_bulk[data_type]):
                e_log_mu = np.einsum('rls -> srl',self.get_e_log_mu(data_type))
                # Here we just need a sum, not einsum
                log_mu_star = e_log_mu[:, np.newaxis, :, :] + log_mu_star
            else:
                e_log_mu = np.einsum('rs -> sr',self.get_e_log_mu(data_type))
                log_mu_star = e_log_mu[:, np.newaxis, :, np.newaxis] + log_mu_star

        # SxKxM, now SxKxRxmaxL
        self.log_mu_star[data_type] = log_space_normalise(log_mu_star, axis=0)

    def _update_mu_star_big_sum(self, data_type):
        # This function will be different in basic_miss_gemm

        IX = self.IX[data_type]  # TxNxM, now TxNxRxmaxL

        # Now it's taking the e_log_epsilon
        e_log_epsilon = self.get_e_log_epsilon(data_type)   # SxT

        # NOTE: einsum seems much faster than the nested for loops!!!

        # the summations below mean
        # log q(G_km=s) = \sum_t \sum_n I(X_nm=t) * E[I(Z_n=k)] * E[log(epsilon_st)] + log(1/2)
        # I am converting m to 1 region of length m, so instead of index m now I have 2 indeces r and l
        # return np.einsum('stnkm, stnkm, stnkm -> skm',
        #                  IX[np.newaxis, :, :, np.newaxis, :],       # IX depends on t (2 groups), n and m
        #                  self.pi_star[np.newaxis, np.newaxis, :, :, np.newaxis],  # pi_star depends on n and k
        #                  e_log_epsilon[:, :, np.newaxis, np.newaxis, np.newaxis])   # e_log_epsilon depends on s and t

        return np.einsum('stnkrl, stnkrl, stnkrl -> skrl',
                         IX[np.newaxis, :, :, np.newaxis, :, :],       # IX depends on t (2 groups), n and m
                         self.pi_star[np.newaxis, np.newaxis, :, :, np.newaxis, np.newaxis],  # pi_star depends on n and k
                         e_log_epsilon[:, :, np.newaxis, np.newaxis, np.newaxis, np.newaxis])   # e_log_epsilon depends on s and t


    ####################

    def _update_gamma_star(self):
        for data_type in self.data_types:
            prior = self.gamma_prior[data_type]
            newterm = self._get_gamma_star_data_term(data_type)
            self.gamma_star[data_type] = prior + newterm
            # print 'Updated gamma_star is {0}'.format(self.gamma_star[data_type])

    def _get_gamma_star_data_term(self, data_type):
        log_mu_star = self.log_mu_star[data_type]

        IX = self.IX[data_type]

        # SxNxKxM   -- converting it to 1 region of length M
        # pi_star_nk + mu_star_skm
        # data_term = np.exp(self.log_pi_star[np.newaxis, :, :, np.newaxis] + log_mu_star[:, np.newaxis, :, :])

        # SxNxKxRxL
        # pi_star_nk + mu_star_skrl
        data_term = self._get_exp_data_term(self.log_pi_star[np.newaxis, :, :, np.newaxis, np.newaxis] + log_mu_star[:, np.newaxis, :, :, :])


        # \sum_n \sum_k \sum_m E[I(Z_n=k)] E[I(G_mk=s)] I(X_nm=t)
        # return np.einsum('stnkm, stnkm -> st',
        #                  data_term[:, np.newaxis, :, :, :],
        #                  IX[np.newaxis, :, :, np.newaxis, :])
        return np.einsum('stnkrl, stnkrl -> st',
                         data_term[:, np.newaxis, :, :, :, :],
                         IX[np.newaxis, :, :, np.newaxis, :, :])

    # faster than np.exp() for 5d array
    @staticmethod
    @njit(parallel=True, cache=True)
    def _get_exp_data_term(log_data):
        shape = log_data.shape
        data_term = np.zeros(shape)
        shape = list(shape)
        for s in prange(shape[0]):
            for c in prange(shape[1]):
                for k in prange(shape[2]):
                    for r in prange(shape[3]):
                        for l in prange(shape[4]):
                            data_term[s][c][k][r][l] = exp(log_data[s][c][k][r][l])

        return data_term

    ####################

    def _update_beta_star(self):
        # MA: added this
        for data_type in self.data_types:
            prior = self.beta_prior[data_type]

            newterm = self._get_beta_star_data_term(data_type)

            # # prior is of size [1,S], newterm is of size [K,S]. I just have to do a simple sum, einsum doesn't do the right thing
            # With bulk, prior is of size [R,maxL,S] and newterm is of size [S,K,R,maxL]
            # self.beta_star[data_type] = np.einsum('krls, krls', prior[np.newaxis,:,:,:] + newterm[:,:,:,:])

            if (self.mu_has_k):
                if (self.include_bulk[data_type]):
                    self.beta_star[data_type] = prior[np.newaxis, :, :, :] + newterm
                    # Above and below seem to be the same
                    #for k in range(self.K):
                    #    self.beta_star[data_type][k] = prior + newterm[k]
                else:
                    self.beta_star[data_type] = prior + newterm
            else:
                if (self.include_bulk[data_type]):
                    self.beta_star[data_type] = prior[np.newaxis, :, :] + newterm
                    # Above and below seem to be the same
                    #for k in range(self.K):
                    #    self.beta_star[data_type][k] = prior + newterm[k]
                else:
                    self.beta_star[data_type] = prior + newterm


    def _get_beta_star_data_term(self, data_type):
        # MA: added this
        mu_star = self.get_mu_star(data_type)        # SxKxR or SxKxRxmaxL

        if (self.mu_has_k):
            # Now beta_star goes over l too
            if (self.include_bulk[data_type]):
                return np.einsum('skrl -> krls', mu_star[:, :, :, :])
            else:
                # \sum_l mu_star
                return np.einsum('skrl -> krs', mu_star[:, :, :, :])
        else:
            # The sum goes over l and also k in this case
            if (self.include_bulk[data_type]):
                return np.einsum('skrl -> rls', mu_star[:, :, :, :])
            else:
                # sum_l sum_k mu_star
                return np.einsum('skrl -> rs', mu_star[:, :, :, :])


    ####################

    def _update_alpha_star(self):
        self.alpha_star = self.alpha_prior + self._get_alpha_star_data_term()

    def _get_alpha_star_data_term(self):
        return self.pi_star.sum(axis=0)

    ####################

    def _update_pi_star(self):
        log_pi_star = np.zeros((self.N, self.K))

        for data_type in self.data_types:
            log_pi_star = log_pi_star + self._get_log_pi_star_d(data_type)

        if (self.Bishop_model_selection is True):
            log_pi_term = np.log(self.pi + 1e-20)
        else:
            log_pi_term = self.get_e_log_pi()

        #cellnum = 12 - 1
        #print 'e_log_pi', e_log_pi[0:4]
        #print 'log_pi_star sum term', log_pi_star[cellnum,0:4]
        log_pi_star = log_pi_term[np.newaxis, :] + log_pi_star

        #print 'log_pi_star before normalize', log_pi_star[cellnum,0:4]
        self.log_pi_star = log_space_normalise(log_pi_star, axis=1)
        self.pi_star = self.set_pi_star()
        #print 'pi_star after normalize', np.exp(self.log_pi_star[cellnum,0:4])

    def _get_log_pi_star_d(self, data_type):
        mu_star = self.get_mu_star(data_type)

        IX = self.IX[data_type]

        e_log_epsilon = self.get_e_log_epsilon(data_type)

        # \sum_s \sum_t \sum_m E[I(G_mk=s)] I(X_nm=t) E[log(\epsilon_st)]
        # log_pi_star = np.einsum('stnkm, stnkm, stnkm -> nk',
        #                   mu_star[:, np.newaxis, np.newaxis, :, :],
        #                   IX[np.newaxis, :, :, np.newaxis, :],
        #                   e_log_epsilon[:, :, np.newaxis, np.newaxis, np.newaxis])

        log_pi_star = np.einsum('stnkrl, stnkrl, stnkrl -> nk',
                          mu_star[:, np.newaxis, np.newaxis, :, :, :],
                          IX[np.newaxis, :, :, np.newaxis, :, :],
                          e_log_epsilon[:, :, np.newaxis, np.newaxis, np.newaxis, np.newaxis])


        return log_pi_star


    ##########################################################

    def _compute_lower_bound(self, skip=0):
        (Elogp, Elogp_X) = self._compute_e_log_p(skip)
        Elogq = self._compute_e_log_q(skip)
        return (Elogp - Elogq, Elogp, Elogp_X, Elogq)

    ####################

    def _compute_e_log_p(self, skip=0):
    # This is the expectation of log of joint for p
    # For Basic-Gemm, this is:
    # = ElogP(X|Z,G,epsilon,pi,mu) + ElogP(Z|pi) + ElogP(G|mu) + ElogP(epsilon) + ElogP(pi) + ElogP(mu)
    # = term 1                     + term 2      + term 3      + term 4         + term 5    + term 6
    # For Basic-Miss-GeMM, there is one more term, term 7
    #       + ElogP(Xmiss|Z,G,epsilon,pi,mu)
        e_log_p_term1 = self._compute_e_log_p_term1()   # term 1
        e_log_p_term2 = self._compute_e_log_p_term2()   # term 2
        e_log_p_term3 = self._compute_e_log_p_term3()   # term 3
        e_log_p_term4 = self._compute_e_log_p_term4()   # term 4

        e_log_p_term5 = 0
        if(self.Bishop_model_selection is False):
            e_log_p_term5 = self._compute_e_log_p_term5()   # term 5

        e_log_p_term6 = self._compute_e_log_p_term6()   # term 6
        e_log_p_term7 = self._compute_e_log_p_term7()   # term 7

        # print('ElogP_epsion: ', e_log_p_term4)
        # print('ElogP_pi: ', e_log_p_term5)
        # print('ElogP_mu: ', e_log_p_term6)
        # print('ElogP_G: ', e_log_p_term3)
        # print('ElogP_Z: ', e_log_p_term2)
        # print('ElogP_X ', e_log_p_term1)


        e_log_p = sum([e_log_p_term1,
                    e_log_p_term2,
                    e_log_p_term3,
                    e_log_p_term4,
                    e_log_p_term5,
                    e_log_p_term6,
                    e_log_p_term7])

        return (e_log_p, e_log_p_term1)


    #################

    def _compute_e_log_p_term1(self):
        # This computes term 1
        # ElogP_T1 = \sum_s \sum_t [(\sum_n \sum_m \sum_k pi*_nk mu*_kms IX_nmt) * ElogEps_st]

        e_log_p = 0
        for data_type in self.data_types:
            # Elog epsilon
            e_log_epsilon = self.get_e_log_epsilon(data_type)

            # below is sum pi_star mu_star IX
            data_term = self._get_gamma_star_data_term(data_type)

            e_log_p += np.sum(safe_multiply(data_term, e_log_epsilon))

        return e_log_p

    #################

    def _compute_e_log_p_term2(self):

        # This below is \sum_n E(I(Zn=k))Elog pi _k     (term 2 in the joint)
        if (self.Bishop_model_selection is False):
            return np.sum(safe_multiply(self.get_e_log_pi(), self._get_alpha_star_data_term()))
        else:
            return np.sum(safe_multiply(self.pi, self._get_alpha_star_data_term()))

    #################

    def _compute_e_log_p_term3(self):
        # term 3 in the joint
        # \sum_k \sum_m \sum_s mu*_kms E_log_mu_kms
        e_log_p = 0

        for data_type in self.data_types:
            mu_star = self._get_beta_star_data_term(data_type)  # size KxS, \sum_m E(I(Gkm=s)), now size KxRxS, with bulk size KxRxmaxLxS
            # this is sum_l mu*_krls, size KxRxS
            e_log_mu = self.get_e_log_mu(data_type)  # size KxS, now KxRxS
            # e_log_mu is shape KxRxS

            # This below is \sum_m E(I(Gkm=s))Elog mu_ks     (term 3 in the joint)
            if (self.mu_has_k):
                if (self.include_bulk[data_type]):
                    # the einsum below does sum_k sum_r sum_l sum_s mu_star[[k,r,l,s] * e_log_mu[k,r,l,s]
                    # I checked this with nested for loops
                    e_log_p += np.einsum('krls,krls', mu_star, e_log_mu)
                    # for k in range(self.K):
                    #     for r in range(self.R[data_type]):
                    #         for l in range(self.maxL[data_type]):
                    #             for s in range(self.S[data_type]):
                    #                 e_log_p += mu_star[k,r,l,s] * e_log_mu[k,r,l,s]

                else:
#                     total = 0
#                     for k in range(self.K):
#                         val = np.einsum('rs,rs', mu_star[k], e_log_mu[k])
#                         print ('ElogP_G for k=', k, ' is: ', val)
#                         total = total + val
#                     print ('Total value: ', total)
                    e_log_p += np.einsum('krs,krs', mu_star, e_log_mu)
            else:
                if (self.include_bulk[data_type]):
                    # the einsum below does sum_k sum_r sum_l sum_s mu_star[[k,r,l,s] * e_log_mu[k,r,l,s]
                    # I checked this with nested for loops
                    e_log_p += np.einsum('krls,krls', mu_star, e_log_mu[np.newaxis,:,:,:])

                else:
                    e_log_p += np.einsum('rs,rs', mu_star, e_log_mu)

        return e_log_p

    #################

    def _compute_e_log_p_term4(self):
        e_log_p = 0

        for data_type in self.data_types:
            # term 4 \sum \sum log Dirichlet(epsilon)
            e_log_p += sum([compute_e_log_p_dirichlet(x, y) for x, y in zip(self.gamma_star[data_type], self.gamma_prior[data_type])])

        return e_log_p

    #################
    def _compute_e_log_p_term5(self):
        # E log P(pi) = E log Dirichlet(alpha_star, alpha_zero)
        # 0 for BMS
        return  compute_e_log_p_dirichlet(self.alpha_star, self.alpha_prior)
        # This above computes E log Dirichlet(alpha_star, alpha_zero)

    #################

    def _compute_e_log_p_term6(self):
        # sum_k sum_m E log P(mu_km) = sum_k sum_m E log Dirichlet(beta_star, beta_zero)
        e_log_p = 0

        for data_type in self.data_types:
            if (self.mu_has_k):
                if (self.include_bulk[data_type]):
                    for k in range(self.K):
                        for r in range(self.R[data_type]):
                            for l in range(self.maxL[data_type]):
                                e_log_p += compute_e_log_p_dirichlet(self.beta_star[data_type][k][r][l], self.beta_prior[data_type][r][l])   # term 6
                         # This above computes E log Dirichlet(beta_star, beta_zero)
                else:
                    for k in range(self.K):
                        for r in range(self.R[data_type]):
                            e_log_p += compute_e_log_p_dirichlet(self.beta_star[data_type][k][r], self.beta_prior[data_type])   # term 6
            else:
                if (self.include_bulk[data_type]):
                    for r in range(self.R[data_type]):
                        for l in range(self.maxL[data_type]):
                            e_log_p += compute_e_log_p_dirichlet(self.beta_star[data_type][r][l], self.beta_prior[data_type][r][l])   # term 6
                         # This above computes E log Dirichlet(beta_star, beta_zero)
                else:
                    for r in range(self.R[data_type]):
                        e_log_p += compute_e_log_p_dirichlet(self.beta_star[data_type][r], self.beta_prior[data_type])   # term 6
        # print ('ElogP is ', e_log_p)
        return e_log_p

    ####################

    def _compute_e_log_p_term7(self):
        # This computes term 7
        return 0

    ####################
    def _compute_e_log_q(self, skip=0):
        # For Basic-GeMM:
        # E(log q) = E(log(eps)) + E(log(pi)) + E(log(mu)) + E(log(G)) + E(log(Z))
        # For Basic-Miss-GeMM:
        # E(log q) = E(log(eps)) + E(log(pi)) + E(log(mu)) + E(log(G)) + E(log(Z)) + E(log(Xmiss))

        e_log_q_epsilon= self._compute_e_log_q_epsilon()

        # NOTE: in Genotyper with samples it goes over the unique samples only!!!
        e_log_q_pi = 0
        if (self.Bishop_model_selection is False):
            e_log_q_pi = compute_e_log_q_dirichlet(self.alpha_star)

        e_log_q_mu = self._compute_e_log_q_mu()

        e_log_q_G = 0

        e_log_q_Xmiss = 0

        for data_type in self.data_types:
#             total = 0
#             for k in range(self.K):
#                 val = compute_e_log_q_discrete(self.log_mu_star[data_type][:,k,:,:])
#                 print ('ElogQ(G) for k=', k, ' is: ', val)
#             total = total + val
#             print ('Total ElogQ(G) is: ', total)
            e_log_q_G += compute_e_log_q_discrete(self.log_mu_star[data_type])
            e_log_q_Xmiss += self._compute_e_log_q_Xmiss(data_type)

        e_log_q_Z = compute_e_log_q_discrete(self.log_pi_star)

        # print('ElogQ_epsion: ', e_log_q_epsilon)
        # print('ElogQ_pi: ', e_log_q_pi)
        # print('ElogQ_mu: ', e_log_q_mu)
        # print('ElogQ_G: ', e_log_q_G)
        # print('ElogQ_Z: ', e_log_q_Z)

        e_log_q = np.sum([e_log_q_epsilon,
                       e_log_q_pi,
                       e_log_q_mu,
                       e_log_q_G,
                       e_log_q_Z,
                       e_log_q_Xmiss])

        return e_log_q


    def _compute_e_log_q_Xmiss(self, data_type):
        return 0

    def _compute_e_log_q_epsilon(self):
        e_log_q = 0

        for data_type in self.data_types:
            for x in self.gamma_star[data_type]:
            # this gives the first row in gamma_star (having 2 elements)
                dir = compute_e_log_q_dirichlet(x)
                # dir is 1 value
                # We need the brackets below for sum to work, although in my case dir is always 1 value
                e_log_q += sum([dir])
            # A shorter way of writing the for loop above
            # e_log_q += sum([compute_e_log_q_dirichlet(x) for x in self.gamma_star[data_type]])

        return e_log_q

    def _compute_e_log_q_mu(self):
        sum = 0

        for data_type in self.data_types:
        # ADDED r regions
            sum = 0
            if (self.mu_has_k):
                if (self.include_bulk[data_type]):
                    for k in range(self.K):
                        for r in range(self.R[data_type]):
                            for l in range(self.maxL[data_type]):
                                x = self.beta_star[data_type][k][r][l]
                                sum += compute_e_log_q_dirichlet(x)
                else:
                    for k in range(self.K):
                        for r in range(self.R[data_type]):
                            x = self.beta_star[data_type][k][r]
                            sum += compute_e_log_q_dirichlet(x)
            else:
                if (self.include_bulk[data_type]):
                    for r in range(self.R[data_type]):
                        for l in range(self.maxL[data_type]):
                            x = self.beta_star[data_type][r][l]
                            sum += compute_e_log_q_dirichlet(x)
                else:
                    for r in range(self.R[data_type]):
                        x = self.beta_star[data_type][r]
                        sum += compute_e_log_q_dirichlet(x)

        return sum

    ###############

    def _diff_lower_bound(self):
        (ELBO, ElogP, ElogP_X, ElogQ) = self._compute_lower_bound()

        self._debug_lower_bound.append(ELBO)

        diff = (self._debug_lower_bound[-1] - self._debug_lower_bound[-2]) / np.abs(self._debug_lower_bound[-1])

        if diff < 0:
            print ('Bound decreased')

        return (ELBO, diff)


    ###############

    # @profile
    def _compute_log_likelihood_times_priors(self):
        logl = self.log_likelihood

        self.whichK = "cluster"     # only the ones used
        # print ('Which K: ', self.whichK)
        logPZ = self._compute_log_P_Z()
        # print("logPZ ", logPZ)
        logPG = self._compute_log_P_G()
        # print("logPG ", logPG)
        logPpi = 0
        if (self.Bishop_model_selection is False):
            logPpi = self._compute_log_P_pi()
        # print("logPpi ", logPpi)
        logPmu = self._compute_log_P_mu()
        # print("logPmu ", logPmu)
        logPeps = self._compute_log_P_epsilon()
        # print("logPeps ", logPeps)
        post1 = logl + logPZ + logPG + logPpi + logPmu + logPeps

        self.whichK = "all"
        print ('Which K: ', self.whichK)
        logPZ = self._compute_log_P_Z()
        print("logPZ ", logPZ)
        logPG = self._compute_log_P_G()
        print("logPG ", logPG)
        logPpi = 0
        if (self.Bishop_model_selection is False):
            logPpi = self._compute_log_P_pi()
        print("logPpi ", logPpi)
        logPmu = self._compute_log_P_mu()
        print("logPmu ", logPmu)
        logPeps = self._compute_log_P_epsilon()
        print("logPeps ", logPeps)
        post2 = logl + logPZ + logPG + logPpi + logPmu + logPeps

        return (post1, post2)

    ###############

    # @profile
    def _compute_log_likelihood(self):
        loglik = 0
        for data_type in self.data_types:
            loglik += self._compute_log_likelihood_helper(self.X[data_type], self.mu_star[data_type], self.gamma_star[data_type], self.N, self.pi_star, self.R[data_type], self.maxL[data_type], self.mean_or_mode)
        return loglik

    @staticmethod
    @njit(parallel=True, cache=True)
    def _compute_log_likelihood_helper(X, mu_star, gamma_star, N, pi_star, R, maxL, mean_or_mode):
        logl = 0
        # first add the likelihood of the methylation data
        # check which is the mode of Z and G
        # log P(X|Z,G,eps) = sum_n sum_r sum_l log  [gamma_star[g_krl,X_nrl]-1 /(gamma_star[g_krl,0]+gamma_star[g_krl,1]-2)]
        for n in prange(N):
            cluster = np.argmax(pi_star[n,])
            # print 'pi_star for n=', n, 'is ', pi_star[n,], ' cluster is ', cluster
            for r in prange(R):
                for l in prange(int(maxL)):
                    X_0 = X[n,r,l]
                    if not np.isnan(X_0):
                        epigeno = int(np.argmax(mu_star[:,cluster,r,l]))
                        # if (epigeno is not int(X_0)):
                        #     print ('IS NOT NA n=', n, ' r=', r, ' l=', l, 'X=', X_0, 'cluster=', cluster)
                        #     print ('Epigeno is ', epigeno)
                        #     print ('gamma_star is ', gamma_star[epigeno, int(X_0)])
                        if (mean_or_mode == "mean"):
                            logl += np.log(gamma_star[epigeno, int(X_0)] / (gamma_star[epigeno,0] + gamma_star[epigeno,1]))
                        #else:
                        #    logl = logl + np.log((gamma_star[epigeno, int(X[n,r,l])]-1) / (gamma_star[epigeno,0] + gamma_star[epigeno,1]-2))
        return logl

    ###############

    def _compute_log_P_Z(self):
    # For every cell, I look at this cluster assignment and add the log prob for that pi_star
        logp = 0
        all_clusters = []
        if (self.whichK is "all"):
            all_clusters = range(self.K)
        else:
            for n in range(self.N):
                cluster = np.argmax(self.pi_star[n,:])
                if cluster not in all_clusters:
                    all_clusters.append(cluster)
        # print('All clusters: ', *all_clusters)
        # print('Length ', len(all_clusters))

        alpha_sum = sum(self.alpha_star[i] for i in all_clusters)

        for data_type in self.data_types:
            for n in range(self.N):
                cluster = np.argmax(self.pi_star[n,:])
                #print(*self.log_pi_star[n,])
                #print("cluster ", cluster, " logpi* ", self.log_pi_star[n,cluster])
                # logp = logp + self.log_pi_star[n,cluster]
                if (self.mean_or_mode is "mean"):
                    if (self.Bishop_model_selection is False):
                        logp = logp + np.log(self.alpha_star[cluster]/alpha_sum)
                    else:
                        logp = logp + np.log(self.pi[cluster])
                #else:
                #    logp = logp + np.log((self.alpha_star[cluster]-1)/(sum(self.alpha_star[i] for i in all_clusters) - len(all_clusters)))
        return logp


    ###############

    def _compute_log_P_G(self):
        logp = 0
        all_clusters = []

        if (self.whichK is "all"):
            all_clusters = range(self.K)
        else:
            # Get the final clusters for each cell
            for n in prange(self.N):
                cluster = np.argmax(self.pi_star[n,:])
                if cluster not in all_clusters:
                   all_clusters.append(cluster)

        all_clusters = np.array(all_clusters)

        for data_type in self.data_types:
            log_mu_star = self.log_mu_star[data_type]
            mu_star = self.get_mu_star(data_type)
            beta_star = self.beta_star[data_type]
            maxL = self.maxL[data_type]
            R = self.R[data_type]

            logp += self._compute_log_P_G_helper(
                log_mu_star,
                mu_star,
                beta_star,
                maxL,
                R,
                self.mu_has_k,
                self.mean_or_mode,
                all_clusters
                )
        return logp


    @staticmethod
    @njit(parallel=True, cache=True)
    def _compute_log_P_G_helper(log_mu_star, mu_star, beta_star, maxL, R, mu_has_k, mean_or_mode, all_clusters):
        logp = 0
    # We use beta_star, not mu_star
    # Here include only the k's that were found - No, that doesn't make sense
        for k in all_clusters:
            for r in prange(R):
                for l in prange(int(maxL)):
                    bests = np.argmax(log_mu_star[:,k,r,l])
                    # print("k ", k)
                    # print(" beta star ", self.beta_star[data_type][k][r][bests])
                    # print(" argmax ", np.argmax(log_mu_star[:,k,r,l]))
                    if (mu_has_k):
                        if (mean_or_mode == "mean"):
                            logp += np.log(beta_star[k][r][bests]/np.sum(beta_star[k][r]))
                        #else:
                        #    logp = logp + np.log((self.beta_star[data_type][k][r][bests]-1)/(sum(self.beta_star[data_type][k][r])-2))
                    else:
                        if (mean_or_mode == "mean"):
                            logp += np.log(np.max(beta_star[r][bests])/np.sum(beta_star[r]))
                        #else:
                        #    logp = logp + np.log((self.beta_star[data_type][r][bests]-1)/(sum(self.beta_star[data_type][r])-2))
        return logp

    ###############

    # Note: If many of the G values are much less than 1 (e.g. 0.72), then logP(G) will be lower, and the log_posterior will be lower
    # That means that the log_posterior score favours the more "certain" results
    # For example, in the following 2 cases log_likelihood is higher for run 17,
    #   but log_posterior is higher for run 39 because logP(G) is much lower for run 17
    # /shahlab/mandronescu/EPI-98_run_synthetic_pipeline/OUTPUT_NREGIONS/RUN/D_NREGIONS_500_9_epiclomal_synthetic/outputs
    # [node0514 outputs]$ grep -i log ../logs/TASK_BASIC_EPICLOMAL__17.o2748085
    # Log likelihood:  2320317.60907
    # logPZ  -0.0424230776958
    # logPG  -22808.0376075
    # logPpi  12.8018274801
    # logPmu  0.0
    # logPeps  8.98091533122
    # Log posterior unnormalized:  2297531.31178
    # [node0514 outputs]$ grep -i log ../logs/TASK_BASIC_EPICLOMAL__39.o2747743
    # Log likelihood:  2319274.73179
    # logPZ  0.0
    # logPG  -3749.627868
    # logPpi  12.8018274801
    # logPmu  0.0
    # logPeps  8.72871796896
    # Log posterior unnormalized:  2315546.63447

    ###############

    def _compute_log_P_pi(self):
        # This is Dirichlet(alpha_star, alpha_zero)
        # Should we include only the found clusters? Or all clusters?
        # Including 2 situations

        logp = 0
        all_clusters = []
        if (self.whichK is "all"):
            all_clusters = range(self.K)
        else:
            for n in range(self.N):
                cluster = np.argmax(self.pi_star[n,:])
                if cluster not in all_clusters:
                    all_clusters.append(cluster)

        this_alpha_star =  [self.alpha_star[i]  for i in all_clusters]
        this_alpha_prior = [self.alpha_prior[i] for i in all_clusters]
        if (self.mean_or_mode is "mean"):
            logp = dirichlet.logpdf (x=this_alpha_star/sum(this_alpha_star), alpha=this_alpha_prior)
        #else:
        #    logp = dirichlet.logpdf (x=(this_alpha_star-np.array(1))/(sum(this_alpha_star)-len(all_clusters)), alpha=this_alpha_prior)
        return logp

    ###############

    def _compute_log_P_mu(self):
        # \sum_k \sum_r Dirichlet(beta_star, beta_zero)
        logp = 0

        for data_type in self.data_types:
            if (self.mu_has_k):
                # TODO: check this
                if (self.include_bulk[data_type]):
                    for k in range(self.K):
                        for r in range(self.R[data_type]):
                            for l in range(self.maxL[data_type]):
                                logp += dirichlet.logpdf(x=self.beta_star[data_type][k][r][l]/sum(self.beta_star[data_type][k][r][l]), alpha=self.beta_prior[data_type][r][l])
                         # This above computes log Dirichlet(beta_star, beta_zero)
                else:
                    for k in range(self.K):
                        for r in range(self.R[data_type]):
                            # print(*self.beta_star[data_type][k][r])
                            # print(*self.beta_prior[data_type])
                            if (self.mean_or_mode is "mean"):
                                logp += dirichlet.logpdf(x=self.beta_star[data_type][k][r]/sum(self.beta_star[data_type][k][r]), alpha=self.beta_prior[data_type])
                            #else:
                            #    logp += dirichlet.logpdf(x=self.beta_star[data_type][k][r]/sum(self.beta_star[data_type][k][r]), alpha=self.beta_prior[data_type])
            else:
                # TODO: check this
                if (self.include_bulk[data_type]):
                    for r in range(self.R[data_type]):
                        for l in range(self.maxL[data_type]):
                            logp += dirichlet.logpdf(x=self.beta_star[data_type][r][l]/sum(self.beta_star[data_type][r][l]), alpha=self.beta_prior[data_type][r][l])
                     # This above computes log Dirichlet(beta_star, beta_zero)
                else:
                    for r in range(self.R[data_type]):
                        # print(*self.beta_star[data_type][k][r])
                        # print(*self.beta_prior[data_type])
                        if (self.mean_or_mode is "mean"):
                            logp += dirichlet.logpdf(x=self.beta_star[data_type][r]/sum(self.beta_star[data_type][r]), alpha=self.beta_prior[data_type])
                        #else:
                        #    logp += dirichlet.logpdf(x=self.beta_star[data_type][k][r]/sum(self.beta_star[data_type][k][r]), alpha=self.beta_prior[data_type])

        return logp

    ###############

    def _compute_log_P_epsilon(self):
        logp = 0

        for data_type in self.data_types:
            # \sum \sum log Dirichlet(epsilon)
            logp += sum([dirichlet.logpdf(x=x/sum(x), alpha=y) for x, y in zip(self.gamma_star[data_type], self.gamma_prior[data_type])])

        return logp


    ###############
    ## Additional functions
    ###############

    def _compute_num_clusters(self):
        all_clusters = []
        # Get the final clusters for each cell
        for n in range(self.N):
            cluster = np.argmax(self.pi_star[n,:])
            if cluster not in all_clusters:
               all_clusters.append(cluster)

        return len(all_clusters)

    ###############

    def _compute_cluster_prevalence(self):
        # Note: this is a bit different from the actual prevalence obtained by counting all the cells in each cluster

        return np.array(self.pi_star.sum(axis=0)/self.N)


    ###############

    def _compute_num_nonzero_pi(self):

        return sum([self.pi[k] > 1e-5 for k in range(self.K)])

    ###############
    ## Additional functions for the SLS bulk
    ###############

    def _compute_different_regions(self, labels_pred, epigenotype):
        different_scores = []
        different_regions = []
        clusters = set(labels_pred)
        # print ("The unique predicted clusters: ", *clusters, " of type ", type(clusters))
        # print ("Epigenotype size: ", epigenotype.shape)
        for data_type in self.data_types:
            # print("There are ", len(self.regions[data_type]['start']), " regions")
            for r in range(len(self.regions[data_type]['start'])):
                # count how many epigenotypes are different in this region
                rscore = 0
                rlen = self.regions[data_type]['end'][r] - self.regions[data_type]['start'][r] + 1
                for cpg in range(self.regions[data_type]['start'][r], self.regions[data_type]['end'][r]+1):
                    # looking only at labels_pred, ignoring the other ones
                    # print ("        CPG ", cpg)
                    if (len(set(epigenotype[list(clusters),cpg])) != 1):       # different
                         rscore = rscore+1
                    # print("        Epigenotypes for this position :", *epigenotype[list(clusters),cpg], " rscore ", rscore)
                rscore = rscore/rlen
                # print("    Region ", r, " with length ", rlen, " score ", rscore)
                if (rscore > 0.5):
                    different_scores.append(rscore)
                    different_regions.append(r)

            #self.regions[data_type]['start']
        # print ("DIfferent regions ", *different_regions)
        return different_regions

    ###############

    def _compute_candidate_cells(self, labels_pred, epigenotype, different_regions):

        candidate_cells = []
        candidate_clusters = defaultdict(list)  # this is a dictionary of lists
        empty_regions = []
        clusters = set(labels_pred)
        for data_type in self.data_types:
            for cell in range(self.N):
                # print("Cell ", cell)
                cell_added = False
                for r in different_regions:
                    # print("    Region ", r, "of length ", int(self.L[data_type][r]))
                    region = self.X[data_type][cell,r,0:int(self.L[data_type][r])]
                    # print(*region)
                    if (all(np.isnan(v) for v in region)):
                        empty_regions.append(str(cell)+"_"+str(r))
                        if not cell_added:
                            candidate_cells.append(cell)
                            cell_added = True
                            # print("IS ALL NAN, adding cell ", cell, " with missing region ", r)
            # print("Cells with missing: ", candidate_cells)
            # print("Empty regions: ", empty_regions)

            c_cells = np.copy(candidate_cells)
            for cell in c_cells:
                candidate_clusters[cell].append(labels_pred[cell])
                similar_clusters = []
                all_clusters_are_candidate = True
                for r in different_regions:
                    if (str(cell)+"_"+str(r) not in empty_regions):
                    # check if this cell belongs to a unique epigenotype for region r
                    # If any non-empty region has unique epigenotype, then the cell is well placed and is not a candidate

                        all_clusters_are_candidate = False
                        rlen = self.regions[data_type]['end'][r] - self.regions[data_type]['start'][r] + 1
                        otherepi = set(labels_pred)
                        # similarity score with the other epigenotypes. The closer to 1 the more similar
                        similarity = defaultdict(int)
                        otherepi.remove(labels_pred[cell])
                        for s in otherepi:
                            similarity[s] = 0
                        # print("Similarity dict ", similarity)
                        for cpg in range(self.regions[data_type]['start'][r], self.regions[data_type]['end'][r]+1):
                            for sim in otherepi:
                                # print("Cell ", cell, " region ", r, " checking other epi ", sim, " my cluster ", labels_pred[cell])
                                if epigenotype[labels_pred[cell],cpg] == epigenotype[sim,cpg]:
                                    similarity[sim]+=1

                        for s in similarity:
                            similarity[s]/=rlen
                            if similarity[s] >= 0.8:
                               similar_clusters.append(s)
                        # print ("Cell ", cell, " region ", r, " similarity score ", similarity, " similar clusters ", *similar_clusters, " len ", len(similar_clusters))
                        if (len(similar_clusters) == 0):     # unique for this region
                            # print ("Removing cell ", cell)
                            candidate_cells.remove(cell)
                            candidate_clusters.pop(cell, None)
                            break
                        else:
                            for x in similar_clusters:
                                if x not in candidate_clusters[cell]:
                                    candidate_clusters[cell].append(x)
                if all_clusters_are_candidate:      # Meaning all the different regions for this cell are empty
                    for l in labels_pred:
                        if l not in candidate_clusters[cell]:
                            candidate_clusters[cell].append(l)


        for key in candidate_clusters:
            print ("Cell ", key, " with candidate clusters ", candidate_clusters[key])
        return candidate_clusters

    ###############

    def _get_relevant_bulk_percentages (self, different_regions):
        # select the different regions from the bulk
        meth_percentages = []
        for data_type in self.data_types:
            for r in different_regions:
                # print ("Region ", r, " from ", self.regions[data_type]['start'][r], " to ", self.regions[data_type]['end'][r])
                for cpg in range (self.regions[data_type]['start'][r], self.regions[data_type]['end'][r]+1):
                    # because in the bulk data it starts from 1
                    # print ("    Bulk in cpg ", cpg+1, " is meth ", self.slsbulk_data['meth_reads'][cpg+1], " and unmeth ", self.slsbulk_data['unmeth_reads'][cpg+1])
                    meth_percentages.append(self.slsbulk_data['meth_reads'][cpg+1]/(self.slsbulk_data['unmeth_reads'][cpg+1]+self.slsbulk_data['meth_reads'][cpg+1]))
        # print(*meth_percentages)
        return meth_percentages


    ###############

    def _get_predicted_percentages(self, labels_pred, epigenotype, different_regions):
    # NOTE: epigenotype is the one from the DIC selection criterion
    #   When I change the prediction
        # select the different regions from the bulk
        meth_percentages = []
        for data_type in self.data_types:
            for r in different_regions:
                # print ("Region ", r, " from ", self.regions[data_type]['start'][r], " to ", self.regions[data_type]['end'][r])
                for cpg in range(self.regions[data_type]['start'][r], self.regions[data_type]['end'][r]+1):
                    # look through all the cells
                    num_meth = 0
                    for cell in range(self.N):
                        # print ("EPigenotype for  cell ", cell)
                        # print (" with label ", labels_pred[cell])
                        # print (" cpg ", cpg)
                        # print (" is ", epigenotype[labels_pred[cell],cpg])
                        num_meth = num_meth + epigenotype[labels_pred[cell],cpg]        # this is 0 for unmethylated and 1 for methylated

                    meth_percentages.append(num_meth/self.N)
        # print("Predicted percentages ", *meth_percentages)
        return meth_percentages

    ###############

    def _get_bulk_score(self, labels_pred, epigenotype, different_regions, bulk_percentages):
        pred_percentages = self._get_predicted_percentages (labels_pred, epigenotype, different_regions)
        bulk_score = mean_absolute_error(bulk_percentages, pred_percentages)
        return bulk_score

    ###############

    def _slsbulk(self, candidate_clusters, labels_pred, epigenotype, different_regions, labels_true):
        # NOTE: Here I am using the true bulk values from labels_true -- BUT this is not correct because it is using the predicted epigenotype
        # Converting labels_true so I can  compare
        # for i in range(len(labels_true)):
        #     if labels_true[i] == 2:
        #         labels_true[i] = 0
        #     if labels_true[i] == 3:
        #         labels_true[i] = 2
        # print("Labels true after conversion")
        # print(*labels_true)
        # bulk_percentages = self._get_predicted_percentages (labels_true, epigenotype, different_regions)
        bulk_percentages = self._get_relevant_bulk_percentages (different_regions)
        # print("Bulk percentages from bulk data")
        # print (*bulk_percentages)
        best_bulk_score = self._get_bulk_score (labels_pred, epigenotype, different_regions, bulk_percentages)
        print("Bulk score at first ", best_bulk_score)
        # MA: I am not print the bulk score of the true any more because sometimes there are more true clusters than K
        # true_bulk_score = self._get_bulk_score (labels_true, epigenotype, different_regions, bulk_percentages)
        # print("Bulk score of the true ", true_bulk_score)
        best_pred = np.copy(labels_pred)
        current_pred = np.copy(labels_pred)

        candidate_keys =  list(candidate_clusters.keys())
        done = False
        for iteration in range(self.slsbulk_iterations):
            if done:
                break
            print("Iteration ", iteration, "\n--------------")

            for key in candidate_keys:
                # print ("Current prediction ", *current_pred)
                print ("    Cell ", key, " with candidate clusters ", candidate_clusters[key])
                for cc in candidate_clusters[key]:
                    if (cc == current_pred[key]):
                        continue
                    current_pred[key] = cc
                    # print ("        Changing cluster for cell ", key, " to ", cc)
                    new_bulk_score = self._get_bulk_score (current_pred, epigenotype, different_regions, bulk_percentages)
                    if (new_bulk_score < best_bulk_score):
                        best_bulk_score = new_bulk_score
                        best_pred = np.copy(current_pred)
                        print ("     Changing cluster for cell ", key, " to ", cc, ", New bulk score ", new_bulk_score, " KEEPING IT")
                    elif  (random() >= 0.8):
                        continue
                        # print ("             Best score ", best_bulk_score, ", Bulk score ", new_bulk_score, " staying with the current prediction anyway")
                    else:
                        current_pred = np.copy(best_pred)
                        # print ("             Best score ", best_bulk_score, ", New bulk score ", new_bulk_score, " THROWING IT")
                if (best_bulk_score == 0):      # this doesn't really happen because of sampling
                    print("Bulk is satisfied, exit!")
                    done = True
                    break
            shuffle(candidate_keys)
        print ("Final best bulk score: ", best_bulk_score)
        return best_pred



