'''
Created on 2017-01-17

@author: Mirela Andronescu
'''

from __future__ import division

import numpy as np
import pandas as pd

from lib.utils import compute_e_log_dirichlet, compute_e_log_q_dirichlet, compute_e_log_p_dirichlet, \
                      compute_e_log_q_discrete, get_indicator_matrix, init_log_EIZ, log_space_normalise, safe_multiply

class VariationalBayesBasicGeMMEpiclomal(object):
    def __init__(self,
                 gamma_prior,
                 alpha_prior,
                 beta_prior,
                 X):
        
        self.K = len(alpha_prior)               # max number of clusters
        
        print self.K, ' max number of clusters'
        
        self.N = X[X.keys()[0]].shape[0]        # number of cells
                
        self.gamma_prior = gamma_prior
        
        self.alpha_prior = alpha_prior
        
        self._init_Elogpi()  
        # MA: this just initializes with 1 for every clone. if 4 clones, it'll be [1 1 1 1]  
        
        self.beta_prior = beta_prior
        
        self.data_types = X.keys()
    
        self.M = {}         
        
        self.S = {}
        
        self.T = {}
        
        self.IX = {}
        
        self.Elogeps = {}
        
        self.Elogmu = {}        
        
        self.log_EIG = {}
        
        for data_type in self.data_types:
            if X[data_type].shape[0] != self.N:
                raise Exception('All data types must have the same number of rows (cells).')
                    
            self.M[data_type] = X[data_type].shape[1]
            # M is the number of loci
            
            self.S[data_type] = gamma_prior[data_type].shape[0]
            # S is 2 for 2 states
            
            self.T[data_type] = gamma_prior[data_type].shape[1]
            # T is 2 when gamma_prior is [99, 1; 1, 99]        
        
            self.IX[data_type] = get_indicator_matrix(range(self.T[data_type]), X[data_type])
            # self.IX is actually  the indicator matrix I(X_nm=t)
            # Now self.IX contains 2 matrices: 
            # the first is I(Xnm=0); this is the opposite of the input matrix
            #   (0 instead of 1 and 1 instead of 0, and with 0 instead of NaNs)
            # the second matrix is I(Xnm=1) and in this case (2 states) it is the same as the input matrix, and with 0 instead of NaNs
            
            self._init_Elogeps(data_type)
            # This just initializes Elogeps with gamma_prior            
            
            self._init_Elogmu(data_type)
            # This just initializes Elogmu for each cluster k with beta_prior, matrix KxS        
        
        self.log_EIZ = init_log_EIZ(self.K, self.N)
        # This function assigns random clusters to the Z variable. 
        # For each cell (row), one of the k values will be 1 and the others 0.
        # Then, return the log of this matrix (log_EIZ) of size N rows x K columns                                                              
        
        self.lower_bound = [float('-inf')]

        self._debug_lower_bound = [float('-inf')]
        
        self.converged = False
    
    def _init_Elogeps(self, data_type):
        self.Elogeps[data_type] = self.gamma_prior[data_type].copy()

    def _init_Elogmu(self, data_type):
    # MA: this has to be a matrix of size KxS
        self.Elogmu[data_type] = np.zeros((self.K, self.S[data_type]))
        for k in range(self.K):
            self.Elogmu[data_type][k] = self.beta_prior[data_type].copy()
            
    def _init_Elogpi(self):
        self.Elogpi = {}
        self.Elogpi = np.ones(self.K)
        
            
    def get_e_log_epsilon(self, data_type):
        # computes expectation(log(epsilon)) using the digamma function
        return compute_e_log_dirichlet(self.Elogeps[data_type])
        
    def get_e_log_mu(self, data_type):
        # computes expectation(log(mu)) using the digamma function
        # I tested this, and it returns K rows, as if I had a for loop
        return compute_e_log_dirichlet(self.Elogmu[data_type])                
        
    def get_e_log_pi(self):
        # computes expectation(log(pi)) using the digamma function
        return compute_e_log_dirichlet(self.Elogpi)
    
    def get_EIG(self, data_type):
        return np.exp(self.log_EIG[data_type])
    
    @property
    def EIG(self):
        EIG = {}
        
        for data_type in self.data_types:
            EIG[data_type] = self.get_EIG(data_type)
        
        return EIG    
    
    @property
    def EIZ(self):
        return np.exp(self.log_EIZ)
    
    def fit(self, convergence_tolerance=1e-4, debug=False, num_iters=100):
        print "Iter  ELBO difference"
        for i in range(num_iters):

            # update E(I(Gmk=s))
            self._update_EIG()
            
            if debug:
                print 'ELBO, diff after update_EIG', self._diff_lower_bound()
            
            # update E (log eps_st)
            self._update_Elogeps()
            
            if debug:
                print 'ELBO, diff after update_Elogeps', self._diff_lower_bound()
            
            # update E(log pi_k)
            self._update_Elogpi()
            
            if debug:
                print 'ELBO, diff after update_Elogpi', self._diff_lower_bound()

            # update E(log mu_k)
            self._update_Elogmu()
            
            if debug:
                print 'ELBO, diff after update_Elogmu', self._diff_lower_bound()
            
            # update E(I(Zn=k))
            self._update_EIZ()
            
            if debug:
                print 'ELBO, diff after update_EIZ', self._diff_lower_bound()
            
            
            # From Blei 2016: computing the ELBO for the full data set may be too expensive, we can compute on a held-out set
            self.lower_bound.append(self._compute_lower_bound())
             
            diff = (self.lower_bound[-1] - self.lower_bound[-2]) / np.abs(self.lower_bound[-1])
             
            print i, self.lower_bound[-1], diff
             
            if abs(diff) < convergence_tolerance:
                print 'Converged'
                
                self.converged = True
                
                break
            
            elif diff < 0:
                print 'Lower bound decreased'
                
                if not debug:
                    self.converged = False
                    
                    break
    
    ####################

    def _update_EIG(self):
        for data_type in self.data_types:
            self._update_EIG_d(data_type)
    
    def _update_EIG_d(self, data_type):        
        
        IX = self.IX[data_type]  # TxNxM

        EIZ = self.EIZ  # NxK
        
        # E[log(epsilon)] = E[log(Dirichlet(gamma))]
        # Now it's taking the e_log_epsilon
        e_log_epsilon = self.get_e_log_epsilon(data_type)   # SxT
        
        # print 'e_log_epsilon: {0}'.format(e_log_epsilon)
        
        # np.newaxis adds an additional unit-length dimension
        # epsilon has indexes st, and we want indexes stnkm
        e_log_epsilon = e_log_epsilon[:, :, np.newaxis, np.newaxis, np.newaxis]
        
        # the summations below mean
        # log q(G_km=s) = \sum_t \sum_n I(X_nm=t) * E[I(Z_n=k)] * E[log(epsilon_st)] + log(1/2)
        log_EIG = np.einsum('stnkm, stnkm, stnkm -> skm',
                          IX[np.newaxis, :, :, np.newaxis, :],       # IX depends on t (2 groups), n and m
                          EIZ[np.newaxis, np.newaxis, :, :, np.newaxis],  # EIZ depends on n and k
                          e_log_epsilon)        # e_log_epsilon depends on s and t
        
        # Now log_EIG has 3 dimensions: s (2 groups), k (total number of clusters), m (num loci)
                          
        # KxS
        # Here I am just transposing the matrix
        e_log_mu = np.einsum('ks -> sk',self.get_e_log_mu(data_type))
        
        # Here we just need a sum, not einsum
        log_EIG = e_log_mu[:, :, np.newaxis] + log_EIG
        
        # SxKxM
        self.log_EIG[data_type] = log_space_normalise(log_EIG, axis=0)
        
    ####################        
        
    def _update_Elogeps(self):
        for data_type in self.data_types:
            prior = self.gamma_prior[data_type]
            newterm = self._get_Elogeps_data_term(data_type)
            self.Elogeps[data_type] = prior + newterm
            # print 'Updated Elogeps is {0}'.format(self.Elogeps[data_type])
            
    def _get_Elogeps_data_term(self, data_type):
        log_EIG = self.log_EIG[data_type]
        
        IX = self.IX[data_type]
        
        # SxNxKxM
        # EIZ_nk + EIG_skm
        data_term = np.exp(self.log_EIZ[np.newaxis, :, :, np.newaxis] + log_EIG[:, np.newaxis, :, :])

        out_dim = 'st'
        
        # \sum_n \sum_k \sum_m E[I(Z_n=k)] E[I(G_mk=s)] I(X_nm=t)
        return np.einsum('stnkm, stnkm -> {0}'.format(out_dim),
                         data_term[:, np.newaxis, :, :, :],
                         IX[np.newaxis, :, :, np.newaxis, :])

    ####################        
        
    def _update_Elogmu(self):
        # MA: added this
        for data_type in self.data_types:
            prior = self.beta_prior[data_type]
            newterm = self._get_Elogmu_data_term(data_type)
            
            # prior is of size [1,S], newterm is of size [K,S]. I just have to do a simple sum, einsum doesn't do the right thing
            self.Elogmu[data_type] = prior + newterm

                    
    def _get_Elogmu_data_term(self, data_type):
        # MA: added this 
        EIG = self.get_EIG(data_type)        # SxKxM
        return np.einsum('skm -> ks', EIG[:, :, :])     
            
    ####################    
    
    def _update_Elogpi(self):
        self.Elogpi = self.alpha_prior + self._get_Elogpi_data_term()
        
    def _get_Elogpi_data_term(self):
        return self.EIZ.sum(axis=0)
                
    ####################                
        
    def _update_EIZ(self):
        log_EIZ = np.zeros((self.N, self.K))
        
        for data_type in self.data_types:
            log_EIZ = log_EIZ + self._get_log_EIZ_d(data_type)
            
        e_log_pi = self.get_e_log_pi()
                        
        log_EIZ = e_log_pi[np.newaxis, :] + log_EIZ
    
        self.log_EIZ = log_space_normalise(log_EIZ, axis=1)
    
    def _get_log_EIZ_d(self, data_type):
        EIG = self.get_EIG(data_type)
        
        IX = self.IX[data_type]
                  
        e_log_epsilon = self.get_e_log_epsilon(data_type)
        
        e_log_epsilon = e_log_epsilon[:, :, np.newaxis, np.newaxis, np.newaxis]
                
        # \sum_s \sum_t \sum_m E[I(G_mk=s)] I(X_nm=t) E[log(\epsilon_st)]
        log_EIZ = np.einsum('stnkm, stnkm, stnkm -> nk',
                          EIG[:, np.newaxis, np.newaxis, :, :],
                          IX[np.newaxis, :, :, np.newaxis, :],
                          e_log_epsilon)
        
        return log_EIZ    
    
    ##########################################################    
    
    def _compute_lower_bound(self):
        Elogp = self._compute_e_log_p()
        Elogq = self._compute_e_log_q()
        return Elogp - Elogq
    
    #################### 
    
    def _compute_e_log_p(self):
    # This is the expectation of log of joint for p
    # = ElogP(X|Z,G,epsilon,pi,mu) + ElogP(Z|pi) + ElogP(G|mu) + ElogP(epsilon) + ElogP(pi) + ElogP(mu)
    # = term 1                     + term 2      + term 3      + term 4         + term 5    + term 6
     
        log_p_Elogeps = self._compute_log_p_Elogeps()   # term 4 + term 1
        
        log_p_Elogpi = self._compute_log_p_Elogpi()     # term 2 + term 5

        log_p_Elogmu = self._compute_log_p_Elogmu()     # term 3 + term 6          
        
        return sum([log_p_Elogeps,
                    log_p_Elogpi,
                    log_p_Elogmu])
                    
    #################                    
                    
    def _compute_log_p_Elogmu(self):
        
        log_p_prior = 0        
        log_p_posterior = 0
        for data_type in self.data_types:
            #log_p_prior += sum([compute_e_log_p_dirichlet(x, y) for x, y in zip(self.Elogeps[data_type], self.gamma_prior[data_type])])        
        
            for k in range(self.K):
                log_p_prior += compute_e_log_p_dirichlet(self.Elogmu[data_type][k], self.beta_prior[data_type])   # term 6
                 # This above computes E log Dirichlet(beta_star, beta_zero)
            log_p_posterior += self._compute_e_log_p_Elogmu_posterior(data_type)  # term 3                
        
        return log_p_prior + log_p_posterior                                        
        
    def _compute_e_log_p_Elogmu_posterior(self, data_type):
        # term 3 in the joint
        e_log_p = 0
        
        EIG = self._get_Elogmu_data_term(data_type)  # size KxS, \sum_m E(I(Gkm=s))
        Elogmu = self.get_e_log_mu(data_type)  # size KxS

        # This below is \sum_m E(I(Gkm=s))Elog mu_ks     (term 3 in the joint)
        e_log_p = np.einsum('ks,ks', EIG, Elogmu)

        return e_log_p        

    #################                    
    
    def _compute_log_p_Elogeps(self):
        log_p_prior = 0
        
        log_p_posterior = 0
        
        for data_type in self.data_types:
            # term 4 \sum \sum log Dirichlet(epsilon)
            log_p_prior += sum([compute_e_log_p_dirichlet(x, y) for x, y in zip(self.Elogeps[data_type], self.gamma_prior[data_type])])
    
            log_p_posterior += self._compute_e_log_p_Elogeps_posterior(data_type)
        
        return log_p_prior + log_p_posterior

    def _compute_e_log_p_Elogeps_posterior(self, data_type):
        # This computes term 1 
        
        # Elog epsilon
        e_log_epsilon = self.get_e_log_epsilon(data_type)
        
        # below is sum EIZ EIG IX
        data_term = self._get_Elogeps_data_term(data_type)
        
        log_p = np.sum(safe_multiply(data_term, e_log_epsilon))
        
        return log_p

    ################# 
    def _compute_log_p_Elogpi(self):
    
        log_p_prior = compute_e_log_p_dirichlet(self.Elogpi, self.alpha_prior)
        # This above computes E log Dirichlet(alpha_star, alpha_zero)
        
        log_p_posterior = self._compute_e_log_p_Elogpi_posterior()
        
        return log_p_prior + log_p_posterior
            
    def _compute_e_log_p_Elogpi_posterior(self):
        e_log_p = 0
        
        # This below is \sum_n E(I(Zn=k))Elog pi _k     (term 2 in the joint)
        e_log_p += np.sum(safe_multiply(self.get_e_log_pi(), self._get_Elogpi_data_term()))
        
        return e_log_p    
            
    ####################                         
    def _compute_e_log_q(self):
        log_q_epsilon = self._compute_log_q_epsilon()
        
        log_q_pi = compute_e_log_q_dirichlet(self.Elogpi)
        
        log_q_mu = self._compute_log_q_mu()
        
        log_q_G = 0
        
        for data_type in self.data_types:
            log_q_G += compute_e_log_q_discrete(self.log_EIG[data_type])
        
        log_q_Z = compute_e_log_q_discrete(self.log_EIZ)

        return np.sum([log_q_epsilon,
                       log_q_pi,
                       log_q_mu,
                       log_q_G,
                       log_q_Z])
    
    def _compute_log_q_epsilon(self):
        log_q = 0
        
        for data_type in self.data_types:
            log_q += sum([compute_e_log_q_dirichlet(x) for x in self.Elogeps[data_type]])
        
        return log_q

    def _compute_log_q_mu(self):
        log_q = 0
        
        for data_type in self.data_types:
            log_q += sum([compute_e_log_q_dirichlet(x) for x in self.Elogmu[data_type]])
        
        return log_q

    ###############    
        
    def _diff_lower_bound(self):
        ELBO = self._compute_lower_bound()
        
        self._debug_lower_bound.append(ELBO)
        
        diff = (self._debug_lower_bound[-1] - self._debug_lower_bound[-2]) / np.abs(self._debug_lower_bound[-1])
        
        if diff < 0:
            print 'Bound decreased',
        
        return (ELBO, diff)
            
if __name__ == '__main__':
    def test_run():
        np.random.seed(0)
    
        model = VariationalBayesBasicEpicaller(gamma_prior, 
                                                 alpha_prior, 
                                                 beta_prior, 
                                                 sim['X'])
 
        model.fit(num_iters=100)
        
        print model.Elogeps['snv'].shape, model.Elogpi.keys()
        print "Results: "
        print model.Z
        
    # We can run the model in some default way or not worry about this main function


