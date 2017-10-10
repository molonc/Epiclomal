'''
Created on 2017-01-17

@author: Mirela Andronescu
'''

# Very nice blog about einsum and how to convert it into nested sums
# https://obilaniu6266h16.wordpress.com/2016/02/04/einstein-summation-in-numpy/

from __future__ import division

import numpy as np
import pandas as pd

from lib.utils import compute_e_log_dirichlet, compute_e_log_q_dirichlet, compute_e_log_p_dirichlet, \
                      compute_e_log_q_discrete, init_log_pi_star, log_space_normalise, safe_multiply

class BasicGeMM(object):
    def __init__(self,
                 gamma_prior,
                 alpha_prior,
                 beta_prior,
                 X,
                 regions,
                 initial_clusters_data):
        
        self.K = len(alpha_prior)               # max number of clusters
        
        print 'Max number of clusters: ', self.K
        
        self.N = X[X.keys()[0]].shape[0]        # number of cells
                
        self.gamma_prior = gamma_prior
        
        self.alpha_prior = alpha_prior
        
        self._init_alpha_star()  
        # MA: this just initializes with 1 for every clone. if 4 clones, it'll be [1 1 1 1]  
                
        self.data_types = X.keys()
    
        # I am actually implementing it as a one region, that way we can extend to regions more easily
        self.M = {}           
        self.R = {}      
        self.maxL = {}    # maxL = max_r(l)  = the maximum region size
        
        self.S = {}
        
        self.T = {}
        
        self.IX = {}
        
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
        
            matrix = self._region_data_matrix(data_type, X[data_type])        
        
            self.IX[data_type] = self._get_indicator_X(data_type, matrix)
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
        
        self.lower_bound = [float('-inf')]

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
        for i, s in enumerate(range(self.S[data_type])):
            Y[i, :, :, :] = (X == s).astype(int)
    
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
            
            # print self.beta_prior[data_type]
            

    def _init_beta_star(self, data_type):
        if (self.include_bulk[data_type]):
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

            
    def _init_alpha_star(self):
        self.alpha_star = {}
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
    
    @property
    def pi_star(self):
        return np.exp(self.log_pi_star)
    
    def fit(self, convergence_tolerance=1e-4, debug=False, num_iters=100):
        print "Iter  ELBO difference"

        # print 'pi_star for cell 14 initially: ', self.pi_star[13,]  
      
        for i in range(num_iters):

            # update alpha_star
            self._update_alpha_star()
            
            if debug:
                print 'ELBO, diff after update_alpha_star', self._diff_lower_bound()

            # update E(I(Gmk=s))
            self._update_mu_star()            
            
            if debug:
                print 'ELBO, diff after update_mu_star', self._diff_lower_bound()
            
            # update beta_star
            self._update_beta_star()
            
            if debug:
                print 'ELBO, diff after update_beta_star', self._diff_lower_bound()            
            
            # update gamma_star
            self._update_gamma_star()
            
            if debug:
                print 'ELBO, diff after update_gamma_star', self._diff_lower_bound()
                        
            # update pi_star
            self._update_pi_star()
            
            if debug:
                print 'ELBO, diff after update_pi_star', self._diff_lower_bound()
                     
            # print 'pi_star for cell 14 after pi_star update: ', self.pi_star[13,]                
                
            
            # update rho_star, but here in the basic model nothing happens
            self._update_rho_star()
            if debug:
                print 'ELBO, diff after update_rho_star', self._diff_lower_bound()            
            
            
            # From Blei 2016: computing the ELBO for the full data set may be too expensive, we can compute on a held-out set
            self.lower_bound.append(self._compute_lower_bound())
             
            diff = (self.lower_bound[-1] - self.lower_bound[-2]) / np.abs(self.lower_bound[-1])
             
            print i, self.lower_bound[-1], diff
             
            if abs(diff) < convergence_tolerance:
                print 'Converged'
                
                self.converged = True
                
                self._unregion_rho_star()
                
                break
            
            elif diff < 0:
                print 'Lower bound decreased'
                
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

        log_mu_star = self._update_mu_star_big_sum (data_type)
               
        # Now log_mu_star has 3 dimensions: s (2 groups), k (total number of clusters), m (num loci)
        # With 1 region, log_mu_star has 4 dimensions: s, k, r, maxl
                         
        # KxS
        # Here I am just transposing the matrix
        if (self.include_bulk[data_type]):
            e_log_mu = np.einsum('krls -> skrl',self.get_e_log_mu(data_type))   
            # Here we just need a sum, not einsum
            log_mu_star = e_log_mu[:, :, :, :] + log_mu_star                
        else:
            e_log_mu = np.einsum('krs -> skr',self.get_e_log_mu(data_type))       
            log_mu_star = e_log_mu[:, :, :, np.newaxis] + log_mu_star      
       

       
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
        data_term = np.exp(self.log_pi_star[np.newaxis, :, :, np.newaxis, np.newaxis] + log_mu_star[:, np.newaxis, :, :, :])
        
        
        # \sum_n \sum_k \sum_m E[I(Z_n=k)] E[I(G_mk=s)] I(X_nm=t)
        # return np.einsum('stnkm, stnkm -> st',
        #                  data_term[:, np.newaxis, :, :, :],
        #                  IX[np.newaxis, :, :, np.newaxis, :])
        return np.einsum('stnkrl, stnkrl -> st',
                         data_term[:, np.newaxis, :, :, :, :],
                         IX[np.newaxis, :, :, np.newaxis, :, :])                         

    ####################        
        
    def _update_beta_star(self):
        # MA: added this
        for data_type in self.data_types:
            prior = self.beta_prior[data_type]                                        
                    
            newterm = self._get_beta_star_data_term(data_type)
            
            # # prior is of size [1,S], newterm is of size [K,S]. I just have to do a simple sum, einsum doesn't do the right thing
            # With bulk, prior is of size [R,maxL,S] and newterm is of size [S,K,R,maxL]
            # self.beta_star[data_type] = np.einsum('krls, krls', prior[np.newaxis,:,:,:] + newterm[:,:,:,:])
            
            if (self.include_bulk[data_type]):
                self.beta_star[data_type] = prior[np.newaxis, :, :, :] + newterm
                # Above and below seem to be the same
                #for k in range(self.K):
                #    self.beta_star[data_type][k] = prior + newterm[k]
            else:
                self.beta_star[data_type] = prior + newterm                    

                    
    def _get_beta_star_data_term(self, data_type):
        # MA: added this 
        mu_star = self.get_mu_star(data_type)        # SxKxR or SxKxRxmaxL
        
        # Now beta_star goes over l too
        if (self.include_bulk[data_type]):
            return np.einsum('skrl -> krls', mu_star[:, :, :, :]) 
        else:
            return np.einsum('skrl -> krs', mu_star[:, :, :, :]) 
                  
                    
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
            
        e_log_pi = self.get_e_log_pi()
                        
        log_pi_star = e_log_pi[np.newaxis, :] + log_pi_star
    
        self.log_pi_star = log_space_normalise(log_pi_star, axis=1)
    
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
    
    def _compute_lower_bound(self):
        Elogp = self._compute_e_log_p()
        Elogq = self._compute_e_log_q()
        return Elogp - Elogq
    
    #################### 
    
    def _compute_e_log_p(self):
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
        e_log_p_term5 = self._compute_e_log_p_term5()   # term 5
        e_log_p_term6 = self._compute_e_log_p_term6()   # term 6
        e_log_p_term7 = self._compute_e_log_p_term7()   # term 7
                    
        e_log_p = sum([e_log_p_term1,
                    e_log_p_term2,
                    e_log_p_term3,                    
                    e_log_p_term4,
                    e_log_p_term5,
                    e_log_p_term6,
                    e_log_p_term7])

        return e_log_p                     
                    
                    
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
        return np.sum(safe_multiply(self.get_e_log_pi(), self._get_alpha_star_data_term()))
                 
    #################
        
    def _compute_e_log_p_term3(self):
        # term 3 in the joint
        # \sum_k \sum_m \sum_s mu*_kms E_log_mu_kms
        e_log_p = 0
        
        for data_type in self.data_types:                    
            mu_star = self._get_beta_star_data_term(data_type)  # size KxS, \sum_m E(I(Gkm=s)), now size KxRxS, with bulk size KxRxmaxLxS
            e_log_mu = self.get_e_log_mu(data_type)  # size KxS, now KxRxmaxLxS

            # This below is \sum_m E(I(Gkm=s))Elog mu_ks     (term 3 in the joint)
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
                e_log_p += np.einsum('krs,krs', mu_star, e_log_mu)    

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
    
        return  compute_e_log_p_dirichlet(self.alpha_star, self.alpha_prior)
        # This above computes E log Dirichlet(alpha_star, alpha_zero)
                  
    #################                   
                    
    def _compute_e_log_p_term6(self):
        # sum_k sum_m E log P(mu_km) = sum_k sum_m E log Dirichlet(beta_star, beta_zero)
        e_log_p = 0        
        
        for data_type in self.data_types:
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
        return e_log_p                                       
                    
    ####################                    
                    
    def _compute_e_log_p_term7(self):
        # This computes term 7
        return 0 
                                                    
    ####################
    def _compute_e_log_q(self):
        # For Basic-GeMM:
        # E(log q) = E(log(eps)) + E(log(pi)) + E(log(mu)) + E(log(G)) + E(log(Z))
        # For Basic-Miss-GeMM:
        # E(log q) = E(log(eps)) + E(log(pi)) + E(log(mu)) + E(log(G)) + E(log(Z)) + E(log(Xmiss))
        
        e_log_q_epsilon = self._compute_e_log_q_epsilon()
        
        # NOTE: in Genotyper with samples it goes over the unique samples only!!!
        e_log_q_pi = compute_e_log_q_dirichlet(self.alpha_star)
        
        e_log_q_mu = self._compute_e_log_q_mu()
        
        e_log_q_G = 0
        
        e_log_q_Xmiss = 0
        
        for data_type in self.data_types:
            e_log_q_G += compute_e_log_q_discrete(self.log_mu_star[data_type])
            e_log_q_Xmiss += self._compute_e_log_q_Xmiss(data_type)
        
        e_log_q_Z = compute_e_log_q_discrete(self.log_pi_star)                
                       
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
        
        return sum

    ###############    
        
    def _diff_lower_bound(self):
        ELBO = self._compute_lower_bound()
        
        self._debug_lower_bound.append(ELBO)
        
        diff = (self._debug_lower_bound[-1] - self._debug_lower_bound[-2]) / np.abs(self._debug_lower_bound[-1])
        
        if diff < 0:
            print 'Bound decreased',
        
        return (ELBO, diff)
            


