'''
Created on 2017-02-23

@author: Mirela Andronescu
'''

from __future__ import division

import numpy as np
import pandas as pd

from lib.utils import compute_e_log_dirichlet, compute_e_log_q_dirichlet, compute_e_log_p_dirichlet, \
                      compute_e_log_q_discrete, init_log_pi_star, log_space_normalise, safe_multiply

from lib.basic_gemm import BasicGeMM

class BasicMissGeMM(BasicGeMM):
    def __init__(self,
                 gamma_prior,
                 alpha_prior,
                 beta_prior,
                 X,
                 regions):
                         
        BasicGeMM.__init__(self, gamma_prior, alpha_prior, beta_prior, X, regions)               

        self.log_rho_star = {}      
        # we keep just the log for numerical reasons, not rho_star
        
        for data_type in self.data_types:  
            matrix = self._region_data_matrix(data_type, X[data_type])      
            self.log_rho_star[data_type] = self._init_log_rho_star (data_type, matrix)
            # print 'rho_star ', self.get_rho_star(data_type)
        
    ###################        
    def _init_log_rho_star (self, data_type, X):
    # Initialize all the nan values with 1/|S|
        states = range(self.T[data_type])
        # rho_star = np.zeros((len(states), self.N, self.M[data_type]))        
        rho_star = np.zeros((len(states), self.N, self.R[data_type], self.maxL[data_type]))

        # initialize with 0.5 for each missing position, 0 otherwise
        for s in range(len(states)):
            rho_star[s, :, :, :] = np.isnan(X)/len(states)            
    
        log_rho_star = np.log(rho_star + 1e-10)       # adding 1e-10 to avoid undetermined log(0)
    
        # Now make sure they sum up to 1 over columns
        return log_space_normalise(log_rho_star, axis=1)
        
    ###################  
            
    def get_rho_star(self, data_type):
        return np.exp(self.log_rho_star[data_type])    
        
    def unregion_rho_star(self, data_type):
    # transform data back from regioned (index rl) to unregioned (index m)
    # this will be done differently in the regions class    
        rho_star = self.get_rho_star(data_type)
        states = range(self.T[data_type])
        return np.argmax(rho_star.reshape((len(states), self.N, self.M[data_type])), axis=0)
        
    ###################             
    def _update_rho_star(self):
        for data_type in self.data_types:
            self._update_rho_star_d(data_type)
        
    
    def _update_rho_star_d(self, data_type):                                
        # Now it's taking the e_log_epsilon
        e_log_epsilon = self.get_e_log_epsilon(data_type)   # SxT   
        
        mu_star = self.get_mu_star(data_type)                 

        # the summations below mean
        # log q(G_km=s) = \sum_k \sum_s mu*_kms pi*_nk Elogeps_st
        # log_rho_star = np.einsum('stnkm, stnkm, stnkm -> tnm',
        #                   mu_star[:, np.newaxis, np.newaxis, :, :],
        #                   self.pi_star[np.newaxis, np.newaxis, :, :, np.newaxis],  # pi_star depends on n and k
        #                   e_log_epsilon[:, :, np.newaxis, np.newaxis, np.newaxis])   # e_log_epsilon depends on s and t                          
               
        # Now turning m into rl               
        log_rho_star = np.einsum('stnkrl, stnkrl, stnkrl -> tnrl',
                          mu_star[:, np.newaxis, np.newaxis, :, :, :],
                          self.pi_star[np.newaxis, np.newaxis, :, :, np.newaxis, np.newaxis],  # pi_star depends on n and k
                          e_log_epsilon[:, :, np.newaxis, np.newaxis, np.newaxis, np.newaxis])   # e_log_epsilon depends on s and t               
                          
                                  
        # SxKxM, now SxKxRxmaxL
        self.log_rho_star[data_type] = log_space_normalise(log_rho_star, axis=0)
        
        
                    
# NOTE: If I don't have _update_mu_star here, then the function called from _update_mu_star 
# named _update_mu_star_d of the parent class is the function in this class if it exists                   
    
    ###################    
    def _update_mu_star_big_sum(self, data_type):    
        # This function is different in basic_miss_gemm
         
        IX = self.IX[data_type]  # TxNxM, now TxNxRxmaxL
        
        # Now it's taking the e_log_epsilon
        e_log_epsilon = self.get_e_log_epsilon(data_type)   # SxT            
         
        rho_star = self.get_rho_star(data_type)
        # IX depends on t (2 groups), n and m (now rl), same for rho_star
        X_sum = IX[np.newaxis, :, :, np.newaxis, :, :] + rho_star[np.newaxis, :, :, np.newaxis, :, :]
        
        # the summations below mean
        # log q(G_km=s) = \sum_t \sum_n [I(X_nm=t) + rho_star_nmt] * E[I(Z_n=k)] * E[log(epsilon_st)] + log(1/2)
        # return np.einsum('stnkm, stnkm, stnkm -> skm',
        #                   X_sum,      
        #                   self.pi_star[np.newaxis, np.newaxis, :, :, np.newaxis],  # pi_star depends on n and k
        #                  e_log_epsilon[:, :, np.newaxis, np.newaxis, np.newaxis])   # e_log_epsilon depends on s and t
         
        # Now turning m into rl         
        return np.einsum('stnkrl, stnkrl, stnkrl -> skrl',
                          X_sum,      
                          self.pi_star[np.newaxis, np.newaxis, :, :, np.newaxis, np.newaxis],  # pi_star depends on n and k
                          e_log_epsilon[:, :, np.newaxis, np.newaxis, np.newaxis, np.newaxis])   # e_log_epsilon depends on s and t
                                           
    ###################      
    def _get_gamma_star_data_term(self, data_type):
        log_mu_star = self.log_mu_star[data_type]
        
        IX = self.IX[data_type]
        
        # SxNxKxM, now SxNxKxRxL
        # pi_star_nk + mu_star_skm
        data_term = np.exp(self.log_pi_star[np.newaxis, :, :, np.newaxis, np.newaxis] + log_mu_star[:, np.newaxis, :, :, :])
        
        rho_star = self.get_rho_star(data_type)
        X_sum = IX[np.newaxis, :, :, np.newaxis, :, :] + rho_star[np.newaxis, :, :, np.newaxis, :, :]              
        
        # \sum_n \sum_k \sum_m E[I(Z_n=k)] E[I(G_mk=s)] I(X_nm=t)
        # return np.einsum('stnkm, stnkm -> st',
        #                  data_term[:, np.newaxis, :, :, :],
        #                  X_sum)
                         
        return np.einsum('stnkrl, stnkrl -> st',
                         data_term[:, np.newaxis, :, :, :, :],
                         X_sum)                         
      
    ###################      
    def _get_log_pi_star_d(self, data_type):
        mu_star = self.get_mu_star(data_type)
        
        IX = self.IX[data_type]
                  
        e_log_epsilon = self.get_e_log_epsilon(data_type)

        rho_star = self.get_rho_star(data_type) 
        X_sum = IX[np.newaxis, :, :, np.newaxis, :, :] + rho_star[np.newaxis, :, :, np.newaxis, :, :]
                
        # \sum_s \sum_t \sum_m E[I(G_mk=s)] I(X_nm=t) E[log(\epsilon_st)]
        # log_pi_star = np.einsum('stnkm, stnkm, stnkm -> nk',
        #                   mu_star[:, np.newaxis, np.newaxis, :, :],
        #                   X_sum,
        #                   e_log_epsilon[:, :, np.newaxis, np.newaxis, np.newaxis])
                          
        log_pi_star = np.einsum('stnkrl, stnkrl, stnkrl -> nk',
                          mu_star[:, np.newaxis, np.newaxis, :, :, :],
                          X_sum,
                          e_log_epsilon[:, :, np.newaxis, np.newaxis, np.newaxis, np.newaxis])
        
        return log_pi_star        
        
        
        
    ### compute the ELBO       
    
    #################

    def _compute_e_log_p_term7(self):
        # This computes term 7 of the joint ElogP(Xmiss|Z,G,epsilon,pi,mu) 
        
        ### To test!!
        
        e_log_p = 0
        for data_type in self.data_types:        
            # Elog epsilon
            e_log_epsilon = self.get_e_log_epsilon(data_type)
                
            log_mu_star = self.log_mu_star[data_type]        
        
            # SxNxKxM, now SxNxKxRxmaxL
            # pi_star_nk + mu_star_skm
            data_term = np.exp(self.log_pi_star[np.newaxis, :, :, np.newaxis, np.newaxis] + log_mu_star[:, np.newaxis, :, :, :])
        
            rho_star = self.get_rho_star(data_type)
            X_sum = rho_star[np.newaxis, :, :, np.newaxis, :, :]              
        
            # \sum_n \sum_k \sum_m E[I(Z_n=k)] E[I(G_mk=s)] E(I(Xmiss_nm=t))
            # data_term = np.einsum('stnkm, stnkm -> st',
            #                      data_term[:, np.newaxis, :, :, :],
            #                      X_sum)
                                 
            data_term = np.einsum('stnkrl, stnkrl -> st',
                                 data_term[:, np.newaxis, :, :, :, :],
                                 X_sum)                                 
        
            e_log_p += np.sum(safe_multiply(data_term, e_log_epsilon))
        
        return e_log_p    
    
    #################    
    
    def _compute_e_log_q_Xmiss(self, data_type):
        return compute_e_log_q_discrete(self.log_rho_star[data_type])
        
        
        
        