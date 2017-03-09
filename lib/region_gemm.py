'''
Created on 2017-03-07

@author: Mirela Andronescu
'''

from __future__ import division

import numpy as np
import pandas as pd

from lib.utils import compute_e_log_dirichlet, compute_e_log_q_dirichlet, compute_e_log_p_dirichlet, \
                      compute_e_log_q_discrete, get_indicator_matrix, init_log_pi_star, log_space_normalise, safe_multiply

from lib.basic_gemm import BasicGeMM

class RegionGeMM(BasicGeMM):
    def __init__(self,
                 gamma_prior,
                 alpha_prior,
                 beta_prior,
                 X,
                 regions):

        self.R = {}
        self.L = {}
        
        for data_type in X.keys():        
            # find the number of regions and the size of each region
            self.R[data_type] = regions[data_type].shape[0]       # the number of regions   
            self.L[data_type] = np.zeros(self.R[data_type])
            for index, row in regions[data_type].iterrows():
                self.L[data_type][index] = row["end"] - row["start"] + 1
            
            print 'R, ', self.R[data_type]
            print 'L, ', self.L[data_type] 
                         
        BasicGeMM.__init__(self, gamma_prior, alpha_prior, beta_prior, X)   
        

    ######################    
    def _init_beta_star(self, data_type):
    # MA: this has to be a matrix of size KxRxS
        self.beta_star[data_type] = np.zeros((self.K, self.R[data_type], self.S[data_type]))
        for k in range(self.K):
            for r in range(self.R[data_type]):
                self.beta_star[data_type][k,r] = self.beta_prior[data_type].copy()    
                
        print 'BETA STAR initial\n', self.beta_star[data_type]                    
        
        
    ######################        
    def _get_beta_star_data_term(self, data_type):
        mu_star = self.get_mu_star(data_type)        # SxKxM

        # Using for loops intead of einsum
        # return np.einsum('skm -> ks', mu_star[:, :, :])   
          
        term = np.zeros((self.K, self.R[data_type], self.S[data_type]))
        for s in range(self.S[data_type]):
            for k in range(self.K):
                for r in range(self.R[data_type]):
                    for l in range(self.L[data_type][r]):
                        term[k,r,s] += mu_star[s,k,r,l]
        return term   
        
    ###################### 
    # override the IX[], same for mu_star[], keep the same data structure, just access the element IX[t,n,r,l] instead of IX[t,n,m] 
    def _get_IX (data_type, t, n, r, l):
        IX = self.IX[data_type]
        m = self.L[r]["start"]+l
        return IX [t, n, m]
    
        
        
        
        