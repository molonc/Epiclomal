'''
Created on 2017-03-07

@author: Mirela Andronescu
'''


# Maybe useful ideas on how to create an efficient manual einsum with loops
# http://stackoverflow.com/questions/35500925/is-numpy-einsum-efficient-compared-to-fortran-or-c

from __future__ import division

import numpy as np
import pandas as pd

from lib.utils import compute_e_log_dirichlet, compute_e_log_q_dirichlet, compute_e_log_p_dirichlet, \
                      compute_e_log_q_discrete, init_log_pi_star, log_space_normalise, safe_multiply

from lib.basic_gemm import BasicGeMM

class RegionGeMM(BasicGeMM):
    def __init__(self,
                 gamma_prior,
                 alpha_prior,
                 beta_prior,
                 X,
                 regions,
                 initial_clusters_data,
                 mu_has_k,
                 Bishop_model_selection,
                 slsbulk_data=None):
                 
        self.Rstart = {}
        self.Rend = {}
        self.L = {}          
                 
        BasicGeMM.__init__(self, gamma_prior, alpha_prior, beta_prior, X, regions, initial_clusters_data, mu_has_k, Bishop_model_selection, slsbulk_data)   
                         
    ###################### 
    
    def _set_region_params(self, data_type, regions):
    # find the number of regions and the size of each region
        self.R[data_type] = regions[data_type].shape[0]       # the number of regions   
        self.L[data_type] = np.zeros(self.R[data_type])
        self.Rstart[data_type] = np.zeros(self.R[data_type])
        self.Rend[data_type] = np.zeros(self.R[data_type])
        for index, row in regions[data_type].iterrows():
            self.Rstart[data_type][index] = row["start"]
            self.Rend[data_type][index] = row["end"]
            self.L[data_type][index] = row["end"] - row["start"] + 1
        self.maxL[data_type] = np.max(self.L[data_type])
        print ('R, ', self.R[data_type])
        # print 'L, ', self.L[data_type] 
        # print 'Rstart ,', self.Rstart[data_type]        
        # print 'Rend ,', self.Rend[data_type]        
        print ('maxL, ', self.maxL[data_type])
       
    ######################       
          
    def _region_data_matrix(self, data_type, X):
    # Here in the regions class I will fill up to 0 the remaining values for each region
    # TO DO: maybe this reshaping can be done more efficiently 
        matrix = np.empty((self.N, self.R[data_type], int(self.maxL[data_type])))
        matrix[:] = np.NAN
        for n in range(self.N):
            for r in range(self.R[data_type]):
                for l in range(int(self.Rstart[data_type][r]), int(self.Rend[data_type][r]+1)):
                    matrix[n,r, int(l - self.Rstart[data_type][r])] = X.values[n,l]                  
        return matrix
        
    ###################### 
            
    def _unregion_data_matrix(self, data_type, X):        
    # TO DO
    # this is necessary when we want to compare the imputed values
        return


    ######################         
        
    def unregion_mu_star(self, data_type):
    # transform data back from regioned (index rl) to unregioned (index m)        
        mu_star = self.get_mu_star(data_type)
        states = range(self.S[data_type])
        unregioned_mu_star = np.zeros((len(states), self.K, self.M[data_type]))
        
        # this is probably very slow, there may be a more efficient solution
        for s in range(len(states)):
            for k in range (self.K):
                location = 0
                for r in range(self.R[data_type]):
                    length = int(self.L[data_type][r])                    
                    unregioned_mu_star[s,k,location:location+length] = mu_star[s,k,r,range(length)]
                    location = location + length
                    
        return unregioned_mu_star                            
         
