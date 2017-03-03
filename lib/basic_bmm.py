'''
Created on 2015-02-03

@author: aroth
'''
from __future__ import division

import numpy as np

from scg.utils import compute_e_log_dirichlet, compute_e_log_q_dirichlet, compute_e_log_p_dirichlet, \
                      compute_e_log_q_discrete, get_indicator_matrix, init_Z, log_space_normalise, safe_multiply

class VariationalBayesDirichletMixtureModel(object):
    def __init__(self, gamma_prior, kappa_prior, X, init_labels=None):    
        self.K = len(kappa_prior)
        
        self.gamma_prior = gamma_prior
        
        self.kappa_prior = kappa_prior
            
        self.kappa = np.ones(self.K)
        
        self.data_types = X.keys()
    
        self.N = X[X.keys()[0]].shape[0]
        
        self.M = {}
        
        self.T = {}
        
        self.X = {}
        
        self.gamma = {}
        
        for data_type in self.data_types:
            if X[data_type].shape[0] != self.N:
                raise Exception('All data types must have the same number of rows (cells).')
                    
            self.M[data_type] = X[data_type].shape[1]
            
            self.T[data_type] = len(gamma_prior[data_type])
        
            self.X[data_type] = get_indicator_matrix(range(self.T[data_type]), X[data_type])

        self.log_Z = init_Z(self.K, self.N, init_labels)
        
        self.lower_bound = [float('-inf')]

        self._debug_lower_bound = [float('-inf')]
        
        self.converged = False

    def get_e_log_mu(self, data_type):
        e_log_mu = np.empty((self.T[data_type], self.K, self.M[data_type]))

        for k in range(self.K):
            for m in range(self.M[data_type]):
                e_log_mu[:, k, m] = compute_e_log_dirichlet(self.gamma[data_type][:, k, m])

        return e_log_mu

    @property
    def e_log_pi(self):
        return compute_e_log_dirichlet(self.kappa)
    
    @property
    def Z(self):
        return np.exp(self.log_Z)
    
    def fit(self, convergence_tolerance=1e-4, debug=False, num_iters=100):
        for i in range(num_iters):
            
            self._update_gamma()
            
            if debug:
                print 'gamma', self._diff_lower_bound()
            
            self._update_kappa()
            
            if debug:
                print 'kappa', self._diff_lower_bound()

            self._update_Z()
            
            if debug:
                print 'Z', self._diff_lower_bound()
            
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
                    
    ##########
   
    # q(mu)    
    def _update_gamma(self):
        for data_type in self.data_types:
            self.gamma[data_type] = self.gamma_prior[data_type][:, np.newaxis, np.newaxis] + self._get_gamma_data_term(data_type)
   
    def _get_gamma_data_term(self, data_type):
        X = self.X[data_type]
        
        # \sum_n E[I(Z_n=k)] I(Xnm=t)
        return np.einsum('tnkm, tnkm -> tkm',
                         self.Z[np.newaxis, :, :, np.newaxis],
                         X[:, :, np.newaxis, :])    
                         
    ##########
             
    # q(pi)                             
    def _update_kappa(self):
        self.kappa = self.kappa_prior + self._get_kappa_data_term()
    
    def _get_kappa_data_term(self):
        return self.Z.sum(axis=0)    
    
    ##########
        
    # q(Z)        
    def _update_Z(self):
        log_Z = self.e_log_pi[np.newaxis, :]
        
        for data_type in self.data_types:
            log_Z = log_Z + self._get_log_Z_d(data_type)
    
        self.log_Z = log_space_normalise(log_Z, axis=1)
    
    def _get_log_Z_d(self, data_type):
        e_log_mu = self.get_e_log_mu(data_type)
   
        X = self.X[data_type]
        
        # \sum_t \sum_m I(X_nm=t) E[log(mu_kmt)]
        log_Z = np.einsum('tnkm, tnkm -> nk',
                          e_log_mu[:, np.newaxis, :, :],
                          X[:, :, np.newaxis, :])
      
        return log_Z
    
    ##########    
    
    def _compute_lower_bound(self):
        return self._compute_e_log_p() - self._compute_e_log_q()
    
    def _compute_e_log_p(self):
        gamma_prior = 0
        
        gamma_posterior = 0
        
        for data_type in self.data_types:
            for k in range(self.K):
                for m in range(self.M[data_type]):
                    gamma_prior += compute_e_log_p_dirichlet(self.gamma[data_type][:, k, m],
                                                             self.gamma_prior[data_type])

            gamma_posterior += self._compute_e_log_p_gamma_posterior(data_type)
        
        kappa_prior = compute_e_log_p_dirichlet(self.kappa, self.kappa_prior)
        
        kappa_posterior = self._compute_e_log_p_kappa_posterior()
      
        return sum([gamma_prior,
                    gamma_posterior,
                    kappa_prior,
                    kappa_posterior])
    
    def _compute_e_log_p_kappa_posterior(self):
        return np.sum(safe_multiply(self.e_log_pi, self._get_kappa_data_term()))
    
    def _compute_e_log_p_gamma_posterior(self, data_type):
        return np.sum(safe_multiply(self.get_e_log_mu(data_type), self._get_gamma_data_term(data_type)))
    
    def _compute_e_log_q(self):
        log_q_mu = 0
        
        for data_type in self.data_types:
            for k in range(self.K):
                for m in range(self.M[data_type]):
                    log_q_mu += compute_e_log_q_dirichlet(self.gamma[data_type][:, k, m])

        log_q_pi = compute_e_log_q_dirichlet(self.kappa)
        
        log_q_z = compute_e_log_q_discrete(self.log_Z)

        return np.sum([log_q_mu,
                       log_q_pi,
                       log_q_z])
    
    
    def _get_kappa_data_term(self):
        return self.Z.sum(axis=0)
    
    def _diff_lower_bound(self):
        self._debug_lower_bound.append(self._compute_lower_bound())
        
        diff = (self._debug_lower_bound[-1] - self._debug_lower_bound[-2]) / np.abs(self._debug_lower_bound[-1])
        
        if diff < 0:
            print 'Bound decreased',
        
        return diff
    
if __name__ == '__main__':
    from sklearn.metrics import v_measure_score
    
    from simulate import get_default_dirichlet_mixture_sim, get_default_genotyper_sim
    
    kappa_prior = np.ones(20)
    
    gamma_prior = {'snv' : np.ones(3), 'breakpoint' : np.ones(2)}
    
    np.random.seed(0)
    
    sim = get_default_dirichlet_mixture_sim()

    v_dmm = []
     
    for i in range(10):
        np.random.seed(i)
        
        model = VariationalBayesDirichletMixtureModel(gamma_prior, kappa_prior, sim['X'])
     
        model.fit(num_iters=100)
     
        Z = model.Z.argmax(axis=1)
        
        v_dmm.append(v_measure_score(sim['Z'], Z))    
    
    np.random.seed(0)
    
    sim = get_default_genotyper_sim()

    v_gen = []
     
    for i in range(10):
        np.random.seed(i)
        
        model = VariationalBayesDirichletMixtureModel(gamma_prior, kappa_prior, sim['X'])
     
        model.fit(num_iters=100)
     
        Z = model.Z.argmax(axis=1)
        
        v_gen.append(v_measure_score(sim['Z'][0][sim['Y'] == 0], Z[sim['Y'] == 0]))
         
    print max(v_dmm), max(v_gen)
