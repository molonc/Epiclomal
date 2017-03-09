'''
Created on 2017-01-17

@author: Mirela Andronescu
'''
from __future__ import division

import gzip
import numpy as np
import os
import pandas as pd
import yaml
import time

from lib.basic_gemm import BasicGeMM
from lib.basic_miss_gemm import BasicMissGeMM
from lib.region_gemm import RegionGeMM
from lib.utils import load_labels

##############################

def run_basic_gemm_model(args):

    t0 = time.clock()
    if args.seed is not None:
        np.random.seed(args.seed)
       
    cell_ids, data, event_ids, priors, regions = load_data(args.config_file, args.data_file)
   
    model = BasicGeMM(priors['gamma'],
                      priors['alpha'],
                      priors['beta'],
                      data)
        
    model.fit(convergence_tolerance=args.convergence_tolerance, num_iters=args.max_num_iters, debug=False)
    
    if args.lower_bound_file is not None:
        write_lower_bound(model, args.lower_bound_file)

    if args.out_dir is not None:
        write_cluster_posteriors(cell_ids, model.pi_star, args.out_dir)
        
        write_genotype_posteriors(event_ids, model.mu_star, args.out_dir)
            
        write_params(model, args.out_dir, event_ids, time.clock() - t0)
    
##############################
    
def run_basic_miss_gemm_model(args):

    t0 = time.clock()
    if args.seed is not None:
        np.random.seed(args.seed)
       
    cell_ids, data, event_ids, priors, regions = load_data(args.config_file)
   
    model = BasicMissGeMM(priors['gamma'],
                          priors['alpha'],
                          priors['beta'],
                          data)
        
    model.fit(convergence_tolerance=args.convergence_tolerance, num_iters=args.max_num_iters, debug=False)
    
    if args.lower_bound_file is not None:
        write_lower_bound(model, args.lower_bound_file)

    if args.out_dir is not None:
        write_cluster_posteriors(cell_ids, model.pi_star, args.out_dir)
        
        write_genotype_posteriors(event_ids, model.mu_star, args.out_dir)
            
        write_params(model, args.out_dir, event_ids, time.clock() - t0)    
        
        write_imputed_data_MAP (cell_ids, model, args.out_dir)    

##############################    

def run_region_gemm_model(args):

    t0 = time.clock()
    if args.seed is not None:
        np.random.seed(args.seed)
       
    cell_ids, data, event_ids, priors, regions = load_data(args.config_file, include_regions=True)
   
    model = RegionGeMM(priors['gamma'],
                      priors['alpha'],
                      priors['beta'],
                      data,
                      regions)
        
    model.fit(convergence_tolerance=args.convergence_tolerance, num_iters=args.max_num_iters, debug=False)
    
    if args.lower_bound_file is not None:
        write_lower_bound(model, args.lower_bound_file)

    if args.out_dir is not None:
        write_cluster_posteriors(cell_ids, model.pi_star, args.out_dir)
        
        write_genotype_posteriors(event_ids, model.mu_star, args.out_dir)
            
        write_params(model, args.out_dir, event_ids, time.clock() - t0)
        
##############################        
    
def load_data(yaml_filename, data_filename, include_regions=False):
# TODO: separate this into different loads: one for config, one for data, one for regions
    with open(yaml_filename) as fh:
        config = yaml.load(fh)
    
    cell_ids = []
    
    data = {}
    
    event_ids = {}
    
    priors = {'gamma' : {}, 'beta' : {}, 'alpha' : {} }
    
    # Used only by the region models
    regions = {}

    for data_type in config['data']:
        # data[data_type] = _load_data_frame(config['data'][data_type]['file'])
        data[data_type] = _load_data_frame(data_filename)
        if (include_regions):
            regions[data_type] = _load_regions_frame(config['data'][data_type]['region_file'])
        
        priors['gamma'][data_type] = np.array(config['data'][data_type]['gamma_prior'])
        
        priors['beta'][data_type] = np.array(config['data'][data_type]['beta_prior'])
        # normalize
        priors['beta'][data_type] = priors['beta'][data_type] / priors['beta'][data_type].sum()
        
        cell_ids.append(data[data_type].index)
        
        event_ids[data_type] = list(data[data_type].columns)

    cell_ids = sorted(set.intersection(*[set(x) for x in cell_ids]))
    
    for data_type in data:
        data[data_type] = data[data_type].loc[cell_ids]
    
    priors['alpha'] = np.ones(config['num_clusters']) * config['alpha_prior']
    
    # I'm giving just one number in the yaml file
    #if 'alpha_prior' in config:
    #    priors['alpha'] = np.array(config['alpha_prior'])
    
    print 'Number of cells: {0}'.format(len(cell_ids))
    
    print 'Number of data types: {0}'.format(len(event_ids))
    
    print 'Alpha priors: {0}'.format(priors['alpha'])            
    
    for data_type in event_ids:
        print 'Number of {0} events: {1}'.format(data_type, len(event_ids[data_type]))
        print 'Beta prior of {0} events: {1}'.format(data_type, priors['beta'][data_type])
        print 'Gamma prior of {0} events: {1}'.format(data_type, priors['gamma'][data_type])
        #print 'Data: \n', data[data_type]
        #print 'Regions: \n', regions[data_type] 
        
        # This is how to iterate through the regions:
        #for index, row in regions[data_type].iterrows():
        #    print index, row["start"], row["end"]        
    
    return  cell_ids, data, event_ids, priors, regions

def _load_data_frame(file_name):
    print 'Loading data file {0}.'.format(file_name)
    df = pd.read_csv(file_name, compression='gzip', index_col='cell_id', sep='\t')
    
    return df

def _load_regions_frame(file_name):
    # Assume the data file contains only regions, and the regions start and end are the column index from the data
    print 'Loading regions file {0}.'.format(file_name)
    df = pd.read_csv(file_name, compression='gzip', index_col='region', sep='\t')
    
    return df


def load_samples(cell_ids, file_name):
    if file_name is not None:
        samples = pd.read_csv(file_name, compression='gzip', index_col='cell_id', sep='\t')
        
        if not set(samples.index) == set(cell_ids):
            raise Exception('Samples file must contain all entries from the data files.')
        
        samples = samples.loc[cell_ids, 'sample']
        
        print 'Number of samples: {0}'.format(samples.nunique())
    
    else:
        samples = None
        
        print 'No samples file supplied. All cells assumed to come from same sample.'
    
    return samples

def write_cluster_posteriors(cell_ids, Z, out_dir):
    file_name = os.path.join(out_dir, 'cluster_posteriors.tsv.gz')

    df = pd.DataFrame(Z, index=cell_ids)
    
    with gzip.GzipFile(file_name, 'w') as fh:
        df.to_csv(fh, index_label='cell_id', sep='\t')
        
        
def write_imputed_data_MAP(cell_ids, model, out_dir):
    # note: here I don't write the posterior probabilities, but the MAP values
    file_name = os.path.join(out_dir, 'imputed_data_MAP.tsv.gz')

    for data_type in model.data_types:
        # first create a single matrix of size NxM from rho_star
        rho_star = model.get_rho_star(data_type)        

        imputed_matrix = np.argmax(rho_star, axis=0)

        df = pd.DataFrame(imputed_matrix, index=cell_ids)
    
        with gzip.GzipFile(file_name, 'w') as fh:
            df.to_csv(fh, index_label='cell_id', sep='\t')                
        
           

def write_double_cluster_posteriors(cell_ids, model, out_dir):
    file_name = os.path.join(out_dir, 'double_cluster_posteriors.tsv.gz')
    
    df = pd.DataFrame(model.Z_1.reshape(model.N, model.K * model.K), index=cell_ids)
    
    with gzip.GzipFile(file_name, 'w') as fh:
        df.to_csv(fh, index_label='cell_id', sep='\t')    

def write_genotype_posteriors(event_ids, G, out_dir):
    file_name = os.path.join(out_dir, 'genotype_posteriors.tsv.gz')
    
    G_out = []
    
    for data_type in event_ids:
        for i, df in enumerate(G[data_type]):
            df = pd.DataFrame(df, columns=event_ids[data_type])
            
            df = df.stack().reset_index()
            
            df.columns = 'cluster_id', 'event_id', 'probability'
            
            df.insert(1, 'event_type', data_type)
            
            df.insert(3, 'event_value', i)
            
            G_out.append(df)
    
    G_out = pd.concat(G_out)
  
    with gzip.GzipFile(file_name, 'w') as fh:
        G_out.to_csv(fh, index=False, sep='\t')

def write_params(model, out_dir, event_ids, cpu_time):
    file_name = os.path.join(out_dir, 'params.yaml')

              #'lower_bound' : float(model.lower_bound[-1]),
    
    params = {
              'lower_bound' : int(model.lower_bound[-1]),
              'converged' : model.converged
              }
    
    params['alpha'] = {}
    
    params['alpha'] = [float(x) for x in model.alpha_star]
    
    params['gamma'] = {}
    
    params['CPU_time'] = int(cpu_time)
    
    for data_type in model.gamma_star:
        params['gamma'][data_type] = [[float(x) for x in row] for row in model.gamma_star[data_type]]
    
    
    # Do I have to change here??
    if hasattr(model, 'alpha'):
        params['alpha'] = [float(x) for x in model.alpha]

    with open(file_name, 'w') as fh:
        yaml.dump(params, fh, default_flow_style=False)
        
def write_lower_bound(model, out_file):
    with open(out_file, 'w') as fh:
        result = {'lower_bound' : float(model.lower_bound[-1]), 'converged' : model.converged}
        
        yaml.dump(result, fh, default_flow_style=False)
