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
import operator

from lib.basic_gemm import BasicGeMM
from lib.basic_miss_gemm import BasicMissGeMM
from lib.region_gemm import RegionGeMM
from lib.region_miss_gemm import RegionMissGeMM
from lib.utils import load_labels


##############################

def run_basic_gemm_model(args):
    run_model('BasicGeMM', args)
    
##############################
    
def run_basic_miss_gemm_model(args):
    run_model('BasicMissGeMM', args)

##############################    

def run_region_gemm_model(args):
    run_model('RegionGeMM', args)
        
##############################    

def run_region_miss_gemm_model(args):
    run_model('RegionMissGeMM', args)      
        
##############################    

def run_model(mtype, args):

    t0 = time.clock()
    if args.seed is not None:
        np.random.seed(args.seed)
       
    include_regions = False 
    impute = False
       
    if (mtype == 'BasicGeMM'):
        className = BasicGeMM
        
    elif (mtype == 'BasicMissGeMM'):
        className = BasicMissGeMM
        impute = True
        
    elif (mtype == 'RegionGeMM'):
        include_regions=True
        className = RegionGeMM
        
    elif (mtype == 'RegionMissGeMM'):
        include_regions=True
        className = RegionMissGeMM
        impute = True
                       
    cell_ids, data, event_ids, priors, regions = load_data (args.config_file, 
                        args.data_file, 
                        args.regions_file, 
                        include_regions)
    model = className(priors['gamma'],
                      priors['alpha'],
                      priors['beta'],
                      data,
                      regions)                                                                                                       
        
    model.fit(convergence_tolerance=args.convergence_tolerance, num_iters=args.max_num_iters, debug=False)
    
    if args.lower_bound_file is not None:
        write_lower_bound(model, args.lower_bound_file)

    if args.out_dir is not None:
        write_cluster_posteriors(cell_ids, model.pi_star, args.out_dir)        
        write_cluster_MAP(cell_ids, model.pi_star, args.out_dir)
        if (impute):
            write_imputed_data_MAP(cell_ids, model, args.out_dir)
            
            
        # to correct this function    
        # write_genotype_MAP(event_ids, model, args.out_dir) 
            
        write_params(model, args.out_dir, event_ids, time.clock() - t0)          
        
##############################        
    
def load_data(yaml_filename, data_filename, regions_filename=None, include_regions=False):
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
            regions[data_type] = _load_regions_frame(regions_filename)
        else:
            regions = None            
        
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

##########################

def write_cluster_posteriors(cell_ids, pi_star, out_dir):
    file_name = os.path.join(out_dir, 'cluster_posteriors.tsv.gz')

    df = pd.DataFrame(pi_star, index=cell_ids)
    
    with gzip.GzipFile(file_name, 'w') as fh:
        df.to_csv(fh, index_label='cell_id', sep='\t')
    
##########################    
       
def write_cluster_MAP(cell_ids, pi_star, out_dir):        
    file_name = os.path.join(out_dir, 'cluster_MAP.tsv.gz') 
        
    labels_pred = []
    # traverse every row, and find out which cluster has the maximum value
    df = pd.DataFrame(pi_star)
    for index, row in df.iterrows():
        max_index, max_value = max(enumerate(row), key=operator.itemgetter(1))
        labels_pred.append(max_index)
    df = pd.DataFrame(np.transpose(labels_pred), index=cell_ids)
    with gzip.GzipFile(file_name, 'w') as fh:
        df.to_csv(fh, index_label='cell_id', sep='\t')          
        
##########################        
        
def write_imputed_data_MAP(cell_ids, model, out_dir):
    # note: here I don't write the posterior probabilities, but the MAP values
    file_name = os.path.join(out_dir, 'imputed_data_MAP.tsv.gz')

    for data_type in model.data_types:
        # first create a single matrix of size NxM from rho_star
        imputed_matrix = model.unregion_rho_star(data_type)        

        df = pd.DataFrame(imputed_matrix, index=cell_ids)
    
        with gzip.GzipFile(file_name, 'w') as fh:
            df.to_csv(fh, index_label='cell_id', sep='\t')                
        
##########################           

# NOT SURE what this would be useful for
def write_double_cluster_posteriors(cell_ids, model, out_dir):
    file_name = os.path.join(out_dir, 'double_cluster_posteriors.tsv.gz')
    
    df = pd.DataFrame(model.Z_1.reshape(model.N, model.K * model.K), index=cell_ids)
    
    with gzip.GzipFile(file_name, 'w') as fh:
        df.to_csv(fh, index_label='cell_id', sep='\t')    

##########################

def write_genotype_MAP(event_ids, model, out_dir):
# TO DO: this is not working properly, to rewrite!!
# this has to be rewritten in a matrix format
    file_name = os.path.join(out_dir, 'genotype_MAP.tsv.gz')
    
    for data_type in model.data_types:
        print 'mu star\n', model.get_mu_star(data_type)
        geno_matrix = model.unregion_mu_star(data_type)
        
        print 'Geno matrix\n', geno_matrix
        df = pd.DataFrame(geno_matrix, index=event_ids)        
  
        with gzip.GzipFile(file_name, 'w') as fh:
            df.to_csv(fh, index=False, sep='\t')

##########################

def write_genotype_posteriors(event_ids, G, out_dir):
# G is actually mu_star, the fitted parameters of G
# this has to be rewritten in a matrix format
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

##########################

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
        
##########################        
        
def write_lower_bound(model, out_file):
    with open(out_file, 'w') as fh:
        result = {'lower_bound' : float(model.lower_bound[-1]), 'converged' : model.converged}
        
        yaml.dump(result, fh, default_flow_style=False)
