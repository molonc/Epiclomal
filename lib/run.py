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
import csv
import os.path
from sklearn.metrics.cluster import v_measure_score, homogeneity_score, completeness_score

from memory_profiler import memory_usage

from lib.basic_gemm import BasicGeMM
from lib.basic_miss_gemm import BasicMissGeMM
from lib.region_gemm import RegionGeMM
from lib.region_miss_gemm import RegionMissGeMM
from lib.utils import load_labels


##############################

def run_basic_gemm_model(args):
    run_model('BasicGeMM', args)    
    
##############################

def run_basic_bayespy_model(args):
    run_model('BasicBayesPy', args)
        
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


    np.set_printoptions(threshold='nan')

    t0 = time.clock()
    if args.seed is not None:
        np.random.seed(args.seed)
       
    include_regions = False 
    impute = False
       
    if (mtype == 'BasicGeMM'):
        className = BasicGeMM
        
    elif (mtype == 'BasicBayesPy'):
        className = BasicBayesPy    
        
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
                       
    cell_ids, data, event_ids, priors, regions, initial_clusters = load_data (args, include_regions)
    # event_ids is "meth" etc.
    model = className(priors['gamma'],
                      priors['alpha'],
                      priors['beta'],
                      data,
                      regions,
                      initial_clusters)                                                                                                       
        
    # mem_usage = memory_usage((model.fit(convergence_tolerance=args.convergence_tolerance, num_iters=args.max_num_iters, debug=False)))
    # maxmem = max(mem_usage)
    maxmem = None
    model.fit(convergence_tolerance=args.convergence_tolerance, num_iters=args.max_num_iters, debug=False)

    # Measuring the memory
    maxmem = max(memory_usage(-1))
    # print(maxmem)


    if args.out_dir is not None:
        if not os.path.exists(args.out_dir):
            os.makedirs(args.out_dir)
            
        write_cluster_posteriors(cell_ids, model.pi_star, args.out_dir)        
        labels_pred = write_cluster_MAP(cell_ids, model.pi_star, args.out_dir)
                      
        write_genotype_posteriors(model, args.out_dir) 
        write_genotype_MAP(model, args.out_dir)         
        
        # function inherited from SCG which doesn't work right now
        # write_genotype_posteriors(event_ids, model.get_mu_star(), args.out_dir)         
         
        vmeasure = None    
        if args.true_clusters_file is not None:    	
            true_clusters = pd.read_csv(args.true_clusters_file, compression='gzip', index_col='cell_id', sep='\t')
            labels_true = np.array(true_clusters['epigenotype_id'])
            print ("True labels:")
            print (*labels_true)
            print ("Predicted labels:")
            print (*labels_pred)  
            vmeasure = v_measure_score(labels_true, labels_pred)
            print ('Vmeasure: {0:.5g}'.format(vmeasure))

        print('Maxmem: %s MB' % maxmem)
            
        write_params(model, args.out_dir, event_ids, time.clock() - t0, vmeasure, maxmem)
        
##############################        
    
def load_data(args, include_regions=False):
    yaml_filename = args.config_file
    meth_filename = args.methylation_file
    regions_filename = args.regions_file
    copynumber_filename = args.copynumber_file
    initial_clusters_data = None
    with open(yaml_filename) as fh:
        config = yaml.load(fh)
    
    if args.K is not None:
        num_clusters = int(args.K)
    else:
        num_clusters = config['num_clusters']        
    
    print ('Num clusters: ', num_clusters)
    
    cell_ids = []
    
    data = {}
    
    event_ids = {}
    
    priors = {'gamma' : {}, 'beta' : {}, 'alpha' : [] }
    
    # Used only by the region models
    regions = {}
    
    if (args.initial_clusters_file != None):
        initial_clusters_data = _load_initial_clusters_frame(args.initial_clusters_file, args.repeat_id)

    # print initial_clusters_data

    for data_type in config['data']:
        # data[data_type] = _load_data_frame(config['data'][data_type]['file'])
        if (data_type == 'meth'):
            data[data_type] = _load_data_frame(meth_filename)
        if (data_type == 'cn'):
            data[data_type] = _load_data_frame(copynumber_filename) 
                               
        if (include_regions):
            regions[data_type] = _load_regions_frame(regions_filename)
        else:
            regions = None             
        
        priors['gamma'][data_type] = np.array(config['data'][data_type]['gamma_prior'])
                
        if (args.bulk_file != None):
            # If the bulk reads file is given, there will be a different parameter for each locus
            priors['beta'][data_type] = _load_bulk_frame(args.bulk_file)
        else:
            priors['beta'][data_type] = np.array(config['data'][data_type]['beta_prior'])
        
        # NOTE: we don't need to normalize because these are parameters of a Beta/Dirichlet distribution        
        # priors['beta'][data_type] = priors['beta'][data_type] / priors['beta'][data_type].sum()
        
        cell_ids.append(data[data_type].index)
        
        event_ids[data_type] = list(data[data_type].columns)

    cell_ids = sorted(set.intersection(*[set(x) for x in cell_ids]))
    
    # not sure what this below does
    for data_type in data.keys():
        data[data_type] = data[data_type].loc[cell_ids]
    
    priors['alpha'] = np.ones(num_clusters) * config['alpha_prior']    
    
    # I'm giving just one number in the yaml file
    #if 'alpha_prior' in config:
    #    priors['alpha'] = np.array(config['alpha_prior'])
    
    print ('Number of cells: {}'.format(len(cell_ids)))
    
    print ('Number of data types: {}'.format(len(event_ids)))
    
    print ('Alpha priors: ', *priors['alpha'])
    
    # print(priors['beta'])
    
    for data_type in data.keys():
        print ('Number of {0} events: {1}'.format(data_type, len(event_ids[data_type])))
        # print 'Beta prior of {0} events: {1}'.format(data_type, priors['beta'][data_type])
        # TODO In python3 I am not able to print the gamma and beta priors            
        # TODO
        #print ('Gamma prior of main events:')
        #print(priors['gamma'][data_type])
        # TODO
        #print ('Beta prior of main events: {0}'.format(priors['beta'][data_type]))
        #print 'Data: \n', data[data_type]
        #print 'Regions: \n', regions[data_type] 
        
        # This is how to iterate through the regions:
        #for index, row in regions[data_type].iterrows():
        #    print index, row["start"], row["end"]        
    
    return  cell_ids, data, event_ids, priors, regions, initial_clusters_data

def _load_data_frame(file_name):
    print ('Loading data file {0}.'.format(file_name))
    df = pd.read_csv(file_name, compression='gzip', index_col='cell_id', sep='\t')    
    return df
    
def _load_initial_clusters_frame(file_name, repeat_id):
    # IF the file doesn't exist, return None
    if (os.path.isfile(file_name) == False):
        print ("Initial clusters file was given (", file_name, "), but it doesn't exist, ignore.")
        return None
    # First check if the number of columns in the file is > repeat_id
    # with gzip.open(file_name, 'r') as file:
    #     reader = list(csv.reader(file, delimiter='\t'))
    #     numcols = len(reader[0])
    df = pd.read_csv(file_name, compression='gzip', sep='\t') 
    numcols = df.shape[1]        
    print ('Num columns in initial_clusters_file is ', numcols, ', repeat_id is ', repeat_id)
    # Make repeat_id bigger by 1, because in kronos repeat_id starts from 0
    repeat_id = repeat_id+1
    if (repeat_id >= numcols):
        print ('Using random initialization')
        return None
    print ('Loading initial clusters file {0}.'.format(file_name))
    df = pd.read_csv(file_name, compression='gzip', index_col=0, sep='\t', usecols=[0,int(repeat_id)])        
    return df    

def _load_regions_frame(file_name):
    # Assume the data file contains only regions, and the regions start and end are the column index from the data
    print ('Loading regions file {0}.'.format(file_name))
    df = pd.read_csv(file_name, compression='gzip', index_col='region_id', sep='\t')
    return df

def _load_bulk_frame(file_name):
    # Assume the header of the bulk file is:
    # position        meth_reads      unmeth_reads
    print ('Loading bulk file {0}.'.format(file_name))
    df = pd.read_csv(file_name, compression='gzip', index_col='position', sep='\t')
    return df

def load_samples(cell_ids, file_name):
    if file_name is not None:
        samples = pd.read_csv(file_name, compression='gzip', index_col='cell_id', sep='\t')
        
        if not set(samples.index) == set(cell_ids):
            raise Exception('Samples file must contain all entries from the data files.')
        
        samples = samples.loc[cell_ids, 'sample']
        
        print ('Number of samples: {0}'.format(samples.nunique()))
    
    else:
        samples = None
        
        print ('No samples file supplied. All cells assumed to come from same sample.')
    
    return samples

##########################

def write_cluster_posteriors(cell_ids, pi_star, out_dir):
    file_name = os.path.join(out_dir, 'cluster_posteriors.tsv.gz')

    df = pd.DataFrame(pi_star, index=cell_ids)
    
#    with gzip.GzipFile(file_name, 'w') as fh:
#        df.to_csv(fh, index_label='cell_id', sep='\t')
    df.to_csv(file_name, index_label='cell_id', sep='\t', mode='w', compression='gzip')
    # For some reason the with gzip.GzipFile doesn't work in python3    

##########################    
       
def write_cluster_MAP(cell_ids, pi_star, out_dir):        
    file_name = os.path.join(out_dir, 'cluster_MAP.tsv.gz') 
     
    just_labels = []    # has only the cluster labels, we need this to compute vmeasure   
    labels_pred = []    # has cluster labels and probabilities
    # traverse every row, and find out which cluster has the maximum value
    df = pd.DataFrame(pi_star)
    for index, row in df.iterrows():
        max_index, max_value = max(enumerate(row), key=operator.itemgetter(1))
        labels_pred.append([max_index, max_value])
        just_labels.append(max_index)
        
    # print labels_pred
    df = pd.DataFrame(labels_pred, index=cell_ids)
    #with gzip.GzipFile(file_name, 'w') as fh:
    #    df.to_csv(fh, index_label='cell_id', sep='\t')          
    df.to_csv(file_name, index_label='cell_id', sep='\t', mode='w', compression='gzip')
    return just_labels
       
        
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
    
    # with gzip.GzipFile(file_name, 'w') as fh:
    #     df.to_csv(fh, index_label='cell_id', sep='\t')    
    df.to_csv(file_name, index_label='cell_id', sep='\t', mode='w', compression='gzip')

##########################

def write_genotype_posteriors(model, out_dir):
# Updated on 20 July 2017
    file_name = os.path.join(out_dir, 'genotype_posteriors.tsv.gz')
    
    for data_type in model.data_types:
        u_mu_star = model.unregion_mu_star(data_type)             
        df = pd.DataFrame(u_mu_star[0], index=range(u_mu_star.shape[1]))
  
        # with gzip.GzipFile(file_name, 'w') as fh:
        #     df.to_csv(fh, index_label='cluster_id', sep='\t')
        df.to_csv(file_name, index_label='cluster_id', sep='\t', mode='w', compression='gzip')    

##########################

def write_genotype_MAP(model, out_dir):
# Updated on 20 July 2017
    file_name = os.path.join(out_dir, 'genotype_MAP.tsv.gz')
    
    for data_type in model.data_types:
        u_mu_star = model.unregion_mu_star(data_type)             
        map_mu_star = np.argmax(u_mu_star, axis=0)
        df = pd.DataFrame(map_mu_star, index=range(u_mu_star.shape[1]))
  
        # with gzip.GzipFile(file_name, 'w') as fh:
        #     df.to_csv(fh, index_label='cluster_id', sep='\t')
        df.to_csv(file_name, index_label='cluster_id', sep='\t', mode='w', compression='gzip')    

##########################

# def write_genotype_posteriors(event_ids, G, out_dir):
# # G is actually mu_star, the fitted parameters of G
# # this has to be rewritten in a matrix format
#     file_name = os.path.join(out_dir, 'genotype_posteriors.tsv.gz')
#     
#     G_out = []
#     
#     for data_type in event_ids:
#         for i, df in enumerate(G[data_type]):
#             df = pd.DataFrame(df, columns=event_ids[data_type])
#             
#             df = df.stack().reset_index()
#             
#             df.columns = 'cluster_id', 'event_id', 'probability'
#             
#             df.insert(1, 'event_type', data_type)
#             
#             df.insert(3, 'event_value', i)
#             
#             G_out.append(df)
#     
#     G_out = pd.concat(G_out)
#   
#     with gzip.GzipFile(file_name, 'w') as fh:
#         G_out.to_csv(fh, index=False, sep='\t')

##########################

def write_params(model, out_dir, event_ids, cpu_time, vmeasure, maxmem):
    file_name = os.path.join(out_dir, 'params.yaml')
    
    params = {
              'lower_bound' : float(model.lower_bound[-1]),
              'log_likelihood' : float(model.log_likelihood),
              'log_posterior_clusterK' : float(model.log_posterior_clusterK),
              'log_posterior_allK' : float(model.log_posterior_allK),              
              'converged' : model.converged
              }
    
    params['alpha_star'] = {}
    
    #params['alpha_star'] = [float(x) for x in model.alpha_star]
    params['alpha_star'] = model.alpha_star.tolist()
    
    params['gamma_star'] = {}
    params['beta_star'] = {}    
    
    params['CPU_time_seconds'] = int(cpu_time)

    if vmeasure is not None:
        params['Vmeasure'] = float(vmeasure)
    
    if maxmem is not None:
        params['Max_memory_MB'] = int(maxmem)
    
    for data_type in model.gamma_star:
        # params['gamma_star'][data_type] = [[float(x) for x in row] for row in model.gamma_star[data_type]]
        params['gamma_star'][data_type] = model.gamma_star[data_type].tolist()
    for data_type in model.beta_star:        
        # params['beta_star'][data_type] = [float(x) for x in model.beta_star[data_type]]
        params['beta_star'][data_type] = model.beta_star[data_type].tolist()        
        
        # I'm not reshaping any more, bulk and region have larger beta_stars
        #if model.include_bulk[data_type]:
        #    params['beta_star'][data_type] = model.beta_star[data_type].reshape(model.K, model.S[data_type]).tolist()
        #else:
        #    params['beta_star'][data_type] = model.beta_star[data_type].reshape(model.K, model.S[data_type]).tolist()

    with open(file_name, 'w') as fh:
        yaml.dump(params, fh, default_flow_style=False)
        
##########################        
        
