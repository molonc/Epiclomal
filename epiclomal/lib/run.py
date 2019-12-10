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
import sys
from sklearn.metrics.cluster import v_measure_score, homogeneity_score, completeness_score
from sklearn.metrics import mean_absolute_error, mean_squared_error


from memory_profiler import memory_usage

from epiclomal.lib.basic_gemm import BasicGeMM
from epiclomal.lib.region_gemm import RegionGeMM


##############################

def run_basic_gemm_model(args):
    if (args.Bishop_model_selection is True):
        print ("Running Bishop model selection")
        run_Bishop_model_selection('BasicGeMM', args)
    else:
        run_model('BasicGeMM', args)

##############################

def run_basic_bayespy_model(args):
    run_model('BasicBayesPy', args)

##############################

def run_region_gemm_model(args):
    run_model('RegionGeMM', args)

##############################

def run_Bishop_model_selection(mtype, args):


    np.set_printoptions(threshold='nan')

    t0 = time.clock()
    if args.seed is not None:
        np.random.seed(args.seed)

    include_regions = False

    if (mtype == 'BasicGeMM'):
        className = BasicGeMM

    elif (mtype == 'BasicBayesPy'):
        className = BasicBayesPy

    elif (mtype == 'RegionGeMM'):
        include_regions=True
        className = RegionGeMM

    cell_ids, data, event_ids, priors, regions, initial_clusters, slsbulk_data = load_data (args, include_regions)
    maxmem = None

    # make a loop to start at K=10 clusters and then go to fewer clusters

    currentK = len(priors['alpha']) + 1
    newK = currentK - 1   # first it should be 10


    while (newK < currentK):
        currentK = newK
        print ('==========================')
        print ('Running with K ', currentK)
        print ('==========================')
        args.repeat_id = 9 + currentK      # this is PBAL 10 in the initial_clusters_file, which has to be given

        # reload the initial data with the new repeat_id
        if (args.initial_clusters_file != None):
            initial_clusters_data = _load_initial_clusters_frame(args.initial_clusters_file, args.repeat_id)

        priors['alpha'] = priors['alpha'][range(currentK)]

        # event_ids is "meth" etc.
        model = className(priors['gamma'],
                          priors['alpha'],
                          priors['beta'],
                          data,
                          regions,
                          initial_clusters,
                          args.mu_has_k,
                          Bishop_model_selection = args.Bishop_model_selection)

        # mem_usage = memory_usage((model.fit(convergence_tolerance=args.convergence_tolerance, num_iters=args.max_num_iters, debug=False)))
        # maxmem = max(mem_usage)

        model.fit(convergence_tolerance=args.convergence_tolerance, num_iters=args.max_num_iters, debug=False)

        newK = model._compute_num_nonzero_pi()


    # Measuring the memory
    maxmem = max(memory_usage(-1))
    # print(maxmem)


    if args.out_dir is not None:
        if not os.path.exists(args.out_dir):
            os.makedirs(args.out_dir)

        write_cluster_posteriors(cell_ids, model.pi_star, args.out_dir)
        labels_pred, labels_prob = write_cluster_MAP(cell_ids, model.pi_star, args.out_dir)

        write_genotype_posteriors(model, args.out_dir)
        write_genotype_MAP(model, args.out_dir)

        # function inherited from SCG which doesn't work right now
        # write_genotype_posteriors(event_ids, model.get_mu_star(), args.out_dir)

        vmeasure = None
        if args.true_clusters_file is not None:
            true_clusters = pd.read_csv(args.true_clusters_file, compression='gzip', sep='\t')
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

def run_model(mtype, args):

    np.set_printoptions()

    t0 = time.clock()

    if args.seed is None:
        args.seed = int(time.time())
    np.random.seed(args.seed)
    print("Using seed: ", args.seed)

    include_regions = False

    if (mtype == 'BasicGeMM'):
        className = BasicGeMM

    elif (mtype == 'BasicBayesPy'):
        className = BasicBayesPy

    elif (mtype == 'RegionGeMM'):
        include_regions=True
        className = RegionGeMM

    cell_ids, data, event_ids, priors, regions, initial_clusters, slsbulk = load_data(args, include_regions)
    # event_ids is "meth" etc.
    model = className(priors['gamma'],
                      priors['alpha'],
                      priors['beta'],
                      data,
                      regions,
                      initial_clusters,
                      args.mu_has_k,
                      Bishop_model_selection = args.Bishop_model_selection,
                      slsbulk_data = slsbulk)

    model.seed = args.seed
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
        labels_pred, labels_prob = write_cluster_MAP(cell_ids, model.pi_star, args.out_dir)

        write_genotype_posteriors(model, args.out_dir)
        epigenotype = write_genotype_MAP(model, args.out_dir)

        # Not doing the correction here, but in the hamming_distance.R script when I check the hamming distance
        # print("Doing genotype correction")
        # epigenotype = write_genotype_MAP_corrected(model, args.out_dir, labels_pred)

        # function inherited from SCG which doesn't work right now
        # write_genotype_posteriors(event_ids, model.get_mu_star(), args.out_dir)

        print ("Predicted labels:")
        print (*labels_pred)

        vmeasure = None
        model.clone_prev_MAE = None
        model.clone_prev_MSE = None
        if args.true_clusters_file is not None:
            true_clusters = pd.read_csv(args.true_clusters_file, compression='gzip', sep='\t')
            # Order the true clusters by the predicted ones
            true_clusters['cell_id'] = pd.Categorical(true_clusters['cell_id'],cell_ids)
            true_clusters.sort_values(by=['cell_id'], inplace=True)
            # Sometimes pred_clusters has fewer cells than true_clusters, so taking only those
            true_clusters = true_clusters[true_clusters['cell_id'].isin(cell_ids)]
            # checking that they are in the same order
            # print("Len of cell_ids")
            # print(len(cell_ids))
            for i in range(len(cell_ids)):
                # print(i)
                # print(i, " ", cell_ids[i])
                if (cell_ids[i] != true_clusters.iloc[i]['cell_id']):
                    print ("Cell ids are different for cell ", i, " predicted ", cell_ids[i], " true ", true_clusters.iloc[i]['cell_id'])
                    sys.exit(2)
            labels_true = np.array(true_clusters['epigenotype_id'])
            print ("True labels:")
            print (*labels_true)
            vmeasure = v_measure_score(labels_true, labels_pred)
            print ('Vmeasure: {0:.5g}'.format(vmeasure))

            # Now calculating the clone prevalence error
            # First calculating the true clone prevalence

            # NOTE: for the calculated prev_true, I am also not including +1 in prev_pred
            prev_true = []
            prev_pred = []

            if args.true_prevalences is None:
                print("True prevalences None, calculating")
                for n in range(model.N):
                    prev_true.append(sum(np.equal(labels_true,labels_true[n]))/model.N)
                    prev_pred.append((sum(np.equal(labels_pred,labels_pred[n])))/(model.N))
            else:
                print("Using given true prevalences")
                prevs = args.true_prevalences.split("_")
                for n in range(model.N):
                    # labels_true start from 1, so I have to subtract 1 to have indexes starting from 0
                    prev_true.append(float(prevs[labels_true[n]-1]))
                    prev_pred.append((sum(np.equal(labels_pred,labels_pred[n]))+1)/(model.N+model.nclusters))

            # print("Prev true ", *prev_true)

            model.clone_prev_MAE = mean_absolute_error(prev_true, prev_pred)
            model.clone_prev_MSE = mean_squared_error(prev_true, prev_pred)


#######
            #prev_true = []
            #prevs = args.true_prevalences.split("_")
            #prev_pred = []
            #for n in range(model.N):
                # For the true, it's better not to count the true cells, but to use a given true prevalence
                # prev_true.append(sum(np.equal(labels_true,labels_true[n]))/model.N)
                # labels_true start from 1, so I have to subtract 1 to have indexes starting from 0
            #    prev_true.append(float(prevs[labels_true[n]-1]))
                # For the prediction, adding 1 for the prior, this will match alpha_star
            #    prev_pred.append((sum(np.equal(labels_pred,labels_pred[n]))+1)/(model.N+model.nclusters))
            #print ("True clone prevalences:")
            #print(*prev_true)
            #print ("Pred clone prevalences:")
            #print(*prev_pred)
            #model.clone_prev_MAE = mean_absolute_error(prev_true, prev_pred)
            #model.clone_prev_MSE = mean_squared_error(prev_true, prev_pred)
            print("Clone prevalence mean absolute error: ", model.clone_prev_MAE)
            print("Clone prevalence mean squared error: ", model.clone_prev_MSE)

        print('Maxmem: %s MB' % maxmem)


        ######################################################
        #
        # SLS BULK
        # First find the regions that are different between the found clusters

        model.slsbulk_vmeasure = None
        model.slsbulk_clone_prev_MAE = None
        model.slsbulk_clone_prev_MSE = None
        model.slsbulk_iterations = int(args.slsbulk_iterations)
        model.uncertainty_tpr = None
        if (args.slsbulk_file != None or args.check_uncertainty):
            different_regions = model._compute_different_regions(labels_pred, epigenotype)
            print('The different regions are ', *different_regions)
            if len(different_regions) > 0:
            # DO this if there is at list one different region, otherwise there's no point
                # Second find the cells that have completely missing data in the different regions
                candidate_cells = model._compute_candidate_cells(labels_pred, epigenotype, different_regions)
                # This is dictionary where the keys are the candidate cells and the values are the possible clusters for each candidate cell
                print('Candidate cells and clusters are ', *candidate_cells)
                if (args.slsbulk_file != None):

                    new_pred = model._slsbulk(candidate_cells, labels_pred, epigenotype, different_regions, labels_true)
                    print ("True labels:")
                    print (*labels_true)
                    print ("New labels:")
                    print (*new_pred)
                    model.slsbulk_vmeasure = v_measure_score(labels_true, new_pred)
                    print ('NEW Vmeasure: {0:.5g}'.format(model.slsbulk_vmeasure))

                    prev_pred = []
                    for n in range(model.N):
                        prev_pred.append((sum(np.equal(new_pred,new_pred[n]))+1) / (model.N + model.nclusters))
                    model.slsbulk_clone_prev_MAE = mean_absolute_error(prev_true, prev_pred)
                    model.slsbulk_clone_prev_MSE = mean_squared_error(prev_true, prev_pred)
                    print ("New clone prev MAE ", model.slsbulk_clone_prev_MAE)
                    write_slsbulk_cluster_MAP(cell_ids, args.out_dir, new_pred)

                if (args.check_uncertainty):
                    true_positive_rate = 0
                    candidate_keys =  list(candidate_cells.keys())
                    for key in candidate_keys:
                        print ("Current prediction ", labels_prob[key])
                        print ("    Cell ", key, " with candidate clusters ", candidate_cells[key], "of size", len(candidate_cells[key]))
                        expected_prob = 1./len(candidate_cells[key])
                        if (labels_prob[key] >= expected_prob - 0.2 and labels_prob[key] <= expected_prob + 0.2):
                            true_positive_rate += 1
                    if len(candidate_cells) > 0:
                        true_positive_rate /= len(candidate_cells)
                        model.uncertainty_tpr = true_positive_rate
                        print("Uncertainty true positive rate: ", true_positive_rate)
                    else:
                        print("Uncertainty tpr is None")

        write_params(model, args.out_dir, event_ids, time.clock() - t0, vmeasure, maxmem)

##############################

def load_data(args, include_regions=False):
    yaml_filename = args.config_file
    meth_filename = args.methylation_file
    regions_filename = args.regions_file
    copynumber_filename = args.copynumber_file
    initial_clusters_data = None
    slsbulk_data = None
    with open(yaml_filename) as fh:
        config = yaml.safe_load(fh)

    if args.K is not None:
        num_clusters = int(args.K)
    else:
        num_clusters = config['num_clusters']

    cell_ids = []

    data = {}

    event_ids = {}

    priors = {'gamma' : {}, 'beta' : {}, 'alpha' : [] }

    # Used only by the region models
    regions = {}

    if (args.initial_clusters_file != None):
        initial_clusters_data = _load_initial_clusters_frame(args.initial_clusters_file, args.repeat_id)

    # 12 Apr 2018: If it's a random initialization, also make num_clusters be random from 1 to K
    if (initial_clusters_data is None):
        print ('Selecting a random number of clusters from 1 to ', num_clusters)
        num_clusters = np.random.randint(low=1, high=num_clusters)
    print ('Num clusters: ', num_clusters)

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

        if (args.slsbulk_file != None):
            # If the bulk reads file is given, we'll calculate a score
            slsbulk_data = _load_bulk_frame(args.slsbulk_file)


        # NOTE: we don't need to normalize because these are parameters of a Beta/Dirichlet distribution
        # In SCG they were normalized because they were probabilities of the G binomial distribution
        # priors['beta'][data_type] = priors['beta'][data_type] / priors['beta'][data_type].sum()

        cell_ids.append(data[data_type].index)
        # print(cell_ids)

        event_ids[data_type] = list(data[data_type].columns)

    cell_ids = sorted(set.intersection(*[set(x) for x in cell_ids]))

    # not sure what this below does
    # PY: me neither, data['meth'].loc[cell_ids].equals(data['meth']) == True
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

    return  cell_ids, data, event_ids, priors, regions, initial_clusters_data, slsbulk_data

def _load_data_frame(file_name):
    print ('Loading data file {0}.'.format(file_name))
    df = pd.read_csv(file_name, compression='gzip', index_col='cell_id', sep='\t', low_memory=False)
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
    df = pd.read_csv(file_name, compression='gzip', index_col='region_id', sep='\t', dtype='int')
    return df

def _load_bulk_frame(file_name):
    # Assume the header of the bulk file is:
    # position        meth_reads      unmeth_reads
    print ('Loading slsbulk file {0}.'.format(file_name))
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
    just_prob = []      # has only the cluster probabilities
    labels_pred = []    # has cluster labels and probabilities
    # traverse every row, and find out which cluster has the maximum value
    df = pd.DataFrame(pi_star)
    for index, row in df.iterrows():
        max_index, max_value = max(enumerate(row), key=operator.itemgetter(1))
        labels_pred.append([max_index, max_value])
        just_labels.append(max_index)
        just_prob.append(max_value)

    # print labels_pred
    #labels_pred.sort_values(by=['cell_id'], inplace=True)
    df = pd.DataFrame(labels_pred, index=cell_ids)
    #print(df)
    #with gzip.GzipFile(file_name, 'w') as fh:
    #    df.to_csv(fh, index_label='cell_id', sep='\t')
    df.to_csv(file_name, index_label='cell_id', sep='\t', mode='w', compression='gzip')
    return just_labels, just_prob


##########################

def write_slsbulk_cluster_MAP(cell_ids, out_dir, predicted_labels):
    file_name = os.path.join(out_dir, 'slsbulk_cluster_MAP.tsv.gz')

    labels_pred = []    # has cluster labels and probabilities
    # traverse every row, and find out which cluster has the maximum value
    for i in range(len(predicted_labels)):
        labels_pred.append([predicted_labels[i], 1.0])

    # print labels_pred
    df = pd.DataFrame(labels_pred, index=cell_ids)

    df.to_csv(file_name, index_label='cell_id', sep='\t', mode='w', compression='gzip')


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
    return map_mu_star


##########################

# not finished
# def correct_genotype (genotype, labels_pred):
# 	print("In correct genotype")
#	print(genotype.shape)
#	print(genotype)
#	print("Labels pred")
#	print(labels_pred)
#	return(genotype)

##########################

# not finished
# def write_genotype_MAP_corrected(model, out_dir, labels_pred):
#     file_name = os.path.join(out_dir, 'genotype_MAP_corrected.tsv.gz')
#
#     for data_type in model.data_types:
#         u_mu_star = model.unregion_mu_star(data_type)
#         map_mu_star = np.argmax(u_mu_star, axis=0)
#         df = pd.DataFrame(map_mu_star, index=range(u_mu_star.shape[1]))
#
#         # 11 Mar 2019: correcting the genotypes
#         df = correct_genotype (df, labels_pred)
#
#         # with gzip.GzipFile(file_name, 'w') as fh:
#         #     df.to_csv(fh, index_label='cluster_id', sep='\t')
#         df.to_csv(file_name, index_label='cluster_id', sep='\t', mode='w', compression='gzip')
#     return map_mu_star

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
                # TODO: doing it all with mean for now, mode was very very similar
              'lower_bound' : float(model.lower_bound[-1]),
              'lower_bound_ElogP' : float(model.ElogP),
              'lower_bound_ElogQ' : float(model.ElogQ),
              'log_likelihood' : float(model.log_likelihood_mean),
              'log_posterior_clusterK' : float(model.log_posterior_mean_clusterK),
              'log_posterior_allK' : float(model.log_posterior_mean_allK),
              'DIC_pd' : float(model.DIC_pd),
              'DIC_measure' : float(model.DIC_measure),
              #'log_likelihood_mode' : float(model.log_likelihood_mode),
              #'log_posterior_mode_clusterK' : float(model.log_posterior_mode_clusterK),
              #'log_posterior_mode_allK' : float(model.log_posterior_mode_allK),
              'converged': model.converged,
              'seed': model.seed,
              'nclusters': int(model.nclusters)
              #'prevalences': model.prevalences
            }
    params['prevalences'] = [float(x) for x in np.array(model.prevalences)]

    if(model.Bishop_model_selection is False):
        params['alpha_star'] = {}

        #params['alpha_star'] = [float(x) for x in model.alpha_star]
        params['alpha_star'] = model.alpha_star.tolist()

    params['gamma_star'] = {}
    params['beta_star'] = {}

    params['CPU_time_seconds'] = int(cpu_time)

    if vmeasure is not None:
        params['Vmeasure'] = float(vmeasure)

    if model.clone_prev_MAE != None:
        params['clone_prev_MAE'] = float(model.clone_prev_MAE)
        params['clone_prev_MSE'] = float(model.clone_prev_MSE)

    if maxmem is not None:
        params['Max_memory_MB'] = int(maxmem)

    for data_type in model.gamma_star:
        # params['gamma_star'][data_type] = [[float(x) for x in row] for row in model.gamma_star[data_type]]
        params['gamma_star'][data_type] = model.gamma_star[data_type].tolist()
    for data_type in model.beta_star:
        # params['beta_star'][data_type] = [float(x) for x in model.beta_star[data_type]]
        #print('Beta star type: ', type(model.beta_star[data_type]))
        #print(*model.beta_star[data_type])
        params['beta_star'][data_type] = model.beta_star[data_type].tolist()

        # I'm not reshaping any more, bulk and region have larger beta_stars
        #if model.include_bulk[data_type]:
        #    params['beta_star'][data_type] = model.beta_star[data_type].reshape(model.K, model.S[data_type]).tolist()
        #else:
        #    params['beta_star'][data_type] = model.beta_star[data_type].reshape(model.K, model.S[data_type]).tolist()

    if model.slsbulk_vmeasure != None:
        params['slsbulk_vmeasure'] = float(model.slsbulk_vmeasure)
        params['slsbulk_clone_prev_MAE'] = float(model.slsbulk_clone_prev_MAE)
        params['slsbulk_clone_prev_MSE'] = float(model.slsbulk_clone_prev_MSE)

    if model.uncertainty_tpr != None:
        params['uncertainty_true_positive_rate'] = float(model.uncertainty_tpr)

    with open(file_name, 'w') as fh:
        yaml.dump(params, fh, default_flow_style=False)

##########################

