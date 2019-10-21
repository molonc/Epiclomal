#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import operator
import sys
from sklearn.metrics.cluster import v_measure_score, homogeneity_score, completeness_score
from sklearn.metrics import mean_absolute_error, mean_squared_error

def main():
    parser = argparse.ArgumentParser(description='Evaluate the clustering v-measure score.')
    parser.add_argument('--true_clusters_file', required=True,
                        help='The file with the true clustering assignment')
    parser.add_argument('--true_prevalences', default="None",
                        help='The file with the true prevalences, give "None" if not available.')
    parser.add_argument('--predicted_clusters_file', required=True,
                        help='The file with the epiclomal predicted clusters')
    parser.add_argument('--clusters_are_probabilities', default=True,
                        help='Set to True if the epiclomal_clusters_file is a posterior, False if it is a point estimate')
    parser.add_argument('--results_file', required=True,
                        help='The output file where to write the results')


    # parser.add_argument('--column', default=2, type=int,
    #                     help='Default is the second column, or set which column to consider')



    args = parser.parse_args()
    # print(args.true_clusters_file)
    # print(args.epiclomal_clusters_file)


    # first read the predicted file and figure out how many clusters there are
    pred_clusters = pd.read_csv(args.predicted_clusters_file, compression='gzip', sep='\t')
    pred_clusters.sort_values(by=['cell_id'], inplace=True)
    #print(pred_clusters)

    labels_pred = []

    if args.clusters_are_probabilities is True:
        # traverse every row, and find out which cluster has the maximum value
        # print "Calculating the MAP value"
        for index, row in pred_clusters.iterrows():
            max_index, max_value = max(enumerate(row), key=operator.itemgetter(1))
            labels_pred.append(max_index)
    else:
        # select the last column using -1, so not using --column any more
        for value in pred_clusters.iloc[:,-1]:
            labels_pred.append(value)

    print (*labels_pred)

    # now read the true clusters file
    true_clusters = pd.read_csv(args.true_clusters_file, compression='gzip', sep='\t')
    true_clusters.sort_values(by=['cell_id'], inplace=True)

    # Sometimes pred_clusters has fewer cells than true_clusters, so taking only those
    true_clusters = true_clusters[true_clusters['cell_id'].isin(pred_clusters['cell_id'])]
    #print(true_clusters)
    labels_true = np.array(true_clusters['epigenotype_id'])
    # print (labels_true.shape)

    # Checking they are in the same order
    for i in range(len(pred_clusters)):
        if (pred_clusters.iloc[i]['cell_id'] != true_clusters.iloc[i]['cell_id']):
            print ("Cell ids are different for cell ", i, " predicted ", pred_clusters.iloc[i]['cell_id'], " true ", true_clusters.iloc[i]['cell_id'])
            sys.exit(2)

    print (*labels_true)
    Vmeasure = v_measure_score(labels_true, labels_pred)
    print ("Vmeasure: ", Vmeasure)

    # Now printing the number of clusters and the clone prevalence errors

    ncells = pred_clusters.shape[0]
    # print ("Number of cells: ", ncells)

    prev_true = []

    if args.true_prevalences == "None":
        print("True prevalences None, calculating")
        for n in range(ncells):
            prev_true.append(sum(np.equal(labels_true,labels_true[n]))/ncells)

    else:
        print("Using given true prevalences")
        prevs = args.true_prevalences.split("_")
        for n in range(ncells):
            prev_true.append(float(prevs[labels_true[n]-1]))

    # print(*prev_true)

    prev_pred = []
    for n in range(ncells):
        prev_pred.append(sum(np.equal(labels_pred,labels_pred[n]))/ncells)

    clone_prev_MAE = mean_absolute_error(prev_true, prev_pred)
    clone_prev_MSE = mean_squared_error(prev_true, prev_pred)


    #print ("MAE: ", clone_prev_MAE)
    #print ("MSE: ", clone_prev_MSE)

    num_true_clusters = len(set(labels_true))
    num_pred_clusters = len(set(labels_pred))

    file = open(args.results_file, 'w+')

    print("best_vmeasure\tclone_prev_MAE\tclone_prev_MSE\tnclusters_true\tnclusters_pred", file=file)
    print(Vmeasure,"\t", clone_prev_MAE, "\t", clone_prev_MSE, "\t", num_true_clusters, "\t", num_pred_clusters, file=file)
    file.close()

if __name__ == '__main__':
    main()
