#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import operator
from sklearn.metrics.cluster import v_measure_score, homogeneity_score, completeness_score

parser = argparse.ArgumentParser(description='Evaluate the clustering v-measure score.')
parser.add_argument('--true_clusters_file', required=True,
                    help='The file with the true clustering assignment')
parser.add_argument('--epiclomal_clusters_file', required=True,
                    help='The file with the epiclomal predicted clusters')

args = parser.parse_args()
# print(args.true_clusters_file)
# print(args.epiclomal_clusters_file)


# first read the predicted file and figure out how many clusters there are
posterior_clusters = pd.read_csv(args.epiclomal_clusters_file, compression='gzip', index_col='cell_id', sep='\t')

labels_pred = []
# traverse every row, and find out which cluster has the maximum value
for index, row in posterior_clusters.iterrows():
    max_index, max_value = max(enumerate(row), key=operator.itemgetter(1))
    labels_pred.append(max_index)
# print labels_pred

# now read the true clusters file
true_clusters = pd.read_csv(args.true_clusters_file, compression='gzip', index_col='cell_id', sep='\t')
labels_true = np.array(true_clusters['epigenotype_id'])
# print labels_true

print '{0:.5g}'.format(v_measure_score(labels_true, labels_pred))
