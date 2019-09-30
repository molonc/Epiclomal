#!/usr/bin/env python

#=======================================================================================================================
# Epiclomal, probabilistic model for clustering and imputing single-cell methylation data
# Author : Mirela Andronescu
#=======================================================================================================================
from lib.run import run_basic_gemm_model, run_basic_bayespy_model, run_region_gemm_model

import argparse

def main():
    parser = argparse.ArgumentParser(prog='Epiclomal')

    parser.add_argument('--version', action='version', version='0.0.1')

    subparsers = parser.add_subparsers()


    def str2bool(v):
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')

    #=======================================================================================================================
    # Shared analysis parser
    #=======================================================================================================================
    analysis_parser = argparse.ArgumentParser(add_help=False)

    analysis_parser.add_argument('--K', default=None,
                                help='''Max number of clusters. If given, use this one instead of the number from the config file.''')

    analysis_parser.add_argument('--config_file', required=True,
                                help='''Path to YAML format configuration file.''')

    analysis_parser.add_argument('--methylation_file', required=True,
                                help='''Path to methylation data input file.''')

    analysis_parser.add_argument('--copynumber_file', default=None,
                                help='''Path to copy number input file.''')

    analysis_parser.add_argument('--regions_file', default=None,
                                help='''Path to regions input file.''')

    analysis_parser.add_argument('--initial_clusters_file', default=None,
                                help='''Start from these clusters instead of random clusters.''')

    analysis_parser.add_argument('--true_clusters_file', default=None,
                                help='''Path to the true_clusters_file, if known. If given, params.yaml will contain the V-measure for this prediction.''')

    analysis_parser.add_argument('--true_prevalences', default=None,
                                help=''' A string with the true prevalences for all the clusters, e.g. 0.33_0.33_0.34''')

    analysis_parser.add_argument('--repeat_id', default=1, type=int,
                                help='''A number >= 0. If there is a column with this number (excluding the first column and starting from 0) in the initial_clusters_file, use that column as initial clusters, else use random initialization.''')

    # Not used any more
    analysis_parser.add_argument('--bulk_file', default=None,
                                help='''A file with 3 columns: locus, #methylated reads, #unmethylated reads. The beta prior will be initialized with these values''')

    analysis_parser.add_argument('--slsbulk_file', default=None,
                                help='''A file with 3 columns: locus, #methylated reads, #unmethylated reads. This will be used to perform an SLS search that optimizes a bulk satisfaction score based on perturbation of uncertain cells.''')

    analysis_parser.add_argument('--slsbulk_iterations', default=10,
                                help='''The number of iterations for the SLSbulk procedure.''')

    analysis_parser.add_argument('--out_dir', default=None,
                                help='''Path where output files will be written.''')

    analysis_parser.add_argument('--mu_has_k', type=str2bool, default=True,
                                help='''True or False depending on whether we want mu to depend on k or not''')

    analysis_parser.add_argument('--convergence_tolerance', default=1e-4, type=float)

    analysis_parser.add_argument('--max_num_iters', default=1000, type=int)

    analysis_parser.add_argument('--seed', default=None, type=int,
                                 help='''Set random seed so results can be reproduced. By default a random seed is
                                 chosen.''')

    analysis_parser.add_argument('--labels_file', default=None,
                                 help='''Path of file with initial labels to use.''')

    analysis_parser.add_argument('--Bishop_model_selection', type=str2bool, default=False,
                                help='''True or False depending on whether we want to apply Corduneanu_Bishop model selection''')

    analysis_parser.add_argument('--check_uncertainty', type=str2bool, default=False,
                                help='''True or False depending on whether we want to check whether the uncertainty is estimated correctly''')


    #----------------------------------------------------------------------------------------------------------------------
    basic_parser = subparsers.add_parser('Basic-GeMM', parents=[analysis_parser],
                                           help='''Analyse single cell methylation data using the Basic-GeMM model.''')

    basic_parser.set_defaults(func=run_basic_gemm_model)

    #----------------------------------------------------------------------------------------------------------------------
    basic_parser = subparsers.add_parser('Basic-BayesPy', parents=[analysis_parser],
                                           help='''Analyse single cell methylation data using the Basic model with BayesPy.''')

    basic_parser.set_defaults(func=run_basic_bayespy_model)

    #----------------------------------------------------------------------------------------------------------------------
    region_parser = subparsers.add_parser('Region-GeMM', parents=[analysis_parser],
                                           help='''Analyse single cell methylation data using the Region-GeMM model.''')

    region_parser.set_defaults(func=run_region_gemm_model)


    #----------------------------------------------------------------------------------------------------------------------
    # the regions parser will add more input arguments
    # doublet_parser = subparsers.add_parser('run_doublet_model', parents=[analysis_parser, genotyper_parser],
    #                                        help='''Analyse single cell data using the single cell genotyper model accounting for doublets.''')
    # doublet_parser.add_argument('--state_map_file', required=True)
    # doublet_parser.set_defaults(func=run_doublet_analysis)


    args = parser.parse_args()

    args.func(args)

if __name__ == '__main__':
    print("in main")
    main()
