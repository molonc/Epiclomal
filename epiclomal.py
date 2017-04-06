#!/usr/bin/env python

#=======================================================================================================================
# Epicaller, probabilistic model for clustering and imputing single-cell methylation data
# Author : Mirela Andronescu
#=======================================================================================================================
from lib.run import run_basic_gemm_model, run_basic_miss_gemm_model, run_region_gemm_model, run_region_miss_gemm_model

import argparse

parser = argparse.ArgumentParser(prog='Epiclomal')

parser.add_argument('--version', action='version', version='0.0.1')

subparsers = parser.add_subparsers()

#=======================================================================================================================
# Shared analysis parser
#=======================================================================================================================
analysis_parser = argparse.ArgumentParser(add_help=False)

analysis_parser.add_argument('--K', default=None,
                            help='''Max number of clusters. If given, use this one instead of the number from the config file.''')

analysis_parser.add_argument('--config_file', required=True,
                            help='''Path to YAML format configuration file.''')

analysis_parser.add_argument('--data_file', required=True,
                            help='''Path to methylation data input file.''')
                            
analysis_parser.add_argument('--copynumber_file', default=None,
                            help='''Path to copy number input file.''')                            
                            
analysis_parser.add_argument('--regions_file', default=None,
                            help='''Path to regions input file.''')                            

analysis_parser.add_argument('--initial_clusters_file', default=None,
                            help='''Start from these clusters instead of random clusters.''')                            

analysis_parser.add_argument('--lower_bound_file', default=None,
                             help='''Path of file where lower bound convergence will be written.''')

analysis_parser.add_argument('--out_dir', default=None,
                            help='''Path where output files will be written.''')

analysis_parser.add_argument('--convergence_tolerance', default=1e-4, type=float)

analysis_parser.add_argument('--max_num_iters', default=1000, type=int)

analysis_parser.add_argument('--seed', default=None, type=int,
                             help='''Set random seed so results can be reproduced. By default a random seed is
                             chosen.''')

analysis_parser.add_argument('--labels_file', default=None,
                             help='''Path of file with initial labels to use.''')

#---------------------------------------------------------------------------------------------------------------------- 
basic_parser = subparsers.add_parser('Basic-GeMM', parents=[analysis_parser],
                                       help='''Analyse single cell methylation data using the Basic-GeMM model.''')

basic_parser.set_defaults(func=run_basic_gemm_model)

#---------------------------------------------------------------------------------------------------------------------- 

basic_miss_parser = subparsers.add_parser('Basic-Miss-GeMM', parents=[analysis_parser],
                                       help='''Analyse single cell methylation data using the Basic-Miss-GeMM model.''')

basic_miss_parser.set_defaults(func=run_basic_miss_gemm_model)

#---------------------------------------------------------------------------------------------------------------------- 
region_parser = subparsers.add_parser('Region-GeMM', parents=[analysis_parser],
                                       help='''Analyse single cell methylation data using the Region-GeMM model.''')

region_parser.set_defaults(func=run_region_gemm_model)

#---------------------------------------------------------------------------------------------------------------------- 
region_parser = subparsers.add_parser('Region-Miss-GeMM', parents=[analysis_parser],
                                       help='''Analyse single cell methylation data using the Region-Miss-GeMM model.''')

region_parser.set_defaults(func=run_region_miss_gemm_model)


#---------------------------------------------------------------------------------------------------------------------- 
# the regions parser will add more input arguments
# doublet_parser = subparsers.add_parser('run_doublet_model', parents=[analysis_parser, genotyper_parser],
#                                        help='''Analyse single cell data using the single cell genotyper model accounting for doublets.''')
# doublet_parser.add_argument('--state_map_file', required=True)
# doublet_parser.set_defaults(func=run_doublet_analysis)


args = parser.parse_args()

args.func(args)
