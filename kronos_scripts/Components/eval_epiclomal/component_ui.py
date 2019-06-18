# -*- coding: utf-8 -*-
"""
Created on Mar 7, 2017

@author: Mirela Andronescu (mandronescu@bccrc.ca)
"""
import argparse


#==============================================================================
# make a UI
#==============================================================================
parser = argparse.ArgumentParser(prog='eval_epiclomal',
                                 description=''' Evaluate epiclomal results''')


parser.add_argument('--input_dir', required=True,
                    help='''Directory with epiclomal results''')

parser.add_argument('--output_dir', required=True,
                    help='''Directory with the output results file''')

parser.add_argument('--model_name', required=True,
                    help='''A description of the model evaluated, such as basic or region''')

parser.add_argument('--calc_vmeasure', 
                    help='''Full path of the software that calculates vmeasure''')                    

parser.add_argument('--calc_hdist', 
                    help='''Full path of software that calculates hamming distance''') 

parser.add_argument('--true_clusters_file', 
                    help='''File with the true clusters, if known''')                     
  
parser.add_argument('--true_epigenotypes_file', 
                    help='''File with the true epigenotypes, if known''')                 
                    
args, unknown = parser.parse_known_args()
