#!/bin/python3

import pyarrow
import pandas as pd
import argparse
import os


sample_name = 'Br16_C'
all_tuples = [["Br16_C33","Br16_C34"]]




data = pd.read_csv(os.path.join('..','input_folder',sample_name,sample_name + '.txt'), delimiter = '\t', header = None)
descriptions = pd.read_csv(os.path.join('..','input_folder', sample_name, sample_name + '_samples_nodeDescription.tsv'), delimiter = '\t', header = None)



color_palette =['orchid', 'orchid1', 'orchid2', 'orchid3', 'orchid4'] 


indices = []
for iterator,sample_list in enumerate(all_tuples):
    indices_temp = []
    wildtype_reads = [0]*data.shape[0]
    mutated_reads = wildtype_reads
    for it in sample_list:
        index = descriptions.iloc[:,0].str.contains(it)
        matched_row = descriptions[index]
        index = matched_row.index.to_list()[0]
        indices_temp.append(index)
        wildtype_reads = wildtype_reads + data.iloc[:,3+2*(index+1)-1]
        mutated_reads  = mutated_reads  + data.iloc[:,3+2*(index+1)  ]   
        
    column_name_wildtype = 'simulated' + str(iterator) + '_wildtype'
    column_name_mutated = 'simulated' + str(iterator) + '_mutated'
    data[column_name_wildtype] = wildtype_reads
    data[column_name_mutated] = mutated_reads
    
    indices.append(indices_temp)
    cluster_name = sample_name + '_sim' + str(iterator)
    newline = [cluster_name, \
               str(len(sample_list)), str(len(sample_list)), 0, \
               '[color=' + color_palette[iterator] + ',label="' + cluster_name + \
               '",fillcolor=' + color_palette[iterator] + ',image="../CTC-cluster-icons/cluster_' \
               + str(len(sample_list)) + '-0.png"]']
    descriptions.loc[len(descriptions)] = newline


descriptions.to_csv(path_or_buf = os.path.join('..','simulations', sample_name, sample_name + '_samples_nodeDescription.tsv'), sep = '\t', header = False, index = False)
data.to_csv(path_or_buf = os.path.join('..','simulations',sample_name,sample_name + '.txt'), sep = '\t', header = False, index = False)


