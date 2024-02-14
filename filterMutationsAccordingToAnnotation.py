#!/bin/python3

import pandas as pd
import os
import argparse


parser = argparse.ArgumentParser(description='Only keep meeningful mutations for tree inference')
parser.add_argument('i', type =str, help='Name of the tree for which to filter mutations')

args = parser.parse_args()


sampleName = args.i


data = pd.read_csv(os.path.join('..','input_folder',sampleName, sampleName + '.txt'), delimiter = '\t', header= None)


data = data.rename(columns={0:'#CHROM', 1:'POS', 2:'REF',3: 'ALT'})


with open(os.path.join('..','input_folder','filtered', 'vcf_files_annotated', sampleName + '.ann.vcf')) as file:
    for line in file:
        if line.startswith('#CHROM'):
            vcf_names = [x for  x in line.split('\t')]
            break
file.close()

vcf = pd.read_csv(os.path.join('..','input_folder','filtered', 'vcf_files_annotated', sampleName + '.ann.vcf'), comment = '#', sep=r'\s+', header=None, names=vcf_names)




include = []
for it, line in enumerate(vcf['INFO']):
    functionalAnnotation = line.split(';')[1].split(',')
    impact = False
    for ann in functionalAnnotation:
        if(ann.split('|')[2] in ['MODERATE','HIGH']):
            impact = True
    include.append(impact)
print(include)

data.head()



includeFunctionalAnnotation = []
for index, line in data.iterrows():
    if(vcf.loc[(vcf['#CHROM'] == line['#CHROM']) & (vcf['POS'] == line['POS'])].shape[0] != 1):
        print(vcf.loc[(vcf['#CHROM'] == line[0]) & (vcf['POS'] == line[1])])
        print('More than one hit in the annotation file. ERROR')
    else:
        includeFunctionalAnnotation.append(include[vcf.loc[(vcf['#CHROM'] == line['#CHROM']) & (vcf['POS'] == line['POS'])].index[0]])
print(includeFunctionalAnnotation)


cgi = pd.read_csv('../input_folder/filtered/CGI/LM2_cgi/alterations.tsv', sep='\t')



print(cgi['CGI-Oncogenic Summary'].unique())
print(cgi['CGI-External oncogenic annotation'].unique())
print(cgi['CGI-Oncogenic Prediction'].unique())


cgiFiltered = cgi.loc[cgi['CGI-Oncogenic Prediction'].isin(['driver (oncodriveMUT)', 'oncogenic (predicted)']),:]



includeDriverAnnotation = [False]*data.shape[0]
for index, cgiLine in cgiFiltered.iterrows():
    if(data.loc[(('chr' + cgiLine['CHROMOSOME']) == data['#CHROM']) & (cgiLine['POSITION'] == data['POS'])].shape[0] != 1):
        print('Annotated genomic position not found in the dataset!')
    else:
        includeDriverAnnotation[data.loc[(('chr' + cgiLine['CHROMOSOME']) == data['#CHROM']) & (cgiLine['POSITION'] == data['POS'])].index[0]] = True
print(includeDriverAnnotation)



#This checks for each mutation sites whether at least one of the samples is mutated. This implies that at least one of the samples has reads for that mutation.
includeIsMutated = (data.iloc[:, 5::2] > 0).sum(axis = 1) > 0


includeFinal = [(a or b) and c for a, b,c in zip(includeDriverAnnotation, includeFunctionalAnnotation, includeIsMutated)]



data.loc[includeFinal,:].to_csv(os.path.join('..', 'input_folder','filtered', sampleName, sampleName + '.txt'), sep = '\t', header = False, index = False)

