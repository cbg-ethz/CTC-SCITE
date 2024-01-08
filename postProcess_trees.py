import sys
import os
import argparse


parser = argparse.ArgumentParser(
                    prog='postProcess_trees',
                    description='Posterior samplings are concatenated and splitting statistics computed')
parser.add_argument('sample')  
args = parser.parse_args()

PATH_TO_CPP_PROGRAM  =  "CTC_parseSampling/CTC_parseSampling" # path to C++ parser



os.system("cat ../input_folder/" + args.sample + "/" + "*postSampling.tsv > " + "../input_folder/" + args.sample + "/" +  "_postSampling.tsv")
path = '../input_folder/' + args.sample + "/"
command = "sbatch --time=24:00:00 "
command += PATH_TO_CPP_PROGRAM
command += ' -i ../input_folder/' + args.input  + "/" +  elem + '.txt'
command += ' -description ../input_folder/' + args.input + "/" +  elem + '_samples_nodeDescription.tsv'
command += ' -samples ../input_folder/' + args.input + '_postSampling.tsv'
os.system(command)
print command
os.system("rm " + elem + "/temp.tsv")
