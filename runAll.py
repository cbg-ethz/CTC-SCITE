
import sys
import os
import glob
import random
import argparse



PATH_TO_EXECUTABLE  =  "./CTC_SCITE_katharina/CTC_SCITE" # path to CTC-SCITE executable



lengths =[1000000] # [500000] #, 1000000]
lengthNames = {'500000': '500K', '1000000': '1M'}
pathName = "../input_folder/"
caseList = ['Br7', 'Pr9', 'Br61', 'LM2']


baseList = [pathName + x + "/" + x for x in caseList]
#print baseList



for i in [1,2,3,4,5,6,7,8,9,10]:
    for l in lengths:
        for base in baseList:
            seed = random.randint(0, 32767)
            run_command ="sbatch --time=24:00:00 --mem-per-cpu=20G --wrap='"
            run_command += PATH_TO_EXECUTABLE + " "
            run_command += "-i "   + base   + ".txt" + " "
            run_command += "-r 1 "
            run_command += "-l "  + str(l)  + " "
            run_command += "-g 1 "
            run_command += "-seed "   + str(seed) + " "
            run_command += "-e 0.2 "
            run_command += "-p " + str(l/5000)+ " " + str(int(0.2*l)) + " "
            run_command += "-samples " + base + "_samples_nodeDescription.tsv" + " "
            run_command += "-o " + base + "_" + lengthNames[str(l)] + "_" + str(i) + "_seed" + str(seed) + "'"
            os.system(run_command)
            print(run_command)
