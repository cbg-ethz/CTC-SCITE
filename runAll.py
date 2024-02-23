
import sys
import os
import glob
import random
import argparse



PATH_TO_EXECUTABLE  =  "./CTC_SCITE_katharina/CTC_SCITE" # path to CTC-SCITE executable



lengths =[1000000] # [500000] #, 1000000]
lengthNames = {'500000': '500K', '1000000': '1M'}
pathName = "../simulations"
caseList = ['Br11']#, 'Br16_AC_max2',  'Br16_AC_max3',  'Br16_AC_max4', 'Br16_B_max2', 'Br16_B_max3',  'Br16_B_max4', 'Br16_C_max1',  'Br16_C_max2', 'Br16_C_max3', 'Br23',  'Br26',  'Br30',  'Br38',  'Br39',  'Br57', 'Brx50', 'Lu2',  'Lu7', 'Ov8', 'Pr6', 'Br16_AC', 'Br16_B', 'Br16_C', 'Pr9', 'LM2', 'Br61', 'Br7']
#caseList = ['Ov8']


baseList = [os.path.join(pathName, x, x) for x in caseList]
#print baseList



for i in [1,2,3,4,5,6,7,8,9,10]:
    for l in lengths:
        for base in baseList:
            seed = random.randint(0, 32767)
            run_command ="sbatch --time=96:00:00  --mem-per-cpu=20G --wrap='"
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
