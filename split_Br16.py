import sys
import os

##############################################################################
#     This splits the Br16 files into separate files based on the cells'     #
#     orgin (different mice). The following files are created:               #
#       * Br16_AC/Br16_AC.txt                                                #
#       * Br16_AC/Br16_AC_samples_nodeDescription.tsv                        #
#       * Br16_B/Br16_B.txt                                                  #
#       * Br16_B/Br16_B_samples_nodeDescription.tsv                          #
#       * Br16_C/Br16_C.txt                                                  #
#       * Br16_C/Br16_C_samples_nodeDescription.tsv                          #
#     The original files are also kept.                                      #
##############################################################################


origCase = sys.argv[1]
componentList = ["_AC", "_B", "_C"]

# original files
nodeDescriptionFileName = origCase + "/" + origCase + "_samples_nodeDescription.tsv"
nodeDescriptionFile = [line.rstrip('\n') for line in open(nodeDescriptionFileName)]
readCountFileName = origCase + "/" + origCase + ".txt"
readCountFile = [line.rstrip('\n') for line in open(readCountFileName)]


# write files for the individual components
for component in componentList:
    outFileName1 = origCase + component + "/" + origCase + component+"_samples_nodeDescription.tsv"
    outFile1 = open(outFileName1, 'w')
    outFileName2 = origCase + component + "/" + origCase + component+".txt"
    outFile2 = open(outFileName2, 'w')
    allCols = []
    
    # Track which lines in the node description file should be used because they belong to the right component
    usedLineIds = []
    for id, line in enumerate(nodeDescriptionFile):
        allCols.append(id)
        if line.startswith(origCase + component):
            outFile1.write(line + "\n");
            #print line
            usedLineIds.append(id)


    # Find out which columns from the read count file should be used
    usedCols = [0, 1, 2, 3]            # always use first four columns (mutation descriptor)
    for id in usedLineIds:
        usedCols.append(4+(2*id))         # get ids of columns corresponding to samples from this component
        usedCols.append(4+(2*id)+1)


    # write the read count file using only columns relevant for the component
    for line in readCountFile:
        outLine = ""
        columnList = line.strip().split("\t")
        for index, col in enumerate(columnList):
            if index in usedCols:
                outLine += col + "\t"
        outFile2.write(outLine.strip() + "\n")
#    print "\nColumns used for " + outFileName2 + ":"
#    print usedCols

    outFile1.close()
    outFile2.close()



# Br16_AC1	3	3	0	[color=lightcoral,label="Br16_AC1",fillcolor=lightcoral,image="../CTC-cluster-icons/cluster_3-0.png"]