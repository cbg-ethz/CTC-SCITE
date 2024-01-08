import sys
import os

##############################################################################
#     This creates a new version of a case with components whose size        #
#     exceeds a given threshold are removed. For example if given            #
#     parameters:   "Br16/"   and   "3"   a new directory "Br16_max3"        #
#     is created where components with more than three cells are removed.    #
#     The original files are also kept.                                      #
##############################################################################


#os.system(cut -f2 Br16/Br16_samples_nodeDescription.tsv | sort | tail -1)



oldDir    = sys.argv[1]
origCase  = oldDir.strip().split("/")[0]
sizeLimit = int(sys.argv[2])

print("\n" + oldDir + "  sizeLimit = " + str(sizeLimit))

# original files
nodeDescriptionFileName = oldDir + origCase + "_samples_nodeDescription.tsv"
nodeDescriptionTable = [line.rstrip('\n') for line in open(nodeDescriptionFileName)]
readCountFileName = oldDir + origCase + ".txt"
readCountFile = [line.rstrip('\n') for line in open(readCountFileName)]

maxSize = 0
for line in nodeDescriptionTable:
    size = int(line.split()[2])
    maxSize = max(maxSize,size)
print("Max component size is " + str(maxSize) + "   in " + nodeDescriptionFileName)


# Test if this case needs to be processed for having large components.
# Exit the program if that is not the case
if maxSize <= sizeLimit:
    # return because there's nothing to do for this case
    print(origCase + " has no components larger than " + str(sizeLimit))
    sys.exit()

# output files
newDir = origCase + "_max" + str(sizeLimit)
os.system("mkdir " + newDir)
outFileName1 = newDir + "/" + newDir + "_samples_nodeDescription.tsv"
outFile1 = open(outFileName1, 'w')
outFileName2 = newDir + "/" + newDir + ".txt"
outFile2 = open(outFileName2, 'w')
print(outFileName1)
print(outFileName2)

# Track which lines in the node description file should be used because they belong to the right component
# and write the new node description file
usedLineIds = []
for id, line in enumerate(nodeDescriptionTable):
    size = int(line.split()[2])
    if size <= sizeLimit:
        outFile1.write(line + "\n");
        usedLineIds.append(id)



# Calculate which columns from the read count file should be used
usedCols = [0, 1, 2, 3]            # always use first four columns (mutation descriptor)
for id in usedLineIds:
    usedCols.append(4+(2*id))         # get ids of columns corresponding to samples from this component
    usedCols.append(4+(2*id)+1)
print(usedCols)

# write the read count file using only columns relevant for the component
for line in readCountFile:
    outLine = ""
    columnList = line.strip().split("\t")
    for index, col in enumerate(columnList):
        if index in usedCols:
            outLine += col + "\t"
    outFile2.write(outLine.strip() + "\n")

outFile1.close()
outFile2.close()



# Br16_AC1	3	3	0	[color=lightcoral,label="Br16_AC1",fillcolor=lightcoral,image="../CTC-cluster-icons/cluster_3-0.png"]
