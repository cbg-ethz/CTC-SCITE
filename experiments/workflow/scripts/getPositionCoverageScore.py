import sys
from collections import Counter


inFileName = sys.argv[1]            # the data table
inFile = open(inFileName, 'r')

outFileName = inFileName.replace(".txt", "_covScore.txt")
outFile = open(outFileName, 'w')

mutCounter = 0
threshold = 3                         # need to have at least this coverage to be counted as callable
cellCount = 0

headerLine = "\tcovScore\n"
outFile.write(headerLine)

for i, line in enumerate(inFile):          # iterate over data file line by line
    lineSplit = line.strip().split("\t")   # split line by tab
    cellCount = 0.5*(len(lineSplit)-4)
    mutCounter += 1
    mut = lineSplit[0] + "_" + lineSplit[1] #+ "_" + lineSplit[2] + lineSplit[3]
    
    covCounter = 0
    
    for index in range(4,len(lineSplit)):    # remaining columns are observed states of this line's mutation
        
        if index%2 != 0:                     # use only every other column, as we want the total read count
            continue
        
        if int(lineSplit[index]) >= threshold:    # if coverage is high enough (callable site), increase counter
            covCounter+=1

        #print lineSplit[index] + "\t" + lineSplit[index+1]

    #print str(covCounter) + "/" + str(cellCount)
    score = float(covCounter)/float(cellCount)
    outputLine = mut + "\t" + str(score) + "\n"
    print(outputLine)
    outFile.write(outputLine)


print(mutCounter)
print(cellCount)


inFile.close()
outFile.close()
