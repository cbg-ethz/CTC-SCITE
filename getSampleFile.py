import sys
import os

##########################################################################
#     This creates the file describing the samples of a patient.         #
#     Parameters:                                                        #
#        * patient id                                                    #
#        * path to the master table                                      #
#     Requirements:                                                      #
#        * file listing patient samples exists in the directory          #
#        * master table is located in parent directory                   #
##########################################################################

### command line parameters
patient = sys.argv[1]
masterTableFileName = sys.argv[2]    # with local path (should be in parent directory)
qualitySamplesFileName = sys.argv[3]
notAll_NA_samplesFileNAme = sys.argv[4]
print("##############       " + patient + "      ############")

### input/output files
sampleFileName=patient+"_samples.tsv"
outFileName=patient+"_samples_nodeDescription.tsv"
sampleFile = open(sampleFileName, 'r')
outFile = open(outFileName, 'w')
masterTable = [line.rstrip('\n') for line in open(masterTableFileName)]


### tracking low-quality samples and samples with too many cells

qualitySamples = [line.rstrip('\n') for line in open(qualitySamplesFileName)]    # list of high-quality samples according to Francesco
notAll_NA_Samples = [line.rstrip('\n') for line in open(notAll_NA_samplesFileNAme)]    # list of high-quality samples according to Francesco
sampleSizeThreshold = 5                                                        # larger samples should be discarded, maybe make this a parameter
removeSampleList = []                          # store here samples that should be removed for having too many cells

### color scheme to represent single cells and samples from the same cluster
w_color = "ghostwhite"
t_color = "gray93"
colors = ["lightcoral", "sandybrown", "skyblue3", "thistle", "lemonchiffon", "violetred3", "lightslateblue", "paleturquoise3", "khaki3", "darkseagreen4", "gold", "plum", "yellowgreen", "navajowhite2", "crimson", "cadetblue", "darkslategray", "deeppink2", "lightpink2", "orangered4", "peachpuff1", "yellow4", "sienna2", "rosybrown4", "springgreen", "palegreen3", "orangered", "mediumorchid4", "burlywood4", "tan", "brown4", "yellow", "wheat", "lawngreen", "indianred", "turquoise4", "tomato", "red", "blue", "midnightblue", "forestgreen", "deeppink", "deepskyblue", "darkorange", "magenta", "deeppink", "mediumaquamarine", "mistyrose", "powderblue", "steelblue", "greenyellow", "honeydew1", "palegoldenrod", "hotpink", "goldenrod", "firebrick", "dodgerblue4", "limegreen", "khaki4", "mistyrose4", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black"]


### dictionary: clusters to color

#print "Split cluster colors:"
clusterColor = {}
colorCounter = 0
for line in masterTable:
    lineSplit = line.strip().split("\t")
    if lineSplit[1]==patient:                     # line (sample) belongs to patient (1: Donor)
        if lineSplit[13] != 'NA':                    # line (sample) comes from cluster (13: clusterID)
            if not lineSplit[13] in clusterColor:       # cluster not yet in dictionary
                clusterColor[lineSplit[13]] = colors[colorCounter]
                colorCounter += 1
#            print lineSplit[13] + ": " + clusterColor[lineSplit[13]]
print( str(colorCounter) + " split cluster(s)")


### create the file describing the samples
sampleCounter = 0
usedSampleCounter = 0
countErrorCounter = 0
typeErrorCounter = 0
largeComponents = 0
tooMuch_NA_Counter = 0
notFoundCounter = 0

print("processing samples...")
for sample in sampleFile:
    sample = sample.strip()
    sampleCounter += 1
    tooManyCells = False       # tag samples with too many cells for removal; they slow down computation and have little information
    tooMuch_NA = False
    sampleQualityIcon = ""       # samples with low quality according to Francesco get marked with a '*'
    
    if sample not in qualitySamples:
        sampleQualityIcon = "*"
    
    
    found = False
    for line in masterTable:
        lineSplit = line.strip().split("\t")
        if lineSplit[0]==sample:
            found = True
    
            # check if sample has too many NA or is large component, print warning
            if sample not in notAll_NA_Samples:
                tooMuch_NA = True
                tooMuch_NA_Counter += 1
                removeSampleList.append(sampleCounter)
            elif int(lineSplit[9]) > 5:
                tooManyCells = True
                removeSampleList.append(sampleCounter)
                print("WARNING LARGE COMPONENT:      " + lineSplit[9] + " cells")
                print(line)
                largeComponents += 1
            
            
            # correct wbc and tumour cell counts in case numbers don't add up
            if int(lineSplit[9]) != int(lineSplit[10]) + int(lineSplit[11]):
                #print "CELL COUNT Correction: "
                #print line
                countErrorCounter += 1
                if lineSplit[8]=="WBC":
                    lineSplit[10] = "0"
                    lineSplit[11] = lineSplit[9]
                else:
                    lineSplit[10] = lineSplit[9]
                    lineSplit[11] = "0"

            # correct cell types based on last column of master table
            if lineSplit[8] != lineSplit[20]:
                #print "CELL TYPE CORRECTION:"
                #print line
                typeErrorCounter += 1
                if lineSplit[20]=="WBC":
                    lineSplit[10] = "0"
                    lineSplit[11] = lineSplit[9]
                else:
                    lineSplit[10] = lineSplit[9]
                    lineSplit[11] = "0"
        
            # write line describing sample for output file
            picture = "cluster_" + lineSplit[10] + "-" + lineSplit[11] + ".png"
            nodeColor = "gray93"
            if lineSplit[13] in clusterColor:
                nodeColor = clusterColor[lineSplit[13]]
            elif int(lineSplit[9])>1:
                nodeColor = colors[colorCounter]
                print(colors[colorCounter])
                colorCounter += 1
            elif int(lineSplit[9])==1 and int(lineSplit[11]) == 1:
                nodeColor = "ghostwhite"
            else:
                "UNCAUGHT CASE: " + line

            # Example node description for a single tumour cell
            # '[color=gray93,label="Br57_CTC_10",fillcolor=gray93,image="/Users/jahnka/Desktop/AcetoData/CTC-cluster-icons/cluster_1-0.png"
            node =  "[color="+ nodeColor + ",label=\"" + sample + sampleQualityIcon + "\",fillcolor=" + nodeColor
            node += ",image=\"../CTC-cluster-icons/" + picture + "\"" + "]"
            #print node

            newLine = sample + '\t' + lineSplit[9] + '\t' + lineSplit[10] + '\t' + lineSplit[11] + '\t' + node
            if not tooManyCells and not tooMuch_NA:                  # only write line if the samples is not discarded for having too many cells
                usedSampleCounter += 1
                outFile.write(newLine + "\n")
    if not found:
        notFoundCounter += 1

# remove columns from the data file that represent discarded samples (e.g. Br16.txt)
if len(removeSampleList) != 0:
    dataFileName = patient + ".txt"
    dataFileNameEdited = patient + "_edited.txt"
    f = open(dataFileName, "r")
    g = open(dataFileNameEdited, "w")
    
    offset = 4
    removeCols = []
    for x in removeSampleList:
        removeCols.append(offset+2*(x-1))
        removeCols.append(offset+2*(x-1)+1)
    for line in f:
        lineSplit = line.strip().split("\t")
        newLine = ""
        for index, item in enumerate(lineSplit):
            if index not in removeCols:
                if index != 0:
                    newLine += "\t"
                newLine += item
        g.write(newLine + "\n")
    os.rename(dataFileNameEdited, dataFileName)
    f.close()
    g.close()

### print summary of processing issues
print( str(sampleCounter) + "\tsamples processed.")
print( str(usedSampleCounter) + "\tsamples written to file.")
print(str(countErrorCounter) + "\tsample(s) needed cell count correction")
print(str(typeErrorCounter) + "\tsample(s) needed cell type correction")
if largeComponents > 0:
    print( str(largeComponents) + "\tsample(s) discarded for having more than " + str(sampleSizeThreshold) + " cells")
if notFoundCounter > 0:
    print( str(notFoundCounter) + "\tsample(s) omitted because not in master table")
if tooMuch_NA_Counter > 0:
    print( str(tooMuch_NA_Counter) + "\tsample(s) omitted because they have too many NAs")
print(removeSampleList)

