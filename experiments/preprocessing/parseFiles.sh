#!/bin/sh

#  parseFiles.sh
#  
#
#  Created by Jahn  Katharina on 08/11/18.
#
#
#  This program takes the list of patients, the master table and a directory name
#  as command line parameters.
#  It copies, respecitvely creates, all input files necessary to run CTC_SCITE
#  The files are stored in the given directory which is created if it not already exists


##################################################################################
###  Function call:
###  ./parseFiles.sh patientList.txt gDNA_sample_annotation_final.tsv testFolder/ dna_samples_passing_francesco_filters_181112.txt
##################################################################################

### Command line parameters
patientListFile=$1
masterTable=$2
newDataDir=$3
sampleQualityFile=$4
mySampleQualityFile=$5
patientList=$(cut -f1 $patientListFile)   # list of patient ids



###  Step 1: create new directory and copy all data files from backup folder
if [ ! -d "$newDataDir" ]; then
mkdir -p $newDataDir
else
echo "Error: directory $newDataDir exists already" 1>&2;
#exit 1
fi

cp "../AcetoData_jahnka/dataBackup/"* "$newDataDir/"
cd $newDataDir
for f in $patientList ; do         # create subdirectory for each patient
mkdir -p $f
mv $f".txt" $f"/"$f".txt"   # move the files into the patient's directory
mv $f".vcf" $f"/"$f".vcf"
done
rm *.txt
rm *.vcf
cd ..
echo "Copied files from backup folder."


###  Step 2: remove the lines from data file that do not have read count information
###  Keep original file as follows: Br11.txt ->  Br11_orig.txt
cd $newDataDir
for f in $patientList ; do
cd $f
#echo $f
mv $f".txt" $f"_orig.txt"               # rename old file
grep "^chr" $f"_orig.txt" > $f".txt"    # store only lines startig with "chr" in original file
cd ..
done
cd ..
echo "Removed non read count lines, kept original files as Br16_orig.txt"


### Step 3: create a file for each patient that lists the samples
### (for cross-reference with the master table)
### sample names are taken from the vcf files

cd $newDataDir
for f in $patientList ; do
cd $f
out=$f
out+="_samples.tsv"
grep "#CHROM" $f.vcf | cut -f10- | tr "\t" "\n" | cut -f1 -d"." > $out
cd ..
done
cd ..
echo "Written file with list of samples"

###  Step 4: create the files describing the samples (needed as inputfile to run CTC_SCITE)
###  Also edits data files by removing samples with too many cells (>5). These components are
###  never used. Alternate files where also smaller mulit-cell components are removed are generated in Step 6.

cd $newDataDir
rm fileGeneration.log
for f in $patientList ; do
cd $f
inFile=$f"_samples.tsv"
outFile=$f"_samples_nodeDescription.tsv"
echo $outFile
python ../../getSampleFile.py $f ../../$masterTable ../../$sampleQualityFile  ../../$mySampleQualityFile >> ../fileGeneration.log
cd ..
done
cd ..
echo "created file describing samples *_samples_nodeDescription.tsv"

###  Step 5: create the input files for the three subsets of Br16: Br16_AC, Br16_B, Br16_C

cd $newDataDir
mkdir Br16_AC Br16_B Br16_C
python ../split_Br16.py Br16
cd ..
echo "Input files split based on mouse of origin"
