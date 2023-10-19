/*
 * findBestTrees_noR.cpp
 *
 *  Created on: Mar 27, 2015
 *      Author: jahnka
 */

#include <stdbool.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>

#include "doublets.h"
#include "matrices.h"
#include "treelist.h"
#include "trees.h"
#include "output.h"
#include "mcmc.h"
#include "rand.h"
#include "scoreTree.h"
#include "binTree_output.h"
#include "CTC_treeScoring.h"

using namespace std;

int** getDataMatrix(int n, int m, string fileName);
double* getErrorRatesArray(double fd, double ad1, double ad2, double cc);
int readParameters(int argc, char* argv[]);
string getOutputFilePrefix(string fileName, string outFile);
string getFileName(string prefix, string ending);
string getFileName2(int i, string prefix, string ending, char scoreType);
vector<string> getGeneNames(string fileName, int nOrig);
vector<string> getSampleNames(string fileName, int count);
vector<double> setMoveProbs();
int* getParentVectorFromGVfile(string fileName, int n);
int getMinDist(int* trueVector, std::vector<bool**> optimalTrees, int n);
void printGeneFrequencies(int** dataMatrix, int n, int m, vector<string> geneNames);
vector<string> getGeneNamesRecMut(string fileName, int nOrig, bool recMutAllowed, int recMut);
vector<string> readMutIds(string fileName);

vector<vector<int> > readTableFile(string fileName);
vector<vector<string> > readStringTableFile(string fileName);
vector<vector<int> > getMutReadCount(vector<vector<int> >& readCounts);
vector<vector<int> > getTotalReadCount(vector<vector<int> >& readCounts);
void filterLowReadCountSites(vector<vector<int> >& readCounts, vector<string>& mutIds, int minCellCount, int minCov);
vector<string> getLeafLabels(vector<vector<string> >& sampleInfo);
vector<int> getClusterIdOfLeaf(vector<vector<string> >& sampleInfo);
vector<int> getAlleleCount(vector<vector<string> >& sampleInfo);
vector<bool> getWbcStatus(vector<vector<string> >& sampleInfo);

double defaultMoveProbsBin[] = {0.4, 0.6};      // default tree move probs: prune&re-attach / swap leaf labels
double errorRateMove = 0.0;
vector<double> treeMoves;

string fileName;             // data file
string outFile;              // the name of the outputfile, only the prefix before the dot
int n;                  // TODO this should be renamed    number of genes
int m;                  // TODO this should be renamed     number of samples
char scoreType = 'm';   // TODO there should be no choice, remove parameter

/**** MCMC  ****/
bool useFixedSeed = false;      // use a predefined seed for the random number generator
unsigned int fixedSeed = 1;   // default seed
int rep;                      // number of repetitions of the MCMC
int loops;                    // number of loops within a MCMC
double gamma_ = 1.0;

/**** sampling ****/
bool sample = false;
int sampleStep;
int burnin;

string sampleNameFile;
string geneNameFile;              // file where the gene names are listed.


//char treeType = 'm';        // the default tree is a mutation tree; other option is 't' for (transposed case), where we have a binary leaf-labeled tree
//int maxTreeListSize = -1;  // defines the maximum size of the list of optimal trees, default -1 means no restriction
//bool doubletModel = false;
//double doubletProb = 0.0;
//bool fixedDoubletProb = false;
//int recMut = -1;                // default case no recurrent mutation, set to number of doublet mutation if applicable
//bool recMutAllowed = false;     // no recurrent mutation allowed
//int z = 0;                    //100000000;     // number of positions sequenced with no mutation called in any of the samples
double alphaMoveRatio = 0.0;
int sampleCount = 0;
int treeRankSize = 10; // Remember the 10 best trees

vector<string> mutLabels;

int main(int argc, char* argv[])
{

	initRand();
	/****************   begin timing  *********************/
			clock_t begin=clock();
	/****************************************************/

	//std::vector<struct treeTheta> optimalTrees;        // list of optimal tree/beta combinations found by MCMC

	/**  read parameters and data file  **/
	readParameters(argc, argv);
	vector<vector<int> > readCounts;            // the original read count table
	vector<vector<string> > sampleInfo;         // the original table of the samples (ids, cell count, type,...)

	vector<vector<int> > mutReadCounts;          // entry[i][j] has mutated read count for position i in sample j
	vector<vector<int> > totalReadCounts;        // entry[i][j] has total read count for position i in sample j

	vector<string> leafLabels;
	vector<int> leafClusterId;
	vector<int> alleleCount;
	vector<bool> wbcStatus;               // wbcStatus[i] is 1 if leaf i represents a wbc, 0 otherwise

	// get variant and total read count
	cout << fileName << endl;
	readCounts = readTableFile(fileName);
	mutLabels  = readMutIds(fileName);
	cout << "here" << endl;
	string prefix = getOutputFilePrefix(fileName, outFile);
	cout << prefix << endl;

	for(int i=0; i<3; i++){
		cout << "read count line: ";
		for(int j=0; j<readCounts.at(i).size(); j++){
			cout << readCounts.at(i).at(j) << "\t";
		}
		cout << endl;
	}

	cout << "before: " << readCounts.size() << endl;
	//filterLowReadCountSites(readCounts, mutLabels, 2, 3);             // filter out positions that are likely sequencing errors
	cout << "after:  " << readCounts.size() << endl;
	mutReadCounts = getMutReadCount(readCounts);            // count of mutated reads after filtering per position per sample
	totalReadCounts = getTotalReadCount(readCounts);        // total read count after filtering per position per sample
	cout << "sample count = " << mutReadCounts.at(0).size()<< endl;
	// get leaf labels and sample sizes
	cout << sampleNameFile << endl;
	sampleInfo = readStringTableFile(sampleNameFile);   // the original table of the samples (ids and cell count)
	sampleCount = sampleInfo.size();
	cout << sampleInfo.size() << endl;
	for(int i=0; i<sampleInfo.size(); i++){
		cout << i << ": ";
		for(int j=0; j<2; j++){
			cout << sampleInfo[i][j] << "\t";
		}
		cout << endl;
	}

	cout << "---------------------------" << endl;
	leafLabels = getLeafLabels(sampleInfo);        // list of leaf label names
	for(int i=0; i<leafLabels.size(); i++){
		cout << leafLabels.at(i) << endl;
	}
	cout << endl;
	//getchar();

	leafClusterId = getClusterIdOfLeaf(sampleInfo);   // maps leaf id to cluster id
	for(int i=0; i<leafClusterId.size(); i++){
		cout << leafClusterId.at(i) << endl;
	}
	m = leafClusterId.size();                 // number of leafs in tree
	n = mutReadCounts.size();                 // number of mutations in tree
	cout << m << " samples" << endl;
	cout << n << " muts" << endl;

	wbcStatus = getWbcStatus(sampleInfo);

	for(int i=0; i<wbcStatus.size(); i++){
		cout << wbcStatus.at(i) << "\t";
	}
	cout << endl;

	for(int i=0; i<leafClusterId.size(); i++){
		cout << leafClusterId.at(i) << "\t";
	}
	cout << endl;
	//getchar();

	vector<string> label;                           // the labels of the tree nodes
	for(int i=0; i<leafLabels.size(); i++){
		//cout << leafLabels.at(i) << endl;
		label.push_back(leafLabels.at(i));
	}
	for(int i=leafLabels.size(); i<(2*m)-1; i++){
		label.push_back("x"+to_string(i));
	}

	cout << "tree labels:\n";
	for(int i=0; i<label.size(); i++){
		cout << label.at(i) << endl;
	}
	cout << "end of tree labels:\n";

	alleleCount = getAlleleCount(sampleInfo);
	cout << "end of allele counts:\n";
	int parentVectorSize = (2*m)-2;                       // binary tree: m leafs and m-1 inner nodes, root has no parent

	double dropOut = 0.20;
	double seqErr = 1e-5;
	double tau = 1.0;
	vector<double> moveProbs = setMoveProbs();

	setClusterColors();
	vector<struct treeTheta> bestTrees;

	vector<string> treeFiles;
	string treeFilePrefix = getOutputFilePrefix(fileName, outFile);
	for(int rank=0; rank<treeRankSize; rank++){
		stringstream filename;
		filename << treeFilePrefix << "_rank_" << rank+1 << ".gv";
		cout << filename.str() << endl;
		treeFiles.push_back(filename.str());
	}
	//getchar();
	cout << "before MCMC\n";
	runMCMCnew(rep, loops, gamma_, moveProbs, n, m, sampleCount, alleleCount, leafClusterId, mutReadCounts, totalReadCounts, dropOut, seqErr, tau, label, treeFiles, mutLabels, wbcStatus, prefix, sampleStep, burnin);

//	for(int i = 0; i < 10000; i++){
//		int*  currTreeParentVec  = getRandomBinaryTree(m);     // sample a random binary tree
//		bool** ancMatrix = parentVector2ancMatrix(currTreeParentVec, (2*m)-2);
//		//cout << getGraphVizBinTree(currTreeParentVec, m, label);
//		writeToFile(getGraphVizBinTree(currTreeParentVec, m, label), "/Users/jahnka/Desktop/AcetoData/CTC-Exome_vcf_180105/Br23.gv");






//		for(int i=0; i<expVarReadCounts.size(); i++){
//			cout << sampleInfo.at(i).at(0) << ": ";
//			for(int j=0; j<expVarReadCounts.at(i).size(); j++){
//				cout << expVarReadCounts[i][j] << " ";
//			}
//			cout << " = " << sampleInfo.at(i).at(1);
//			cout << endl;
//		}

//		computeMaxAlleleCount(sampleInfo);
//
//		cout << "computing score\n";
//		double score = scoreTree(m, n, sampleCount, ancMatrix, alleleCount, leafClusterId, mutReadCounts, totalReadCounts, dropOut, seqErr, tau);
//		cout << score << endl;
//		cout << "...\n";
//		free_boolMatrix(ancMatrix);
//		delete [] currTreeParentVec;
//	}
	//getchar();


	//double* errorRates = getErrorRatesArray(fd, ad1, ad2, cc);



	/* initialize the random number generator, either with a user defined seed, or a random number */
	useFixedSeed? srand(fixedSeed) : initRand();


	/**  Run MCMC  **/

	//print_intMatrix(dataMatrix, m, n, ' ');
	//sampleOutput = runMCMCnew(optimalTrees, errorRates, rep, loops, gamma_, moveProbs, n, m, dataMatrix, scoreType, trueParentVec, sampleStep, sample, chi, priorSd_beta, priorSd_alpha, useTreeList, treeType, doubletModel, doubletProb, recMut, z, alphaMoveRatio, fixedDoubletProb);

//	/***  output results  ***/
//	//string prefix = getOutputFilePrefix(fileName, outFile);
//
//
//
//	/* output the optimal trees individually */
//	//double** logScores = getLogScores(fd, ad1, ad2, cc);
//	//int parentVectorSize = n;
//	if(treeType=='t'){parentVectorSize = (2*m)-2;}   // transposed case: binary tree, m leafs and m-1 inner nodes, root has no parent
//
//	int outputSize = optimalTrees.size();
//
//
//	//vector<string> geneNamesList = getGeneNamesRecMut(geneNameFile, n, recMutAllowed, recMut);
//	vector<string> sampleNamesList = getSampleNames(sampleNameFile, m);
//	//for(int i=0; i<geneNamesList.size(); i++){
//	//	cout << geneNamesList.at(i) << endl;
//	//}
//	for(int i=0; i<sampleNamesList.size(); i++){
//		cout << sampleNamesList.at(i) << endl;
//	}
//
//	for(int i=0; i<outputSize; i++){
//
//		int* parentVector = optimalTrees.at(i).tree;
//		bool** ancMatrix = parentVector2ancMatrix(parentVector, parentVectorSize);
//		vector<vector<int> > childLists = getChildListFromParentVector(parentVector, parentVectorSize);
//
//		// print newick presentation
//		stringstream newick;
//		string outputFile = getFileName2(i, prefix, ".newick", scoreType);
//		newick << getNewickCode(childLists, parentVectorSize) << "\n";
//		//newick <<  getNewickCodeMutNames(childLists, parentVectorSize, nodeLabels);
//		writeToFile(newick.str(), outputFile);
//		//cout << newick.str() << "\n";
//
//		// print GraphViz representation
//		outputFile = getFileName2(i, prefix, ".gv", scoreType);
////		if(errorRateMove != 0.0){
////			updateLogScoresAlphaBeta(logScores, optimalTrees[i].beta, optimalTrees[i].beta);
////		}
//
//		if(treeType == 'm'){
//			string output;
//
//			//print_intArray(parentVector, parentVectorSize);
//	//		output = getMutTreeGraphVizString(parentVector, parentVectorSize, geneNamesList, attachSamples, sampleNamesList, dataMatrix, logScores, recMutAllowed, recMut);
//
//			//output = getGraphVizFileContentNames(parentVector, parentVectorSize, geneNamesList, attachSamples, ancMatrix, m, logScores, dataMatrix);
//			writeToFile(output, outputFile);
//		}
//		else{
//	//		int* bestPlacement = getHighestOptPlacementVector(dataMatrix, n, m, logScores, ancMatrix);
//
//	//		vector<string> bestBinTreeLabels = getBinTreeNodeLabels((2*m)-1, bestPlacement, n, geneNamesList, sampleNamesList);
//			//getBinTreeGraphVizString
//			//cout << bestBinTreeLabels.size() << endl;
////			for(int i=0; i<bestBinTreeLabels.size(); i++){
////				cout << bestBinTreeLabels.at(i) << endl;
////			}
////			//string temp =  getGraphVizBinTree(optimalTrees.at(0).tree, (2*m)-1, m, bestBinTreeLabels);
////			string temp = getBinTreeGraphVizString(optimalTrees.at(i).tree, parentVectorSize, bestBinTreeLabels, sampleNamesList);
////			writeToFile(temp, outputFile);
////			//cout << getNewickCodeGeneNames(getChildListFromParentVector(parentVector, parentVectorSize), (2*m)-1, bestBinTreeLabels);
//		}
////		free_boolMatrix(ancMatrix);
//	}

	stringstream treeFileName;
	if(scoreType == 'm'){
		treeFileName << prefix << "_ml<INT>.gv";
	}
	else{
		treeFileName << prefix << "_map<INT>.gv";
	}
	cout << "optimal trees written to files:    " << treeFileName.str() << "\n";

	//delete [] logScores[0];
	//delete [] logScores;
//	free_intMatrix(dataMatrix);
	//cout << optimalTrees.size() << " opt trees \n";
	//emptyVectorFast(optimalTrees, n);


	/****************   end timing  *********************/
  		clock_t end=clock();
  		double diffticks=end-begin;
  		double diffms=(diffticks*1000)/CLOCKS_PER_SEC;
  		cout << "Time elapsed: " << diffms << " ms"<< endl;

  	/****************************************************/
}


vector<string> getSampleNames(string fileName, int count){

	vector<string> list;
	ifstream in(fileName.c_str());

	if (!in) {
		cout << "Can't find sample names file " << fileName << ", ";    // in case no labels are given, use sample ids instead
	    cout << "using ids instead.\n";
	    vector<string> list;
	    for(int i=0; i<count; i++){
	    	stringstream id;
	    	id << "s_";
	    	id << i+1;
	    	list.push_back(id.str());
	    }
	    in.close();
	    return list;
	}
	else{
		for (int i = 0; i < count; i++) {          // read labels from gene names file
			string temp;
			in >> temp;
			list.push_back(temp);
		}
		in.close();
		return list;
	}
}


vector<string> getGeneNamesRecMut(string fileName, int nOrig, bool recMutAllowed, int recMut){

	vector<string> v;
	ifstream in(fileName.c_str());

	if(recMutAllowed){
		n = nOrig+1;         // additional node name for recurrent mutation
	}
	else{
		n = nOrig;
	}
	if (!in) {
		cout << "Cannot open gene names file " << fileName << ", ";    // in case no labels are given, use node ids instead
	    cout << "using ids instead.\n";
	    vector<string> empty;
	    for(int i=0; i<nOrig; i++){
	    	stringstream id;
	    	id << i+1;
	    	empty.push_back(id.str());
	    }
	    if(recMutAllowed){              // add label for recurrent mutation
	    	stringstream temp;
	    	temp << recMut+1;
	    	temp << "_copy";
	    	empty.push_back(temp.str());
	    }
	    empty.push_back("Root");          // add label for root
	    return empty;
	}

	for (int i = 0; i < nOrig; i++) {          // read labels from gene names file
		string temp;
	    in >> temp;
	    v.push_back(temp);
	}

	if(recMutAllowed){            // add label for recurrent mutation
		stringstream temp;
		temp << v.at(recMut);
		temp << "_copy";
		v.push_back(temp.str());
	}
	v.push_back("Root");            // add label for the root

	in.close();
	return v;
}

void printGeneFrequencies(int** dataMatrix, int n, int m, vector<string> geneNames){
	for(int i=0; i<n; i++){
		int freq = 0;
		for(int j=0; j<m; j++){
			if(dataMatrix[j][i]==1 || dataMatrix[j][i]==2){
				freq++;
			}
		}
		cout << freq << "\t" << geneNames.at(i) << "\n";
	}
}



int* getParentVectorFromGVfile(string fileName, int n){
	int* parentVector = new int[n];
	std::vector<std::string> lines;
	std::ifstream file(fileName.c_str());
	std::string line;
	while ( std::getline(file, line) ) {
	    if ( !line.empty() )
	        lines.push_back(line);
	}
	for(int i=0; i < lines.size(); i++){

		std::size_t found = lines[i].find(" -> ");
		if (found!=std::string::npos){
			int parent = atoi(lines[i].substr(0, found).c_str());
			int child = atoi(lines[i].substr(found+3).c_str());
			parentVector[child-1] = parent-1;
	   }
	}
	return parentVector;
}



int getMinDist(int* trueVector, std::vector<bool**> optimalTrees, int n){
	int minDist = n+1;
	for(int i=0; i<optimalTrees.size(); i++){
		int dist = getSimpleDistance(trueVector, ancMatrixToParVector(optimalTrees.at(i), n), n);
		minDist = min(minDist, dist);
	}
	return minDist;
}


string getOutputFilePrefix(string fileName, string outFile){
	if(outFile.empty()){
		int lastIndex = fileName.find_last_of(".");
		return fileName.substr(0, lastIndex);
	}
	return outFile;
}


string getFileName(string prefix, string ending){
	stringstream fileName;
	fileName << prefix << ending;
	return fileName.str();
}

string getFileName2(int i, string prefix, string ending, char scoreType){
	stringstream fileName;
	if(scoreType == 'm'){
		fileName << prefix << "_ml" << i << ending;
	}
	else{
		fileName << prefix << "_map" << i << ending;
	}
	return fileName.str();
}

int readParameters(int argc, char* argv[]){
	for (int i = 1; i < argc; ++i) {

		if (strcmp(argv[i], "-i") == 0) {                      // name of data file
			if (i + 1 < argc) { fileName = argv[++i];}
		}else if (strcmp(argv[i], "-samples")==0) {           // name of file describing the samples
			if (i + 1 < argc) { sampleNameFile = argv[++i];}
		} else if(strcmp(argv[i], "-o")==0) {                  // (optional) output file if not generic name created from input file name
			if (i + 1 < argc) { outFile = argv[++i];}
		} else if(strcmp(argv[i], "-r") == 0) {           // # restarts of the MCMC
			if (i + 1 < argc) { rep = atoi(argv[++i]);}
		} else if(strcmp(argv[i], "-l")==0) {               // length of each MCMC chain
			if (i + 1 < argc) { loops = atoi(argv[++i]);}
		} else if(strcmp(argv[i], "-g")==0) {                // gamma, default 1, only adjust when NOT sampling from posterior
			if (i + 1 < argc) { gamma_ = atof(argv[++i]);}
		} else if(strcmp(argv[i],"-e")==0) {                    // probability of MCMC to pick move for proposing new error rate
					if (i + 1 < argc) { errorRateMove = atof(argv[++i]);}
		} else if(strcmp(argv[i], "-p")==0) {       // set the step size for sampling from posterior and invoke sampling
			if (i + 1 < argc) {
				sampleStep = atoi(argv[++i]);
				sample = true;
				if (i + 1 < argc){
					string next = argv[i+1];
					if(next.compare(0, 1, "-") != 0){
						burnin = atof(argv[++i]);
					}
				}
			}
		}else if (strcmp(argv[i], "-move_probs")==0) {        // sets the probabilities of different tree moves in MCMC
			vector<double> newMoveProbs;
			if (i + 1 < argc) { treeMoves.push_back(atof(argv[++i]));}
			if (i + 1 < argc) { treeMoves.push_back(atof(argv[++i]));}
			if (i + 1 < argc){
				string next = argv[i+1];
				if(next.compare(0, 1, "-") != 0){
					treeMoves.push_back(atof(argv[++i]));
				}
			}
			//cout << move1_prob << " " << move2_prob << " " << move3_prob << "\n";
//		} else if(strcmp(argv[i],"-rec")==0) {
//					//std::cerr << "recurrent mut " << std::endl;
//					if (i + 1 < argc && strncmp(argv[i+1], "-", 1)!=0) {
//						recMutAllowed = true;
//						//std::cerr << "recurrent mut " << std::endl;
//						recMut = atoi(argv[++i]);
//					}
		}else if (strcmp(argv[i], "-seed")==0) {
			useFixedSeed = true;
			if (i + 1 < argc) { fixedSeed = atoi(argv[++i]);}
		} else if (strcmp(argv[i],"-s")==0) {
			scoreType = 's';
		} else {
			std::cerr << "unknown parameter " << argv[i] << std::endl;
			getchar();
			return 1;
		}
	}
	return 0;
}


vector<double> setMoveProbs(){
	vector<double> moveProbs;

	moveProbs.push_back(errorRateMove);

	if(treeMoves.size()==0){                                       // use default probabilities
		moveProbs.push_back(defaultMoveProbsBin[0]);
		moveProbs.push_back(defaultMoveProbsBin[1]);
	}
	else{                                                                            // use probabilities from command line
		double sum = 0.0;
		for(int i=0; i< treeMoves.size(); i++){ sum += treeMoves[i]; }
		if(sum != 1.0){
			cerr << "move probabilities do not sum to 1.0, recalculating probabilities\n";     // normalize to sum to one
			for(int i=0; i< treeMoves.size(); i++){
				treeMoves[i] = treeMoves[i]/sum;
			}
			cout << "new move probabilities:";
			for(int i=0; i< treeMoves.size(); i++){ cout << " " << treeMoves[i];}
			cout << "\n";
		}
		for(int i=0; i< treeMoves.size(); i++){
			moveProbs.push_back(treeMoves[i]);
		}
	}
	treeMoves.clear();
	return moveProbs;
}


int** getDataMatrix(int n, int m, string fileName){

    int** dataMatrix = init_intMatrix(n, m, -1);

    ifstream in(fileName.c_str());

    cout << fileName << endl;
    if (!in) {
    	cout << "2 Cannot open file " << fileName << "\n";
      cout << fileName << endl;
      return NULL;
    }

    //cout << fileName << endl;
    for (int i = 0; i < n; i++) {
    	//cout << i << ": ";
        for (int j = 0; j < m; j++) {
            in >> dataMatrix[i][j];
           // cout << dataMatrix[i][j] << " ";
        }
        //cout << endl;
    }

    in.close();
    int** transposedMatrix = transposeMatrix(dataMatrix, n, m);
    free_intMatrix(dataMatrix);

    return transposedMatrix;
}


vector<string> getGeneNames(string fileName, int nOrig){

	vector<string> v;
	ifstream in(fileName.c_str());


	n = nOrig;

	if (!in) {
		//cout << "Cannot open gene names file " << fileName << ", ";
	    //cout << "using ids instead.\n";
	    vector<string> empty;
	    for(int i=0; i<=n; i++){
	    	stringstream id;
	    	id << i+1;
	    	empty.push_back(id.str());
	    }
	    return empty;
	}

	for (int i = 0; i < nOrig; i++) {
		string temp;
	    in >> temp;
	    v.push_back(temp);
	}
	v.push_back("Root"); // the root
	return v;
}



double* getErrorRatesArray(double fd, double ad1, double ad2, double cc){
	double* array = new double[4];
	array[0] = fd;
	array[1] = ad1;
	array[2] = ad2;
	array[3] = cc;
	return array;
}

vector<string> readMutIds(string fileName){
    fstream in;
    in.open(fileName);
    string line;
    vector<string> v;
    int i = 0;

    if (in.is_open()) {
    	/* ok, proceed with output */
    	cout << fileName << " open" << endl;
    }
    else{
    	cout << "not open" << endl;
    }

    if(!getline(in, line)){
    	cout << "no line read" << endl;
    }
    bool not_done = true;

    while (not_done)
    {
    	//cout << line << endl;
    	istringstream iss(line);
    	string token_1;
    	string token_2;

    	getline(iss, token_1, '\t');
    	getline(iss, token_2, '\t');

		stringstream id;
    	id << token_1 << "_" << token_2;
    	v.push_back(id.str());
    	if(!getline(in, line)){
    		not_done = false;
    	}
    }


    for(int i=0; i<v.size(); i++){
    	cout << i << "\t";
    	for(int j=0; j<v.at(i).size(); j++){
    		cout << v.at(i).at(j) << "\t";
    	}
    	cout << endl;
    }
    cout << i << " lines read" << endl;
    //getchar();
    in.close();
    return v;
}

/****   reads the read count file into a vector     ****/
vector<vector<int> > readTableFile(string fileName){
    fstream in;
    in.open(fileName);
    string line;
    vector<std::vector<int> > v;
    int i = 0;
    cout << "opening read count file ..." << endl;
    if (in.is_open()) {
    	/* ok, proceed with output */
    	//cout << fileName << " open" << endl;
    }
    else{
    	cout << "not open" << endl;
    }

    if(!getline(in, line)){
    	cout << "no line read" << endl;
    }
    bool not_done = true;

    while (not_done)
    {
    	//cout << "line " << i << endl;

    	istringstream iss(line);
    	string token;

        v.push_back(vector<int>());
        cout << line << endl;
        int count = 0;
        while(getline(iss, token, '\t'))
        {

        	if(count>3 && token!="\r"){
        		//cout << token << endl;
        		v[i].push_back(stoi(token));
        	}
        	count++;
        }
        ++i;
        if(!getline(in, line)){
        	not_done = false;
        }
    }


    for(int i=0; i<v.size(); i++){
//    	cout << i << "\t";
//    	for(int j=0; j<v.at(i).size(); j++){
//    		cout << v.at(i).at(j) << "\t";
//    	}
//    	cout << endl;
    }
    cout << i << " lines read from " << fileName << endl;
    //getchar();
    in.close();
    return v;
}

vector<vector<string> > readStringTableFile(string fileName){
    std::fstream in(fileName);
    std::string line;
    std::vector<std::vector<string> > v;
    int i = 0;

    while (std::getline(in, line))
    {
        string value;

        //cout << "line      : " << line << endl;

        std::string delimiter = "\t";

        size_t pos = 0;
        std::string token;
        v.push_back(vector<string>());
        while ((pos = line.find(delimiter)) != std::string::npos) {
            token = line.substr(0, pos);
            value = token;
            v[i].push_back(value);
            //std::cout << ":: " << token << std::endl;
            line.erase(0, pos + delimiter.length());
        }
        v[i].push_back(line);
       //std::cout << "read: " << line << std::endl;
        ++i;
       // cout << "line size: " << v[i].size() << endl;
    }

    return v;
}

//vector<vector<string> > readStringTableFile(string fileName){
//    fstream in(fileName);
//    string line;
//    vector<std::vector<string> > v;
//
//    int i = 0;
//
//    while (std::getline(in, line))
//    {
//        string entry;
//        stringstream ss(line);
//
//       // cout << "read: " << ss.str() << endl;
//        v.push_back(vector<string>());
//
//        while (ss >> entry)
//        {
//        	//cout << "entry: " << entry << endl;
//            v[i].push_back(entry);
//        }
//        ++i;
//        cout << "line size: " << v[i].size() << endl;
//    }
//    cout << endl;
//    return v;
//}

void filterLowReadCountSites(vector<vector<int> >& readCounts, vector<string>& mutIds, int minCellCount, int minCov){

	vector<int> remove;
	for(int i=0; i<readCounts.size(); i++){
		int count = 0;
		for(int j=0; j<readCounts.at(i).size(); j+=2){
			if(readCounts[i][j]>=minCov){ count++; }
		}
		if(count < minCellCount){
			remove.push_back(i);
		}
	}
	for(int i=remove.size()-1; i>=0; i--){
//		cout << "erasing: ";
//		for(int j=0; j<readCounts.at(remove.at(i)).size(); j++){
//			cout << readCounts[remove.at(i)][j] << "\t";
//		}
		//cout << endl;
		readCounts.erase(readCounts.begin() + remove.at(i));
		mutIds.erase(mutIds.begin() + remove.at(i));
	}
}

vector<vector<int> > getMutReadCount(vector<vector<int> >& readCounts){
	vector<vector<int> > mutCounts;
	for(int i=0; i<readCounts.size(); i++){             // get for each mutation
		mutCounts.push_back(vector<int>());
		for(int j=1; j<readCounts.at(i).size(); j+=2){      // the count of mutated reads per sample
			mutCounts.at(i).push_back(readCounts[i][j]);
		}
	}
	return mutCounts;
}

vector<vector<int> > getTotalReadCount(vector<vector<int> >& readCounts){
	vector<vector<int> > mutCounts;
	for(int i=0; i<readCounts.size(); i++){
		mutCounts.push_back(vector<int>());
		for(int j=0; j<readCounts.at(i).size(); j+=2){
			mutCounts.at(i).push_back(readCounts[i][j]);

//			if(readCounts[i][j] > readCounts[i][j+1]){
//				cout << "tot = " << readCounts[i][j] << "    mut = " << readCounts[i][j+1] <<  "   line = " << i << endl;
//			}
//			if(readCounts[i][j] < readCounts[i][j+1]){
//				cout << "tot = " << readCounts[i][j] << "    mut = " << readCounts[i][j+1] <<  "   line = " << i << endl;
//				getchar();
//			}
		}
	}
	return mutCounts;
}

vector<string> getLeafLabels(vector<vector<string> >& sampleInfo){
	vector<string> leafLabels;
	cout << "getting leaf labels\n";
	for(int i=0; i<sampleInfo.size(); i++){

		cout << i <<  ": element count: " << sampleInfo.at(i).size() << endl;
//		for(int p=0; p<sampleInfo.at(i).size(); p++){
//			cout << "sample info: " << sampleInfo.at(i).at(p) << "\t";
//		}
		cout << endl;
		int count = stoi(sampleInfo.at(i).at(1));
		cout << "count is " << count << endl;
		for(int j=1; j<=count; j++){
			stringstream temp;
			//label << "\"";
			//label << sampleInfo.at(i).at(0);
			//if(count>1){
			//	label << "_" << char (j+64);
			//}
			//label << "\n" << sampleInfo.at(i).at(4);
			//label << "\"";
			temp << sampleInfo.at(i).at(4);
			string label = temp.str();
			if(count>1){
				std::string delimiter = "\"";
				std::string name;
				size_t left  = label.find(delimiter);
				size_t right = label.find(delimiter,left+1);
				name = label.substr(left, right-left);  // sample name from sample name file
				stringstream extended;
				extended << name <<  "_" << char (j+64);
				//cout << "token: " <<  extended.str() << endl;
				size_t start_pos = label.find(name);
				//cout << start_pos << endl;
				if(start_pos != std::string::npos){
					label.replace(start_pos, name.length(), extended.str());
				}
			}
			cout<< "label: " << label << endl;
			leafLabels.push_back(label);

		}
	}
	return leafLabels;
}

vector<int> getClusterIdOfLeaf(vector<vector<string> >& sampleInfo){
	vector<int> clusterId;

	for(int i=0; i<sampleInfo.size(); i++){

		int count = stoi(sampleInfo.at(i).at(1));  // position 1 has total cell count
		for(int j=0; j<count; j++){

			clusterId.push_back(i);
		}
	}
	return clusterId;
}

/** Goes through all cells in all clusters and sets the cell state to be wbc or tumour. **/
/** This is stored in a boolean vector: 0 for tumour and 1 for wbc  **/
vector<bool> getWbcStatus(vector<vector<string> >& sampleInfo){
	vector<bool> wbcStatus;

	for(int i=0; i<sampleInfo.size(); i++){
		int cellCount = stoi(sampleInfo.at(i).at(1));  // position 1 has total cell count in cluster
		int wbcCount  = stoi(sampleInfo.at(i).at(3));  // position 3 has total wbc count in cluster
		for(int j=0; j<cellCount; j++){
			if(j < wbcCount){
				wbcStatus.push_back(1);
			}
			else{
				wbcStatus.push_back(0);
			}
		}
	}
	return wbcStatus;
}


/***     computes the  number of alleles (normal and mutated) of a cell (cluster) ***/
vector<int> getAlleleCount(vector<vector<string> >& sampleInfo){
	vector<int> alleleCount;
	for(int i=0; i<sampleInfo.size(); i++){
		alleleCount.push_back(stoi(sampleInfo.at(i).at(1))*2);
		cout << "sample " << i << ": " << alleleCount.at(i) << " alleles" << endl;
	}
	return alleleCount;
}

int** getAlpha(int clusterCount, int* parentVec, int parentVectorSize, vector<int>& leafClusterId){

	int** alpha = init_intMatrix(clusterCount,parentVectorSize+1, 0);  // init counts with zero
	int leafCount = leafClusterId.size();

	for(int id=0; id<clusterCount; id++){           // compute counts for each cluster

		for(int node=0; node<leafCount; node++){

			if(leafClusterId.at(node)==id){       // init counts for leafs
				alpha[id][node]++;
			}
		}

		for(int node=0; node<parentVectorSize; node++){

			alpha[id][parentVec[node]] += alpha[id][node];
		}
	}
	return alpha;
}
