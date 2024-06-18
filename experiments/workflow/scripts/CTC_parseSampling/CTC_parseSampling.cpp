//============================================================================
// Name        : CTC_parseSampling.cpp
// Author      : Katharina Jahn
// Version     :
// Copyright   : Your copyright notice
// Description : Parse the sampling output of CTC_SCITE. Get information such as
//               separating mutations for leaf pairs, frequency of tree splits,...
//============================================================================

#include <stdbool.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <time.h>
#include <string>
#include <map>
#include <set>
#include <utility>

using namespace std;

int readParameters(int argc, char* argv[]);
vector<vector<int> > readTableFile(string fileName);
vector<string> readMutIds(string fileName);
vector<vector<int> > getMutReadCount(vector<vector<int> >& readCounts);
vector<vector<int> > getTotalReadCount(vector<vector<int> >& readCounts);
vector<vector<string> > readStringTableFile(string fileName);
vector<string> getLeafLabels(vector<vector<string> >& sampleInfo);
vector<int> getClusterIdOfLeaf(vector<vector<string> >& sampleInfo);
vector<bool> getWbcStatus(vector<vector<string> >& sampleInfo);
vector<int> getAlleleCount(vector<vector<string> >& sampleInfo);
void computeMaxAlleleCount(vector<vector<string> >& sampleInfo);
vector<vector<int> > getExpVarAlleleCount(int m, int sampleCount, bool** ancMatrix, vector<int>& leafClusterId);
double logBinomCoeff(int n, int k);
double logBetaBinom(int k, int r, int alpha, int beta, int p, int q, double tau);
double computeLogSampleScore(int expCount, int alleleCount, int obsMutCount, int obsTotalCount, double dropoutRate, double seqErr, double tau);
int sampleEvent(std::vector<double> prob);
vector<int> getMutPlacement(int m, int n, int sampleCount, bool** ancMatrix, vector<int>& alleleCount, vector<int>& leafClusterId, vector<vector<int> >& mutReadCounts, vector<vector<int> >& totalReadCounts, double dropoutRate, double seqErr, double tau, vector<bool>& wbcStatus);
bool** parentVector2ancMatrix(int* parent, int n);
bool** init_boolMatrix(int n, int m, bool value);
bool** allocate_boolMatrix(int n, int m);
string getActualLabel(string s);
int getLCA(int v, int w, int* parentVector, vector<vector<bool > >& ancMatrix, int leafCount);
vector<vector<bool > > parentVector2ancMatrixWithRoot(int* parent, int n);
bool isThreeCellSample(string str);
string getNodeColor(string s);
int getHPA(int v, vector<vector<bool > >& sameOrigCluster, int* parentVector, vector<vector<bool > >& ancMatrix, int leafCount);
void updateSplitCounts(std::map<string, int> &splitCounts, vector<vector<bool> > &ancMatrix, vector<string> &leafLabels, int nodeCount, int* parent);
string dataFileName;            // file with mutations/read counts
string descriptionFileName;     // file describing the CTC/small cell clusters
string sampleFileName;          // file containing the MCMC samples

int m;  // number of sequenced sampled
int n;  // number of mutations

double chi_wbc = 1.0; // penalty for placing a mutation above a wbc

int	maxAlleleCount = 0;

int swapCount = 0;
int noSwapCount = 0;

int main(int argc, char* argv[]) {

	vector<string> mutLabels;
	vector<vector<int> > readCounts;            // the original read count table
	vector<vector<string> > sampleInfo;         // the original table of the samples (ids, cell count, type,...)
	int sampleCount;

	vector<vector<int> > mutReadCounts;          // entry[i][j] has mutated read count for position i in sample j
	vector<vector<int> > totalReadCounts;        // entry[i][j] has total read count for position i in sample j

	vector<string> leafLabels;
	vector<int> leafClusterId;
	vector<int> alleleCount;
	vector<bool> wbcStatus;               // wbcStatus[i] is 1 if leaf i represents a white blood cell, 0 otherwise

	// get variant and total read count
	readParameters(argc, argv);
	readCounts = readTableFile(dataFileName);
	mutLabels  = readMutIds(dataFileName);

	int lastIndex = sampleFileName.find_last_of(".");
	string prefix =  sampleFileName.substr(0, lastIndex);

	mutReadCounts = getMutReadCount(readCounts);            // count of mutated reads after filtering per position per sample
	totalReadCounts = getTotalReadCount(readCounts);        // total read count after filtering per position per sample
	cout << "sample count = " << mutReadCounts.at(0).size()<< endl;

	// get leaf labels and sample sizes
	cout << descriptionFileName << endl;
	sampleInfo = readStringTableFile(descriptionFileName);   // the original table of the samples (ids and cell count)
	sampleCount = sampleInfo.size();
	cout << sampleInfo.size() << endl;
	leafLabels = getLeafLabels(sampleInfo);        // list of leaf label names


 	leafClusterId = getClusterIdOfLeaf(sampleInfo);   // maps leaf id to cluster id
	for(int i=0; i<leafClusterId.size(); i++){
		cout << leafClusterId.at(i) << endl;
	}
	m = leafClusterId.size();                 // number of leafs in tree
	n = mutReadCounts.size();                 // number of mutations in tree
	cout << m << " samples" << endl;
	cout << n << " muts" << endl;

	wbcStatus = getWbcStatus(sampleInfo);

	vector<string> label;                           // the labels of the tree nodes
	for(int i=0; i<leafLabels.size(); i++){
		//cout << leafLabels.at(i) << endl;
		label.push_back(leafLabels.at(i));
	}
	for(int i=leafLabels.size(); i<(2*m)-1; i++){
		label.push_back("x"+to_string(i));
	}

//	for(int i=0; i<leafLabels.size(); i++){
//		cout << i<< "\t" << leafLabels.at(i) << endl;
//	}
//	cout << endl;
//	getchar();
	computeMaxAlleleCount(sampleInfo);
	alleleCount = getAlleleCount(sampleInfo);
	cout << "end of allele counts:\n";
	int parentVectorSize = (2*m)-2;                       // binary tree: m leafs and m-1 inner nodes, root has no parent

	vector<vector<int> > mutPlacementCounts;     // counts for each mutation and each tree node how often the mutation
	for(int i=0; i<n; i++){                      // is placed at this node in the posterior sampling
		vector<int> vect((2*m)-1, 0);
		mutPlacementCounts.push_back(vect);
	}


	/******************************************************************************************************/
	/*  This is the new initialisation where we count for each mutation pair and each original cluster    */
	/*  how often the pair co-occurs in one of the clusters's private branches.                           */
	/*  Here a private branch is defined as the path up from a leaf just before the lca of the leaf       */
	/*	and another leaf of the same cluster. If the original cluster was broken into multiple samples,   */
	/*	we can also define private branches for each unbroken (multi-cell) sample. That means a branch is */
	/*	private until an lca with a cell from one of the other (multi-cell) samples is reached.           */
	/******************************************************************************************************/

	vector<vector<vector<int > > > coOccCount;            // for each leaf and each mutation pair count how often
	for(int leaf=0; leaf<leafClusterId.size(); leaf++){      //  the pair co-occurs in the private lineage of the leaf
		vector<vector<int> > temp1;
		for(int mut=0; mut<n; mut++){
			vector<int> temp2(n,0);
			temp1.push_back(temp2);
		}
		coOccCount.push_back(temp1);
	}

	vector<vector<bool > > sameOrigCluster;    // Is '1' if cell i and cell j are from same original cluster (same color)
	for(int leaf=0; leaf<leafClusterId.size(); leaf++){             // Initialize default to '0'
		vector<bool> temp(leafClusterId.size(),false);
		sameOrigCluster.push_back(temp);
	}
	for(int cell_1=0; cell_1<leafClusterId.size(); cell_1++){
		string color1 = getNodeColor(leafLabels.at(cell_1));
		if(color1 == "gray93" || color1 =="ghostwhite"){   // Skip sample if it has no proper color,
			continue;                                      // sample is a separately sequenced single cell (tumor or wbc),
		}
		for(int cell_2=0; cell_2<leafClusterId.size(); cell_2++){
			string color2 = getNodeColor(leafLabels.at(cell_2));
			if(color2 == "gray93" || color1 == "ghostwhite"){      // Skip if second sample is separately sequenced single cell
				continue;
			}
			if(getNodeColor(leafLabels.at(cell_1)) == getNodeColor(leafLabels.at(cell_2))){    // both samples have same color
				sameOrigCluster.at(cell_1).at(cell_2) = true;                               // so they're from same orig cluster
			}
		}
	}

	/* counts for each leaf pair the average number of separating mutations (minimum of both branches)*/
	vector<vector<int > > numberOfMinSeparatingMuts;
	for(int i=0; i<leafLabels.size(); i++){
		vector<int> temp(leafLabels.size(),0);
		for(int j=0; j<leafLabels.size(); j++){
			numberOfMinSeparatingMuts.push_back(temp);
		}
	}

	vector<int> sepMutCountHist;

	/* gets the split statistics that Francesc uses */
	/* Lu7_6, Lu7_8, Lu7_10, Lu7_12 | Lu7_2, Lu7_3, Lu7_5 | all other samples */
	map<string,int> splitCounts;

	/******************************************************************************************/
	/*    Now we read the file containing the posterior samples, each line is one sample,     */
	/*    and extract all the needed information                                              */
	/******************************************************************************************/
	string line;
	ifstream f (sampleFileName);
	int counter = 0;
	while (getline(f, line))
	{
	    stringstream instream(line);
	    vector<string> temp;
	    string token;
	    cout << "line " << counter  << "    : " << line << endl;
	    counter++;
	    if(counter % 100 == 0){
	    	cout << counter << endl;
	    }
	    while (std::getline(instream, token, '\t'))
	    {
	        temp.push_back(token);
	    }
	    if(temp.size()-1 !=4){
	    	cout << "wrong number of tabs in line: " << temp.size() << endl;
	    	getchar();
	   	}

	    double currSeqErrRate = stod(temp.at(1));
	    double currDropOutRate = stod(temp.at(2));
	    double currLogTau = stod(temp.at(3));

	    std::string buf;                        // Have a buffer string
	    std::stringstream ss(temp.at(4));       // Insert the string into a stream
	    std::vector<std::string> elems;        // Create vector to hold node elems


	    while (ss >> buf)
	    	elems.push_back(buf);

	    int* currParentVector = new int[elems.size()];
	    if(parentVectorSize != elems.size()){
	    	cout << "size error" << endl;
	    	getchar();
	    }

	    for(int i=0; i<elems.size(); i++){
	    	currParentVector[i] = stoi(elems.at(i));
	    }


	    bool** currAncMatrix = parentVector2ancMatrix(currParentVector, parentVectorSize);
	    vector<vector<bool > > currAncMatrixWithRoot = parentVector2ancMatrixWithRoot(currParentVector, parentVectorSize);
	    vector<int> mutPlacement = getMutPlacement(m, n, sampleCount, currAncMatrix, alleleCount, leafClusterId, mutReadCounts, totalReadCounts, currDropOutRate, currSeqErrRate, currLogTau, wbcStatus);

	    /* Update the Split Counts for Francesco */
	    updateSplitCounts(splitCounts, currAncMatrixWithRoot, leafLabels, (2*m)-2, currParentVector);

	    for(int mut=0; mut<n; mut++){
	    	mutPlacementCounts.at(mut).at(mutPlacement.at(mut)) += 1;  // update counts how often a mutation is placed at
	    }                                                              // at a tree node in the posterior sampling

	    vector<int> leafHPA;                        // get for each leaf the highest ancestor in private branch
	    for(int v=0; v<leafLabels.size(); v++){     //  (right below the first lca with other leaf from same original cluster)
	    	leafHPA.push_back(getHPA(v, sameOrigCluster, currParentVector, currAncMatrixWithRoot, leafLabels.size()));
	    }
	    //cout << "HPA done" << endl;
		for(int i=0; i<leafLabels.size(); i++){                   // go through all leafs and mutation pairs
			//cout << getActualLabel(leafLabels.at(i)) << endl;	  // update co-occurrence counter if the pair
			for(int mut_1=0; mut_1<n; mut_1++){			          // co-occurs in the leafs private lineage
				int mutNode_1 = mutPlacement[mut_1];

				// skip if first mutation is not in private lineage of leaf
				if(!(currAncMatrixWithRoot[mutNode_1][i] && currAncMatrixWithRoot[leafHPA.at(i)][mutNode_1])){
					continue;
				}

				for(int mut_2=0; mut_2<n; mut_2++){
					int mutNode_2 = mutPlacement[mut_2];

					// skip if second mutation is not in private lineage of leaf
					if(!(currAncMatrixWithRoot[mutNode_2][i] && currAncMatrixWithRoot[leafHPA.at(i)][mutNode_2])){
						continue;
					}
					coOccCount.at(i).at(mut_1).at(mut_2)++;
				}
			}
		}
		//cout << "co-occ init done" << endl;

		// count the number of separating mutations for all leaf pairs, i.e. the number of mutation s
		for(int i=0; i<leafLabels.size(); i++){
			for(int j=i+1; j<leafLabels.size(); j++){
				int counter1 = 0;
				int counter2 = 0;
				int lca = getLCA(i, j, currParentVector, currAncMatrixWithRoot, leafLabels.size());
				for(int mut=0; mut<n; mut++){
					int mutNode = mutPlacement[mut];
					if(currAncMatrixWithRoot[lca][mutNode] && currAncMatrixWithRoot[mutNode][i] && mutNode!=lca){
						counter1++;
					}
					if(currAncMatrixWithRoot[lca][mutNode] && currAncMatrixWithRoot[mutNode][j] && mutNode!=lca){
						counter2++;
					}
				}
				numberOfMinSeparatingMuts[i][j] += min(counter1, counter2);

				// For getting the histogram of min separation counts for ambiguous cluster in Pr9
//				if(getActualLabel(leafLabels.at(i))=="\"Pr9_CTC_26_A\"" && getActualLabel(leafLabels.at(j))=="\"Pr9_CTC_26_B\""){
//					sepMutCountHist.push_back(min(counter1, counter2));
//				}
			}
		}
	    /************************************************************************************************************/
//	    int root = (2*m)-2;
//	    for(int mut=0; mut<n; mut++){
//	    	for(int i=0; i<leafLabels.size(); i++){
//	    		for(int j=i+1; j<leafLabels.size(); j++){
//	    			//cout << mut << "\t" << i << "\t" << j << endl;
//	    			int cand = mutPlacement.at(mut);
//	    			if(cand==root){
//	    				otherMuts.at(i).at(j) += 1;
//	    				continue;
//	    			}
//	    			if(currAncMatrix[cand][i]){
//	    				if(currAncMatrix[cand][j]){
//	    					otherMuts.at(i).at(j) += 1;
//	    				}
//	    				else{
//	    					privMuts_1.at(i).at(j) += 1;
//	    					if(swap){
//	    						separatingMuts_2.at(i).at(j).at(mut) += 1;
//	    					}
//	    					else{
//	    						separatingMuts_1.at(i).at(j).at(mut) += 1;
//	    					}
//	    				}
//	    			}
//	    			else if(currAncMatrix[cand][j]){
//	    				privMuts_2.at(i).at(j) += 1;
//	    				if(swap){
//	    					separatingMuts_1.at(i).at(j).at(mut) += 1;
//	    				}
//	    				else{
//	    					separatingMuts_2.at(i).at(j).at(mut) += 1;
//	    				}
//	    			}
//	    			else{
//	    				otherMuts.at(i).at(j) += 1;
//	    			}
//	    		}
//	    	}
//	    }

	}        // file with posterior samples read
	cout << "file read" << endl;

	/* output minimum number of separating mutations per cluster (average over all pairs) */
	stringstream dummyOutFileName;
	dummyOutFileName << prefix << "_sepMinStats.txt";
	ofstream sepMinStatsOutput (dummyOutFileName.str());     // dummy output for average number of separating mutations (sum branches)


	for(int leaf=0; leaf<leafClusterId.size(); leaf++){
		vector<int> clusterList;
		bool isFirst = true;
		for(int leaf2=0; leaf2<leafClusterId.size(); leaf2++){
			if(sameOrigCluster[leaf][leaf2]){
				if(leaf2 < leaf){
					isFirst=false;
					break;
				}
				else{
					clusterList.push_back(leaf2);
				}
			}
		}
		if(clusterList.size()>0){
			sepMinStatsOutput << getNodeColor(leafLabels.at(leaf)) << endl;
		}
		for(int i=0; i<clusterList.size(); i++){
			for(int j=i+1; j<clusterList.size(); j++){
				sepMinStatsOutput << getActualLabel(leafLabels.at(clusterList.at(i))) << "   ";
				sepMinStatsOutput << getActualLabel(leafLabels.at(clusterList.at(j))) << "   ";
				sepMinStatsOutput << (1.0*numberOfMinSeparatingMuts[clusterList.at(i)][clusterList.at(j)])/counter << endl;
			}
		}
		sepMinStatsOutput << endl;
	}

	stringstream scatterOutFileName;
	scatterOutFileName << prefix << "_sepMinStatsScatter.txt";
	ofstream scatterSepMinStatsOutput (scatterOutFileName.str());

	for(int leaf=0; leaf<leafClusterId.size(); leaf++){
		for(int leaf2=leaf+1; leaf2<leafClusterId.size(); leaf2++){
			scatterSepMinStatsOutput << getActualLabel(leafLabels.at(leaf)) << "_";
			scatterSepMinStatsOutput << getActualLabel(leafLabels.at(leaf2)) << "\t";
			scatterSepMinStatsOutput << (1.0*numberOfMinSeparatingMuts[leaf][leaf2])/counter << endl;
		}
	}


	// For getting the histogram of min separation counts for ambiguous cluster in Pr9
//	for (int i=0; i<sepMutCountHist.size(); i++){
//		sepMinStatsOutput << sepMutCountHist.at(i) << endl;
//	}

	/* output top separating mutations */
	for(int leaf=0; leaf<leafClusterId.size(); leaf++){
		vector<int> relatedList;
		bool isFirst = true;
		for(int leaf2=0; leaf2<leafClusterId.size(); leaf2++){
			if(sameOrigCluster[leaf][leaf2]){
				if(leaf2 < leaf){
					isFirst=false;
					break;
				}
				else{
					relatedList.push_back(leaf2);
				}
			}
		}
		if(!isFirst || relatedList.size()==0){
			continue;
		}
		cout << getNodeColor(leafLabels.at(leaf));
		for(int i=0; i<relatedList.size(); i++){
			cout << "\t" << getActualLabel(leafLabels.at(relatedList.at(i)));
		}
		cout << endl;

		stringstream dummyOutFileName2;
		dummyOutFileName2 << prefix << "_" << getNodeColor(leafLabels.at(leaf)) << ".txt";
		ofstream distOutput (dummyOutFileName2.str());     // outfile for table with pairwise mutation distances

		for(int mut_2=0; mut_2<n; mut_2++){                   // header
			distOutput << "\t" << mutLabels.at(mut_2);
		}
		distOutput << endl;

		for(int mut_1=0; mut_1<n; mut_1++){
			distOutput << mutLabels.at(mut_1);
			for(int mut_2=0; mut_2<n; mut_2++){
				int coOccSum = 0;
				for(int i=0; i<relatedList.size(); i++){
					coOccSum += coOccCount.at(relatedList.at(i)).at(mut_1).at(mut_2);
				}
				if(mut_1==mut_2){
					distOutput << "\t" << 0.0;
				}
				else{
					distOutput << "\t" << 1.0 - (1.0*coOccSum/counter);
				}
		//					//tableCoOccOutput << "\t" << coOccCount.at(i).at(j).at(mut_1).at(mut_2);
			}
			distOutput << endl;
		}
		distOutput << endl;

	}


//	stringstream dummyOutFileName;
//	dummyOutFileName << prefix << "_coOccProbs.txt";
//	ofstream tableCoOccOutput (dummyOutFileName.str());             // dummy output for posterior co-occurrence probability
//	if (!tableCoOccOutput.is_open()){
//		cout << "Unable to open file ";
//	}
//	cout << dummyOutFileName.str() << endl;
//
//
//
//	for(int i=0; i<leafLabels.size(); i++){
//		for(int j=i+1; j<leafLabels.size(); j++){
//			stringstream temp;
//			temp << getActualLabel(leafLabels.at(i)) << "-" << getActualLabel(leafLabels.at(j));
//			string leafPair = temp.str();
//			leafPair.erase(std::remove(leafPair.begin(), leafPair.end(), '\"'), leafPair.end());
//			tableCoOccOutput << leafPair << endl;
//			for(int mut_2=0; mut_2<n; mut_2++){
//				tableCoOccOutput << "\t" << mutLabels.at(mut_2);
//			}
//			tableCoOccOutput << endl;
//			//getchar();
//			for(int mut_1=0; mut_1<n; mut_1++){
//				tableCoOccOutput << mutLabels.at(mut_1);
//				for(int mut_2=0; mut_2<n; mut_2++){
//					if(mut_1==mut_2){
//						tableCoOccOutput << "\t" << 0.0;
//					}
//					else{
//						tableCoOccOutput << "\t" << 1.0 - (1.0*coOccCount.at(i).at(j).at(mut_1).at(mut_2)/counter);
//					}
//					//tableCoOccOutput << "\t" << coOccCount.at(i).at(j).at(mut_1).at(mut_2);
//				}
//				tableCoOccOutput << endl;
//			}
//			tableCoOccOutput << endl;
//		}
//	}
//	cout << "max count = " << counter << endl;
//
//	string tableFile = prefix + "_table.tsv";     // output file containing for each mutation and each tree node how
//	ofstream tableOutput (tableFile);             // often the mutation was placed at that node in the posterior sample
//	if (!tableOutput.is_open()){
//		cout << "Unable to open file";
//	}
//	cout << tableFile << endl;
//	//getchar();
//	for(int mut=0; mut<n; mut++){
//		for(int node=0; node<(2*m)-1;node++){
//			if(node!=0){
//				tableOutput << "\t";
//			}
//			tableOutput << mutPlacementCounts.at(mut).at(node);
//		}
//		tableOutput << endl;
//	}



//	string privFile = prefix + "_privMuts.tsv";  // output file containing for each leaf pair the average number of
//	ofstream privOutput (privFile);              // private mutation in each of the two leaf lineages and average
//	if (!privOutput.is_open()){                  // number of other mutations
//		cout << "Unable to open file";
//	}
//	cout << privFile << endl;
//	privOutput << "leaf_1" << "\t" << "leaf_2" << "\t" << "privateMuts_1" << "\t" <<  "privateMuts_1" << "\t" << "otherMuts" << endl;
//	for(int i=0; i<leafLabels.size(); i++){
//		for(int j=i+1; j<leafLabels.size(); j++){
//			privOutput << getActualLabel(leafLabels.at(i)) << "\t" << getActualLabel(leafLabels.at(j)) << "\t" << (1.0*privMuts_1.at(i).at(j))/counter << "\t";
//			privOutput << (1.0*privMuts_2.at(i).at(j))/counter << "\t" << 1.0*(otherMuts.at(i).at(j))/counter << "\t";
//			privOutput << (1.0*(privMuts_1.at(i).at(j)+privMuts_2.at(i).at(j)+otherMuts.at(i).at(j)))/counter << "\t" << counter << endl;
//		}
//	}

//	string sepFile = prefix + "_separationStat.tsv";       // output file listing for each leaf pair and each mutation the
//	ofstream sepOutput (sepFile);                          // probability that the mutation is in the left or right
//	if (!sepOutput.is_open()){                             // private lineage and the max and min of these values
//		cout << "Unable to open file";
//	}
//	cout << sepFile << endl;
//	for(int i=0; i<leafLabels.size(); i++){
//		for(int j=i+1; j<leafLabels.size(); j++){
//			sepOutput << getActualLabel(leafLabels.at(i)) << "\t" << getActualLabel(leafLabels.at(j)) << endl;
//			for(int mut=0; mut<n; mut++){
//				sepOutput << mutLabels.at(mut) << ":\t";
//				sepOutput << 1.0*separatingMuts_1.at(i).at(j).at(mut)/counter << "\t";
//				sepOutput << 1.0*separatingMuts_2.at(i).at(j).at(mut)/counter << "\t";
//				sepOutput << max(1.0*separatingMuts_1.at(i).at(j).at(mut)/counter, 1.0*separatingMuts_2.at(i).at(j).at(mut)/counter) << "\t";
//				sepOutput << min(1.0*separatingMuts_1.at(i).at(j).at(mut)/counter, 1.0*separatingMuts_2.at(i).at(j).at(mut)/counter) << endl;
//			}
//		}
//	}

	//getchar();
//	string topSepFile = prefix + "_topSeparators.tsv";
//	ofstream topSepOutput (topSepFile);
//	if (!topSepOutput.is_open()){
//		cout << "Unable to open file";
//	}
//	for(int i=0; i<leafLabels.size(); i++){
//		for(int j=i+1; j<leafLabels.size(); j++){
//			//cout << "\nTop separators for " << getActualLabel(leafLabels.at(i)) << " and " << getActualLabel(leafLabels.at(j)) << ":" << endl;
//			topSepOutput << "\nTop separators for " << getActualLabel(leafLabels.at(i)) << " and " << getActualLabel(leafLabels.at(j)) << ":" << endl;
//			vector<int> copy_1 = separatingMuts_1.at(i).at(j);
//			vector<int> copy_2 = separatingMuts_2.at(i).at(j);
//			vector<string> topSepMuts_1;
//			vector<string> topSepMuts_2;
//			for(int top=0; top<10; top++){
//				int maxSepIndex=0;
//				int maxSep=0;
//				for(int mut=0; mut<n; mut++){
//					//cout << abs(copy_1.at(mut)-copy_2.at(mut)) << endl;
//					int sep = copy_1.at(mut);
//					if(sep>maxSep){
//						maxSepIndex = mut;
//						maxSep = sep;
//						//cout << "top: " << mut << endl;
//					}
//				}
//				copy_1.at(maxSepIndex) = 0;
//				stringstream s;
//				s << mutLabels.at(maxSepIndex) << "\t" << 1.0*maxSep/counter;
//				topSepMuts_1.push_back(s.str());
//			}
//			for(int top=0; top<10; top++){
//				int maxSepIndex=0;
//				int maxSep=0;
//				for(int mut=0; mut<n; mut++){
//					//cout << abs(copy_1.at(mut)-copy_2.at(mut)) << endl;
//					int sep = copy_2.at(mut);
//					if(sep>maxSep){
//						maxSepIndex = mut;
//						maxSep = sep;
//						//cout << "top: " << mut << endl;
//					}
//				}
//				copy_2.at(maxSepIndex) = 0;
//				stringstream s;
//				s << mutLabels.at(maxSepIndex) << "\t" << 1.0*maxSep/counter;
//				topSepMuts_2.push_back(s.str());
//			}
//			topSepOutput << getActualLabel(leafLabels.at(i)) << ":" << endl;
//			for(int t=0; t<10; t++){
//				topSepOutput << topSepMuts_1.at(t) << endl;
//			}
//			topSepOutput << getActualLabel(leafLabels.at(j)) << ":" << endl;
//			for(int t=0; t<10; t++){
//				topSepOutput << topSepMuts_2.at(t) << endl;
//			}
//			//cout << endl;
//		}
//	}

	stringstream splitsOutFileName;
	splitsOutFileName << prefix << "_splits.txt";
	ofstream splitsOutput (splitsOutFileName.str());

	for(auto elem : splitCounts)
	{
		splitsOutput << 1.0*elem.second/counter << "\t" << elem.first << "\n";
	}

//	auto cmp = [](const auto &p1, const auto &p2)
//	{
//		return p2.second < p1.second || !(p1.second < p2.second) && p1.first < p2.first;
//	};
//
//	std::set < std::pair<string, int>, decltype( cmp )> s(splitCounts.begin(), splitCounts.end(), cmp);
//
//	for (const auto &p : s)
//	{
//		std::cout << "{ " << p.first << ", " << p.second << " }\n";
//	}
//	std::cout << std::endl;

}

/*  Compute the highest private ancestor of a leaf, i.e. the earliest ancestor of a sample that is not an ancestor  */
/*  of another sample from the same original cluster (same color)  */
int getHPA(int v, vector<vector<bool > >& sameOrigCluster, int* parentVector, vector<vector<bool > >& ancMatrix, int leafCount){
	int root = 2*leafCount-2;
	int curr = v;
	int anc = parentVector[curr];

	bool done = false;
	while(anc!=root && !done){
		for(int w=0; w<leafCount; w++){
			if(v==w){
				continue;
			}
			if(sameOrigCluster.at(v).at(w)){
				if(ancMatrix.at(anc).at(w)){
					done = true;
					break;
				}
			}
		}
		if(!done){
			curr = anc;
			anc = parentVector[anc];
		}
	}
	return curr;
}

int getLCA(int v, int w, int* parentVector, vector<vector<bool > >& ancMatrix, int leafCount){
	int root = 2*leafCount-2;
	int anc = v;

	while(anc!=root){
		if(ancMatrix.at(anc).at(w)){
			return anc;
		}
		anc = parentVector[anc];
	}
	return root;
}

string getActualLabel(string s){
	string temp = s;
	string left = ",label=";
	string right = ",fillcolor=";
	int pos = temp.find(left);
	temp.erase(0, pos + left.length());
	pos = temp.find(right);
	string token = temp.substr(0, pos);
	return token;
}

string getNodeColor(string s){
	string temp = s;
	string left = ",fillcolor=";
	string right = ",image=";
	int pos = temp.find(left);
	temp.erase(0, pos + left.length());
	pos = temp.find(right);
	string token = temp.substr(0, pos);
	return token;
}

	//print_intArray(currTreeParentVec, 2*m-2);
	    			//getchar();

	    			/* Sample from the posterior if required and past the burn-in phase */
	    			//        		if(sample && it>=burnIn && it % step == 0){
	    			//        			// get mutation matrix and update counts
	    			//        			for(int mut=0; mut<n; mut++){
	    			//        				for(int leaf=0; leaf<m; leaf++){
	    			//        					if(currTreeAncMatrix[bestPlacementPoint[mut]][leaf]){
	    			//        						postSamplingMatrix[mut][sample]++;
	    			//        					}
	    			//        				}
	    			//        			}
	    			//        		}


//	    			/* Get for each sample/attachment point pair the probability of placement */
//	    			vector<vector<double> > currTable = getPlacementScoreTable(m, n, sampleCount, currTreeAncMatrix, alleleCount, leafClusterId, mutReadCounts, totalReadCounts, dropoutRate, seqErr, tau, bestPlacementPoint, wbcStatus);
//	    			//cout << "table size: " << currTable.size() << endl;
//
//	    			/* compute for all leaf pairs, the separation score*/
//	    			string line = getSeparationScores(currTable, currTreeAncMatrix, n, m, currTreeParentVec, leafPairs);
//	    			sampleOutput << line;
//	    			//cout << line << endl;
//
//
//	    			string partition = getPartition(currTreeAncMatrix, m, label);
//	    			treeSampleOutput << partition;
//}





//string getSeparationScores(vector<vector<double> > table, bool** ancMatrix, int n, int m, int* parentVector, vector<string> leafPairs){
//	vector<vector<double> > normScoreTable;
//	for(int mut=0; mut<table.size(); mut++){
//		vector<double> normScores;
//		double bestLogScore = table.at(mut).at(0);           // best placement score of mut in the precomputed table
//
//		for(int node=0; node<table.at(mut).size(); node++){             // find the best placement score for mut
//			bestLogScore = max(bestLogScore, table.at(mut).at(node));  // among all possible placement points
//		}
//
//		double tempScore = 0.0;
//		for(int node=0; node<table.at(mut).size(); node++){              // sum of placement scores for mutation
//			tempScore += exp(table.at(mut).at(node)-bestLogScore);   // subtraction of best score to calculate with score differences (smaller values)
//		}
//		double logSumScore = log(tempScore)+bestLogScore;              // transform back to log scores and change from score differences to actual scores
//
//		//cout << exp(logSumScore) << endl;
//
//
//		for(int node=0; node<table.at(mut).size(); node++){
//			double normScore = table.at(mut).at(node)-logSumScore;    // set table with normalized scores
//			normScores.push_back(normScore);
//			//cout << "norm score mut " << mut << ": " << exp(normScore) << endl;
//		}
//		//getchar();
//		normScoreTable.push_back(normScores);
//	}



//	for(int mut=0; mut<normScoreTable.size(); mut++){
//		for(int node=0; node<normScoreTable.at(mut).size(); node++){
//			cout << "mut " << mut << ": " << exp(normScoreTable.at(mut).at(node)) << " ";
//		}
//		cout << "\t";
//	}
//
//	cout << "done" << endl;

	//vector<double> separation;
//	int leafCount = m;
//
//	int count = 0;
//	stringstream sepLine;
//	for(int v=0; v<leafCount-1; v++){                 // all leaf pairs
//		for(int w=v+1; w<leafCount; w++){
//			//cout << "lca of " << v << " and " << w << endl;
//			int lca = getLCA(v, w, ancMatrix, m);    // lca of the two leafs
//			//cout << "lca: " << lca << endl;
//			//getchar();
//			sepLine << leafPairs.at(count++);
//
//			double v_expPathMutCount = getExpPathMutCount(v, lca, parentVector, normScoreTable);
//			double w_expPathMutCount = getExpPathMutCount(w, lca, parentVector, normScoreTable);
//			sepLine << "\t" << max(v_expPathMutCount,w_expPathMutCount) << "\t" << v_expPathMutCount << "\t" << w_expPathMutCount;
//			for(int mut=0; mut<table.size(); mut++){
//
//				double v_sep = mutSeparationOdds(mut, v, lca, parentVector, normScoreTable);
//				double w_sep = mutSeparationOdds(mut, w, lca, parentVector, normScoreTable);
//				double maxSep = max(v_sep,w_sep);
//				sepLine  << "\t" <<  maxSep;
//			}
//			sepLine << "\n";
//		}
//	}
//	//cout << sepLine.str();
//	//getchar();
//	return sepLine.str();
//}
//
//double mutSeparationOdds(int mut, int v, int lca, int* parentVector, vector<vector<double> >& normScoreTable){
//
//	double logSepProb = 0.0;
//	int parent = v;
//
//	while(parent !=lca){
//				//cout << "node: " << parent << endl;
//		logSepProb+= normScoreTable.at(mut).at(parent);
//		parent = parentVector[parent];
//	}
//			//cout << "done" << endl;
//	//		if(lca != 108){
//	//			getchar();
//	//		}
//	return logSepProb/(1-logSepProb);
//}
//
//
///**  Computes the expected number of mutations on a path from a leaf to an **/
///**  inner node (lca) based on the relative placement weights.  **/
//double getExpPathMutCount(int v, int lca, int* parentVector, vector<vector<double> >& normScoreTable){
//
//	double pathExpMutCount = 0.0;
//
//	for(int mut=0; mut<normScoreTable.size(); mut++){
//		double pathMutFrac = 0.0;
//		int parent = v;
//
//		while(parent !=lca){
//			//cout << "node: " << parent << endl;
//			pathMutFrac+= exp(normScoreTable.at(mut).at(parent));
//			parent = parentVector[parent];
//		}
//		//cout << "done" << endl;
////		if(lca != 108){
////			getchar();
////		}
//		pathExpMutCount += pathMutFrac;
//	}
//	return pathExpMutCount;
//}
//
//int getLCA(int v, int w, bool** ancMatrix, int m){
//
//	int root = 2*m-2;            // root is default LCA
//	int currLCA = root;
//	//cout << "initLCA = " << currLCA << "  for " << v << " and " << w << endl;
//	for(int node=0; node<2*m-2; node++){
//		if(ancMatrix[node][v] && ancMatrix[node][w]){   // common ancestor
//			if(currLCA==root || ancMatrix[currLCA][node]){             // now also lowest common ancestor
//				currLCA = node;
//				//cout << "new LCA: " << currLCA << endl;
//				//getchar();
//			}
//		}
//	}
//	return currLCA;
//}



int readParameters(int argc, char* argv[]){
	for (int i = 1; i < argc; ++i) {

		if (strcmp(argv[i], "-i") == 0) {                      // name of data count file
			if (i + 1 < argc) {
				dataFileName = argv[++i];
				cout << dataFileName << endl;
			}
		}else if (strcmp(argv[i], "-description")==0) {           // name of file describing the samples (i.e. sample_nodeDescription.tsv)
			if (i + 1 < argc) { descriptionFileName = argv[++i];}
		}else if (strcmp(argv[i], "-samples")==0) {           // name of file containing the MCMC samples
					if (i + 1 < argc) { sampleFileName = argv[++i];}
		}
	}
	return 0;
}

vector<vector<int> > readTableFile(string fileName){
    fstream in;
    in.open(fileName);
    string line;
    vector<std::vector<int> > v;
    int i = 0;
    cout << "opening read count file: " << fileName << endl;
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
       // cout << line << endl;
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
    	i++;
    }


//    for(int i=0; i<v.size(); i++){
//    	cout << i << "\t";
//    	for(int j=0; j<v.at(i).size(); j++){
//    		cout << v.at(i).at(j) << "\t";
//    	}
//    	cout << endl;
//    }
    cout << i << " lines read" << endl;
    //getchar();
    in.close();
    return v;
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

vector<vector<string> > readStringTableFile(string fileName){
    std::fstream in(fileName);
    std::string line;
    std::vector<std::vector<string> > v;
    int i = 0;

    while (std::getline(in, line))
    {
        string value;

        cout << "line      : " << line << endl;

        std::string delimiter = "\t";

        size_t pos = 0;
        std::string token;
        v.push_back(vector<string>());
        while ((pos = line.find(delimiter)) != std::string::npos) {
            token = line.substr(0, pos);
            value = token;
            v.at(i).push_back(value);
            //std::cout << ":: " << token << std::endl;
            line.erase(0, pos + delimiter.length());
        }
        v.at(i).push_back(line);
       //std::cout << "read: " << line << std::endl;
        ++i;
       // cout << "line size: " << v[i].size() << endl;
    }

    return v;
}

vector<string> getLeafLabels(vector<vector<string> >& sampleInfo){
	vector<string> leafLabels;
	cout << "getting leaf labels\n";
	for(int i=0; i<sampleInfo.size(); i++){

//		cout << i <<  ": element count: " << sampleInfo.at(i).size() << endl;
//		for(int p=0; p<sampleInfo.at(i).size(); p++){
//			cout << "sample info: " << sampleInfo.at(i).at(p) << "\t";
//		}
//		cout << endl;
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

vector<vector<int> > getExpAlleleCountFreqs(int m, int sampleCount, bool** ancMatrix, vector<int>& leafClusterId){

	vector<vector<int> > expVarReadCount = getExpVarAlleleCount(m, sampleCount, ancMatrix, leafClusterId);
	vector<vector<int> > expVarReadCountFreqs;

//	for(int i=0; i< expVarReadCount.size(); i++){
//		for(int j=0; j< expVarReadCount.at(i).size(); j++){
//			cout << expVarReadCount.at(i).at(j) << " ";
//		}
//		cout << endl;
//	}
	for(int sample=0; sample<sampleCount; sample++){

		vector<int> freqCounter;
		int maxExpVarReadCounter = 0;
		for(int i=0; i<expVarReadCount.at(sample).size(); i++){
			maxExpVarReadCounter = max(maxExpVarReadCounter, expVarReadCount.at(sample).at(i));
		}
		for(int i=0; i<=maxExpVarReadCounter; i++){
			freqCounter.push_back(0);                    // set frequency of each expected allele frequency to zero
		}

		for(int node=0; node<(2*m)-1; node++){                        // increase frequency counter of expected frequency
			freqCounter.at(expVarReadCount.at(sample).at(node))+=1;  // found at this node by one
		}
		expVarReadCountFreqs.push_back(freqCounter);

//		int sum = 0;
//		cout << "sample " << sample << ":  ";
//		for(int j=0; j<freqCounter.size(); j++){
//			cout << j << " expected: " << freqCounter.at(j) << " times\t|\t";
//			sum += freqCounter.at(j);
//		}
//		cout << endl;
//		if(sum != (2*m)-1){
//			cout << "Error" << endl;
//			getchar();
//		}
	}
	return expVarReadCountFreqs;
}


/** this computes for each node in the current tree how many wbc's are below it (including node itself)**/
vector<int> getWBC_count_below(int m, int sampleCount, bool** ancMatrix, vector<bool>& wbcStatus){
	vector<int> wbcBelowCount;

	for(int node=0; node<(2*m)-1; node++){  // init vector of wbc counts
		wbcBelowCount.push_back(0);
	}

	for(int leaf=0; leaf<m; leaf++){
		if(wbcStatus.at(leaf)==1){                                // for each wbc leaf
			for(int node=0; node<(2*m)-2; node++){
				if(ancMatrix[node][leaf]==1){
					wbcBelowCount.at(node) += 1;      // add to wbc counter for all nodes
				}
			}
		}
	}
	return wbcBelowCount;
}

/* precompute for each pair (sample, mutation placement) the number of expected mutated alleles in the sample */
vector<vector<int> > getExpVarAlleleCount(int m, int sampleCount, bool** ancMatrix, vector<int>& leafClusterId){

	vector<vector<int> > expVarReadCount;
	for(int sample=0; sample<sampleCount; sample++){
		expVarReadCount.push_back(vector<int>());
		for(int node=0; node<(2*m)-1; node++){
			expVarReadCount.at(sample).push_back(0);
		}
	}

	for(int leaf=0; leaf<m; leaf++){
		int clusterId = leafClusterId[leaf];
		for(int node=0; node<(2*m)-2; node++){
			if(ancMatrix[node][leaf]==1){
				expVarReadCount[clusterId][node] += 1;   // count how many leafs belonging to the sample have
			}                                            // the node where the mutation is placed as ancestor
		}
		expVarReadCount[clusterId][(2*m)-2] += 1;
	}
	return expVarReadCount;
}


/* initializes constant data structures used in the tree scoring */
void computeMaxAlleleCount(vector<vector<string> >& sampleInfo){
	maxAlleleCount = 0;
	for(int i=0; i<sampleInfo.size(); i++){
		if(2*stoi(sampleInfo.at(i).at(1)) > maxAlleleCount){
			maxAlleleCount = 2*stoi(sampleInfo.at(i).at(1));
		}
	}
}



double logBinomCoeff(int n, int k){

//	cout << "n = " << n << "  k = " << k << endl;
//
//	cout << "lgamma(" << n+1 << ") = " << lgamma(n+1.0) << endl;
//	cout << "lgamma(" << k+1 << ") = " << lgamma(k+1.0) << endl;
//	cout << "lgamma(" << n+1-k << ") = " << lgamma(n+1.0-k) << endl;
//	cout << "total binom coeff: " << lgamma(n+1.0)-(lgamma(k+1.0)+lgamma(n+1.0-k)) << endl;
//	if(k>n){
//			getchar();
//		}
	return lgamma(n+1.0)-(lgamma(k+1.0)+lgamma(n+1.0-k));
}

/* computes the probability mass function of the beta-binomial distribution */
double logBetaBinom(int k, int r, int alpha, int beta, int p, int q, double tau){

	double logScore = 0.0;
	double M = ((alpha-p)+(beta-q))*tau;
	double mu = (1.0*alpha-p)/((alpha-p)+(beta-q));

//	cout << "M  = " << M << endl;
//	cout << "mu = " << mu << endl;
//	cout << r << " chose " << k << endl;
//	cout << "add:    " << lgamma(M) << " + " << logBinomCoeff(r,k) << " + " << lgamma(k+(M*mu)) << " + " << lgamma((r-k)+M*(1-mu)) << endl;
//	cout << "substr: " << lgamma(M*mu) << "+" << lgamma(M*(1-mu)) << "+" << lgamma(r+M) << endl;

	logScore += lgamma(M) - (lgamma(M*mu) + lgamma(M*(1-mu)));
	logScore += lgamma(r+1.0)-(lgamma(k+1.0)+lgamma(r+1.0-k));
	logScore += lgamma(k + M*mu) + lgamma((r-k)+ M*(1-mu)) - lgamma(r+M);

	return logScore;
}

double computeLogSampleScore(int expCount, int alleleCount, int obsMutCount, int obsTotalCount, double dropoutRate, double seqErr, double tau){

	int alpha = expCount;
	int beta = alleleCount - alpha;
	int k = obsMutCount;
	int r = obsTotalCount;
	double delta = dropoutRate;


	vector<double> logScores;
	double logSumScore;
	double maxLogScore = -DBL_MAX;

//	cout << "alpha: " << alpha << endl;
//	cout << "beta:  " << beta << endl;
//	if(alpha>20){
//		getchar();
//	}

//	cout << "--------------------------------\n";
//	cout << "allele count: " << alleleCount << endl;
//	cout << "mutated = " << alpha << endl;
//	cout << " normal = " << beta << endl;
//	cout << "    read count = " << r << endl;
//	cout << "var read count = " << k << endl;

	if(r==0){
		return 0;
	}
	// both alleles still present: use beta-binomial for modeling drop-out
	for(int p=0; p<alpha; p++){
		//cout << "p " << p << endl;
		for(int q=0; q<beta; q++){
			//cout << "q " << q << endl;
			double logScore = -DBL_MAX;
			//double test_1 = logBinomCoeff(alpha, p);
			//cout << "part 1: " << test_1 << endl;
			//double test_2 = logBinomCoeff(beta,q);
			//cout << "part 2: " << test_2 << endl;
			//double test_3 = logBetaBinom(k, r, alpha, beta, p, q, tau);
			//cout << "part 3: " << test_3 << endl;
			logScore = logBinomCoeff(alpha, p) + logBinomCoeff(beta,q) + (p+q)*log(delta) + ((alpha-p)+(beta-q)) * log(1-delta) + logBetaBinom(k, r, alpha, beta, p, q, tau);
			//cout << "        (i)" << logScore << endl;
			logScores.push_back(logScore);
			maxLogScore = max(maxLogScore, logScore);
		}
	}
	// only mutated alleles amplified: use binomial for modeling sequencing errors
	for(int p=0; p<alpha; p++){
		double logScore = logBinomCoeff(alpha, p) + (beta+p)*log(delta) + (alpha-p)*log(1-delta) + logBinomCoeff(r, k) + (r-k)*log(seqErr) + k* log(1-seqErr);
		logScores.push_back(logScore);
		//cout << "       (ii)" << logScore << endl;
		maxLogScore = max(maxLogScore, logScore);
	}
	//only normal allele amplified: use binomial for modeling sequencing errors
	for(int q=0; q<beta; q++){
		double logScore = logBinomCoeff(beta, q) + (alpha+q)*log(delta) + (beta-q)*log(1-delta) + logBinomCoeff(r, k) + k*log(seqErr) + (r-k)* log(1-seqErr);
		logScores.push_back(logScore);

//		cout << "logBinomCoeff(beta, q) = " << logBinomCoeff(beta, q) << endl;
//		cout << "(alpha+q)*log(delta) = " << (alpha+q)*log(delta) << endl;
//		cout << "(beta-q)*log(1-delta) = " << (beta-q)*log(1-delta) << endl;
//		cout << "logBinomCoeff(r, k) = " << logBinomCoeff(r, k) << endl;
//		cout << "k*log(seqErr) = " << k*log(seqErr) << endl;
//		cout << "(r-k)* log(1-seqErr) = " << (r-k)* log(1-seqErr) << endl;
//
//
//		cout << "      (iii)" << logScore << endl;
		maxLogScore = max(maxLogScore, logScore);
	}
	// these scores need to be summed -> go from log space to normal space
	double sumScore = 0.0;
	for(int i=0; i<logScores.size(); i++){              // sum over all allele states, exp is necessary as scores are actually log scores
		sumScore += exp(logScores.at(i)-maxLogScore);   // subtraction of best score to calculate with score differences (smaller values)
	}
	logSumScore = log(sumScore)+maxLogScore;      // transform back to log scores and change from score differences to actual scores
	//getchar();
	if(k>r){
		cout << "k > r: " << k << " > " << r<< endl;
		getchar();
	}
	return logSumScore;
}

int sampleEvent(std::vector<double> prob){ // picks an event based on the given the event probabilities, given probs need to sum to 1

    double percent = (rand() % 100)+1;    // between 1 and 100
    double probSum = prob[0];
//    for(int i=0; i<prob.size()-1; i++){
//    	cout << prob.at(i) << "\t";
//    }
//    cout << endl;
    for(int i=0; i<prob.size()-1; i++){
        if(percent <= probSum*100){
          return i;
        }
        probSum += prob[i+1];
    }
    return prob.size()-1;
}

vector<int> getMutPlacement(int m, int n, int sampleCount, bool** ancMatrix, vector<int>& alleleCount, vector<int>& leafClusterId, vector<vector<int> >& mutReadCounts, vector<vector<int> >& totalReadCounts, double dropoutRate, double seqErr, double tau, vector<bool>& wbcStatus){

	vector<vector<int> > expVarAlleleCount = getExpVarAlleleCount(m, sampleCount, ancMatrix, leafClusterId);
	vector<vector<int> > expVarAlleleFreqs = getExpAlleleCountFreqs(m, sampleCount, ancMatrix, leafClusterId);
	vector<int> wbcBelowCount = getWBC_count_below(m, sampleCount, ancMatrix, wbcStatus);

	double logScore = 0.0;

	vector<int> mutPlacement;
	int parentVectorSize = (2*m)-2;
	for(int mut=0; mut< n; mut++){               // compute score separately for each mutation
		//cout << "mut " << mut << " of " << n << endl;
		vector<double> logAttachmentScores;
		double bestLogAttachmentScore = -DBL_MAX;
		//bestPlacementPoint.push_back(-1);

		/***  pre-compute scores per expected number of mutated alleles ***/
		vector<vector<double> > scorePrecomp;
		for(int sample=0; sample<sampleCount; sample++){
			vector<double> sampleScorePrecomp;
			//cout << "sample: " << sample << endl;
			for(int expFreq=0; expFreq<expVarAlleleFreqs.at(sample).size(); expFreq++){
				//cout << "exp freq (" << mut << "): " << expFreq << "   allele count at sample " << alleleCount.at(sample) << " (" << mutReadCounts.at(mut).at(sample) << "," << totalReadCounts.at(mut).at(sample) << ")" << endl;
				double newScore = computeLogSampleScore(expFreq, alleleCount.at(sample), mutReadCounts[mut][sample], totalReadCounts[mut][sample], dropoutRate, seqErr, exp(tau));
				sampleScorePrecomp.push_back(newScore);
				//cout << "done" << endl;
			}
			scorePrecomp.push_back(sampleScorePrecomp);
		}

		//cout << "end 2" << endl;
		for(int att=0; att< (2*m)-1; att++){           // for each attachment point
			//cout << "wbcBelowCount.at(" << att << ") = " << wbcBelowCount.at(att) << endl;

			double logAttachmentScore = 0.0;
			double wbc_penalty = chi_wbc * wbcBelowCount.at(att);
			for(int sample=0; sample<sampleCount; sample++){   // for each sample
				//cout << "sample " << sample << endl;
				int expCount = expVarAlleleCount[sample][att];
				//cout << "expCount: " << expCount << "\n";
//				cout << alleleCount.at(sample) << endl;
//				cout << mut << ", " << sample << endl;
//				cout << mutReadCounts.size() << endl;
//				cout << mutReadCounts.at(0).size() << endl;
//				cout << mutReadCounts[mut][sample] << endl;

				double newScore = scorePrecomp.at(sample).at(expCount);
//				double newScoreOld = computeLogSampleScore(expCount, alleleCount.at(sample), mutReadCounts[mut][sample], totalReadCounts[mut][sample], dropoutRate, seqErr, tau);
//				if(newScore!=newScoreOld){
//					cout << "different scores" << endl;
//					getchar();
//				}
				logAttachmentScore += newScore;
				//cout << newScore << endl;
			}
			//cout << "\n mut score: " << logAttachmentScore << endl;
			logAttachmentScores.push_back(logAttachmentScore - wbc_penalty);
			bestLogAttachmentScore = max(bestLogAttachmentScore, logAttachmentScore);
		}

		double totalScore = 0.0;
		for(int node=0; node<=parentVectorSize;node++){
			totalScore += exp(logAttachmentScores.at(node));
		}
		//cout << "sumScore: " << sumScore << endl;

		vector<double> nodeProb;
		for(int node=0; node<=parentVectorSize;node++){
			//cout << exp(logAttachmentScores.at(node)) << " / " << totalScore << endl;
			nodeProb.push_back(exp(logAttachmentScores.at(node))/totalScore);
		}
		int sampledNode = sampleEvent(nodeProb);
		mutPlacement.push_back(sampledNode);
		//getchar();
	}

//	for(int i=0; i< mutPlacement.size(); i++){
//		cout << mutPlacement.at(i) << "\t";
//	}
//	cout << endl;
	return mutPlacement;
}

/* transforms a parent vector to an ancestor matrix*/
bool** parentVector2ancMatrix(int* parent, int n){
	bool** ancMatrix = init_boolMatrix(n, n, false);
	int root = n;
	for(int i=0; i<n; i++){
		int anc = i;
		int its =0;
		while(anc < root){                              // if the ancestor is the root node, it is not represented in the adjacency matrix
			if(parent[anc]<n){
				ancMatrix[parent[anc]][i] = true;
			}

			anc = parent[anc];
			its++;
		}
	}
	for(int i=0; i<n; i++){
		ancMatrix[i][i] = true;
	}
	return ancMatrix;
}

/* transforms a parent vector to an ancestor matrix*/
vector<vector<bool > > parentVector2ancMatrixWithRoot(int* parent, int n){

	int root = n;                                // the root has the highest id

	vector<vector<bool > > ancMatrix;            // ancestor matrix with root represented  (0,1,...,n-1, n)
	for(int node=0; node<=n; node++){            //  set default to false
		vector<bool> temp(n+1,false);
		ancMatrix.push_back(temp);
	}

	for(int node=0; node<=n; node++){
		int anc = node;
		bool done = false;
		while(!done){                              // follow the line of ancestors of each node
			ancMatrix.at(anc).at(node) = true;     // set ancestor states to true
			if(anc==root){                         // once the root is reached, we are done for this node
				done = true;
			}
			else{
				anc = parent[anc];
			}
		}
	}

//	for(int node=0; node<=n; node++){
//		for(int node2=0; node2<=n; node2++){
//			cout << " " << ancMatrix.at(node).at(node2);
//			if(node==node2 && ancMatrix.at(node).at(node2)==false){
//				cout << endl;
//				getchar();
//			}
//		}
//		cout << endl;
//	}
	return ancMatrix;
}

bool** init_boolMatrix(int n, int m, bool value){

    bool** matrix = allocate_boolMatrix(n, m);     // allocate

    for (int i=0; i<n; ++i)             // initialize
    {
         for (int j=0; j<m; ++j)
      {
        	matrix[i][j] = value;
    	}
    }
    return matrix;
}

bool** allocate_boolMatrix(int n, int m){

    bool** matrix = new bool*[n];
    matrix[0] = new bool[n*m];
    for (int i=1; i<n; ++i)
    {
        matrix[i] = matrix[i-1] + m;
    }
    return matrix;
}

bool isThreeCellSample(string str){
	if (str.find("cluster_3-0.png") != string::npos) {
		return true;
	}
	if (str.find("cluster_2-1.png") != string::npos) {
		return true;
	}
	if (str.find("cluster_1-2.png") != string::npos) {
		return true;
	}
	return false;
}

void updateSplitCounts(std::map<string, int> &splitCounts, vector<vector<bool> > &ancMatrix, vector<string> &leafLabels, int nodeCount, int* parent){
	vector<string> subTreeLabels;
	for(int i=0; i<nodeCount; i++){
		vector<string> subTreeLeafs;
		for(int j=0; j<leafLabels.size(); j++){
			if(ancMatrix.at(i).at(j)){
				string label = getActualLabel(leafLabels.at(j));
				label.erase(0, 1);     // remove leading quotation mark
				label.pop_back();      // remove trailing quotation mark
				int length = label.length();
				char last = label.at(length-1);
				char secondLast = label.at(length-2);
				if(isalpha(last) && secondLast == '_'){      // remove last two characters if ends on "_A", "_B", etc.
					label.pop_back();
					label.pop_back();
				}
				subTreeLeafs.push_back(label);
			}
		}
		sort(subTreeLeafs.begin(), subTreeLeafs.end());
		string subTreeLabel = "";
		for(int j=0; j<subTreeLeafs.size(); j++){
			string sep = ", ";
			if(j>0){
				subTreeLabel += ", ";
			}
			subTreeLabel += subTreeLeafs.at(j);
		}
		subTreeLabels.push_back(subTreeLabel);
	}

	for(int i=m; i<nodeCount; i++){

		vector<int> childNodes;
		for(int j=0; j<nodeCount; j++){
			if(parent[j] == i){
				childNodes.push_back(j);
			}
		}

		if(childNodes.size()!=2){
			cerr << "not a binary tree" << endl;
		}

		string temp1 = subTreeLabels.at(childNodes.at(0)) + " | " + subTreeLabels.at(childNodes.at(1)) + " | all other samples";
		string temp2 = subTreeLabels.at(childNodes.at(1)) + " | " + subTreeLabels.at(childNodes.at(0)) + " | all other samples";

		if(splitCounts.find(temp1) != splitCounts.end()){ // if split is already in list increment count
			splitCounts[temp1] += 1;
		}
		else if(splitCounts.find(temp2) != splitCounts.end()){ // if split is already in list (with swapped branches) increment count
			splitCounts[temp2] += 1;
		}
		else{
			splitCounts.insert ( std::pair<string,int>(temp1,1) );  // else insert into list and increment
		}

	}
}


// This is the version where we look at the subtree leaf labels
//		if(splitCounts.find(subTreeLabel) != splitCounts.end()){ // if split is already in list increment count
//			splitCounts[subTreeLabel] += 1;
//		}
//		else{
//			splitCounts.insert ( std::pair<string,int>(subTreeLabel,1) );  // else insert into list and increment
//		}



