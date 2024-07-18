/*
 * binTree_output.cpp
 *
 *  Created on: Jan 15, 2018
 *      Author: jahnka
 */

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <float.h>
#include <algorithm>

#include "output.h"
#include "scoreTree.h"
#include "matrices.h"
#include "trees.h"
#include "binTree_output.h"


using namespace std;

vector<string> nodeColors;

/* writes the given string to file */
void writeToFile(string content, string fileName){
	ofstream outfile;
	outfile.open (fileName.c_str());
	outfile << content;
	outfile.close();
}

/* Score contribution by a specific mutation when placed at the root, that means all samples should have it */
/* This is the same for all trees and can be precomputed */
double binTreeRootScore(int** obsMutProfiles, int mut, int m, double ** logScores){
	double score = 0.0;
	for(int sample=0; sample<m; sample++){
		score += logScores[obsMutProfiles[sample][mut]][1];
	}
	return score;
}

/* computes the best placement of a mutation, the highest one if multiple co-optimal placements exist*/
int getHighestOptPlacement(int** obsMutProfiles, int mut, int m, double ** logScores, bool** ancMatrix){

	int nodeCount = (2*m)-1;
	int bestPlacement = (2*m)-2;   // root
	double bestPlacementScore = binTreeRootScore(obsMutProfiles, mut, m, logScores);
	//cout << bestPlacementScore << " (root)\n";
	//print_boolMatrix(bool** array, int n, int m);
	for(int p=0; p<nodeCount-1; p++){                           // try all possible placements (nodes in the mutation tree)

		double score = 0.0;                   // score for placing mutation at a specific node
		for(int sample=0; sample<m; sample++){
			//cout << p << " " << sample << "\n";
			if(ancMatrix[p][sample] == 1){
				score += logScores[obsMutProfiles[sample][mut]][1]; // sample should have the mutation
			}
			else{
				score += logScores[obsMutProfiles[sample][mut]][0]; // sample should not have the mutation
			}
		}
		if(score > bestPlacementScore){
			bestPlacement = p;
			bestPlacementScore = score;
			//cout << bestPlacementScore << " (non-root)\n";
		}
		else if (score == bestPlacementScore && ancMatrix[p][bestPlacement] == true){
			bestPlacement = p;
		}
	}

	//if(bestPlacement == (2*m)-2){
	//	cout<< "best placed at root\n";
	//	getchar();
	//}
	return bestPlacement;
}

/* computes the best placement of a mutation, the highest one if multiple co-optimal placements exist*/
int* getHighestOptPlacementVector(int** obsMutProfiles, int n, int m, double ** logScores, bool** ancMatrix){
	int* bestPlacements = init_intArray(n, -1);
	for(int mut=0; mut<n; mut++){
		bestPlacements[mut] = getHighestOptPlacement(obsMutProfiles, mut, m, logScores, ancMatrix);
	 }
	//print_intArray(bestPlacements, n);
	return bestPlacements;
}

vector<string> getBinTreeNodeLabels(int nodeCount, int* optPlacements, int n, vector<string> geneNames, vector<string>& sampleNames){
	vector<string> v;
	int count = 0;
	int sampleCount = (nodeCount+1)/2;
	for(int i = 0; i < nodeCount; i++){
		v.push_back("");
	}

	for(int mut=0; mut<n; mut++){
		string toAppend;
		int lastBreak = v.at(optPlacements[mut]).find_last_of("\n");
		string lastLine = v.at(optPlacements[mut]).substr(lastBreak+1);
		if(v.at(optPlacements[mut]) == ""){
			toAppend = geneNames.at(mut);
			count++;
		}
		else if(lastLine.length()>40){
			toAppend = "\n " + geneNames.at(mut);
			count++;
		}
		else{
			toAppend = " " + geneNames.at(mut);
			count++;
		}
		//cout << "        " << j << "\n";
		//cout << "                     "<< optPlacements[j] << "\n";

		v.at(optPlacements[mut]) += toAppend;
	}
	if(v.at(nodeCount-1) == ""){
		v.at(nodeCount-1) = "root";
	}

	for(int i=0; i< sampleCount; i++){
		if(v.at(i)==""){
			v.at(i) += sampleNames.at(i);
		}
	}
	for(int i = 0; i < nodeCount; i++){
		if(v.at(i).find(" ") != string::npos){
			v.at(i) = "\"" + v.at(i) + "\"";
		}
	}

//	for(int i = 0; i < nodeCount; i++){
//		int mutCount = 1;
//		mutCount += (int)std::count(v.at(i).begin(),v.at(i).end(),',');
//		if(mutCount >=10){
//			v.at(i) = "\"+" + std::to_string(mutCount) + " " + std::to_string(i) + "\"";
//		}
//	}

	//cout << "added mutations " << count << "\n";

	return v;
}

/* returns the lca of a node that has a non-empty label, the root is assumed to always have a label */
int getLcaWithLabel(int node, int* parent, vector<string> label, int nodeCount){
	int root = nodeCount -1;
	int p = parent[node];;
	while(p != root && label[p]==""){
		p = parent[p];
	}
	return p;
}


/* create the content of the GraphViz output file for a mutation tree */
string getBinTreeGraphVizString(int* parentVector, int parentVectorSize, vector<string>& nodeLabels, vector<string>& sampleNames){

	cout << "..............................................\n";
	for(int i=0; i<parentVectorSize; i++){
		cout << nodeLabels[i] << endl;
	}

	stringstream str;
	int leafCount = (parentVectorSize+1)/2;
	str << "digraph g{\n";
	str << "node [fontname = \"helvetica\"]";
	str << "node [color=" << nodeColors.at(0) << ", style=filled, fontcolor=black, shape=circle];  \n";   // sample nodes
	for(int i=0; i<leafCount; i++){
		if(nodeLabels[i]!=sampleNames[i]){
			str << "node [color=deeppink4, style=filled, fontcolor=white, shape=box];	\n";
			str << i << "[label=" << nodeLabels[i] << "];\n";
			str << "node [color=lightgrey, style=filled, fontcolor=black, shape=circle];  \n";   // sample nodes
			str << "L" << i << "[label=" << sampleNames[i] << "];\n";
			str << i << " -> " << "L" << i << ";\n";

		}
		str << i << "[label=" << nodeLabels[i] << "];\n";
	}

	str << "node [color=deeppink4, style=filled, fontcolor=white, shape=box];	\n";   // style of tree nodes
	for(int i=leafCount; i<nodeLabels.size(); i++){
		str << i << "[label=" << nodeLabels[i] << "];\n";
	}

	for(int i=0; i<parentVectorSize; i++){
		str << parentVector[i] << " -> " << i << ";\n";
	}
	str << "}\n";
	return str.str();
}


void updateInnerNodeLabels(vector<string>& label, int n, int m, vector<int>& bestPlacementPoint){

	vector<int> counter;
	for(int i=0; i<2*m-1; i++){
		counter.push_back(0);
	}
	//cout << n << " mutations to be placed" << endl;
	for(int i=0; i<n; i++){
		counter.at(bestPlacementPoint[i])++;
	}
	for(int i=m; i<2*m-1; i++){
		//cout << label.at(i) << ": " << counter.at(i) << endl;
		label.at(i) = "+" + to_string(counter.at(i));
	}

//	for(int i=0; i<2*m-1; i++){
//		cout << label.at(i) << ": " << counter.at(i) << endl;
//	}
}

void updateInnerNodeLabelsMutIds(vector<string>& label, int n, int m, vector<int>& bestPlacementPoint, vector<string>& mutIds){

	for(int i=m; i<2*m-1; i++){
		label.at(i) = "";
	}
	//cout << n << " mutations to be placed" << endl;
	for(int mut=0; mut<n; mut++){
		if(bestPlacementPoint.at(mut)>=m){
			stringstream newLabel;
			int updateNode = bestPlacementPoint.at(mut);
			newLabel << label.at(updateNode) << mutIds.at(mut) << "\n";
			label.at(updateNode) = newLabel.str();
		}
	}
//	for(int i=0; i<2*m-1; i++){
//		cout << label.at(i) << ": " << counter.at(i) << endl;
//	}
}
string getGraphVizBinTree(int* parents, int m, vector<string>& label, double score, double dropoutRate, double seqErrorRate, vector<int>& clusterId, vector<int>& bestPlacementPoint, vector<string>& mutLabels){
	std::stringstream content;
	int currColor = 0;
	content << "digraph G {\n";
	//content << "label= \"Score = " << score << "\";\n";
	//content << "label= \"Score = " << score << "\";\n";

	content << "node [color=black, shape=box, fontcolor=black, fontsize=18, fontname=Verdana],fontcolor=black;\n";
	content << "scoreBox [label= \"Tree score = " << score << "\n Dropout rate = " << dropoutRate << "\n Sequencing error rate = " << seqErrorRate << "\"];\n";
	//updateInnerNodeLabels(label, bestPlacementPoint.size(), m, bestPlacementPoint);
	updateInnerNodeLabelsMutIds(label, bestPlacementPoint.size(), m, bestPlacementPoint, mutLabels);

	content << "node [color=black, style=empty, fontcolor=black, shape=oval];\n";
	for(int i=m; i<(2*m)-1; i++){
		content << i << " [label=" << "\"" << label[i] << "\"" << "];\n";     // inner nodes
	}

	content << "node [color=" << nodeColors.at(0) << ", style=filled, fontcolor=black, shape=oval];  \n";   // sample nodes
	for(int i=0; i<m; i++){
		//cout << "getting color for " << i << "\n";
		string nodeColor = getNodeColor(i, currColor, clusterId);
	//	cout << nodeColor << endl;
	//	cout << i << " [color=" << nodeColor << ", label=" << label[i] << "];\n";
		content << i << " [color=" << nodeColor << ", label=" << label[i] << "];\n";

	}
	for(int i=0; i<(2*m)-2; i++){
		content << parents[i] << " -> " << i << ";\n";
	}

	content <<  "}\n";
	return content.str();
}

string getPartition(bool** ancMatrix, int leafCount, vector<string>& label){
	stringstream partition;
	for(int i=leafCount; i<(2*leafCount)-2; i++){
		stringstream subtree;
		for(int leaf=0; leaf<leafCount; leaf++){
			//cout << i << "/" << leaf << endl;
			if(ancMatrix[i][leaf]==1){
				string s = label.at(leaf);
				string l_delim = ",label=";
				string r_delim = ",fillcolor";
				string token = s.substr(s.find(l_delim)+8, s.find(r_delim)-(s.find(l_delim)+9));
				int extCount = count(token.begin(), token.end(), '_');
				if (extCount>1){
					string newtoken = token.substr(0, token.length()-2);
					token = newtoken;
				}
				subtree << token << " ";
			}
		}
		partition << subtree.str() << endl;
	}
	return partition.str();
}

string getFancyGraphVizBinTree(int* parents, int m, vector<string>& label, double score, double dropoutRate, double seqErrorRate, vector<int>& clusterId, vector<int>& bestPlacementPoint, vector<string>& mutLabels){
	std::stringstream content;
	int currColor = 0;
	content << "digraph G {\n";
	//content << "label= \"Score = " << score << "\";\n";
	//content << "label= \"Score = " << score << "\";\n";

	content << "node [color=gray, shape=box, fontcolor=black, fontsize=20, fontname=Helvetica];\n";
	content << "scoreBox [label= \"Tree score = " << score << "\n Dropout rate = " << dropoutRate << "\n Sequencing error rate = " << seqErrorRate << "\"];\n";
	//updateInnerNodeLabels(label, bestPlacementPoint.size(), m, bestPlacementPoint);
	updateInnerNodeLabelsMutIds(label, bestPlacementPoint.size(), m, bestPlacementPoint, mutLabels);
	content << "edge [penwidth=4];\n";
	content << "node [color=honeydew4, style=filled, fontcolor=white, shape=oval];\n";
	for(int i=m; i<(2*m)-1; i++){
		content << i << " [label=" << "\"" << label[i] << "\"" << "];\n";     // inner nodes
	}

	//content << "node [color=" << nodeColors.at(0) << ", style=filled, fontcolor=black, shape=oval];  \n";   // sample nodes
	content << "node [fontname=helvetica,fontcolor=black,shape=\"box\",penwidth=10,style=\"rounded,filled,bold\",imagepos=\"bc\",imagescale=true, labelloc=b]; \n";
	for(int i=0; i<m; i++){
		//cout << "getting color for " << i << "\n";
		string nodeColor = getNodeColor(i, currColor, clusterId);
	//	cout << nodeColor << endl;
	//	cout << i << " [color=" << nodeColor << ", label=" << label[i] << "];\n";
	//	content << i << " [color=" << nodeColor << ", label=" << label[i] << "];\n";
		content << i << label[i] << ";\n";
		//cout << label[i] << endl;
		//getchar();

	}
	for(int i=0; i<(2*m)-2; i++){
		content << parents[i] << " -> " << i << ";\n";
	}

	content <<  "}\n";
	return content.str();
}


string getNodeColor(int i, int& currColor, vector<int>& clusterId){

	if(i==0){                                     // first node
		if(clusterId[i]==clusterId[i+1]){
			currColor++;                           // new cluster -> new Color
			return nodeColors.at(currColor);
		}
		return nodeColors.at(0);                                 // default color for non-cluster leafs
	}
	else if(i==clusterId.size()-1){                // last node
		if(clusterId[i]==clusterId[i-1]){
			return nodeColors.at(currColor);                    // same cluster as before
		}
		return nodeColors.at(0);
	}
	else{
		if(clusterId[i]==clusterId[i-1]){      // same cluster as before
			return nodeColors.at(currColor);
		}
		else if(clusterId[i]==clusterId[i+1]){   // new cluster
			currColor++;
			return nodeColors.at(currColor);
		}
		else{
			return nodeColors.at(0);
		}
	}
}


void setClusterColors(){
	nodeColors.push_back("gray");
	nodeColors.push_back("lightcoral");
	nodeColors.push_back("skyblue3");
	nodeColors.push_back("sandybrown");
	nodeColors.push_back("paleturquoise3");
	nodeColors.push_back("thistle");
	nodeColors.push_back("lemonchiffon");
	nodeColors.push_back("darkolivegreen3");
	nodeColors.push_back("lightpink");
	nodeColors.push_back("mediumpurple");
	nodeColors.push_back("darkseagreen3");
	nodeColors.push_back("yellowgreen");
	nodeColors.push_back("gold");
	nodeColors.push_back("indianred1");
	nodeColors.push_back("seagreen2");
	nodeColors.push_back("khaki1");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
	nodeColors.push_back("gray");
}
