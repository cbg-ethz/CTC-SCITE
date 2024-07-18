/*
 * binTree_output.h
 *
 *  Created on: Jan 18, 2018
 *      Author: jahnka
 */

#ifndef BINTREE_OUTPUT_H_
#define BINTREE_OUTPUT_H_

using namespace std;

string getFancyGraphVizBinTree(int *parents, int m, vector<string> &label,
                               double score, double dropoutRate,
                               double seqErrorRate, vector<int> &clusterId,
                               vector<int> &bestPlacementPoint,
                               vector<string> &mutLabels);
void writeToFile(string content, string fileName);
string getBinTreeGraphVizString(int *parentVector, int parentVectorSize,
                                vector<string> &nodeLabels,
                                vector<string> &sampleNames);
double binTreeRootScore(int **obsMutProfiles, int mut, int m,
                        double **logScores);
int getHighestOptPlacement(int **obsMutProfiles, int mut, int m,
                           double **logScores, bool **ancMatrix);
int *getHighestOptPlacementVector(int **obsMutProfiles, int n, int m,
                                  double **logScores, bool **ancMatrix);
vector<string> getBinTreeNodeLabels(int nodeCount, int *optPlacements, int n,
                                    vector<string> geneNames,
                                    vector<string> &sampleNames);
int getLcaWithLabel(int node, int *parent, vector<string> label, int nodeCount);
// std::string getGraphVizBinTree(int* parents, int m, vector<string>& label,
// double score);
string getGraphVizBinTree(int *parents, int m, vector<string> &label,
                          double score, double dropoutRate, double seqErrorRate,
                          vector<int> &clusterId,
                          vector<int> &bestPlacementPoint,
                          vector<string> &mutLabels);
// string getGraphVizBinTree(int* parents, int m, vector<string>& label, double
// score, double dropoutRate, double seqErrorRate, vector<int>& clusterId,
// vector<int>& bestPlacementPoint);
string getNodeColor(int i, int &currColor, vector<int> &clusterId);
void setClusterColors();
void updateInnerNodeLabelsMutIds(vector<string> &label, int n, int m,
                                 vector<int> &bestPlacementPoint,
                                 vector<string> &mutIds);
string getPartition(bool **ancMatrix, int leafCount, vector<string> &label);
#endif /* BINTREE_OUTPUT_H_ */
