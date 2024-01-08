/*
 * CTC_treeScoring.h
 *
 *  Created on: Jan 22, 2018
 *      Author: jahnka
 */

#ifndef CTC_TREESCORING_H_
#define CTC_TREESCORING_H_

using namespace std;


//std::string runMCMCnew(vector<struct treeTheta>& bestTrees, int noOfReps, int noOfLoops, double gamma_, vector<double> moveProbs, int n, int m, int sampleCount, vector<int>& alleleCount, vector<int>& leafClusterId, vector<vector<int> >& mutReadCounts, vector<vector<int> >& totalReadCounts, double dropoutRate, double seqErr, double tau, vector<string>& label, string treeFile);
std::string runMCMCnew(int noOfReps, int noOfLoops, double gamma_, vector<double> moveProbs, int n, int m, int sampleCount, vector<int>& alleleCount, vector<int>& leafClusterId, vector<vector<int> >& mutReadCounts, vector<vector<int> >& totalReadCounts, double dropoutRate, double seqErr, double tau, vector<string>& label, vector<string>& treeFiles, vector<string> mutLabels, vector<bool>& wbcStatus, string prefix, int sampleStep, int burnin);
vector<vector<int> > getExpVarAlleleCount(int m, int sampleCount, bool** ancMatrix, vector<int>& leafClusterId);
int** getBinomCoeffTable(int maxAlleleCount);
//double scoreTree(int sampleCount, int attachmentPoint, bool** ancMatrix,  vector<int> alleleCount, vector<int> leafClusterId);
//double scoreTree(int m, int n, int sampleCount, bool** ancMatrix, vector<int>& alleleCount, vector<int>& leafClusterId, vector<vector<int> >& mutReadCounts, vector<vector<int> >& totalReadCounts, double dropoutRate, double seqErr, double tau, vector<int>& bestPlacementPoint, vector<bool>& wbcStatus);
double scoreTree(int m, int n, int sampleCount, bool** ancMatrix, vector<int>& alleleCount, vector<int>& leafClusterId, vector<vector<int> >& mutReadCounts, vector<vector<int> >& totalReadCounts, double dropoutRate, double seqErr, double tau, vector<int>& bestPlacementPoint, vector<bool>& wbcStatus);

void computeMaxAlleleCount(vector<vector<string> >& sampleInfo);
double logBetaBinom(int k, int r, int alpha, int beta, int p, int q, double tau);
double computeLogSampleScore(int expCount, int alleleCount, int obsMutCount, int obsTotalCount, double dropoutRate, double seqErr, double tau);
double logBinomCoeff(int n, int k);
vector<vector<int> > getExpAlleleCountFreqs(int m, int sampleCount, bool** ancMatrix, vector<int>& leafClusterId);
vector<int> getWBC_count_below(int m, int sampleCount, bool** ancMatrix, vector<bool>& wbcStatus);
double getExpPathMutCount(int v, int lca, int* parentVector, vector<vector<double> >& normScoreTable);
string getSeparationScores(vector<vector<double> > table, bool** ancMatrix, int n, int m, int* parentVector, vector<string> leafPairs);
vector<vector<double> > getPlacementScoreTable(int m, int n, int sampleCount, bool** ancMatrix, vector<int>& alleleCount, vector<int>& leafClusterId, vector<vector<int> >& mutReadCounts, vector<vector<int> >& totalReadCounts, double dropoutRate, double seqErr, double tau, vector<int>& bestPlacementPoint, vector<bool>& wbcStatus);
int getLCA(int v, int w, bool** ancMatrix, int n);
double mutSeparationOdds(int mut, int v, int lca, int* parentVector, vector<vector<double> >& normScoreTable);
double proposeNewLogTau(double currLogTau, int currTauCount, double currLogTausSum, double currLogTausSquaredSum);
#endif /* CTC_TREESCORING_H_ */
