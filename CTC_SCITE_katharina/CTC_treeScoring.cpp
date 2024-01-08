/*
 * CTC_treeScoring.cpp
 *
 *  Created on: Jan 22, 2018
 *      Author: jahnka
 */

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <float.h>
#include <math.h>
#include <deque>
#include "binTree_output.h"
#include "scoreTree.h"
#include "matrices.h"
#include "rand.h"
#include "trees.h"
#include "matrices.h"
#include "mcmcBinTreeMove.h"
#include <map>
#include <iomanip>
#include <sstream>
#include <limits>
//#include "binTree_output.h"
//#include "beta_distr.h"
#include "mcmc.h"

#include "CTC_treeScoring.h"

#define UNIFORM_BETA_DISTR_seqErr  false  // if true, use uniform distribution for the prior of seqErr
#define UNIFORM_BETA_DISTR_dropout false  // if true, use uniform distribution for the prior of the dropout rate


using namespace std;

vector<double> logSum;
vector<double> logSquareSum;

int maxAlleleCount = 0;

double chi_wbc = 1.0; // penalty for placing a mutation above a wbc

//int burnIn; // = 300;
//int step;  //= 100;
bool doSampling = true;


/* a new value for the sequencing error is sampled from a normal distribution around the current error rate */
double proposeNewSeqErr(double currRate, double jumpSd){
	double sampledValue = sampleNormal(0, jumpSd);
	double propRate = currRate+sampledValue ;                   //rnorm(1,0,jumpsd)
	if(propRate < 0){
		propRate = fabs(propRate);
	}
	if(propRate > 1){
		propRate = propRate - 2*(propRate-1);
	}
    return propRate;
}

bool changeSeqErr(double prob){
	 double percent = (rand() % 100)+1;    // between 1 and 100
	 if(percent <= prob*100){
		 return true;
	 }
	 return false;
}

char pickParam(){
	double percent = (rand() % 100)+1;
	if(percent <= 50.0){
		return 'd';
	}
//	else if(percent<=66.0){
//		return 's';
//	}
	else{
		return 's';
	}
}

double proposeNewLogTau(double currLogTau, int currTauCount, double currLogTausSum, double currLogTausSquaredSum){
	double first = currLogTausSquaredSum/currTauCount;
	double second = ((currLogTausSum*currLogTausSum)/(currTauCount*currTauCount))/currTauCount;
	double stdev = sqrt(first - second);
	if(first-second < 0){
		cout << "sqrt(" << first << "," << second << ")" << endl;
		//cout << "stdev = sqrt(" << (currLogTausSquaredSum/currTauCount) - ((currLogTausSum*currLogTausSum)/(currTauCount*currTauCount)) << ")" << endl;
		cout << "stdev: " << stdev << "    logTauSum: " << currLogTausSum << "    squared: " << currLogTausSquaredSum << "   count: " << currTauCount << endl;
		getchar();
	}

	cout << "stdev = sqrt(" << (currLogTausSquaredSum/currTauCount) - ((currLogTausSum*currLogTausSum)/(currTauCount*currTauCount)) << ")" << endl;
	cout << "stdev: " << stdev << "    logTauSum: " << currLogTausSum << "    squared: " << currLogTausSquaredSum << "   count: " << currTauCount << endl;
	stdev = max(stdev, 0.1);
	double newLogTau = sampleNormal(currLogTau, stdev*2.38);
	//cout << "proposed tau: " << exp(newLogTau) << endl;
//	for(int i=0; i<logSquareSum.size(); i++){
//			cout << logSquareSum.at(i) << "\t" << logSum.at(i) << endl;
//	}
	//getchar();
	return newLogTau;
}

std::string runMCMCnew(int noOfReps, int noOfLoops, double gamma_, vector<double> moveProbs, int n, int m, int sampleCount, vector<int>& alleleCount, vector<int>& leafClusterId, vector<vector<int> >& mutReadCounts, vector<vector<int> >& totalReadCounts, double dropoutRate, double seqErr, double tau, vector<string>& label, vector<string>& treeFiles, vector<string> mutLabels, vector<bool>& wbcStatus, string prefix, int step, int burnIn){


	//Beta_Distr seqErrPrior(seqErr, 0.1, UNIFORM_BETA_DISTR_seqErr);          // prior distribution sequencing error
	//Beta_Distr dropOutPrior(dropoutRate, 0.1, UNIFORM_BETA_DISTR_dropout);   // prior distribution allelic dropout
	vector<string> leafIds;
	for(int i=0; i<m; i++){
		string copy = label.at(i);
		size_t pos = copy.find("\"");
		string token = copy.substr(0, pos);
		copy.erase(0, pos + 1);
		string name = copy.substr(0, copy.find("\""));
		leafIds.push_back(name);
		cout << label.at(i) << " to " << name << endl;
	}

	//getchar();
	vector<string> leafPairs;
	for(int i=0; i<m-1; i++){
		for(int j=i+1; j<m; j++){
			string pair = leafIds.at(i) + "-" + leafIds.at(j);
			leafPairs.push_back(pair);
			//cout << pair << endl;
		}
	}

	stringstream pairNames;
	for(int i=0; i<leafPairs.size(); i++){
		pairNames << leafPairs.at(i) << "\t";
	}
	pairNames << endl;

	string sampleFile = prefix + "_postSampling.tsv";
	ofstream sampleOutput (sampleFile);
	if (!sampleOutput.is_open()){
		cout << "Unable to open file";
	}


//	for(int i=0; i<leafPairs.size(); i++){
//		cout << i << ": " << leafPairs.at(i) << "\n";
//	}
//	cout << endl;


	if(m != sampleCount){

		cout << "Error: " << m << " !+ " << sampleCount << endl;
	}
	cout << "Note: " << m << " !+ " << sampleCount << endl;
	int parentVectorSize = (2*m)-2;    // transposed case: binary tree, m leafs and m-1 inner nodes, root has no parent
	double bestTreeLogScore = -DBL_MAX;          // log score of T in best (T,beta)
	int*   bestTreeParentVec = nullptr;
	std::map<double,std::string> topTreeList;

	double bestScore = -DBL_MAX;                 // log score of best combination (T, beta)
	double bestSeqErrRate = seqErr;
	double bestDropOutRate = dropoutRate;
	//double bestTau = tau;
	bool sampleStep = false;


	//int** postSamplingMatrix = init_intMatrix(n, m, 0);   // counts for each (mutation,sample) pair how often the
                                                          //  mutation occurs in the sample in the posterior samples
	double currTreeLogScore = -DBL_MAX;
	for(int r=0; r<noOfReps; r++){   // repeat the MCMC, start over with random tree each time, only best score and list of best trees is kept between repetitions

		cout << "MCMC repetition " << r << "\n";
		int*   currTreeParentVec;
		cout << m << " leafs" << endl;
		currTreeParentVec = getRandomBinaryTree(m);                     // transposed case: random binary tree
		double currSeqErrRate = seqErr;
		double currDropOutRate = dropoutRate;
		double currLogTau = log(tau);
		int currTauCount = 1;
		//double currLogTausSum = currLogTau;
		//double currLogTausSquaredSum = currLogTau*currLogTau;
		//logSquareSum.push_back(currLogTausSquaredSum);
		//logSum.push_back(currLogTau);
		bool** currTreeAncMatrix =  parentVector2ancMatrix(currTreeParentVec,parentVectorSize);
		vector<int> bestPlacementPoint;
		for(int mut=0; mut< n; mut++){               // compute score separately for each mutation
			bestPlacementPoint.push_back(-1);
		}
		vector<vector<double> > currTable;
		currTreeLogScore = scoreTree(m, n, sampleCount, currTreeAncMatrix, alleleCount, leafClusterId, mutReadCounts, totalReadCounts, dropoutRate, seqErr, tau, bestPlacementPoint, wbcStatus);
		for(int it=0; it<noOfLoops; it++){                                     // run the iterations of the MCMC
			//cout << "iteration " << it << ": " << endl;
        	if(it % 10000 == 0){
        		cout.precision(16);
        		cout << "At mcmc repetition " << r+1 << "/" << noOfReps << ", step " << it << ": ";
        		cout << "best overall score " << bestScore << ", best tree score " << bestTreeLogScore;
        		cout << endl;
        	}

        	//cout << "it=" << it << "   burnin=" << burnIn << "   step=" << step << endl;
        	if(doSampling && it>=burnIn && it % step == 0){   // set to true if need to sample in this step
        		sampleStep = true;
        	}
        	else{
        		sampleStep=false;
        	}

        	bool isAcceptedMove = false;                               // Is the MCMC move accepted?
        	bool moveChangesParams = changeSeqErr(moveProbs[0]);        // true if this move changes beta, not the tree

        	if(moveChangesParams){                     // new theta is proposed, log scores change, tree is copy of current tree
        		//cout << "change params\n";
        		double propSeqErrRate = currSeqErrRate;
        		double propDropOutRate = currDropOutRate;
        		//double propLogTau = currLogTau;
        		double propLogTau = 0.0;
        		char paramType = pickParam();
        		if(paramType=='s'){
        			propSeqErrRate = proposeNewSeqErr(currSeqErrRate, 0.1*currSeqErrRate);
        		}
        		else if(paramType=='d'){
        			propDropOutRate = proposeNewSeqErr(currDropOutRate, 0.1*currDropOutRate);
        		}
        		//else{
        		//	propLogTau = proposeNewLogTau(currLogTau, currTauCount, currLogTausSum, currLogTausSquaredSum);
        			//cout << "proposed: " << exp(propLogTau) << endl;
        		//}

        		//double propSeqErrScore = beta.logBetaPDF(propSeqErrRate);
        		double propTreeLogScore;

        		propTreeLogScore = scoreTree(m, n, sampleCount, currTreeAncMatrix, alleleCount, leafClusterId, mutReadCounts, totalReadCounts, propDropOutRate, propSeqErrRate, exp(propLogTau), bestPlacementPoint, wbcStatus);

        		//cout << propTreeLogScore << "\t" << currTreeLogScore << "\t" << bestTreeLogScore << endl;

        		isAcceptedMove = sample_0_1() < exp((propTreeLogScore-currTreeLogScore)*gamma_);

        		if (isAcceptedMove){               // the proposed move is accepted
        			currTreeLogScore  = propTreeLogScore;                                       // update score of current tree
        			currSeqErrRate = propSeqErrRate;                                                        // the current AD rate
        			currDropOutRate = propDropOutRate;
        			currLogTau = propLogTau;
        			//logSquareSum.push_back(currLogTausSquaredSum);
        			//logSum.push_back(currLogTau);
        			//currScore = currTreeLogScore;                          // combined score of current tree and current beta
        			//cout << "new tau: " << exp(currLogTau) << endl;
        			//cout << "new param" << endl;
        			//print_intArray(currTreeParentVec, parentVectorSize);
        		}
        		//currTauCount++;
        		//currLogTausSum += currLogTau;
        		//currLogTausSquaredSum += (currLogTau*currLogTau);
        	}
        	else{      // move changed tree

        		//cout << "change tree\n";
        		double nbhcorrection = 1.0;
        		int* propTreeParVec;
        		double propTreeLogScore;
        		propTreeParVec = proposeNextBinTree(moveProbs, m, currTreeParentVec, currTreeAncMatrix);
        		bool** propTreeAncMatrix = parentVector2ancMatrix(propTreeParVec,parentVectorSize);
        		propTreeLogScore = scoreTree(m, n, sampleCount, propTreeAncMatrix, alleleCount, leafClusterId, mutReadCounts, totalReadCounts, currDropOutRate, currSeqErrRate, tau, bestPlacementPoint, wbcStatus);
//        		cout << propTreeLogScore << "\t" << currTreeLogScore << bestTreeLogScore << endl;

        		isAcceptedMove = sample_0_1() < nbhcorrection*exp((propTreeLogScore-currTreeLogScore)*gamma_);

        		if(isAcceptedMove){                                    // the proposed state is accepted
        			free_boolMatrix(currTreeAncMatrix);                                            // discard outdated tree
        			delete[] currTreeParentVec;
        			currTreeAncMatrix = propTreeAncMatrix; //parentVector2ancMatrix(propTreeParVec,parentVectorSize); // update matrix of current tree
        			currTreeParentVec = propTreeParVec;                                         // update parent vector of current tree
        			currTreeLogScore  = propTreeLogScore;                                       // update score of current tree
        		}
        		else{
        			delete [] propTreeParVec;            // discard proposed tree
        			free_boolMatrix(propTreeAncMatrix);
        		}
        		//cout << "new tree" << endl;
        		//print_intArray(currTreeParentVec, parentVectorSize);
        	}

        	/* Update best tree in case we have found a new best one */
			if(currTreeLogScore > bestScore){
				//optStatesAfterBurnIn = 0;                    // new opt state found, discard old count
				bestTreeLogScore = currTreeLogScore;
				delete [] bestTreeParentVec;
				bestTreeParentVec = deepCopy_intArray(currTreeParentVec, parentVectorSize);
				bestScore = bestTreeLogScore;                 // log score of best combination (T, beta)
				bestSeqErrRate = currSeqErrRate;
				bestDropOutRate = currDropOutRate;
				//bestTau = exp(currLogTau);
				cout << it <<  "\tnew best score: " << bestTreeLogScore << endl;
			}

        	if(isAcceptedMove){        // update top tree list if necessary
        		bool insert = false;
        		if(topTreeList.size() < 10){
        			insert = true;
        		}
        		else{
        			auto minScore = topTreeList.begin()->first;
        			if(currTreeLogScore>minScore){
        				topTreeList.erase(minScore);
        				insert = true;
        			}
        		}
        		if(insert){
        			string gv = getFancyGraphVizBinTree(bestTreeParentVec, m, label, bestTreeLogScore,  bestDropOutRate, bestSeqErrRate,leafClusterId, bestPlacementPoint, mutLabels);
        			topTreeList.insert ( std::pair<double,string>(currTreeLogScore,gv) );
        		}
        		//cout << topTreeList.size() << endl;
        		std::map<double,string>::iterator it;// = mymap.begin();
        		int pos = 0;
        		for (it=topTreeList.begin(); it!=topTreeList.end(); ++it){
        			//std::cout << it->first << "  " << pos << '\n';
        			pos++;
        		}
        		int rank = topTreeList.size()-1;
        		for (it=topTreeList.begin(); it!=topTreeList.end(); ++it){ // write all 10 trees to file
        			//writeToFile(topTreeList.at(rank), treeFiles.at(rank));
        			writeToFile(it->second, treeFiles.at(rank));
        			//cout << rank << endl;
        			//std::cout << "here:  " << it->first << "  " << rank+1 << '\n';
        			rank--;
        		}
        		//getchar();


        	}

    		if(sampleStep){
    			//getSeparationScores(vector<vect>)
    			stringstream treeString;
    			for(int node=0; node<2*m-2; node++){
    				if(node!=0){
    					treeString << " ";
    				}
    				treeString << currTreeParentVec[node];

    			}
    			//cout << treeString.str() << endl;
    			//cout << std::setprecision( std::numeric_limits<int>::max() ) << currSeqErrRate << "\t" << currDropOutRate << "\t" << currLogTau << "\t" << treeString.str() << endl;
    			stringstream parameters;
    			parameters << std::setprecision(std::numeric_limits<double>::digits10 + 1 ) << currTreeLogScore << "\t" << currSeqErrRate << "\t" << currDropOutRate << "\t" << currLogTau;
    			sampleOutput << parameters.str() << "\t" << treeString.str() << endl;
    		}
        }
        delete [] currTreeParentVec;
        delete [] bestTreeParentVec;
        //free_doubleMatrix(currLogScores);
        free_boolMatrix(currTreeAncMatrix);
	}                                              // last repetition of MCMC done

//	unsigned int noStepsAfterBurnin = noOfReps*(noOfLoops-burnIn);
//	cout.precision(17);
//	cout << "best log score for tree:\t" << bestTreeLogScore <<  "\n";
//	cout << "#optimal steps after burn-in:\t" << optStatesAfterBurnIn << "\n";
//	cout << "total #steps after burn-in:\t" << noStepsAfterBurnin << "\n";
//	cout << "%optimal steps after burn-in:\t" << (1.0*optStatesAfterBurnIn)/noStepsAfterBurnin << "\n";
//	if(moveProbs[0]!=0.0){
//		cout << "best value for beta:\t" << bestBeta << "\n";
//		cout << "best value for alpha:\t" << bestAlpha << "\n";
//		cout << "best log score for (T, theta):\t" << bestScore << "\n";
//	}
	sampleOutput.close();
	return "";
}




double scoreTree(int m, int n, int sampleCount, bool** ancMatrix, vector<int>& alleleCount, vector<int>& leafClusterId, vector<vector<int> >& mutReadCounts, vector<vector<int> >& totalReadCounts, double dropoutRate, double seqErr, double tau, vector<int>& bestPlacementPoint, vector<bool>& wbcStatus){

	vector<vector<int> > expVarAlleleCount = getExpVarAlleleCount(m, sampleCount, ancMatrix, leafClusterId);
	vector<vector<int> > expVarAlleleFreqs = getExpAlleleCountFreqs(m, sampleCount, ancMatrix, leafClusterId);
	vector<int> wbcBelowCount = getWBC_count_below(m, sampleCount, ancMatrix, wbcStatus);

	double logScore = 0.0;

	for(int mut=0; mut< n; mut++){               // compute score separately for each mutation
		//cout << "mut " << mut << " of " << n << endl;
		vector<double> logAttachmentScores;
		double bestLogAttachmentScore = -DBL_MAX;
		//bestPlacementPoint.push_back(-1);

		/* precompute scores per expected number of mutated alleles */
		vector<vector<double> > scorePrecomp;
		for(int sample=0; sample<sampleCount; sample++){
			vector<double> sampleScorePrecomp;
			//cout << "sample: " << sample << endl;
			for(int expFreq=0; expFreq<expVarAlleleFreqs.at(sample).size(); expFreq++){
				//cout << "exp freq (" << mut << "): " << expFreq << "   allele count at sample " << alleleCount.at(sample) << " (" << mutReadCounts.at(mut).at(sample) << "," << totalReadCounts.at(mut).at(sample) << ")" << endl;
				double newScore = computeLogSampleScore(expFreq, alleleCount.at(sample), mutReadCounts[mut][sample], totalReadCounts[mut][sample], dropoutRate, seqErr, tau);
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


//			cout << "mut placement list size: " << bestPlacementPoint.size() <<  " at mut " << mut << endl;
//			cout << "previous best: " << bestLogAttachmentScore << endl;
//			cout << "          new: " << logAttachmentScore << endl;
			if(bestLogAttachmentScore < logAttachmentScore){      // update best placement point of mut in tree
				//cout << "better placement" << endl;
				bestPlacementPoint.at(mut) = att;
			}
			else if(bestLogAttachmentScore == logAttachmentScore){
//				cout << "equally good placement" << endl;
				//cout << "current best placement: " << bestPlacementPoint.at(mut) << endl;
				//cout << "   candidate placement: " << att << endl;
				if(att==(2*m)-2){
					bestPlacementPoint.at(mut) = att;
				}
				else if(ancMatrix[bestPlacementPoint.at(mut)][att]==1){
					bestPlacementPoint.at(mut) = att;
				}
			}
			//cout << "done" << endl;

			bestLogAttachmentScore = max(bestLogAttachmentScore, logAttachmentScore);
		}

		double sumScore = 0.0;
		for(int i=0; i<logAttachmentScores.size(); i++){              // sum over all allele states, exp is necessary as scores are actually log scores
			sumScore += exp(logAttachmentScores.at(i)-bestLogAttachmentScore);   // subtraction of best score to calculate with score differences (smaller values)
		}

		double logMutScore = log(sumScore)+bestLogAttachmentScore;              // transform back to log scores and change from score differences to actual scores
		logScore += logMutScore;
	}

	return logScore;
}


vector<vector<double> > getPlacementScoreTable(int m, int n, int sampleCount, bool** ancMatrix, vector<int>& alleleCount, vector<int>& leafClusterId, vector<vector<int> >& mutReadCounts, vector<vector<int> >& totalReadCounts, double dropoutRate, double seqErr, double tau, vector<int>& bestPlacementPoint, vector<bool>& wbcStatus){

	vector<vector<int> > expVarAlleleCount = getExpVarAlleleCount(m, sampleCount, ancMatrix, leafClusterId);
	vector<vector<int> > expVarAlleleFreqs = getExpAlleleCountFreqs(m, sampleCount, ancMatrix, leafClusterId);
	vector<int> wbcBelowCount = getWBC_count_below(m, sampleCount, ancMatrix, wbcStatus);

	double logScore = 0.0;

	vector<vector<double> > table;

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
				double newScore = computeLogSampleScore(expFreq, alleleCount.at(sample), mutReadCounts[mut][sample], totalReadCounts[mut][sample], dropoutRate, seqErr, tau);
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


//			cout << "mut placement list size: " << bestPlacementPoint.size() <<  " at mut " << mut << endl;
//			cout << "previous best: " << bestLogAttachmentScore << endl;
//			cout << "          new: " << logAttachmentScore << endl;
			if(bestLogAttachmentScore < logAttachmentScore){      // update best placement point of mut in tree
				//cout << "better placement" << endl;
				bestPlacementPoint.at(mut) = att;
			}
			else if(bestLogAttachmentScore == logAttachmentScore){
//				cout << "equally good placement" << endl;
//				cout << "current best placement: " << bestPlacementPoint.at(mut) << endl;
//				cout << "   candidate placement: " << att << endl;
				if(ancMatrix[bestPlacementPoint.at(mut)][att]==1){
					bestPlacementPoint.at(mut) = att;
				}
			}
			//cout << "done" << endl;

			bestLogAttachmentScore = max(bestLogAttachmentScore, logAttachmentScore);
		}

		table.push_back(logAttachmentScores);
	}

	return table;
}

/* computes the logScore for one sample for a given mutation and placement */
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

/* precompute all binomial coefficients used in the score computation */
int** getBinomCoeffTable(int maxAlleleCount){

	int** table = init_intMatrix(maxAlleleCount+1, maxAlleleCount+1, 0);

	for (int i=0; i<maxAlleleCount+1; i++){
		for (int j=0; j<=i; j++){

			if (j == 0 || j == i){   // base cases
				table[i][j] = 1;
			}
			else{                                         // recursive cases
				table[i][j] = table[i-1][j-1] + table[i-1][j];
	        }
	    }
	}
	return table;
}

double logBinom(int n, int k){
	return lgamma(n+1.0)-(lgamma(k+1.0)+lgamma(n+1.0-k));
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

double* precomputeLogGamma(int maxGamma){

	double* table = init_doubleArray(maxGamma+1, 0.0);
	for(int i=0; i<=maxGamma; i++){
		double value = lgamma(i);
		table[i] = value;
	}
	return table;
}

//double getLogGamma(double )
