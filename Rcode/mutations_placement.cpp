#include <Rcpp.h>

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <float.h>
#include <math.h>
#include <deque>
#include <map>
#include <iomanip>
#include <sstream>
#include <limits>



double chi_wbc = 1.0; // penalty for placing a mutation above a wbc

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



/* precompute for each pair (sample, mutation placement) the number of expected mutated alleles in the sample */

std::vector<std::vector<int> > getExpVarAlleleCount(int m, int sampleCount, std::vector<std::vector<bool> > ancMatrix, std::vector<int> leafClusterId){
  
  std::vector<std::vector<int> > expVarReadCount;
  for(int sample=0; sample<sampleCount; sample++){
    expVarReadCount.push_back(std::vector<int>());
    for(int node=0; node<(2*m)-1; node++){
      expVarReadCount.at(sample).push_back(0);
    }
  }
 
  for(int leaf=0; leaf<m; leaf++){
    int clusterId = leafClusterId[leaf];
    
    for(int node=0; node<(2*m)-2; node++){
      if(ancMatrix[node][leaf]==1){
        expVarReadCount[clusterId][node] += 1;   // count how many leaves belonging to the sample have
      }                                          // the node where the mutation is placed as ancestor
    }
    
  
    expVarReadCount[clusterId][(2*m)-2] += 1; 
  }
   
  return expVarReadCount;
}


std::vector<std::vector<int> > getExpAlleleCountFreqs(int m, int sampleCount, std::vector<std::vector<bool> > ancMatrix, std::vector<int> leafClusterId){
  
  std::vector<std::vector<int> > expVarReadCount = getExpVarAlleleCount(m, sampleCount, ancMatrix, leafClusterId);
  std::vector<std::vector<int> > expVarReadCountFreqs;
  
  
  for(int sample=0; sample<sampleCount; sample++){
    
    std::vector<int> freqCounter;
    int maxExpVarReadCounter = 0;
    for(int i=0; i<expVarReadCount.at(sample).size(); i++){
      maxExpVarReadCounter = std::max(maxExpVarReadCounter, expVarReadCount.at(sample).at(i));
    }
    for(int i=0; i<=maxExpVarReadCounter; i++){
      freqCounter.push_back(0);                    // set frequency of each expected allele frequency to zero
    }
    
    for(int node=0; node<(2*m)-1; node++){                        // increase frequency counter of expected frequency
      freqCounter.at(expVarReadCount.at(sample).at(node))+=1;  // found at this node by one
    }
    expVarReadCountFreqs.push_back(freqCounter);
  }
  return expVarReadCountFreqs;
}


/** this computes for each node in the current tree how many wbc's are below it (including node itself)**/
std::vector<int> getWBC_count_below(int m, int sampleCount, std::vector<std::vector<bool> > ancMatrix, std::vector<bool> wbcStatus){
  std::vector<int> wbcBelowCount;
  
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


/* computes the logScore for one sample for a given mutation and placement */
double computeLogSampleScore(int expCount, int alleleCount, int obsMutCount, int obsTotalCount, double dropoutRate, double seqErr, double tau){
  
  int alpha = expCount;
  int beta = alleleCount - alpha;
  int k = obsMutCount;
  int r = obsTotalCount;
  double delta = dropoutRate;
  
  //std::cout << "alleleCount " << alleleCount << std::endl;
  //std::cout <<"alpha " << alpha << std::endl;
  
  std::vector<double> logScores;
  double logSumScore;
  double maxLogScore = -DBL_MAX;
  
  
  if(r==0){
    return 0;
  }
  // both alleles still present: use beta-binomial for modeling drop-out
  for(int p=0; p<alpha; p++){
    //std::cout << "p " << p << std::endl;
    for(int q=0; q<beta; q++){
      //std::cout << "q " << q << std::endl;
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
      maxLogScore = std::max(maxLogScore, logScore); 
    }
  }

  // only mutated alleles amplified: use binomial for modeling sequencing errors
  for(int p=0; p<alpha; p++){
    double logScore = logBinomCoeff(alpha, p) + (beta+p)*log(delta) + (alpha-p)*log(1-delta) + logBinomCoeff(r, k) + (r-k)*log(seqErr) + k* log(1-seqErr);
    logScores.push_back(logScore);
    //cout << "       (ii)" << logScore << endl;
    maxLogScore = std::max(maxLogScore, logScore);
  }
  //only normal allele amplified: use binomial for modeling sequencing errors
  for(int q=0; q<beta; q++){
    double logScore = logBinomCoeff(beta, q) + (alpha+q)*log(delta) + (beta-q)*log(1-delta) + logBinomCoeff(r, k) + k*log(seqErr) + (r-k)* log(1-seqErr);
    logScores.push_back(logScore);
    
    
    //		cout << "      (iii)" << logScore << endl;
    maxLogScore = std::max(maxLogScore, logScore);
  }
  // these scores need to be summed -> go from log space to normal space
  double sumScore = 0.0;
  for(int i=0; i<logScores.size(); i++){              // sum over all allele states, exp is necessary as scores are actually log scores
    sumScore += exp(logScores.at(i)-maxLogScore);   // subtraction of best score to calculate with score differences (smaller values)
  }
  logSumScore = log(sumScore)+maxLogScore;      // transform back to log scores and change from score differences to actual scores
  //getchar();
  if(k>r){
    std::cout << "k > r: " << k << " > " << r<< std::endl;
    getchar();
  }
  
  return logSumScore;
}





// [[Rcpp::export]]
std::vector<int> getMutationPlacement(int m, int n, int sampleCount,
                                      std::vector<std::vector<bool> > ancMatrix,
                                      std::vector<int> alleleCount,
                                      std::vector<int> leafClusterId,
                                      std::vector<std::vector<int> > mutReadCounts,
                                      std::vector<std::vector<int> > totalReadCounts,
                                      double dropoutRate, 
                                      double seqErr,
                                      double tau,
                                      std::vector<bool> wbcStatus){

  std::vector<int> bestPlacementPoint;
  for(int mut=0; mut< n; mut++){               // compute score separately for each mutation
    bestPlacementPoint.push_back(-1);
  }
  
	std::vector<std::vector<int> > expVarAlleleCount = getExpVarAlleleCount(m, sampleCount, ancMatrix, leafClusterId);
  
	std::vector<std::vector<int> > expVarAlleleFreqs = getExpAlleleCountFreqs(m, sampleCount, ancMatrix, leafClusterId);
	std::vector<int> wbcBelowCount = getWBC_count_below(m, sampleCount, ancMatrix, wbcStatus);
  
	//double logScore = 0.0;
  
	for(int mut=0; mut<n; mut++){               // compute score separately for each mutation
		//cout << "mut " << mut << " of " << n << endl;
		std::vector<double> logAttachmentScores;
		double bestLogAttachmentScore = -DBL_MAX;
		//bestPlacementPoint.push_back(-1);

		// precompute scores per expected number of mutated alleles
		std::vector<std::vector<double> > scorePrecomp;
		
		for(int sample=0; sample<sampleCount; sample++){
			std::vector<double> sampleScorePrecomp;
			//cout << "sample: " << sample << endl;
		 
			for(int expFreq=0; expFreq<expVarAlleleFreqs.at(sample).size(); expFreq++){
				//std::cout << "exp freq (" << mut << "): " << expFreq << "   allele count at sample " << alleleCount.at(sample) << " (" << mutReadCounts.at(mut).at(sample) << "," << totalReadCounts.at(mut).at(sample) << ")" << std::endl;
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

				if(att==(2*m)-2){
					bestPlacementPoint.at(mut) = att;
				}
				else if(ancMatrix[bestPlacementPoint.at(mut)][att]==1){
					bestPlacementPoint.at(mut) = att;
				}
			}

			bestLogAttachmentScore = std::max(bestLogAttachmentScore, logAttachmentScore);
		}

		//double sumScore = 0.0;
		//for(int i=0; i<logAttachmentScores.size(); i++){              // sum over all allele states, exp is necessary as scores are actually log scores
		//	sumScore += exp(logAttachmentScores.at(i)-bestLogAttachmentScore);   // subtraction of best score to calculate with score differences (smaller values)
		//}

		//double logMutScore = log(sumScore)+bestLogAttachmentScore;              // transform back to log scores and change from score differences to actual scores
		//logScore += logMutScore;
	}
  
	return bestPlacementPoint;
}

// [[Rcpp::export]]
std::vector<std::vector<bool>> parentVector2ancMatrix(std::vector<int> parent, int n){
  
  std::vector<std::vector<bool>> ancMatrix(n, std::vector<bool>(n, false));
  
  int root = n;
  
  for(int i = 0; i < n; ++i){
    int anc = i;
    while(anc < root){                              
      if(parent[anc] < n){
        ancMatrix[parent[anc]][i] = true;
      }
      anc = parent[anc];
    }
  }
  
  for(int i = 0; i < n; ++i){
    ancMatrix[i][i] = true;
  }
  
  return ancMatrix;
}
