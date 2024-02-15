#include <Rcpp.h>

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <float.h>
#include <math.h>
#include <deque>
#include <map>
#include <iomanip>
#include <sstream>
#include <limits>
#include <random>




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



std::vector< double > normalizeDistribution(std::vector< double >& probabilityWeights) {
  double sum = 0;
  int lengthOfVector = probabilityWeights.size();
  for (int i = 0; i < lengthOfVector; i++) {
    sum += probabilityWeights[i];
  }
  std::vector< double > normalizedProbabilityWeights;
  
  for (int i = 0; i < lengthOfVector; i++) {
    normalizedProbabilityWeights.push_back(exp(log(probabilityWeights[i])-log(sum)));
  }
  return normalizedProbabilityWeights;
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


std::vector < int > drawWithReplacement(std::vector <double> probabilityWeights, int sampleSize) {
  std::vector <double> probabilityDistribution = normalizeDistribution(probabilityWeights);
  std::discrete_distribution<int> distribution(probabilityDistribution.begin(), probabilityDistribution.end());
  std::random_device rd;
  
  std::vector < int > sampledNumbers;
  for (int i = 0; i < sampleSize; i++) {
    std::mt19937 generator(rd());
    //std::mt19937 generator(seed);
    sampledNumbers.push_back(distribution(generator));  
  }
  return sampledNumbers;
}




// [[Rcpp::export]]
std::vector< std::vector<int> > transposeMatrix(std::vector< std::vector<int> > matrix, int nrows, int ncols) {
  std::vector<std::vector<int>> transposedMatrix(ncols, std::vector<int>(nrows));
  for(int i = 0; i < nrows; i++){
    for(int j = 0; j< ncols; j++){
      transposedMatrix[j][i] = matrix[i][j];
    }
  }
  return transposedMatrix;
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



/** computes the logScore for one sample for a given mutation and placement
 *
 *Since the mutation placements are assumed to be independent, the placement of
 *each mutation is optimised individually.
 *
 * @param expCount The estimated number of mutated alleles in the cell cluster (alpha)
 * @param alleleCount The number alleles in the cell cluster. Equals number of 
 * cells * 2 (alpha + beta)
 * @param obsMutCount The observed mutated allele count (after multiple displacement
 * amplification)
 * @param obsTotalCount The observed total allele count (after MDA)
 * @param dropoutRate (the inferred dropout rate)
 * @param seqErr (the inferred sequencing error rate (which is modelled as the false positive rate))
 * @param tau 
 * 
 * @output The logarithm of the probability of k mutated reads under the parameters given as input
 * and the read count model 
 * 
 */
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




/** Computes the best mutation placement on the tree
 *
 *Since the mutation placements are assumed to be independent, the placement of
 *each mutation is optimised individually.
 *
 * @param m The total number of samples
 * @param n The total number of mutations
 * @param sampleCount the total number of samples 
 * @param ancMatrix tree in ancestor-matrix format
 * @param alleleCount the number of alleles in each of the samples
 * @param leafClusterID indicates, which cluster each of the leaves belongs to.
 *  Format: value k in position i means that the i'th leaf of the tree is a cell
 *  belonging to CTC-cluster k.
 * @param mutReadCounts a matrix of mutated read counts per sample and mutation
 * @param totalReadCounts a matrix of total read counts per sample and mutation
 * @param dropoutRate The estimated allelic dropout rate
 * @param seqErr The estimated sequencing error rate (as in false-positive rate)
 * @param tau
 * @param wbcStatus 1 at position i indicates that leaf i is a white blood cell,
 *  otherwise it is a tumour cell
 * 
 * @output An integer-valued vector which for each mutations indicates the node it
 * is optimally placed on.
 * 
 */
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
		
		/*  Given a fixed mutation, for all possible alphas and samples (determines alpha + beta)
		 compute the probability of the observes read count */
		
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
		 
    
    /* Given a mutation, for each attachment point compute the joint probability of observing the read counts
		and having a particular attachment, given tree and other model parameters.
		Report the optimal attachment point of the mutation.*/
		//cout << "end 2" << endl;
		for(int att=0; att< (2*m)-1; att++){           
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
			if(bestLogAttachmentScore < (logAttachmentScore - wbc_penalty) ){      // update best placement point of mut in tree
				//cout << "better placement" << endl;
				bestPlacementPoint.at(mut) = att;
			}
			else if(bestLogAttachmentScore == (logAttachmentScore - wbc_penalty) ){

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


/** Computes the posterior distribution of mutation placements. Since mutation
 * placements of different mutations are assumed to be independent, it is enough
 * to give a finite distribution for each of the mutations.
 *
 *Since the mutation placements are assumed to be independent, the placement of
 *each mutation is optimised individually.
 *
 * @param m The total number of cells
 * @param n The total number of mutations
 * @param sampleCount the total number of samples, i.e. CTC clusters 
 * @param ancMatrix tree in ancestor-matrix format
 * @param alleleCount the number of alleles in each of the samples
 * @param leafClusterID indicates, which cluster each of the leaves belongs to.
 *  Format: value k in position i means that the i'th leaf of the tree is a cell
 *  belonging to CTC-cluster k.
 * @param mutReadCounts a matrix of mutated read counts per sample and mutation
 * @param totalReadCounts a matrix of total read counts per sample and mutation
 * @param dropoutRate The estimated allelic dropout rate
 * @param seqErr The estimated sequencing error rate (as in false-positive rate)
 * @param tau
 * @param wbcStatus 1 at position i indicates that leaf i is a white blood cell,
 *  otherwise it is a tumour cell
 * 
 * @output An integer-valued matrix which for each mutations (row) gives an unnormalized 
 * probability weight to be placed on a certain node (column). 
 * 
 */
// [[Rcpp::export]]
std::vector< std::vector<double> > computeMutationDistribution(int m, int n, int sampleCount,
                                                               std::vector<std::vector<bool> > ancMatrix,
                                                               std::vector<int> alleleCount,
                                                               std::vector<int> leafClusterId,
                                                               std::vector<std::vector<int> > mutReadCounts,
                                                               std::vector<std::vector<int> > totalReadCounts,
                                                               double dropoutRate, 
                                                               double seqErr,
                                                               double tau,
                                                               std::vector<bool> wbcStatus){

  //std::vector<int> bestPlacementPoint;
  //for(int mut=0; mut< n; mut++){               // compute score separately for each mutation
  //  bestPlacementPoint.push_back(-1);
  //}
  
  std::vector<std::vector<int> > expVarAlleleCount = getExpVarAlleleCount(m, sampleCount, ancMatrix, leafClusterId);
  
  std::vector<std::vector<int> > expVarAlleleFreqs = getExpAlleleCountFreqs(m, sampleCount, ancMatrix, leafClusterId);
  std::vector<int> wbcBelowCount = getWBC_count_below(m, sampleCount, ancMatrix, wbcStatus);
  
  //double logScore = 0.0;
  

  
  std::vector< std::vector<double> > logAttachmentScores;
  
  for(int mut=0; mut<n; mut++){               // compute score separately for each mutation
    //cout << "mut " << mut << " of " << n << endl;
    
    logAttachmentScores.push_back(std::vector<double>());
    
    //double bestLogAttachmentScore = -DBL_MAX;
    //bestPlacementPoint.push_back(-1);
    
    // precompute scores per expected number of mutated alleles
    std::vector<std::vector<double> > scorePrecomp;
    
    /*  Given a fixed mutation, for all possible alphas and samples (determines alpha + beta)
     compute the probability of the observes read count */

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
    

    /* Given a mutation,compute the joint probability of observing the read counts
     and having a particular attachment, given tree and other model parameters.
     Report the optimal attachment point of the mutation.*/
    //cout << "end 2" << endl;
    for(int att=0; att< (2*m)-1; att++){           
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
      logAttachmentScores.back().push_back(logAttachmentScore - wbc_penalty);
      
      
      
      
      
      //			cout << "mut placement list size: " << bestPlacementPoint.size() <<  " at mut " << mut << endl;
      //			cout << "previous best: " << bestLogAttachmentScore << endl;
      //			cout << "          new: " << logAttachmentScore << endl;
      //if(bestLogAttachmentScore < logAttachmentScore){      // update best placement point of mut in tree
      //cout << "better placement" << endl;
      //    bestPlacementPoint.at(mut) = att;
      //}
      //else if(bestLogAttachmentScore == logAttachmentScore){
      
      //  if(att==(2*m)-2){
      //  bestPlacementPoint.at(mut) = att;
      //}
      //else if(ancMatrix[bestPlacementPoint.at(mut)][att]==1){
      //  bestPlacementPoint.at(mut) = att;
      //}
      //}
      
      //  bestLogAttachmentScore = std::max(bestLogAttachmentScore, logAttachmentScore);
    }
    
    //double sumScore = 0.0;
    //for(int i=0; i<logAttachmentScores.size(); i++){              // sum over all allele states, exp is necessary as scores are actually log scores
    //	sumScore += exp(logAttachmentScores.at(i)-bestLogAttachmentScore);   // subtraction of best score to calculate with score differences (smaller values)
    //}
    
    //double logMutScore = log(sumScore)+bestLogAttachmentScore;              // transform back to log scores and change from score differences to actual scores
    //logScore += logMutScore;
  }
  
  return logAttachmentScores;
}









/** Computes the posterior distribution of mutation placements. Since mutation
 * placements of different mutations are assumed to be independent, it is enough
 * to give a finite distribution for each of the mutations.
 *
 *Since the mutation placements are assumed to be independent, the placement of
 *each mutation is optimised individually.
 *
 * @param nCells The total number of cells
 * @param nMutations The total number of mutations
 * @param sampleCount the total number of samples (i.e. clusters)
 * @param ancMatrix tree in ancestor-matrix format
 * @param alleleCount the number of alleles in each of the samples
 * @param leafClusterID indicates, which cluster each of the leaves belongs to.
 *  Format: value k in position i means that the i'th leaf of the tree is a cell
 *  belonging to CTC-cluster k.
 * @param mutReadCounts a matrix of mutated read counts per sample and mutation
 * @param totalReadCounts a matrix of total read counts per sample and mutation
 * @param dropoutRate The estimated allelic dropout rate
 * @param seqErr The estimated sequencing error rate (as in false-positive rate)
 * @param tau
 * @param wbcStatus 1 at position i indicates that leaf i is a white blood cell,
 *  otherwise it is a tumour cell
 * 
 * @output An integer-valued matrix which for each mutations (row) gives an unnormalized 
 * probability weight to be placed on a certain node (column). 
 * 
 */
// [[Rcpp::export]]
std::vector< std::vector<int> > sampleMutationsPlacement(int nSamplingEvents, int nMutations,
                                                         std::vector< std::vector<double> > logMutationPlacementProbabilities){
  

  
  std::vector< std::vector<int> > mutationPlacements;

  for(int m = 0; m < nMutations; m++){
    std::vector<double> mutationPlacementProbabilities;
    for(int it = 0; it < logMutationPlacementProbabilities[m].size(); it++){
      mutationPlacementProbabilities.push_back(exp(logMutationPlacementProbabilities[m][it]));
    }
    mutationPlacements.push_back(drawWithReplacement(mutationPlacementProbabilities, nSamplingEvents));
  }
  
  std::vector< std::vector<int> > mutationPlacementsTransposed = transposeMatrix(mutationPlacements, nMutations, nSamplingEvents);
  
  return mutationPlacementsTransposed;
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

// // [[Rcpp::export]]
//std::vector<std::vector<bool>> ComputeDistanceLeavesPosterior() {
//  std::vector<double> distance;
//  std::vector<int> clusterIdentityofdistance;
//
//  for (int c = 0; c < nClusters; ++c) {
//    std::vector<int> cellsInCluster;
//  
//    for (int i = 0; i < ClusterID.size(); ++i) {
//      if (ClusterID[i] == c) {
//        cellsInCluster.push_back(i);
//      }
//    }
  
//    for (int i = 0; i < cellsInCluster.size(); ++i) {
//      for (int j =0; j< i; j++) {
//        distance.push_back(produce_Distance_Posterior(cellsInCluster[i], cellsInCluster[j], treePosterior, treeName));
//        clusterIdentityofdistance.push_back(c);
//        }
//      }
//    }
//  }
//  return ....;
//}



//Input: a tree in parent vector format, meaning that the i'th entry of the vector
// is te parent node of the entry i. Nodes are counted from zero and the root is 
// length(Tree)
//Output: A list with three entries:
//  - the first entry is a vector of nodes tracing back leaf 1 to the root
//  - the second entry is a vector of nodes tracing back leaf2 to the descendant
//    of the MRCA in the lineage
//  - the thirs entry is the MRCA


/** Computes the most recent common ancestor and the shortest path between two
 * leaves of the tree
 *
 *
 * @param treeParentVectorFormat
 * @param leaf1 A number from 0 to the number of leaves-1 denoting the leaf index
 * @param leaf2 like leaf1
 * 
 * @output A list with 4 entries:
 *  - The first is the lineage of leaf 1 encoded as
 *    a sequence of nodes starting with leaf 1 and ending with the root
 *  - The second is the lineage of leaf2 starting with leaf2 and ending at the 
 *    daughter of the most recent common ancestor (MRCA)
 *  - The third one is the MRCA
 *  - The fourth one is the shortest path between leaf1 and leaf 2. If leaf1==leaf2
 *    this will output leaf1.
 * 
 */
// [[Rcpp::export]]
std::vector<std::vector<int>> findMostRecentCommonAncestor(const std::vector<int>& treeParentVectorFormat, int leaf1, int leaf2) {
  std::vector<int> lineage1{leaf1};
  
  while (true) {
    lineage1.push_back(treeParentVectorFormat[lineage1.back()]);
    // std::cout << lineage1[lineage1.back()] << std::endl;
    if (lineage1.back() == treeParentVectorFormat.size()) break;
  }
  
  std::vector<int> lineage2;
  int nextParent = leaf2;
  
  while (std::find(lineage1.begin(), lineage1.end(), nextParent) == lineage1.end()) {
    lineage2.push_back(nextParent);
    nextParent = treeParentVectorFormat[nextParent];
  }
  
  int MRCA = nextParent;
  
  
  // Calculate path between leaves

  auto it = std::find(lineage1.begin(), lineage1.end(), MRCA);
  
  std::vector<int> pathBetweenLeaves(lineage1.begin(), it + 1);
  

  std::vector<int> reversedLineage2(lineage2);
  std::reverse(reversedLineage2.begin(), reversedLineage2.end());
  pathBetweenLeaves.insert(pathBetweenLeaves.end(),reversedLineage2.begin(), reversedLineage2.end());
  
  
  
  
  return {lineage1, lineage2, {MRCA},{pathBetweenLeaves}};
}

/** Takes a tree and computes the posterior probability
 * distribution of mutation placements in the tree.
 * WORK IN PROGRESS
 *
 * @param tree A tree as ancestor matrix
 * @param nMutations The total number of mutations
 * @param nCells The total number of Cells
 * 
 * @output A matrix giving the probability of mutation i being placed on branch j.
 * 
 */

//std::vector<std::vector<double> > ComputeMutationPlacementPosterior(
//    std::vector<std::vector<int> > tree, int nMutations, int nCells, 
//    std::vector<std::vector<int> > expCount, std::vector<int> alleleCount,
//   std::vector<std::vector<int> > obsMutCount,
//   std::vector<std::vector<int> > obsTotalCount, double dropoutRate,
//   double seqErr) {
// std::vector<double> LogAttachmentScores;
// double normalisingConstant = ;
// for(int i = 0; i<n, ++i){
//   LogAttachmentScores.push_back((std::vector<int>())
//   for(int j = 0;j< number of nodes in the tree(which is a function of nCells); j++){
//     for(int l = 0; l<nCells, ++l){
//       double logSampleScore = computeLogSampleScore(expCount, alleleCount, obsMutCount, obsTotalCount, dropoutRate, seqErr, 1));
//       double attachmentPrior = ;
//       
//     }
//     LogAttachmentScores.back().push_back(LogAttachmentScore)
//   }
// }
//} 



/** For a given tree and a given pair of leaves, the evolutionary distance of
 * the leaves is computed by counting the number of mutations that are on the
 * shortest path from one leaf to the other.
 *
 *
 * @param treeData The DataFrame containing the information from the posterior Sampling
 * @param leaf1 A number from 0 to the number of leaves-1 denoting the leaf index
 * @param leaf2 like leaf1
 * @param nCells the total number of cells in the experiment
 * @param nMutations the total number of mutations in the experiment
 * @param nClusters the total number of CTC-clutsers (inclduing single CTCs)
 * in the experiment
 * @param allelCount A vector that counts the number of alleles for each of the 
 * clusters
 * @param ClusterID A vector of the length corresponding to the number of cells
 * with entries denoting an index for the ClusterID. Corresponds to a number from
 * 0 to number of Clusters -1.
 * @param mutatedReadCounts A matrix with mutations as rows and clusters as columns
 * @param totalReadCounts As mutatedReadCounts
 * @param wbcStatus A binary vector with as many entries as there are cells. 1 stands
 * for being a white blood cell.
 * 
 * @output the total number of mutations separating the genotypes of the two leaves.
 * 
 */
// [[Rcpp::export]]
int computePairwiseDistanceOfLeaves(Rcpp::DataFrame treeData, int leaf1, int leaf2, int nCells,
                                        int nMutations, int nClusters,
                                        std::vector<int> alleleCount, std::vector<int> ClusterID,
                                        std::vector<std::vector<int> > mutatedReadCounts, std::vector<std::vector<int> > totalReadCounts,
                                        std::vector<bool> wbcStatus) {
  
  std::string tree = treeData["Tree"];
  
  //Now need to split the string into single numbers and turn them into integers
  // Using a stringstream to split the string
  std::istringstream iss(tree);
  std::vector<int> treeParentVectorFormat;
  // Iterate over each substring and convert to integer
  int num;
  while (iss >> num) {
    treeParentVectorFormat.push_back(num);
  }
  
  double dropoutRate = Rcpp::as<double>(treeData["DropoutRate"]);
  double seqErrRate = Rcpp::as<double>(treeData["SequencingErrorRate"]);
  
  // Preprocess tree
  std::vector<std::vector<bool> > ancestorMatrix = parentVector2ancMatrix(treeParentVectorFormat, treeParentVectorFormat.size());
  
  // Find best Mutation placement
  std::vector<int> bestMutationPlacement = getMutationPlacement(nCells, nMutations, nClusters,
                                                             ancestorMatrix, alleleCount, ClusterID,
                                                             mutatedReadCounts, totalReadCounts,
                                                             dropoutRate, seqErrRate, 1, wbcStatus);
  
  
  //std::vector<std::vector <double> > mutationDistribution = computeMutationDistribution(nCells, nMutations, nClusters,
  //                                                              ancestorMatrix, alleleCount, ClusterID,
  //                                                              mutatedReadCounts, totalReadCounts,
  //                                                              dropoutRate, seqErrRate, 1, wbcStatus);
  
  
  // Finding most recent common ancestor
  std::vector<std::vector<int> > pairwiseGenealogy = findMostRecentCommonAncestor(treeParentVectorFormat, leaf1, leaf2);
  
  std::vector<int> firstLeafToMRCA = pairwiseGenealogy[0];
  std::vector<int> secondLeafToMRCA = pairwiseGenealogy[1];
  
  // Calculate path between leaves
  std::vector<int> pathBetweenLeaves = firstLeafToMRCA;
  std::reverse(secondLeafToMRCA.begin(), secondLeafToMRCA.end());
  pathBetweenLeaves.insert(pathBetweenLeaves.end(),secondLeafToMRCA.begin(), secondLeafToMRCA.end());
  
  
  // Count mutations in the shortest path between leaves excluding MRCA
  int positionOfMRCA = 0;
  for (int i = 0; i < pathBetweenLeaves.size(); ++i) {
    if (pathBetweenLeaves[i] == firstLeafToMRCA[positionOfMRCA]) {
      positionOfMRCA = i;
      break;
    }
  }
  
  int result = 0;
  for (int i = 0; i < bestMutationPlacement.size(); ++i) {
    if (std::find(pathBetweenLeaves.begin(), pathBetweenLeaves.end(), bestMutationPlacement[i]) !=
        pathBetweenLeaves.end() && bestMutationPlacement[i] != firstLeafToMRCA[positionOfMRCA]) {
      result++;
    }
  }
  
  return result;
}









/** For a given tree, pair of leaves and mutation placement on the tree,
 * the evolutionary distance of the leaves is computed by counting the number of
 *  mutations that are on the shortest path from one leaf to the other.
 *
 *
 * @param treeParentVectorFormat
 * @param leaf1 A number from 0 to the number of leaves-1 denoting the leaf index
 * @param leaf2 like leaf1
 * @param mutationPlacement a vector of tree node indices indicating the placement 
 * of the i'th mutation in its i'th entry
 * @param pairwiseGenealogy: The output of findMostRecentCommonAncestor
 * 
 * @output the total number of mutations separating the genotypes of the two leaves.
 * 
 */
// [[Rcpp::export]]
std::vector< int > computePairwiseDistanceOfLeaves2(std::vector<int> treeParentVectorFormat, int leaf1, int leaf2,
                                     std::vector<int> mutationPlacement, std::vector< std::vector<int> > pairwiseGenealogy) {
  
  //std::string tree = treeData["Tree"];
  
  //Now need to split the string into single numbers and turn them into integers
  // Using a stringstream to split the string
  //std::istringstream iss(tree);
  //std::vector<int> treeParentVectorFormat;
  // Iterate over each substring and convert to integer
  //int num;
  //while (iss >> num) {
  //  treeParentVectorFormat.push_back(num);
  //}
  
  // Find best Mutation placement
  //std::vector<int> bestMutationPlacement = getMutationPlacement(nCells, nMutations, nClusters,
  //                                                              ancestorMatrix, alleleCount, ClusterID,
  //                                                              mutatedReadCounts, totalReadCounts,
  //                                                              dropoutRate, seqErrRate, 1, wbcStatus);
  
  

  
  std::vector<int> firstLeafToMRCA = pairwiseGenealogy[0];
  std::vector<int> secondLeafToMRCA = pairwiseGenealogy[1];
  //const std::vector<int>& nestedVector = pairwiseGenealogy[2];
  int MRCA = pairwiseGenealogy[2][0];
  std::vector<int> pathBetweenLeaves =pairwiseGenealogy[3];
  
  
  int distance = 0;
  int mutationsOnFirstBranch = 0;
  int mutationsOnSecondBranch = 0;
  int split = 0;
  for (int i = 0; i < mutationPlacement.size(); ++i) {
    if ((std::find(pathBetweenLeaves.begin(), pathBetweenLeaves.end(), mutationPlacement[i]) != pathBetweenLeaves.end()) && (mutationPlacement[i] != MRCA)) {
        //pathBetweenLeaves.end() && mutationPlacement[i] != firstLeafToMRCA[positionOfMRCA]) {
      distance++;
      if((std::find(firstLeafToMRCA.begin(), firstLeafToMRCA.end(), mutationPlacement[i]) != firstLeafToMRCA.end())){
        mutationsOnFirstBranch++;
      }
      else{
        mutationsOnSecondBranch++;
      }
    }
  }
  if(mutationsOnFirstBranch != 0 & mutationsOnSecondBranch != 0){
    split = 1;
  }
  
  std::vector<int> output;
  
  output.push_back(distance);
  output.push_back(split);
  return output;
}






/** This takes a tree and the mutation placement probabilites of all mutations.
 * It computes for a pair of leaves (given by their pairwise genealogy) and their two
 * branches connecting them
 * the maximal probability of any mutation to land on either of the branches.
 * It output the minimum of the two branaches. Intuitiely, if there is no branching
 * evolution, then at least one of the branches should be populated with mutations
 * with very low probability.
 *
 *
 * @param pairwiseGenealogy as output by findMostRecentCommonAncestor
 * 
 * @param logMutationPlacementProbabilities As output by computeMutationDistribution
 * @param nMutations The total number of mutations
 * @param nCells The total number of cells.
 * 
 * @output The minimum (over the two branches) of the maximum (over all mutations)
 * of the probability of a mutation to land on one of the branches
 * 
 */
// [[Rcpp::export]]
double ComputePerMutationProbabilityOfPolyclonality(std::vector< std::vector<int> > pairwiseGenealogy,
                                                 std::vector< std::vector<double> > logMutationPlacementProbabilities,
                                                 int nMutations, int nCells){
  
  std::vector<int> lineage1 = pairwiseGenealogy[0];
  auto it = std::find(lineage1.begin(), lineage1.end(), pairwiseGenealogy[2][0]);
  lineage1.erase(it, lineage1.end());
  std::vector<int> lineage2 = pairwiseGenealogy[1];
  

  
  double maxProbabilityOfMutationInLineage1 = 0;
  double maxProbabilityOfMutationInLineage2 = 0;
  
  std::vector<double> mutationPlacementProbabilities;
  for(int mut = 0; mut < nMutations; mut++){
    double ProbabilityOfMutationInLineage1 = 0;
    double ProbabilityOfMutationInLineage2 = 0;
    std::vector<double> mutationPlacementProbabilityWeights;
    
    for(int it = 0; it < logMutationPlacementProbabilities[mut].size(); it++){
      mutationPlacementProbabilityWeights.push_back(exp(logMutationPlacementProbabilities[mut][it]));
    }
    mutationPlacementProbabilities = normalizeDistribution(mutationPlacementProbabilityWeights);
    for(int i = 0; i < 2*nCells-1; i++){
      if(std::find(lineage1.begin(), lineage1.end(), i) != lineage1.end()){
        ProbabilityOfMutationInLineage1 += mutationPlacementProbabilities[i];
      }
      if(std::find(lineage2.begin(), lineage2.end(), i) != lineage2.end()){

        ProbabilityOfMutationInLineage2 += mutationPlacementProbabilities[i];
      }
    }
    maxProbabilityOfMutationInLineage1 = std::max(maxProbabilityOfMutationInLineage1, ProbabilityOfMutationInLineage1);
    maxProbabilityOfMutationInLineage2 = std::max(maxProbabilityOfMutationInLineage2, ProbabilityOfMutationInLineage2);

  }

  //return mutationPlacementProbabilities;
  return std::min(maxProbabilityOfMutationInLineage1, maxProbabilityOfMutationInLineage2);
}





/** For a given tree and pair of leaves find mutation placements on the tree and
 * compute for each placement the evolutionary distance of the leaves by
 * counting the number of mutations that are on the shortest path from one leaf to the other.
 *
 *
 * @param treeData The DataFrame containing the information from the posterior Sampling
 * @param leaf1 A number from 0 to the number of leaves-1 denoting the leaf index
 * @param leaf2 like leaf1
 * @param nCells the total number of cells in the experiment
 * @param nMutations the total number of mutations in the experiment
 * @param nClusters the total number of CTC-clutsers (inclduing single CTCs)
 * in the experiment
 * @param allelCount A vector that counts the number of alleles for each of the 
 * clusters
 * @param ClusterID A vector of the length corresponding to the number of cells
 * with entries denoting an index for the ClusterID. Corresponds to a number from
 * 0 to number of Clusters -1.
 * @param mutatedReadCounts A matrix with mutations as rows and clusters as columns
 * @param totalReadCounts As mutatedReadCounts
 * @param wbcStatus A binary vector with as many entries as there are cells. 1 stands
 * for being a white blood cell.
 * 
 * @output A vector containing the total number of mutations separating the
 * genotypes of the two leaves for each of the mutation placements.
 * 
 */
// [[Rcpp::export]]
std::vector<std::vector<double> > computePairwiseDistanceOfLeavesGivenTree(Rcpp::DataFrame treeData, int leaf1, int leaf2, int nCells,
                                                           int nMutations, int nClusters,
                                                           std::vector<int> alleleCount, std::vector<int> ClusterID,
                                                           std::vector<std::vector<int> > mutatedReadCounts, std::vector<std::vector<int> > totalReadCounts,
                                                           std::vector<bool> wbcStatus, int nSamplingEvents){
  std::string tree = treeData["Tree"];
  
  //Now need to split the string into single numbers and turn them into integers
  // Using a stringstream to split the string
  std::istringstream iss(tree);
  std::vector<int> treeParentVectorFormat;
  // Iterate over each substring and convert to integer
  int num;
  while (iss >> num) {
    treeParentVectorFormat.push_back(num);
  }
  
  double dropoutRate = Rcpp::as<double>(treeData["DropoutRate"]);
  double seqErrRate = Rcpp::as<double>(treeData["SequencingErrorRate"]);
  
  // Preprocess tree
  std::vector<std::vector<bool> > ancestorMatrix = parentVector2ancMatrix(treeParentVectorFormat, treeParentVectorFormat.size());
  
  
  
  
  std::vector< std::vector<double> > logMutationPlacementProbabilities = computeMutationDistribution(nCells, nMutations, nClusters,
                                                                                                     ancestorMatrix,
                                                                                                     alleleCount,
                                                                                                     ClusterID,
                                                                                                     mutatedReadCounts,
                                                                                                     totalReadCounts,
                                                                                                     dropoutRate, 
                                                                                                     seqErrRate,
                                                                                                     1,
                                                                                                     wbcStatus);
  
  
  std::vector< std::vector<int> >mutationPlacements = sampleMutationsPlacement(nSamplingEvents,
                                                                               nMutations,
                                                                               logMutationPlacementProbabilities);
  
  

  //mutationPlacements.push_back(getMutationPlacement(nCells, nMutations, nClusters,
  //                                                  ancestorMatrix, alleleCount, ClusterID,
  //                                                  mutatedReadCounts, totalReadCounts,
  //                                                  dropoutRate, seqErrRate, 1, wbcStatus));
  
  
  // Finding most recent common ancestor
  std::vector< std::vector<int> > pairwiseGenealogy = findMostRecentCommonAncestor(treeParentVectorFormat, leaf1, leaf2);
  
  std::vector<int> distanceVector;
  std::vector<int> splittingStats;
  for(int it = 0; it < mutationPlacements.size(); it++){
    std::vector<int> output = computePairwiseDistanceOfLeaves2(treeParentVectorFormat, leaf1, leaf2,
                                     mutationPlacements[it], pairwiseGenealogy);
    
    distanceVector.push_back(output[0]);
    splittingStats.push_back(output[1]);
  }
  
  std::vector<double> PerMutationProbabilityOfPolyclonality = {ComputePerMutationProbabilityOfPolyclonality(pairwiseGenealogy,
                                               logMutationPlacementProbabilities,
                                               nMutations, nCells)};
  
  std::vector<double> distanceVectorDouble(distanceVector.begin(), distanceVector.end());
  std::vector<double> splittingStatsDouble(splittingStats.begin(), splittingStats.end());
  
  std::vector< std::vector<double> > result;
  result.push_back(distanceVectorDouble);
  result.push_back(splittingStatsDouble);
  result.push_back(PerMutationProbabilityOfPolyclonality);
  
  return result;
}
