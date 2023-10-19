/*
 * Beta_Distr.h
 *
 *  Created on: Feb 26, 2018
 *      Author: jahnka
 */

#ifndef BETA_DISTR_H_
#define BETA_DISTR_H_

using namespace std;

class Beta_Distr {

public:
    double mean, sd, alpha, beta;
    bool uniform;

public:
    Beta_Distr (double priorMean, double priorSD, bool uni){
    	mean = priorMean;
    	sd = priorSD;
    	uniform = uni;
    	if(uniform){
    		alpha = 1.0;                                    // setting alpha and beta of the beta distribution for alpha to 1 gives a uniform prior
    		beta = 1.0;
    	}
    	else{
    		alpha = ((1-mean)*mean*mean/(sd*sd)) - mean;    // <-10.13585344 turn the mean and sd into parameters of the beta distribution
    	  	beta = alpha*((1/mean)-1);                      //<-13.38666556
    	}
    };
    double logBetaPDF (double x) {return lgamma(alpha+beta)+(alpha-1)*log(x)+(beta-1)*log(1-x)-lgamma(alpha)-lgamma(beta);}
};



#endif /* BETA_DISTR_H_ */
