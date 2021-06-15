
write("data{
    		int<lower=1> N; //number of samples
    		int<lower=0> enc[N]; //encounters
    		vector[N] hypox; //hypoxia
    		vector[3] prior_mu;
    		vector[3] prior_sd;
    	}
    	parameters{
    		real logit_occ;
    		real logit_det;
    		real hypox_p; //hypoxia slope
    	}
    	transformed parameters{
    		real<lower=0,upper=1> occ; // occupancy
    		real<lower=0,upper=1> det; // detection
    		occ = inv_logit(logit_occ);
    		det = inv_logit(logit_det);
    	}
    	model{
    	  // storage
    	    vector[N] occ_eff;
    
    	  // weakly informative priors
    		logit_occ ~ normal(prior_mu[1],prior_sd[1]);
    		logit_det ~ normal(prior_mu[2],prior_sd[2]);
    		hypox_p ~ normal(prior_mu[3],prior_sd[3]);
    
    	// likelihoods
    		for(i in 1:N){
    			occ_eff[i] = logit_occ + hypox_p*hypox[i];
    			if(enc[i] > 0){
    				//the site was occupied, and you detected it
    				target += log_inv_logit(occ_eff[i]) + bernoulli_logit_lpmf(1| logit_det );
    				
    			}else{
    				//the site was occupied but you didn't detect it & the site was unoccupied
    				target += log_sum_exp(log_inv_logit(occ_eff[i]) + bernoulli_logit_lpmf(0| logit_det),  log1m_inv_logit(occ_eff[i]));
    			}
    		}
    	}", 
      "occ1.stan")

occ.mod1 <- "occ1.stan"


# target +=  includes all additive constants in the log density 
# is adding whatever is after it to the target density (essentially like the joint density of params)
# log_inv_logit(): natural logarithm of the inverse logit function of x
# bernoulli_logit_lpmf(): The log Bernoulli probability mass of y given chance of success inv_logit(alpha)
## chance-of-success param. is more stable if in logit scale & factors additive terms
# log_sum_exp(): Return the natural logarithm of the sum of the natural exponent of x and the natural exponent of y
# log1m_inv_logit(): natural logarithm of 1 minus the inverse logit function of x
#jacobian adjustment - aabsolute derivative of the inverse of the transformation; not present/detected
### posterior depends on parameterization 
#target += log_inv_logit(logit_p) + bernoulli_logit_lpmf(1| logit_psi + hypox_p*hypox[i]);
