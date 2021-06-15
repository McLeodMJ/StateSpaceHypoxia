write("data{
    		int<lower=1> N; //number of samples
    		int<lower=0> enc[N]; //encounters
    		vector[N] hypox; //hypoxia
    		vector[4] prior_mu;
    		vector[4] prior_sd;
    	}
    	parameters{
    		real logit_occ;
    		real logit_det;
    		real hypox_p; //hypoxia slope
    		real lambda; // poisson param
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
    		lambda ~ normal(prior_mu[4],prior_sd[4]);
    
    	// likelihoods
    		for(i in 1:N){
    			occ_eff[i] = logit_occ + hypox_p*hypox[i];
    			if(enc[i] > 0){
    				//the site was occupied, and you detected it
    				target += log_inv_logit(occ_eff[i]) + bernoulli_logit_lpmf(1| logit_det ) + poisson_lpmf(enc[i] | lambda) - log1m_exp(-lambda);
    				
    			}else{
    				//the site was occupied but you didn't detect it & the site was unoccupied
    				target += log_sum_exp(log_inv_logit(occ_eff[i]) + bernoulli_logit_lpmf(0| logit_det),  log1m_inv_logit(occ_eff[i]));
    			}
    		}
    	}", 
      "occ3.stan")

occ.mod3 <- "occ3.stan"





