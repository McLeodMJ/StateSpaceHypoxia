# May need to write this code so that we have just a vector of enc with 0 or 1 then index the 1s to be vectors of size bins

write("data{
    		int<lower=1> N; //number of samples
    		int<lower=1> mesh; //number of size bins
    		vector[N] enc; //list of encounters
    		vector[N] hypox; //hypoxia
    		vector[3] prior_mu;
    		vector[3] prior_sd;
    		matrix[mesh, N] lambda; //lambda is IPM pop matrix
    		int pop[mesh, N]; // matrix of just size bins that have sp present
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
    			
    			// Hurdle model
    			if(enc[i] == 0){
    			   	//the site was occupied but you didn't detect it & the site was unoccupied
    				target += log_sum_exp(log_inv_logit(occ_eff[i]) + bernoulli_logit_lpmf(0| logit_det),  log1m_inv_logit(occ_eff[i]));
    				// lambda model prediction - enc = vector for each bin
    				
    			}else{
    				//the site was occupied, and you detected it
    				  target += log_inv_logit(occ_eff[i]) + bernoulli_logit_lpmf(1| logit_det ) + sum(poisson_lpmf(pop[,i] | lambda[,i]) - log1m_exp(-lambda[,i]));
    			  // poisson_lpmf(enc[[i]] | lambda[,i])
    		}
        }
    	}", 
      "occ4.stan")

occ.mod4 <- "occ4.stan"
