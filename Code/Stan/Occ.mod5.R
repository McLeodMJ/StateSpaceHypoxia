
write("data{
    		int<lower=1> N; //number of samples
    		int<lower=1> mesh; //number of size bins
    		vector[N] enc; //list of encounters
    		vector[N] hypox; //hypoxia
    		vector[4] prior_mu;
    		vector[4] prior_sd;
    		// p.filter variables
    			int nact[mesh, N];
    			matrix[mesh, N] Dat;
    			vector[4] fish;
    			int<lower=1> selc;
    			int<lower=1> q;
    			int<lower=1> Time;
    			real<lower=0> Cv;
    		  real<lower=0> Cv_q;
    		  real<lower=0> Sigma_p;
    		  matrix[N,q] likel;
    		  
    	}
    	parameters{
    		real logit_det;
    		real hypox_p; //hypoxia slope
    	  real<lower=0,upper=1> fi; // not sure how to include fishing rate fi
    	}
    	transformed parameters{
    		real<lower=0,upper=1> det; // detection
    		det = inv_logit(logit_det);
    		matrix[N,q] likel;
    		likel = p.filter(Dat, nact, fish[1], hypox, fi, Selc, mesh, q, N, Cv, Cv_q, Sigma_p)$likelihood;
    	}
    	// rbinom()hypox.p
    	
    	model{
    	  // storage
    	    vector[N] occ_eff;
    
    	  // weakly informative priors
    		logit_occ ~ normal(prior_mu[1],prior_sd[1]);
    		logit_det ~ normal(prior_mu[2],prior_sd[2]);
    		hypox_p ~ normal(prior_mu[3],prior_sd[3]);
    		fi ~ lognormal(prior_mu[4],prior_sd[4]);
    
    	// likelihoods
    		for(i in 1:N){
    			occ_eff[i] = hypox_p*hypox[i];
    			
    			// Hurdle model
    			if(enc[i] == 0){
    			   	//the site was occupied but you didn't detect it & the site was unoccupied
    				target += log_sum_exp(occ_eff[i] + bernoulli_logit_lpmf(0| logit_det),  occ_eff[i]);
    				// lambda model prediction - enc = vector for each bin
    				
    			}else{
    				//the site was occupied, and you detected it
    				  target += occ_eff[i] + bernoulli_logit_lpmf(1| logit_det ) + likel[i] ;
    		// add likelihood from p.filter to binomial calc but need to add fi - fishing rate
    		// 
    		  }
    		}
    	}
      generated quantiles{
      
      
      }",  
      "occ5.stan")

occ.mod5 <- "occ5.stan"
model <- stan_model(occ.mod5) # remove when working
