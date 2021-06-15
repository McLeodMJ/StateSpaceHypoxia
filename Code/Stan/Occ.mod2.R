
write("
	data{
		int<lower=1> N; //number of samples
		int<lower=0> enc[N]; //encounters
		vector[2] prior_mu;
		vector[2] prior_sd;
	}
	parameters{
		real logit_occ;
		real logit_det;
	}
	transformed parameters{
		real<lower=0,upper=1> occ; // occupancy
		real<lower=0,upper=1> det; // detection
		occ = inv_logit(logit_occ);
		det = inv_logit(logit_det);
	}
	model{

	  // weakly informative priors
		logit_occ ~ normal(prior_mu[1],prior_sd[1]);
		logit_det ~ normal(prior_mu[2],prior_sd[2]);

	// likelihoods
		for(i in 1:N){
			if(enc[i] > 0){
				//the site was occupied, and you detected it
				target += log_inv_logit(logit_occ) + bernoulli_logit_lpmf(1| logit_det );
				
			}else{
				//the site was occupied but you didn't detect it & the site was unoccupied
				target += log_sum_exp(log_inv_logit(logit_occ) + bernoulli_logit_lpmf(0| logit_det), log1m_inv_logit(logit_occ));
			}
		}
	}", 
      "occ2.stan")

occ.mod2 <- "occ2.stan"
