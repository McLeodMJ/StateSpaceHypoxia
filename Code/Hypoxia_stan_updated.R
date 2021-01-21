################################################################################################################################ 
#occ = prob. space is occupied by thst sp.
#p1 = prob. of detection w/ covariate
#hypox.p = variation around hypoxia
##########################################
setwd("~/Box Sync/McLeod_thesis/code")
library(rstan)
library(bayesplot)
library(shinystan)
library('dplR')
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

############################################## Simulated Data ######################################################
load("DO_sims.Rdata")

good_yr <- DsD[DsD$Year == 2015,5]
bad_yr <- DsD[DsD$Year == 2018,5]


############################################ Functions for transformations ##########################################
logit <- function(x){
  log(x/(1-x))
}
inv_logit <- function(x){
  exp(x)/(exp(x)+1)
}

############################################  Initial variables # OK NOW NEW VALUES FOR THESE THINGS? ######################
N = length(good_yr)
occ = 0.75 #occ.
det = 0.25 #det.
hypox.p = 1 # not sure what to do for this value which is used below in the inv-logit eqn but 1+ works
#subtract true value from value/ true value --> SD

# calculate the bias and precis.
############################################ factoring in hypoxia as a covariate ############################################ 
occ.full <- inv_logit(logit(occ) + hypox.p * good_yr)
## Why inverset then logit fxn?
# encounters [0 or 1]
enc <- rbinom(N, 1, occ.full) * rbinom(N, 1, det)

###############################################  STAN MODEL  ############################################ 

occ.mod2 <- "
	data{
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
				target += log_sum_exp(log_inv_logit(occ_eff[i]) + bernoulli_logit_lpmf(0| logit_det), log1m_inv_logit(occ_eff[i]));
			}
		}
	}"
# target as a vector then take sum of all values for occupied and non-occupied. 

###############################  Creating data in list format necessary for Stan ####################################
occ.stan.dat <- list(N = N,
                     enc = enc,
                     hypox = good_yr,
                     prior_mu = c(0, logit(det), 0),
                     prior_sd = c(1,1,1)) #flatter- 1
# Run Stan model 
occ.stan <- stan(model_code = occ.mod2,
                 pars = c('occ','det','hypox_p','logit_det','logit_occ'),
                 data = occ.stan.dat,
                 chains = 4,
                 warmup = 2000,
                 iter = 6000,
                 control = list(adapt_delta = 0.95)) #target acceptance prob.


traceplot(occ.stan) 
sum <-summary(occ.stan)
pairs(occ.stan)
 

# recompile to avoid crashing 
# m <- stan_model(occ.mod2)
# fit <- sampling(m, data= occ.stan.dat)
#traceplot(fit)
