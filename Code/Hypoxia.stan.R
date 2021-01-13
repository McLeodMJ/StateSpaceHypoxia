# Hypoxia STAN Model

################################################################################################################################ 
#psi = prob. space is occupied by thst sp.
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
occ = 0.5 #occupied
det = 0.2 #det.
hypox.p = 1 # not sure what to do for this value which is used below in the inv-logit eqn but 1+ works

par_gs <- c(occ, det, hypox.p)

############################################ factoring in hypoxia as a covariate ############################################ 
occ.full <- inv_logit(logit(occ) + hypox.p * good_yr)
## Why inverset then logit fxn?
# encounters [0 or 1]
enc <- rbinom(N, 1, occ.full) * rbinom(N, 1, det)
table(enc) #good year: 1(observ.) makes up 9% and bad yr= 11%

###############################################  STAN MODEL  ############################################ 
occ.mod <- "
data{
	int<lower=1> N; //number of samples
	int<lower=0> enc[N]; //encounters
	vector[N] hypox; //hypoxia
	vector[3] par_gs; // parameter priors
	real<lower=0> cv_guess;
}
parameters{
  real<lower=0, upper=1> occ;
	real<lower=0, upper=1> det;
	real hypox_p; //hypoxia slope
}
transformed parameters{
	real logit_occ;
	real logit_det;
	
	logit_occ = logit(occ);
	logit_det = logit(det);
	
}
model{
  // uninformative priors
  occ ~ normal(par_gs[1], par_gs[1]* cv_guess);
	det ~ normal(par_gs[2], par_gs[2]* cv_guess);
	hypox_p ~ normal(par_gs[3], par_gs[3]* cv_guess);

// likelihoods
	for(i in 1:N){
		if(enc[i] > 0){
			//the site was occupied, and you detected it
			target += log_inv_logit(logit_det) + bernoulli_logit_lpmf(1| logit_occ + hypox_p*hypox[i]);
		}else{
			//the site was occupied but you didn't detect it & the site was unoccupied
			target += log_sum_exp(log_inv_logit(logit_det) + bernoulli_logit_lpmf(0| logit_occ + hypox_p*hypox[i]), log1m_inv_logit(logit_det));
		}
	}
}"
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

###############################  Creating data in list format necessary for Stan ####################################
occ.stan.dat <- list(N = N,
                     enc = enc,
                     hypox = good_yr,
                     par_gs = par_gs[1:3],
                     cv_guess = 0.1)
 


# Run Stan model 
occ.stan <- stan(model_code = occ.mod,
                 data = occ.stan.dat,
                 pars = c('occ', 'det', 'hypox_p'),
                 chains = 2,
                 warmup = 1000,
                 iter = 3000,
                 init = par_gs[1:3], #may need to be species specific 
                 control = list(adapt_delta = 0.95)) #target acceptance prob.
             

# Check Stan model runs
pairs(occ.stan)
# below the diagnal (below medial acceptance rates) shows red then it is likely an issue with model
# above the diagnal (above median acccpetance rates) red[divergence] - likely an issure with the adapt delta 
traceplot(occ.stan)
summary(occ.stan)$summary %>% head() # use to see the posterior estimates
launch_shinystan(occ.stan)


