# Hypoxia STAN Model

################################################################################################################################ 
#psi = prob. space is occupied by thst sp.
#p1 = prob. of detection w/ covariate
#hypox.p = variation around hypoxia
##########################################
setwd("~/Box Sync/McLeod_thesis/code")
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#library("devtools")
#devtools::install_github("andybunn/dplR")
library('dplR')



############################################## Simulated Data ######################################################
#load('Groundfish.Hyp.Rdata')
#RUN CODE BELOW TO REPLACE - detrend(Dsole$o2_at_gear_ml_per_l_der) 
#hypox <- lm(Lingcod$o2_at_gear_ml_per_l_der ~ Lingcod$year)$resid #Take the residuals of a stright line
### Should we be taking the absolute value so we do not have negative DO concentrations?

#z <- fft(hypox) #fourier transform
### DO WE RUN A sample() FROM HERE?

# Convert time to a more useful scale (proportions of a year):
#Time.fft = (1:length(z)) / 12
#plot(Time.fft, na.omit(z)^2, type='l', xlab='Frequency', ylab='Variance' ) #frequency is alrgely swayed to the left and right
load("DO_sims.Rdata")

good_yr <- DsD[DsD$Year == 2015,4]
bad_yr <- DsD[DsD$Year == 2018,4]
############################################ Functions for transformations ##########################################
logit <- function(x){
  log(x/(1-x))
}
inv_logit <- function(x){
  exp(x)/(exp(x)+1)
}


############################################  Initial variables # OK NOW NEW VALUES FOR THESE THINGS? ######################
N = length(bad_yr)
psi = 0.5 #occ.
p = 0.2 #det.
hypox.p = 1 # not sure what to do for this value which is used below in the inv-logit eqn but 1+ works



############################################ factoring in hypoxia as a covariate ############################################ 
psi.full <- inv_logit(logit(psi) + hypox.p * bad_yr)
## Why inverset then logit fxn?
# encounters [0 or 1]
enc <- rbinom(N, 1, psi.full) * rbinom(N, 1, p)


###############################################  STAN MODEL  ############################################ 
occ.mod <- "
data{
	int<lower=1> N; //number of samples
	int<lower=0> enc[N]; //encounters
	vector[N] hypox; //hypoxia
}
parameters{
	real logit_psi; //logit occupancy
	real logit_p; //logit detection
	real hypox_p; //hypoxia slope
}
transformed parameters{
	real<lower=0, upper=1> psi;
	real<lower=0, upper=1> p;

	psi = inv_logit(logit_psi);
	p = inv_logit(logit_p);
}
model{
  // uninformative priors
	logit_psi ~ normal(0,10);
	logit_p ~ normal(0,10);
	hypox_p ~ normal(0,10);

// likelihoods
	for(i in 1:N){
		if(enc[i] > 0){
			//the site was occupied, and you detected it
			target += log_inv_logit(logit_psi) + bernoulli_logit_lpmf(1| logit_p + hypox_p*hypox[i]);
		}else{
			//the site was occupied but you didn't detect it & the site was unoccupied
			target += log_sum_exp(log_inv_logit(logit_psi) + bernoulli_logit_lpmf(0| logit_p + hypox_p*hypox[i]), log1m_inv_logit(logit_psi));
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
                     hypox = bad_yr)

# Run Stan model 
occ.stan <- stan(model_code = occ.mod,
                 pars = c('psi','p','hypox_p'),
                 data = occ.stan.dat,
                 chains = 4,
                 warmup = 1000,
                 iter = 2000,
                 #init = inits, #may need to be species specific 
                 control = list(adapt_delta = 0.95)) #target acceptance prob.
             

# Check Stan model runs
pairs(occ.stan)
# below the diagnal (below medial acceptance rates) shows red then it is likely an issue with model
# above the diagnal (above median acccpetance rates) red[divergence] - likely an issure with the adapt delta 
traceplot(occ.stan)
summary(occ.stan)$summary %>% head() # use to see the posterior estimates
launch_shinystan(occ.stan)


