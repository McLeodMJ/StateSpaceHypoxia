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

#time series simulation from spectrum analysis
load("DO_sims.Rdata")

good_yr <- good_yr$Sim_DO_matrix[1,] 
bad_yr <- bad_yr$Sim_DO_matrix[1,]

############################################ Functions for transformations ##########################################
source("./Library/logit.R")
source("./Library/inv.logit.R")

############################################  Initial variables # OK NOW NEW VALUES FOR THESE THINGS? ######################
N = length(good_yr)
occ = 0.5 #occ.
det = 0.5 #det.
hypox.p = 1 # not sure what to do for this value which is used below in the inv-logit eqn but 1+ works

############################################ factoring in hypoxia as a covariate ############################################ 
occ.full <- inv_logit(logit(occ) + hypox.p * good_yr)

# encounters [0 or 1]
enc <- rbinom(N, 1, occ.full) * rbinom(N, 1, det)

###############################################  STAN MODEL  ############################################ 
source("Occ.mod1.R")
model <- stan_model(occ.mod1) # recompile to avoid crashing 

###############################  Creating data in list format necessary for Stan ####################################
occ.stan.dat <- list(N = N,
                     enc = enc,
                     prior_mu = c(0, logit(det)),
                     prior_sd = c(1,1)) #flatter- 1
# Run Stan model 

fit <- sampling(object = model, 
                data = occ.stan.dat,
                chains = 4, 
                cores = 4,
                iter = 8000, 
                warmup = 2000,
                control = list(adapt_delta = 0.99))


sum <- as.data.frame( summary(fit)$summary )

traceplot(fit) 
sum <-summary(fit)
pairs(fit)
 
# ANother possible way to write the likelihood
################################################################################################################
// likelihoods
for(i in 1:N){
  if(enc[i] > 0){
    //the site was occupied, and you detected it
    1 ~ bernoulli_logit(logit_psi + hypox_p*hypox[i]);
    enc[i] ~ bernoulli_logit(logit_p);
  }else{
    //the site was occupied but you didn't detect it & the site was unoccupied
				target += log_sum_exp(bernoulli_logit_lpmf(1 | logit_psi + hypox_p*hypox[i]) + bernoulli_logit_lpmf(0 | logit_p), bernoulli_logit_lpmf(0 | logit_psi + hypox_p*hypox[i]));
			}
		}
		
