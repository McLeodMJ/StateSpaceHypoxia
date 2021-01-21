### Testing initial values of MCMC

setwd("~/Box Sync/McLeod_thesis/code")
library(rstan)
library(bayesplot)
library(shinystan)
library('dplyr')
library(purrr)
library(readxl)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

############################################## Simulated Data ######################################################
load("DO_sims.Rdata")

good_yr <- DsD[DsD$Year == 2015,5]
bad_yr <- DsD[DsD$Year == 2018,5]


######################################### initual values to test #####################################################
inits <- read_excel("~/Box Sync/McLeod_thesis/Data/inits.csv")
inits$occ_bias <- rep(NA, nrow(inits))
inits$det_bias <- rep(NA, nrow(inits))
inits$hypox.p_bias <- rep(NA, nrow(inits))


############################################ Functions for transformations ##########################################
logit <- function(x){
  log(x/(1-x))
}
inv_logit <- function(x){
  exp(x)/(exp(x)+1)
}

############################################  START OF FXN: inital variables ##########################################################
# dat = good_yr or bad_yr
test.init <- function(dat) {
  N = length(dat)
  for (j in 1:nrow(inits)){
   
    occ = inits$occ[j]
    det = inits$det[j]
    hypox.p = inits$hypox.p[j]
    
    ############################################ factoring in hypoxia as a covariate ############################################ 
    occ.full <- inv_logit(logit(occ) + hypox.p * dat)
    
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
                         hypox = dat,
                         prior_mu = c(0, logit(det), 0),
                         prior_sd = c(1,1,1)) #flatter- 1

occ.stan <- stan(model_code = occ.mod2,
             pars = c('occ','det','hypox_p','logit_det','logit_occ'),
             data = occ.stan.dat,
             chains = 4,
             warmup = 2000,
             iter = 6000,
             control = list(adapt_delta = 0.95)) 
 
sum <- summary(occ.stan)

  # Create if else statement of convergence using rhat > 1.05 to eliminate posteriors that do not converge
  x <- rep(0, nrow(sum$summary))
    for(s in 1: nrow(sum$summary)){
      # insufficient rhat >=1.05
      if(sum$summary[s, 10] >= 1.05){
        x[s] = NA
        x <- x[is.na(x)]
      } 
      #insufficient effect size <= 100
      if(sum$summary[s, 9] <= 100){
        x[s] = NA
        x <- x[is.na(x)]
      }
    } #end small for loop
  
  # when insufficient --> make NA
    if(is.na(x)){
      inits$occ_bias[j] <- NA;
      inits$det_bias[j] <- NA;
      inits$hypox.p_bias[j] <- NA;
   }else{
     # when sufficient --> report posterior values
      inits$occ_bias[j] <-abs(sum$summary[1, 1] - inits$occ[j]) / sum$summary[1, 1] * 100;
      inits$det_bias[j] <- abs(sum$summary[2, 1] - inits$det[j]) / sum$summary[2, 1] * 100;
      inits$hypox.p_bias[j] <- abs(sum$summary[3, 1] - inits$hypox.p[j]) / sum$summary[3, 1] * 100;
    }
  } #end fxn for loop
  return(inits)
} # end fxn

### run the function for good and bad years
inits_good <- test.init(good_yr)
inits_bad <- test.init(bad_yr)


### check for smallest values of bias
inits_sm_g <- inits_good %>% 
  filter(occ_bias < 50) %>% 
  filter(det_bias < 50) %>% 
  filter(hypox.p_bias < 50) #0.75/.25/1 -- .75 /.5 /1 -- 0.75/.75/1 [1st - 2nd - 3rd]

inits_sm_b <- inits_bad %>% 
  filter(occ_bias < 50) %>% 
  filter(det_bias < 50) %>% 
  filter(hypox.p_bias < 50) #0.75/.75/1 -- .25 /.25 /1 -- 0.75/.25/1 [1st - 2nd - 3rd] then .75/.5/1


write.csv(zzz, "w_o_warnings_inits.csv")


########################### failed checks for warnings ######################################################
# check for warnings
tryCatch({
  warning("Examine the pairs() plot to diagnose sampling problems")
}, warning=function(w) {
  ## do something about the warning, maybe return 'NA'
  message("handling warning: ", conditionMessage(w))
  NA }) 

