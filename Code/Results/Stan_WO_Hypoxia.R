################################### Run the stan model W/O hypoxia 

setwd("~/Box Sync/McLeod_thesis/code")
library(rstan)
library(bayesplot)
library(shinystan)
library('dplR')
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source("Occ.model2.R")

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
# hypox.p = 1 # not sure what to do for this value which is used below in the inv-logit eqn but 1+ works

############################################ factoring in hypoxia as a covariate ############################################ 
#occ.full <- inv_logit(logit(occ) + hypox.p * good_yr)

# encounters [0 or 1]
enc <- rbinom(N, 1, occ) * rbinom(N, 1, det)

###############################  Creating data in list format necessary for Stan ####################################
occ.stan.dat <- list(N = N,
                     enc = enc,
                     prior_mu = c(0, logit(det)),
                     prior_sd = c(1,1)) #flatter- 1
############################################### Run the Stan Model ###############################################
model <- stan_model(occ.mod2)

fit <- sampling(object = model, 
                data = occ.stan.dat,
                chains = 4, 
                cores = 4,
                iter = 8000, 
                warmup = 2000,
                control = list(adapt_delta = 0.99))


sum <- as.data.frame( summary(fit)$summary )
summ <- sum[1:2, c(1,3)]
save(summ, file = "NOHypox_post.Rdata")
