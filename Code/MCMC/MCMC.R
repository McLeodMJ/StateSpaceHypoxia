require(rstan)
require(bayesplot)
require(shinystan)
require('dplyr')
require(purrr)
require(readxl)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

setwd("~/Documents/THESIS/StateSpaceHypoxia/Code")
load("../Data/DO_sims.Rdata")
load("./Results/Priors_fromNO_Hypox.Rdata")
dat <- bad_yr$Sim_DO_matrix[1,]
priors <- summ

source("./Library/logit.R")
source("./Library/inv.logit.R")

occ = 0.5
det = 0.5
hypox.p = 4
N <- length(dat)
lambda = 4

  #################################### factoring in hypoxia as a covariate ############################################ 
# operating model wtih IPM  

occ.full <- inv_logit(logit(occ) + hypox.p * dat)
  
  # encounters [0 or 1]
  enc <- rbinom(N, 1, occ.full) * rbinom(N, 1, det)
   Pr <- as.data.frame(table(enc))$Freq[2]  # stores total number of sp. present interactions
   enc <- ifelse( enc == 1, rpois(Pr, lambda), enc) # gives present values a count number
  # if 1 then take data from sim IPM
   
  model <- stan_model(occ.mod3)
  ###############################  Creating data in list format necessary for Stan ####################################
  occ.stan.dat <- list(N = N,
                       enc = enc,
                       hypox = dat,
                       prior_mu = c(priors$mean[1], priors$mean[2], 0, 0),
                       prior_sd = c(priors$sd[1], priors$sd[2], 10, 10)) #flatter- 1 
  
  fit <- sampling(object = model, 
                  data = occ.stan.dat,
                  chains = 1, 
                  cores = 4,
                  iter = 8000, 
                  warmup = 2000,
                  control = list(adapt_delta = 0.99))

  pairs(fit)
  