require(rstan)
require(bayesplot)
require(shinystan)
require('dplyr')
require(purrr)
require(readxl)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

setwd("~/Documents/THESIS/StateSpaceHypoxia/Code")

#Load MCMC files/ fxns
  load("../Data/DO_sims.Rdata")
  load("./Results/Priors_fromNO_Hypox.Rdata") # priors from MCMC bias runs [occ & det]?
  source("./Library/logit.R")
  source("./Library/inv.logit.R")

# Load IPM Files/ fxns
  load("./Results/Expected_Fishing_rate.Rdata") #from the SPR vs fishing values ## F priors
  load("../Data/Hist_spp.Rdata") 
  source("./Library/IPM.R")
  source("./Library/params.R")
  source("./Library/kernmat.R")
  source("./Library/fecmat.R")
  
# Format data
  dat <- bad_yr$Sim_DO_matrix[1,]
  priors <- summ
  Dsole_f <- Sp_depl[Sp_depl$Species == "Dsole", 3] #expected Dsole fishing rate based on SPR analysis
 
# Set true parameter values 
  occ = 0.5
  det = 0.5
  hypox.p = 4
  N <- length(dat)

#################################### factoring in hypoxia as a covariate ############################################ 
# encounters [0 or 1]
  occ.full <- inv_logit(logit(occ) + hypox.p * dat)
  enc <- rbinom(N, 1, occ.full) * rbinom(N, 1, det)  
  #Pr <- as.data.frame(table(enc))$Freq[2] # counts how many are recorded as present
  
# operating model wtih IPM 
  Pop <- IPM("Dsole", Dsole_f, 100, 200, NA, 0) # run IPM once to get to equilibrium 
  Pop.var <- IPM("Dsole", Dsole_f,length(enc), 200, Pop$SAD, 0.01)
  
  Pop.list <- Pop.var$Pop.matrix

# Poisson sums
  Pop.count <- matrix(rpois(c(length(enc)*200), Pop.var$Pop.matrix), nrow= 200) # to run through poisson
  
  model <- stan_model(occ.mod4)
###############################  Creating data in list format necessary for Stan ####################################
occ.stan.dat <- list(N = N,
                     mesh = 200,
                     enc = enc,
                     hypox = dat,
                     pop = Pop.count,
                     lambda = Pop.list,
                     prior_mu = c(priors$mean[1], priors$mean[2], 0),
                     prior_sd = c(priors$sd[1], priors$sd[2], 10)) 

  #stan does not support NAs in data
  
fit <- sampling(object = model, 
                data = occ.stan.dat,
                chains = 3, 
                cores = 3,
                iter = 3000, 
                warmup = 1000,
                control = list(adapt_delta = 0.99))

pairs(fit)
traceplot(fit)
summary(fit)$summary
