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
source("./Library/p.filter.R")

# Format data
data <- bad_yr$Sim_DO_matrix[1,]
priors <- summ
Dsole_f <- Sp_depl[Sp_depl$Species == "Dsole", 3] #expected Dsole fishing rate based on SPR analysis
N <- length(data)
mesh = 200

# Set true parameter values 
#occ = 0.5
det = 0.5
hypox.p = 4

###############################################################################################

# Operating model wtih IPM 
Pop.eq <- IPM("Dsole", Dsole_f, N, mesh, NA, 0) # run IPM once to get to equilibrium 
Pop.var <- IPM("Dsole", Dsole_f, N, mesh, Pop.eq$SAD, 0.01) #IPM w/ stochasticity
Pop <- round(Pop.var$Pop.matrix)

# Take sample of population 
Pop.samp <- Pop* 0.001
Pop.sim <- rpois(length(Pop.samp), Pop.samp) # sample of pop serves as mean
dim(Pop.sim) <- dim(Pop.samp)

#' sampling - fraction of pop 
#' rpois - mean of sampling draw -> pop 
# encounters [0 or 1]
occ.full <- inv_logit(hypox.p * data) #linear dep. on hypoxia
enc <- rbinom(N, 1, occ.full) * rbinom(N, 1, det)  # gets 0 or 1 for detection


# Particle Filter Sim Pop & Likelihood
#Pop.count <- p.filter(Pop.eq$Pop.matrix, Pop.sim, "Dsole", 4, Dsole_f, NA, mesh, 100, N, 0.01, 0.1, 0.1) # uses the IPM sims w/ variation as data

# Accounting for false F's
##' cutpointr() library has rate calclulators



model <- stan_model(occ.mod5)
###############################  Creating data in list format necessary for Stan ####################################
occ.stan.dat <- list(N = N,
                     mesh = 200,
                     enc = enc,
                     hypox = dat,
                     nact = Pop.list,
                     Dat = Pop.eq$Pop.matrix,
                     fish = fish,
                     selc = 1,
                     q = 100,
                     Cv = 0.01,
                     Cv_q = 0.1,
                     Sigma_p = 0.1,
                     prior_mu = c(priors$mean[1], priors$mean[2], 0, Dsole_f),
                     prior_sd = c(priors$sd[1], priors$sd[2], 10, 0.01))


fit <- sampling(object = model, 
                data = occ.stan.dat,
                chains = 1, 
                cores = 3,
                iter = 2000, 
                warmup = 1000,
                control = list(adapt_delta = 0.95))

pairs(fit)
traceplot(fit)
summary(fit)$summary
