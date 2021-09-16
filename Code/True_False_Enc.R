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
Sp_selex <- read.csv("../Data/Sp.selex.csv") 
###############################################################################################
# Format data
data = spline(bad_yr$Sim_DO_matrix[1,],n=200)$y
priors = summ
Dsole_f = Sp_depl[Sp_depl$Species == "Dsole", 3] #expected Dsole fishing rate based on SPR analysis
time = length(data)
mesh = 200
Q = 100
cv_q = 0.1

# Set true parameter values 
#occ = 0.5
det = 0.5
hypox_p = 4

###############################################################################################

# Operating model wtih IPM 
Pop.eq <- IPM("Dsole", Dsole_f, time, mesh, NA, 0) # run IPM once to get to equilibrium 
Pop.var <- IPM("Dsole", Dsole_f, time, mesh, Pop.eq$SAD, 0.01) #IPM w/ stochasticity
Pop <- round(Pop.var$Pop.matrix)

# Take sample of population 
Pop.samp <- Pop* 0.001
Pop.sim <- rpois(length(Pop.samp), Pop.samp) # sample of pop serves as mean
dim(Pop.sim) <- dim(Pop.samp)

###############################################################################################
#' sampling - fraction of pop 
# encounters [0 or 1]
occ.full <- inv_logit(hypox.p * data) #linear dep. on hypoxia
enc <- rbinom(time, 1, occ.full) * rbinom(time, 1, det)  # gets 0 or 1 for detection
###############################################################################################


# Particle filter: (following Knape & deValpine 2012)
N <- matrix(-999, nrow = nrow(Pop.eq$Pop.matrix), ncol = time)

# Generate Q particles (independent simulations of N):
Nf = array(rep(-999, nrow(N) * time * Q), c(nrow(N), time, Q)) #pop.dist. x time x particles
Nf[,1,] <- Pop.eq$Pop.matrix[,time] #saves last column of equil. data

# Simulate random variation for first particle
Nf[,1,] <- Nf[,1,] + rnorm(nrow(N), 0, cv_q * Pop.eq$Pop.matrix[,time])
Nf[Nf<0] = 0 # check for all values greater than zero 

# Step 1: Initialize the model 
# resample for accuracy
ftmp <- matrix(0, nrow= time, ncol = Q)
Avg.likel <- rep(-999, time) #average of the ftmp
ftmp[1,] <- rep(1, Q)
Wgt = cumsum(ftmp[1,]/ sum(ftmp[1,]))  
Rnd = runif(Q)
Wgt = matrix(Wgt, nrow=Q, ncol=Q) # same across row?
Rnd =  t(matrix(Rnd, nrow=Q, ncol=Q)) # same across column?
Pass <- matrix(Rnd < Wgt, nrow=Q, ncol= Q)
Ind = Q - colSums(Pass) + 1

Nf[,1,] <- Nf[ ,1, Ind] # replace with resampled values -- based on Ind variable that pulls the best representable samples 

Ind2 = sample(1:Q, 1) # index so must be less than 100
N[,1] = Nf[,1,Ind2] # pick one randomly to be *the* distribution to carry forward to the next step

## Step 2: simulate encounters for Null or absences based on hypoxia parameter  - DONE in p.filter model
## Step 3: run model through time T 
# Create the kernel:
pars <- params("Dsole")
WClen = pars$selc 

meshmax = pars$Linf * 2
x <- seq(from = 1, to = meshmax, length.out = mesh)
dx = diff(x)[1]

## size distribution of recruits
Rvec <- dnorm(x, pars$Rlen, pars$Rlen.sd)  
ones <- rep(1, Q)

K <- kernmat(x, pars, Dsole_f)
Fe <- fecmat(x, pars)
###############################################################################################
Parts <- rnorm(nrow(N), 0, cv_q * dat[,time])
iters = 2000
Unif <- array(rep(NA, iters  * Q), c(iters, Q)) 
Nrand = array(rep(NA, iters  * mesh), c(iters, mesh)) 
Rndmizer <- rep(NA, iters)
for(s in 1:iters){
  Rndmizer[s] <- sample(1:Q, 1)
  Unif[s,] <- runif(Q)
  Nrand[s,] = rnorm(mesh, 0, sigma.p)
}

model <- stan_model(occ.mod5)
###############################  Creating data in list format necessary for Stan ####################################
occ.stan.dat <- list(time = time,
                     mesh = mesh,
                     enc = enc,
                     hypox = dat,
                     Nact = Pop.sim,
                     Dat = Pop.eq$Pop.matrix,
                     fish = fish,
                     part = Parts,
                     rndmizer = Rndmizer,
                     unif = Unif,
                     Q = 100,
                     Cv = 0.01,
                     Cv_q = 0.1,
                     Sigma_p = 0.1,
                     prior_mu = c(priors$mean[1], priors$mean[2], 0, Dsole_f, 0),
                     prior_sd = c(priors$sd[1], priors$sd[2], 10, 0.01, 0.35))# rec deve is based on what I saw in the lit but need to check other species


fit <- sampling(object = model, 
                data = occ.stan.dat,
                chains = 1, 
                cores = 3,
                iter = iters, 
                warmup = 1000,
                control = list(adapt_delta = 0.95))

pairs(fit)
traceplot(fit)
summary(fit)$summary
