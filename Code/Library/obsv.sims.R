## OBSERVATION MODEL
## Run Model selection on P.filter b4 running MCMC
#setwd("~/Box Sync/McLeod_thesis/Code")

#Load MCMC files/ fxns
# load("../Data/DO_sims.Rdata")
#load("./Results/Priors_fromNO_Hypox.Rdata") # summ - priors from MCMC bias runs [occ & det]?
#source("./Library/inv.logit.R")

# Load IPM Files/ fxns
#load("./Results/Expected_Fishing_rate.Rdata") #Sp_depl - from the SPR vs fishing values ## F priors
#load("../Data/Hist_spp.Rdata") 
#source("./Library/IPM.R")
#source("./Library/params.R")
#source("./Library/kernmat.R")
#source("./Library/fecmat.R")
#source("./Library/p.filter.R")put

####################################  Format data ####################################
# OBSERVATION MODEL 
#' DO_data : DO data from newport moorings or DO from wcgbts for real data
#' fish : fish of interest - Dsole, Ling, Yrock, Grock
#' hypox_a: dep. hypoxia parameter - based on sim. results that give at least 4 sim. absences - hypox_a= 3
#' hypox_b: intercept hypoxia parameter - if DO = 0 & prob is very low [0.001] logit(0.001/.999) --> hypox_b = -7
#' fi_spp : fishing rate dependent on spp. - if we want to fix fi, then make "NA" and it will call fi from the sp_depl.csv 
#' selec : selectivity of the trawl - F when simulating data & T when calling selec from param.R 


obsv.sims <- function(DO_data, fish, hypox_a, hypox_b, fi_spp, selec){

  #format DO_data by preferred length
DO = spline(DO_data$Sim_DO_matrix[9,],n=20)$y
time = length(DO)
mesh = 200
Q = 100

# simulations where fi rate is fixed 
if(is.na(fi_spp)){
# fishing rate by species
fi_spp = Sp_depl[Sp_depl$Species == fish, 3] #expected Dsole fishing rate based on SPR analysis
fi_spp = log(fi_spp)
}

### Operating model wtih IPM 
Pop.eq <- IPM(fish, fi_spp, time, mesh, NA, rep(0, time)) # run IPM once to get to equilibrium 
Pop.var.R <- IPM(fish, fi_spp, time, mesh, Pop.eq$SAD, rep(log(0.1), time)) #IPM w/ recruitment variation
Pop.samp <- round(Pop.var.R$Pop.matrix)

## hypoxia dependence
det <- inv_logit(DO * hypox_a + hypox_b) # detection is depen. on hypoxia and hypox param.
pres <- rbinom(time, 1, det) # pres/abs sims.
pres
#logit plot of probability of detection
 # obs <- cbind(DO, det)
 # plot(obs, ylim = c(0,1))

# Take poisson sample dist. of population 
Pop.sim <- rpois(length(Pop.samp), Pop.samp) # sample of pop serves as mean
dim(Pop.sim) <- dim(Pop.samp)
Pop.sim[ ,which(pres == 0)] <- rep(0, mesh)

## WCGBTS selectivity 
if(selec == T){ # indicate true if we want to include selectivity
  selec <- pars$selc 
}else{
  selec = 1
}

Pop.sim <- Pop.sim * selec


### Format data for estimation model
Data_log <- list(Q = Q,
     Mesh = mesh,
     time = time, 
     fish = fish,
     DO = DO,
     Nint = Pop.eq$SAD,
     N.act = Pop.sim,
     hypox_p = hypox_p, #normal dist.
     fi = fi_spp, #log-normal
     rec_var = rep(log(0.1), time), #log-normal
     sigma.p = log(0.5), #log-normal
     prior_mu = c(0, 0.1, 0, 0.1), #hypox, fi, rec_var, sigma.p
     prior_sd = c(10, 0.1, 0.1, NA))

return(Data_log)
}


