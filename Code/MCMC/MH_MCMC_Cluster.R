### MCMC - Metroplosis Hastings with Particle Filter likelihood 

#' Chain indicates which MCMC chain is currently being run (if not given, or given as NaN, 
#' then it will run the number of chains specified and run them in series (rather than parallel)

source("logit.R")
source("inv.logit.R")
source("pois.likel.R")

# Load IPM Files/ fxns
source("IPM.R")
source("params.R")
source("kernmat.R")
source("fecmat.R")
source("p.filter.R")
source("MH.MCMC.R")
Sp_selex <- read.csv("Sp.selex.csv") 
load("Expected_Fishing_rate.Rdata") #from the SPR vs fishing values ## F priors
load("MH.Sim.Dsole.Rdata")

library(matrixStats)
library(parallel)
library(foreach)
library(doParallel)
library(invgamma)
library(tidyverse)
library(dplyr)

# MCMC standards
chains =2
iters = 5000
  # 100,000
# burn = iters/4
# thin = 1 / length(parm_vec)
pars <- params("Dsole")

# test parameters over MCMC
test_params <- cbind( c(MH.Sim.Dsole$fi, seq(0.001, pars$M, length.out = 5), rep(pars$f, 5) ), # fi
c(MH.Sim.Dsole$hypox_p, rep(MH.Sim.Dsole$hypox_p, 5), seq(1, 10, length.out = 5)) ) #hypox_p

# multiResultClass <- function(result1=NULL,result2=NULL)
# {me <- list(
#     result1 = result1,
#     result2 = result2)
#   
#   ## Set the name for the class
#   class(me) <- append(class(me),"multiResultClass")
#   return(me)
# }


# for each to run MCMC fxn over multiple parameters
registerDoParallel(10)
MH.Sim.Dsole_Results <- foreach(p = 1:nrow(test_params)) %dopar% {
 list(index = p, MH.MCMC(MH.Sim.Dsole, p, test_params[p,], chains, iters, NA, NA, 100))
}
# registerDoParallel(4)
# MH.Sim.Dsole_Results <- foreach(p = 1:nrow(test_params[1:5,])) %dopar% {
#   list(index = p, MH.MCMC(MH.Sim.Dsole, p, test_params[p,], 2, 5, NA, NA, 100))
# }

#save(MH.Sim.Dsole_Results, file= paste('~/Documents/MCMC/Library/MH.Sim.Dsole_Results', name, '.Rda', sep=''))
########################################## Run on Machine ###############################################################


# # MCMC fxn
# MH.Sim.Dsole_Results <- MH.MCMC(MH.Sim.Dsole, test_params[1,], 1, 3, NA, NA)
#  Mh.data = MH.Sim.Dsole
# # parallel fxn
# #mclapply(zzz, MH.MCMC, mc.cores= 2) , chains, iters, NA, NA, MH.Sim.Dsole$DO, MH.Sim.Dsole$Nint, MH.Sim.Dsole$N.act, MH.Sim.Dsole$rec_var, MH.Sim.Dsole$prior_mu, MH.Sim.Dsole$prior_sd)
# #mclapply(norm_DO, test.init, mc.cores= 16, inits)
# 
# setwd("~/Documents/THESIS/StateSpaceHypoxia/Code")
# #setwd("~/Box Sync/McLeod_thesis/Code")#Load MCMC files/ fxns
# load("../Data/DO_sims.Rdata")
# load("./Results/Priors_fromNO_Hypox.Rdata") # priors from MCMC bias runs [occ & det]?
# source("./Library/logit.R")
# source("./Library/inv.logit.R")
# source("./Library/pois.likel.R")
# 
# # Load IPM Files/ fxns
# load("./Results/Expected_Fishing_rate.Rdata") #from the SPR vs fishing values ## F priors
# #load("../Data/Hist_spp.Rdata") 
# source("./Library/IPM.R")
# source("./Library/params.R")
# source("./Library/kernmat.R")
# source("./Library/fecmat.R")
# source("./Library/p.filter.R")
# source("./Library/MH.MCMC.R")
# Sp_selex <- read.csv("../Data/Sp.selex.csv") 
# 
# library(matrixStats)
# library(parallel)
# library(foreach)
# library(doParallel)
# library(invgamma)
# library(tidyverse)
# library(dplyr)
# 
# load("../Data/MH.Sim.Dsole.Rdata")

