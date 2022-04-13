# SS MLE estimates current stage
# 4/1/2022

require(broom)
require(httr)
require(invgamma)
require(bbmle)
require(dplyr)
require(tidyverse)
require(optimx)


load("../Data/DO_sims.Rdata")
# load("./Results/Priors_fromNO_Hypox.Rdata") # priors from MCMC bias runs [occ & det]?
source("./Library/logit.R")
source("./Library/inv.logit.R")

# Load IPM Files/ fxns
load("./Results/Expected_Fishing_rate.Rdata") #from the SPR vs fishing values ## F priors
#load("../Data/Hist_spp.Rdata") 
#load("../Data/MH.Sim.Dsole.Rdata")
source("./Library/IPM.R")
source("./Library/params.R")
source("./Library/kernmat.R")
source("./Library/fecmat.R")
source("./Library/p.filter.4.timesteps.R")
source("./Library/p.filter.R")
source("./Library/obsv.sims.R")
Sp_selex <- read.csv("../Data/Sp.selex.csv") # 

#load("../Results/Test_true&0.2_fi_updated/Rep_Dsole_det_15.Rdata")

pars <- params("Dsole")

#data Nact
#set.seed(1342)
Dsole_data1 <- obsv.sims(bad_yr, 4, "Dsole", 3, -7, log(pars$f), F, 0.1)

#bounds
#low <-c("hypox_a" =2,"hypox_b" = -14, "fi" = -1e2, "sigma_p"= -1e2, "rec1"= -1e6, "rec2"= -1e6,"rec3"= -1e6,"rec4"= -1e6)
#, "rec5"= -1e6,"rec6"= -1e6,"rec7"= -1e6, "rec8"= -1e6,"rec9"= -1e6,"rec10"= -1e6, "rec11"= -1e6, "rec12"= -1e6,"rec13"= -1e6,"rec14"= -1e6, "rec15"= -1e6,"rec16"= -1e6,"rec17"= -1e6, "rec18"= -1e6,"rec19"= -1e6,"rec20"= -1e6)
#upp <-c("hypox_a" = 4, "hypox_b" = 14, "fi" = 1e2, "sigma_p"= 1e2, "rec1"= 1e6, "rec2"= 1e6, "rec3"= 1e6, "rec4"= 1e6)
#,  "rec5"= 1e6,"rec6"= 1e6,"rec7"= 1e6, "rec8"= 1e6,"rec9"= 1e6,"rec10"= 1e6, "rec11"= 1e6, "rec12"= 1e6,"rec13"= 1e6,"rec14"= 1e6, "rec15"= 1e6,"rec16"= 1e6,"rec17"= 1e6, "rec18"= 1e6,"rec19"= 1e6,"rec20"= 1e6)


low <-c("hypox_a" =2,"hypox_b" = -14, "fi" = -6, "sigma_p"= -1e1, "rec1"= -1e2, "rec2"= -1e2,"rec3"= -1e2,"rec4"= -1e2) #, "rec5"= -1e2,"rec6"= -1e2,"rec7"= -1e2, "rec8"= -1e2,"rec9"= -1e2,"rec10"= -1e2, "rec11"= -1e2, "rec12"= -1e2,"rec13"= -1e2,"rec14"= -1e2, "rec15"= -1e2,"rec16"= -1e2,"rec17"= -1e2, "rec18"= -1e2,"rec19"= -1e2,"rec20"= -1e2)
upp <-c("hypox_a" = 8, "hypox_b" = 14, "fi" = 1e0, "sigma_p"= 1e1, "rec1"= 1e2, "rec2"= 1e2, "rec3"= 1e2, "rec4"= 1e2) #,  "rec5"= 1e2,"rec6"= 1e2,"rec7"= 1e2, "rec8"= 1e2,"rec9"= 1e2,"rec10"= 1e2, "rec11"= 1e2, "rec12"= 1e2,"rec13"= 1e2,"rec14"= 1e2, "rec15"= 1e2,"rec16"= 1e2,"rec17"= 1e2, "rec18"= 1e2,"rec19"= 1e2,"rec20"= 1e2)


#starting values
log_start <- list(hypox_a = 4,
                  hypox_b= 1,
                  fi = log(0.8),
                  cv_q = 0.5,
                  sigma_p = log(0.1),
                  rec1 =  1,
                  rec2 =  1,
                  rec3 =  1,
                  rec4 =  1) #,
                  rec5 =  1,
                  rec6 =  1,
                  rec7 =  1,
                  rec8 =  1,
                  rec9 =  1,
                  rec10 =  1,
                  rec11 =  1,
                  rec12 =  1,
                  rec13 =  1,
                  rec14 =  1,
                  rec15 =  1,
                  rec16 =  1,
                  rec17 =  1,
                  rec18 =  1,
                  rec19 =  1,
                  rec20 =  1)

start_time <- Sys.time()
start_time

#MLE - typically runs for < 2 minutes
SS_MLE.check7 <- mle2(minuslogl =  p.filter.4, #p.filter.wcbgts
                       start = log_start,
                       fixed = list(cv_q = 0.5),
                       data = list(dat = Dsole_data1,
                                   scale = 1),
                      optimizer = "optimx",
                      # method = "bobyqa",
                       lower = low, 
                       upper = upp,
                       control=list(maxit=2000, trace=5), skip.hessian = TRUE)
print(SS_MLE.check7)
end_time <- Sys.time()
start_time - end_time


# try just running optmix w/o bbmle
test <- mle2(minuslogl =  p.filter.4, #p.filter.wcbgts
             start = log_start,
             fixed = list(cv_q = 0.5),
             data = list(dat = Dsole_data1,
                         scale = 1),
             optimizer = "optimx",
             # method = "bobyqa",
             lower = low, 
             upper = upp,
             control=list(maxit=2000, trace=5, all.methods = T), skip.hessian = TRUE)




