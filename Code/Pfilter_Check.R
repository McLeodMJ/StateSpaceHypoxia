## Run Model selection on P.filter b4 running MCMC
setwd("~/Box Sync/McLeod_thesis/Code")

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
#Sp_selex <- read.csv("../Data/Sp.selex.csv") 

library(ggplot2)
library(tidyverse)
library(stats)

####################################  Format data ####################################
DO = spline(bad_yr$Sim_DO_matrix[9,],n=20)$y
Dsole_f = Sp_depl[Sp_depl$Species == "Dsole", 3] #expected Dsole fishing rate based on SPR analysis
time = length(DO)
mesh = 200
Q = 100
#cv_q = 0.1
#rec_dev <- rlnorm(time, 0.01, 0.1)
#f_rate <- rlnorm(1, Dsole_f, 0.01)


### Operating model wtih IPM 
Pop.eq <- IPM("Dsole", log(Dsole_f), time, mesh, NA, rep(0, time)) # run IPM once to get to equilibrium 
Pop.var.R <- IPM("Dsole", log(Dsole_f), time, mesh, Pop.eq$SAD, rep(log(0.1), time)) #IPM w/ recruitment variation

Pop.samp <- round(Pop.var.R$Pop.matrix)


# Take sample of population 
Pop.sim <- rpois(length(Pop.samp), Pop.samp) # sample of pop serves as mean
dim(Pop.sim) <- dim(Pop.samp)
Pop.sim[ ,c(sample(time, c(length(DO)/4)))] <- rep(0, mesh)


#actual parameters changing for each run
Hypox_p = c(4, 1, 4, 4, 4)   # 2nd run is different hypox_p
Fi = c(Dsole_f, Dsole_f, 0.5, Dsole_f, Dsole_f) # 3rd run is different fi
Rec_var <- matrix(rep(rec_dev, 5), nrow= 5, ncol = time) 
Rec_var[4,] <- 1/rec_dev # 4th run is different rec_var
Sigma_p = c(0.1, 0.1, 0.1, 0.1, 10) # 5th run is different sigma_p

#****p.filter() is expecting a vector of rec devs, not a scalar
################################################################################################################################################
## Nact data ???
#N.act <- Pop.sim
# 
# Pfilter_Check <- data.frame(Incorrect_Param = c(rep('sigma_p', time),
#                                                 rep('rec_var', time),
#                                                 rep('fi', time),
#                                                 rep('hypox_p', time),
#                                                 rep('TRUE', time)),
#                             Time = rep(1:time, times = 5),
#                             Avg_likel = rep(NA, time * 5))
# 
# check <- NULL
# for(c in 1:length(Sigma_p)){
#   z <- p.filter(Pop.eq$Pop.matrix, Pop.sim, DO, "Dsole", Hypox_p[c], Fi[c], mesh, Q, time, Rec_var[c,], 0.1, Sigma_p[c])$likelihood
#   check <- c(z, check)
# }
# Pfilter_Check$Avg_likel <- check
#  Pfilter_Check$Avg_likel <- ifelse(is.na(Pfilter_Check$Avg_likel), 0, Pfilter_Check$Avg_likel)
# 
#  
# # plot of differences between p.filter processing incoorrect parameters 
# ggplot(Pfilter_Check, aes(x=Time, y= Avg_likel))+
#   geom_line(aes(color = Incorrect_Param))
# 
# # summary of model selection
# Pfilter_check <- Pfilter_Check %>% 
#   group_by(Incorrect_Param) %>% 
#   summarise(Mean = mean(Avg_likel)) %>% 
#   mutate(AIC = -2 * Mean + 2 * 5) %>% 
#   arrange(AIC)
# 
# # plot of populations
# ggplot(Pfilter_Check, aes(x=Time, y=Pop))+
#   geom_line(aes(color = Incorrect_Param)) +
#   geom_line(aes(Time, True_data))
# 
# 
# 
###############################################################################################################


Pfilter_Check <- data.frame(Incorrect_Param = c(rep('sigma_p', time),
                                                rep('rec_var', time),
                                                rep('fi', time),
                                                rep('hypox_p', time),
                                                rep('TRUE', time)),
                            Time = rep(1:time, times = 5),
                            Avg_likel = rep(NA, time * 5))

check <- NULL
check_N <- NULL
for(c in 1:length(Sigma_p)){
  z <- p.filter(Pop.eq$Pop.matrix* 0.00001, round(Pop.sim  * 0.00001), DO, "Dsole", Hypox_p[c], Fi[c], mesh, Q, time, Rec_var[c,], 0.1, Sigma_p[c])$likelihood
  y <-p.filter(Pop.eq$Pop.matrix*0.00001, round(Pop.sim*0.00001), DO, "Dsole", Hypox_p[c], Fi[c], mesh, Q, time, Rec_var[c,], 0.1, Sigma_p[c])$Pop
  check <- c(z, check)
  check_N <- c(colSums(y), check_N)
}

Pfilter_Check$Avg_likel <- check
Pfilter_Check$Pop <- check_N
Pfilter_Check$Avg_likel <- ifelse(is.na(Pfilter_Check$Avg_likel), 0, Pfilter_Check$Avg_likel)
Pfilter_Check$Pop <- ifelse(is.na(Pfilter_Check$Pop), 0, Pfilter_Check$Pop)

Pfilter_Check$True_data <- rep(colSums(round(Pop.sim*0.00001)), 5)

# plot of differences between p.filter processing incoorrect parameters 
ggplot(Pfilter_Check, aes(x=Time, y= Avg_likel))+
  geom_line(aes(color = Incorrect_Param))

#Pfilter_Check[Pfilter_Check$Incorrect_Param != "rec_var" & Pfilter_Check$Incorrect_Param != "sigma_p",]
# plot of populations
ggplot(Pfilter_Check, aes(x=Time, y=Pop))+
  geom_line(aes(color = Incorrect_Param)) +
  geom_line(aes(Time, True_data))

# summary of model selection
Pfilter_check <- Pfilter_Check %>% 
  group_by(Incorrect_Param) %>% 
  summarise(Mean = mean(Avg_likel)) %>% 
  mutate(AIC = -2 * Mean + 2 * 5) %>% 
  arrange(AIC)

Pfilter_check









################################################################ Just zeros ################################################################
# changing parameters 
Hypox_p = c(4, 1, 4, 4, 4)   # 2nd run is different hypox_p
Fi = c(0, 0, 0.5, 0, 0) # 3rd run is different fi
Rec_var <- matrix( rep(1, time), nrow= 5, ncol = time) 
Rec_var[4,] <- rep(0.1, time) # 4th run is different rec_var
Sigma_p = c(0, 0, 0, 0, 0.1) # 5th run is different sigma_p

## simulated population data at fi = 0 & no variation
# Operating model wtih IPM 
Pop.eq <- IPM("Dsole", 0, time, mesh, NA, 0) # run IPM once to get to equilibrium 
Pop.var.R.F <- IPM("Dsole", 0, time, mesh, Pop.eq$SAD, 0.01) #IPM w/ recruitment variation 
Pop.samp2 <- round(Pop.var.R.F$Pop.matrix * 0.001)


Pop.sim2 <- rpois(length(Pop.samp2), Pop.samp2) # sample of pop serves as mean
dim(Pop.sim2) <- dim(Pop.samp2)
Pop.sim2[ ,c(sample(time, 30))] <- rep(0, mesh)

popp <- round(Pop.eq$Pop.matrix * 0.001)
Pop.sim3 <- rpois(length(popp), popp) # sample of pop serves as mean
dim(Pop.sim3) <- dim(popp)
Pop.sim3[ ,c(sample(time, 30))] <- rep(0, mesh)

pop <- list( pop = Pop.sim3,
             pop1 =  Pop.sim3,
            pop_var = Pop.sim2,
             pop_var = Pop.sim2)

Pfilter_Check2 <- data.frame(Param = c(rep('fi_var', time), 
                                       rep('TRUE_var', time),
                                       rep('fi', time),
                                       rep('TRUE', time)),
                             Time = rep(1:time, times = 4),
                             Avg_likel = rep(NA, time * 4))

check<- NULL
for(c in 1:4){
  z <- p.filter(Pop.eq$Pop.matrix, pop[[c]], DO, "Dsole", 4, c(0,0.5, 0, 0.5)[c], mesh, Q, time, rep(1, time), 0.1, 0)
  check <- c(z, check)
}
Pfilter_Check2$Avg_likel <- check
Pfilter_Check2$Avg_likel <- ifelse(is.na(Pfilter_Check2$Avg_likel), 0, Pfilter_Check2$Avg_likel)


# plot of differences between p.filter processing incoorrect parameters 
ggplot(Pfilter_Check2, aes(x=Time, y= Avg_likel))+
  geom_line(aes(color = Param))

# summary of model selection
Pfilter_check <- Pfilter_Check2 %>% 
  group_by(Param) %>% 
  summarise(Mean = mean(Avg_likel)) %>% 
  mutate(AIC = -2 * Mean + 2 * 5) %>% 
  arrange(AIC)

