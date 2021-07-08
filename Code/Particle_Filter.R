
#require(ggplot2)
#require(patchwork)
setwd("~/Box Sync/McLeod_thesis/Code")
load("./Results/Expected_Fishing_rate.Rdata") #from the SPR vs fishing values 
load("../Data/DO_sims.Rdata")
load("../Data/Hist_spp.Rdata") # data made into histogram
Sp_selex <- read.csv("../Data/Sp.selex.csv")

Dsole_f <- Sp_depl[Sp_depl$Species == "Dsole", 3] #expected Dsole fishing rate based on SPR analysis
DO_data <- bad_yr$Sim_DO_matrix[1,]


source("./Library/params.R")
source("./Library/kernmat.R")
source("./Library/fecmat.R")

source("./Library/IPM.R")
source("./Library/p.filter.R")
source("./Library/inv.logit.R")


Pop_IPM_eq <- IPM("Dsole", Dsole_f, 100, 200, NA, 0) # run IPM once to get to equilibrium 
Pop_IPM_var <- IPM("Dsole",Dsole_f, 100, 200, Pop_IPM_eq$SAD, 0.01) # run IPM to include recruitment variation and use N0 from first IPM run [N0 is N[,time] from first IPM run]
Pop_var <- round(Pop_IPM_var$Pop.matrix)

Pfil_Out <- p.filter(Pop_IPM_eq$Pop.matrix, Pop_var, DO_data, 1, 4, Dsole_f, 200, 100, 100, 0.01, 0.1, 0.1) # particle filter using 1st IPM runs data to initialize the model and 2nd IPM run as the Nact data

par(mfrow=c(1,2))
plot(colSums(Pfil_Out$Pop), type='l'); plot(colSums(Pop_var), type='l') 




