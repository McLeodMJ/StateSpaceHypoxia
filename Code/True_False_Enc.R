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
dat <- bad_yr$Sim_DO_matrix[1,]
priors <- summ
Dsole_f <- Sp_depl[Sp_depl$Species == "Dsole", 3] #expected Dsole fishing rate based on SPR analysis
N <- length(dat)
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
Pop.temp <- rpois(length(Pop.samp), Pop.samp) # sample of pop serves as mean
dim(Pop.temp) <- dim(pop.temp)

#' sampling - fraction of pop 
#' rpois - mean of sampling draw -> pop 
# encounters [0 or 1]
occ.full <- hypox.p * dat #linear dep. on hypoxia
occ.prop <- occ.full/ occ.full[which.max(occ.full)] # new occ is >1 so this settles them proportional to max value
enc <- rbinom(N, 1, occ.prop) * rbinom(N, 1, det)  # gets 0 or 1 for detection

# predicted absences  
Pop.temp <- Pop 
for( i in 1:N){
  if(enc[i] == 0){
    Pop.temp[,i] <- rep(-999, 200) 
  }} #removes zeros to be 'NAs'

# Particle Filter Sim Pop & Likelihood
#Pop.count <- p.filter(Pop.eq$Pop.matrix, Pop.temp, "Dsole", 4, Dsole_f, NA, mesh, 100, N, 0.01, 0.1, 0.1) # uses the IPM sims w/ variation as data

# Accounting for false F's
##' cutpointr() library has rate calclulators
##' 