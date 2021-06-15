require(ggplot2)
require(patchwork)

setwd("/Users/montanamcleod/Documents/THESIS/StateSpaceHypoxia/Code")
load("./Results/Expected_Fishing_rate.Rdata") #from the SPR vs fishing values ## F priors


source("./Library/params.R")
source("./Library/kernmat.R")
source("./Library/fecmat.R")

fish <- c("Lingcod", "Yrock", "Grock", "Dsole")
time=100

df <- data.frame( Time = rep(c(50:time), 4), 
                  Type = rep(fish, each = (length(50:time))),
                  Pop.Size =  rep(NA, length(c(50:time)*4 )))

for (f in fish){
    
    pars <- params(f)
    # IPM integration parameters:
    meshsize = 200
    meshmin = 0
    meshmax = pars$Linf * 2
    
    x <- seq(from = meshmin, to = meshmax, length.out = meshsize) # different size classes
    dx = diff(x)[1] #width of the 'rectangle' to do the midpoint rule *cough cough* left rule
    
    ## size distribution of recruits
    Rvec <- dnorm(x, pars$Rlen, pars$Rlen.sd)  
    
    ## Kernel functions 
    K <- kernmat(x, pars, Sp_depl[Sp_depl$Species == f, 3], 1)
    Fe<- fecmat(x, pars)
    
    ### Initialize the model:
    N = matrix(0, nrow = meshsize, ncol = time) #pop size w/ growth/mortal kernel
    E <- NULL
    Recruits <- NULL
    
      N[,1] <- Rvec * exp(pars$R0) # SAD - not known - initialize with one pulse of new recruits
   
    ### Run the model
    for (t in 2:time){
      N[,t] <- K %*% N[,t-1]  * dx  # midpoint rule integration
      E <- Fe %*% N[,t-1] * dx
      
      ### BH Recruitment 
      # BH eqn: (Methot and Tylor 2011 - eqn A.7)
      # Ry = 4h*R0*Eggs / S0(1-h) + Eggs(5h -1)
      Recruits <- as.vector((4 * pars$steep * exp(pars$R0) * E) / ((pars$S0 * (1 - pars$steep)) + (E * (5 * pars$steep - 1))))
      N[,t] = N[,t] + Recruits * Rvec  # + rnorm()*process noise 
    } #end of model
    
    N0 <- N[,time] # save SAD

  
  ### Run the model WITH variation
  Nv = matrix(0, nrow = meshsize, ncol = time)
  Nv[,1] <- N0 #initalize with SAD
  E <- NULL
  Recruits <- NULL
  RR <- NULL
  cv = 0.01 # coef. of variation to introduce noise around recruitment 
  
  for (t in 2:time){
    Nv[,t] <- K %*% Nv[,t-1]  * dx  # midpoint rule integration
    E <- Fe %*% Nv[,t-1] * dx 
    
    Recruits <- as.vector((4 * pars$steep * exp(pars$R0) * E) / ((pars$S0 * (1 - pars$steep)) + (E * (5 * pars$steep - 1))))
    RR <- exp(rnorm(1, mean = cv *log(Recruits) - ((cv * log(Recruits))^2)/2, sd= cv * log(Recruits) ))# change cv to 0 for no variation
    # should there be another cv * in front of mean = log(recruits) ??????
    Nv[,t] = Nv[,t] + Recruits * Rvec * RR  
  }
  #df$Pop.Size <- scale(df$Pop.Size) #noramlize data

pop <- colSums(Nv[ ,50:time])
df[df$Type == f, 3] <- pop
} #end of f loop 



d <- ggplot(df[df$Type =="Dsole", ], aes(Time, Pop.Size))+
  geom_line(colour = "steelblue")+
  labs(title= "Dover Sole Model Simulations", y= "Dover Sole Biomass")

l <- ggplot(df[df$Type =="Lingcod", ], aes(Time, Pop.Size))+
  geom_line(color = "salmon1")+
  labs(title = "Lingcod Model Simulations", y="Lingcod Biomass")

g <- ggplot(df[df$Type =="Grock", ], aes(Time, Pop.Size))+
  geom_line(colour = "darkolivegreen")+
  labs(title="Greenstriped Rockfish Model Simulations",  y="Greenstriped Rockfish Biomass")

y <- ggplot(df[df$Type =="Yrock", ], aes(Time, Pop.Size))+
  geom_line(colour = "yellow3")+
  labs(title="Yellowtail Rockfish Model Simulations",  y="Yellowtail Rockfish Biomass")

(l | g) / (y | d) #organize plots

