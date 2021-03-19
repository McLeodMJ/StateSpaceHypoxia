# Greenstriped rockfish 
# Fishing rate is referred to as exploitation rate
require(ggplot2)
require(patchwork)


source("Code/Library/params.R")
source("./Library/kernmat.R")
source("./Library/fecmat.R")


fish <- c("Lingcod", "Yrock", "Grock", "Dsole")

time = 1e2 #model runs in time
df <- data.frame( Time = rep(c(50:time), 4), 
                  Type = rep(fish, each = (length(50:time))),
                  Pop.Size =  rep(NA, length(c(50:time)*4 )))


for ( f in fish){
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
    K1 <- kernmat(x, pars, 1)
    Fe1<- fecmat(x, pars)
    
  ### Initialize the model:
    N1 = matrix(0, nrow = meshsize, ncol = time) #pop size w/ growth/mortal kernel
    N1[,1] <- Rvec * exp(pars$R0) # initialize with one pulse of new recruits
    E1 <- NULL
    Recruits1 <- NULL
    
  ### Run the model
  for (t in 2:time){
    N1[,t] <- K1 %*% N1[,t-1]  * dx  # midpoint rule integration
    E1[t] <- Fe1 %*% N1[,t-1] * dx
    
  ### BH Recruitment 
    # BH eqn: (Methot and Tylor 2011 - eqn A.7)
      # Ry = 4h*R0*Eggs / S0(1-h) + Eggs(5h -1)
    Recruits1[t] <- (4 * pars$steep * exp(pars$R0) * E1[t]) / ((pars$S0 * (1 - pars$steep)) + (E1[t] * (5 * pars$steep - 1)))
    N1[,t] = N1[,t] + Recruits1[t] * Rvec # + rnorm()*process noise 
    # Y[,t] = rnorm()*observation error
  } #end of model
  
   N01 <- N1[,time] # save SAD
   
   ### Run the model WITH variation
   Nv1 = matrix(0, nrow = meshsize, ncol = time)
   Nv1[,1] <- N01 #initalize with SAD
   E1 <- NULL
   Recruits1 <- NULL
   cv = 0.1 # coef. of variation to introduce noise around recruitment 
   
   for (t in 2:time){
     Nv1[,t] <- K1 %*% Nv1[,t-1]  * dx  # midpoint rule integration
     E1[t] <- Fe1 %*% Nv1[,t-1] * dx 

     Recruits1[t] <- (4 * pars$steep * exp(pars$R0) * E1[t]) / ((pars$S0 * (1 - pars$steep)) + (E1[t] * (5 * pars$steep - 1)))
     RR <- exp(rnorm(1, mean = log(Recruits1[t]) - ((cv * log(Recruits1[t]))^2)/2, sd= cv * log(Recruits1[t]) ))# change cv to 0 for no variation
     Nv1[,t] = Nv1[,t] + Recruits1[t] * Rvec * RR  
   }
   # lognormal variation - 
   
   pop <- colSums(Nv1[ ,50:time])
   df[df$Type == f, 3] <- pop
   #plot(pop, type='l')

} #end of f loop
#df$Pop.Size <- scale(df$Pop.Size) #noramlize data

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

