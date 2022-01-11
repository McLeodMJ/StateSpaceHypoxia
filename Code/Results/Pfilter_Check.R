## Run Model selection on P.filter b4 running MCMC

Pop.sim[ ,c(sample(time, 350))] <- rep(0, mesh)
#fish = "Dsole"
DO = spline(bad_yr$Sim_DO_matrix[1,],n=200)$y
Hypox_p = c(4, 1, 4, 4, 4)   # 2nd run is different hypox_p
Fi = c(Dsole_f, Dsole_f, 0.3, Dsole_f, Dsole_f) # 3rd run is different fi
#mesh = 200
#Q = 100
#time = 100
Rec_var <- matrix( rep(0.1, time), nrow= 5, ncol = time) 
Rec_var[4,] <- rep(0.6, time) # 4th run is different rec_var
#cv_q = 0.1
Sigma_p = c(0.1, 0.1, 0.1, 0.1, 100) # 5th run is different sigma_p

#source("./Library/p.filter.R")
library(ggplot2)
library(tidyverse)

Pfilter_Check <- data.frame(Param = c(rep('TRUE', time), 
                                      rep('hypox_p', time),
                                      rep('fi', time),
                                      rep('rec_var', time),
                                      rep('sigma_p', time)),
                            Time = rep(1:time, times = 5),
                            Avg_likel = rep(NA, time * 5))

check <- NULL
for( x in 1:length(Sigma_p)){
  z <- p.filter(Pop.eq$Pop.matrix, Pop.sim, DO, "Dsole", Hypox_p[x], Fi[x], 200, Q, time, Rec_var[x,], 0.01, Sigma_p[x])$likelihood
  check <- c(z, check)
}
Pfilter_Check$Avg_likel <- check

ggplot(Pfilter_Check, aes(x=Time, y= Avg_likel))+
  geom_line(aes(color = Param))

Pfilter_check <- Pfilter_Check %>% 
  group_by(Param) %>% 
  summarise(Mean = mean(Avg_likel)) %>% 
  mutate(AIC = -2 * Mean + 2 * 5) %>% 
  arrange(AIC)



