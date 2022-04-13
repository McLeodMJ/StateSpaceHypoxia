load("DO_sims.Rdata")
DO = spline(DO_data$Sim_DO_matrix[9,],n=20)$y
hypox_a = 3
hypox_b = -7

## hypoxia dependence
det <- inv_logit(DO * hypox_a + hypox_b) # detection is depen. on hypoxia and hypox param.

#logit plot of probability of detection
obs <- cbind(DO, det)
obs <- as.data.frame(obs)
plot(obs, ylim = c(0,1))

#plot of the detections with DO levels
ggplot(obs, aes(DO, det))+
  geom_point(position = "jitter")+
  geom_smooth(method="loess", linetype="dashed")+
  theme_classic()+
  xlab("DO Concentrations ml/l")+
  ylab("Probability of Detection")

Obs <- NULL
  for(t in 1:6){
    DO <- Dsole_dat$DO[[t]]$mean_hyp
    det <- inv_logit(DO * hypox_a + hypox_b) # detection is depen. on hypoxia and hypox param.
    
    #logit plot of probability of detection
    obs <- as.data.frame(cbind(DO, det))
    obs$Year <- rep(names(Dsole_dat$DO)[t], nrow(obs))
    Obs <- rbind(Obs, obs)
    
  }

#plot of the detections with DO levels
ggplot(Obs, aes(DO, det, color = Year))+
  geom_point(position = "jitter", aes(color= Year, height = 0.5))+
  geom_smooth(method="loess", linetype="dashed")+
  facet_wrap(~Year)+
  theme_classic()
