## Plot of logit detectabilties based on NLL plot values

det_est <- coef(summary(Rep_Dsole_det_1[[1]][[1]]))[,1]
#det_pop <- Dsole_tidyr_det_4[[1]][["MLE"]]@data[["dat"]][["N.act"]]
det_DO <- Rep_Dsole_det_1[[1]][["MLE"]]@data[["dat"]][["DO"]]

P.det <- inv_logit(det_DO * det_est[1] + det_est[2])
p.det.high <- inv_logit(det_DO * 0.5 + -12)
p.det.low <- inv_logit(det_DO * 8 + 5)


Obs <-data.frame(DO = rep(det_DO, 3),
           Scenario = rep(c("True", "High NLL", "Low NLL"), each = 20),
           Detectability = c(P.det, p.det.high, p.det.low))

ggplot(Obs, aes(DO, Detectability, color= Scenario))+
  geom_line()
