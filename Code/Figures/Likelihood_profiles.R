

det_est <- coef(summary(Rep_Dsole_det_1[[1]][[1]]))[,1]
det_pop <- Rep_Dsole_det_1[[1]][[2]]

#Dsole_data <- obsv.sims(bad_yr, 20, "Dsole", 3, -7, test_params[tt,1], F, 0.1)

fi_sims = seq(0.001, 0.4, length = 10) 
hy_a_sims = seq(0.1, 10, length = 10) 
hy_b_sims = seq(-15, 10, length = 10) 

fi_likel = NULL 
hy.a_likel = NULL
hy.b_likel = NULL
det_0.2_likel <- NULL
for (i in 1:10) { 
  fi_likel[i] <- p.filter(det_pop, det_est[1], det_est[2], log(fi_sims[i]), 10, det_est[4], 1, det_est[5], det_est[6], det_est[7], det_est[8], det_est[9], det_est[10], det_est[11], det_est[12], det_est[13], det_est[14], det_est[15], det_est[16], det_est[17], det_est[18], det_est[19], det_est[20], det_est[21], det_est[22], det_est[23], det_est[24])
  hy.a_likel[i] <- p.filter(det_pop, hy_a_sims[i], det_est[2], det_est[3], 10, det_est[4], 1, det_est[5], det_est[6], det_est[7], det_est[8], det_est[9], det_est[10], det_est[11], det_est[12], det_est[13], det_est[14], det_est[15], det_est[16], det_est[17], det_est[18], det_est[19], det_est[20], det_est[21], det_est[22], det_est[23], det_est[24])
  hy.b_likel[i] <- p.filter(det_pop, det_est[1], hy_b_sims[i], det_est[3], 10, det_est[4], 1, det_est[5], det_est[6], det_est[7], det_est[8], det_est[9], det_est[10], det_est[11], det_est[12], det_est[13], det_est[14], det_est[15], det_est[16], det_est[17], det_est[18], det_est[19], det_est[20], det_est[21], det_est[22], det_est[23], det_est[24])
  det_0.2_likel[i] <- p.filter(det_pop.2, det_est.2[1], det_est.2[2], log(fi_lg_sims[i]), 10, det_est.2[4], 1, det_est.2[5], det_est.2[6], det_est.2[7], det_est.2[8], det_est.2[9], det_est.2[10], det_est.2[11], det_est.2[12], det_est.2[13], det_est.2[14], det_est.2[15], det_est.2[16], det_est.2[17], det_est.2[18], det_est.2[19], det_est.2[20], det_est.2[21], det_est.2[22], det_est.2[23], det_est.2[24])
  
  }

par(mfrow=c(1,3))
plot(fi_sims, fi_likel, type='l', xlab= "Fishing Rate", ylab= "-LL", main ="Fishing Rate 0.008 Likelihood Slice")
abline(v=exp(det_est[3]), col="blue")

plot(hy_a_sims, hy.a_likel, type='l', xlab= "Hypoxia Parameter", ylab= "-LL", main ="Hypoxia_a Likelihood Slice")
abline(v=det_est[1], col="blue")

plot(hy_b_sims, hy.b_likel, type='l', xlab= "Hypoxia Parameter", ylab= "-LL", main ="Hypoxia_b Likelihood Slice")
abline(v=det_est[2], col="blue")

######################### No detection sims ##############


Nodet_est <- coef(summary(Rep_Dsole_NOdet_1[[1]][[1]]))[,1]
NOdet_pop <- Rep_Dsole_NOdet_1[[1]][[2]]
Nodet.2_est <- coef(summary(Rep_Dsole_NOdet_1[[2]][[1]]))[,1]
NOdet.2_pop <- Rep_Dsole_NOdet_1[[2]][[2]]

det_est.2 <- coef(summary(Rep_Dsole_det_1[[2]][[1]]))[,1]
det_pop.2 <- Rep_Dsole_det_1[[2]][[2]]

fi_lg_sims <- seq(0.05, 0.95, length = 10) 
10

Nodet_likel <- NULL
Nodet.2_likel <- NULL
for (i in 1:10){
  Nodet_likel[i] <- p.filter.NO.det(NOdet_pop, log(fi_sims[i]), 10, Nodet_est[2], 1, Nodet_est[3], Nodet_est[4], Nodet_est[5], Nodet_est[6], Nodet_est[7], Nodet_est[8], Nodet_est[9], Nodet_est[10], Nodet_est[11], Nodet_est[12], Nodet_est[13], Nodet_est[14], Nodet_est[15], Nodet_est[16], Nodet_est[17], Nodet_est[18], Nodet_est[19], Nodet_est[20], Nodet_est[21], Nodet_est[22])
  Nodet.2_likel[i] <- p.filter.NO.det(NOdet.2_pop, log(fi_lg_sims[i]), 10, Nodet.2_est[2], 1, Nodet.2_est[3], Nodet.2_est[4], Nodet.2_est[5], Nodet.2_est[6], Nodet.2_est[7], Nodet.2_est[8], Nodet.2_est[9], Nodet.2_est[10], Nodet.2_est[11], Nodet.2_est[12], Nodet.2_est[13], Nodet.2_est[14], Nodet.2_est[15], Nodet.2_est[16], Nodet.2_est[17], Nodet.2_est[18], Nodet.2_est[19], Nodet.2_est[20], Nodet.2_est[21], Nodet.2_est[22])
  det_0.2_likel[i] <- p.filter(det_pop.2, det_est.2[1], det_est.2[2], log(fi_lg_sims[i]), 10, det_est.2[4], 1, det_est.2[5], det_est.2[6], det_est.2[7], det_est.2[8], det_est.2[9], det_est.2[10], det_est.2[11], det_est.2[12], det_est.2[13], det_est.2[14], det_est.2[15], det_est.2[16], det_est.2[17], det_est.2[18], det_est.2[19], det_est.2[20], det_est.2[21], det_est.2[22], det_est.2[23], det_est.2[24])
}

# p.filter with hypoxia parameters 
par(mfrow=c(1,2))
plot(fi_sims, fi_likel, type='l', xlab= "Fishing Rate", ylab= "-LL", main ="Fishing Rate 0.008 Likelihood Slice")
abline(v=exp(det_est[3]), col="blue")

plot(fi_lg_sims, det_0.2_likel, type='l', xlab= "Fishing Rate", ylab= "-LL", main ="Fishing Rate 0.2 Likelihood Slice")
abline(v=exp(det_est.2[3]), col="blue")


# p.filter without hypoxia parameters 
par(mfrow=c(1,2))
plot(fi_sims, Nodet_likel, type='l', xlab= "Fishing Rate", ylab= "-LL", main ="Fishing Rate 0.008 Likelihood Slice NO Detection")
abline(v=exp(Nodet_est[1]), col="blue")

plot(fi_lg_sims, Nodet.2_likel, type='l', xlab= "Fishing Rate", ylab= "-LL", main ="Fishing Rate 0.2 Likelihood Slice NO Detection")
abline(v=exp(Nodet.2_est[1]), col="blue")

