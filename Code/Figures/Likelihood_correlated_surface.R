# likelihood profile for correlated hypoxia parameters
fi_est <- seq(0.001, 0.5, length= 10)
hy.a_sim10 <- seq(1, 10, length = 10) # col
hy.b_sim10 <- seq(-12, 8, length= 10) # row
det_estt <- coef(summary(Rep_Dsole_det_1[[1]][[1]]))[,1]
dattaa <- Rep_Dsole_det_1[[1]][["MLE"]]@data[["dat"]]

#NLL <- matrix(NA, 10, 10)
#row.names(NLL) <- round(hy.b_sim10,1)
#colnames(NLL) <- round(hy.a_sim10,1)
NLL_corr <- array(rep(NA, 10*10*10), c(10,10,10))
# make a 3D array

for(i in 1:10){
  for(j in 1:10){
    for(k in 1:10){
    NLL_corr[i,j,k] <- p.filter(dattaa, hy.a_sim10[j], hy.b_sim10[i], fi_est[k], 10, det_estt[4], 1, det_estt[5], det_estt[6], det_estt[7], det_estt[8], det_estt[9], det_estt[10], det_estt[11], det_estt[12], det_estt[13], det_estt[14], det_estt[15], det_estt[16], det_estt[17], det_estt[18], det_estt[19], det_estt[20], det_estt[21], det_estt[22], det_estt[23], det_estt[24])
    }
  }
}

plot_ly(M, type="surface")
plot_ly(x=hy.b_sim10, y=hy.a_sim10, z=NLL_corr[,,4], type="surface")

# make a 3d plot from 3D array
library(reshape2)
library(rgl)
M=melt(NLL_corr)
points3d(M$Var1,M$Var2,M$Var3)

points3d(M$Var1,M$Var2,M$Var3,size=10,color=rainbow(10)[M$value*10])


contour(hy.b_sim10, hy.a_sim10, NLL_corr[,,1]/100000000, col = "pink", method = "edge", vfont = c("sans serif", "plain"), xlim=c(-12, -10), ylim=c(1,2))
