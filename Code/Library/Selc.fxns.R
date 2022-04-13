


Dbl.norm <- function(selc.50, peak, asc, dsc){
  ### selectivity eqn A.1.30 [appendix]
  
  Len <- 1:100
  j1 <- 1/ 1 + exp(-20 *( (Len - selc.50) / ( 1 + abs(Len - selc.50))) )
  j2 <- 1/ 1 + exp(-20 *( (Len - peak) / ( 1 + abs(Len - peak))) )
  
  Selc <- asc *(1 - j1) + j1 *((1 - j2) + j2 * dsc)  
  return(Selc)
}
z <- Dbl.norm(par$selc.50, par$peak, par$asc, par$dsc )

### selectivity eqn A.1.29 [appendix]
#S = 1 /1 + e(-log(19)(L-B1) / B2)
## B1: 50% female B2: diff btw 95%-50% l: length at time t
Selc <- 1/ (1 + exp(-log(19) * (Size_c - pars$selc.50 ) / (pars$selc.100 - pars$selc.50)))

