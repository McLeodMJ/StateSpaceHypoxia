fecmat <- function (x, pars){
  
  ### Maturity - size when females become mature
  # eqn A.1.16
  ## 1/ (1 +exp(slope of maturity * ( length - B1)))
  # b1= length at 50% maturity
  Matr <- 1 / (1 + exp(pars$Mat.slp * (x - pars$Mat.len))) 
  
  ### calculate fecundity by species instructions
  if(pars$name == "Grock" ){
    Fec <- x * (pars$fec.const + pars$fec.exp * x)
  }else{
    Fec <- pars$fec.int + pars$fec.const * pars$fec.slp * x^(pars$fec.exp)  #incldes both constant and eggs per weight params.
  }
  Fe <- Fec * Matr 
  return(Fe)
}
