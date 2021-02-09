##################
# inv -logit function
##################
inv_logit <- function(x){
  exp(x)/(exp(x)+1)
}
