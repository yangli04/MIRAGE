MLE_for_lambda1 <- function(N_fixed_treatment, N_total_treatment, delta = 0.001, initial_lambda1 = 0.2){

  num.sites <- length(N_fixed_treatment)
  Likelihood <- function(lambda1) {
    # Target function
    aa <- 1 - delta/3
    bb <- lambda1*aa + (1-lambda1)*delta
    #nt <- 0
    #for(i in 1:num.sites){
    #  nt <- nt + (N_total_treatment[i]- N_fixed_treatment[i])*log(bb) + N_fixed_treatment[i]*log(1 - bb)
    #}
    nt <- sum(apply(cbind(N_fixed_treatment,N_total_treatment),1,function(x) {(x[2]-x[1])*log(bb)+x[1]*log(1-bb)}))
    -1*nt
  }
  res <- optim(initial_lambda1,  Likelihood, method = "Brent", lower = 0, upper = 1)
  return(res$par)

}
