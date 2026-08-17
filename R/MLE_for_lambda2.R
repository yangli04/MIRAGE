MLE_for_lambda2 <- function(N_fixed_treatment, N_total_treatment, num.sites = 50000, delta = 0.001, initial_lambda2 = 0.1, seed = NULL){

  if (!is.null(seed))
    set.seed(seed)

  total.len <- length(N_fixed_treatment)
  if(total.len <= num.sites){
    num.sites <- total.len
  }
  selected <- sample(1:total.len, num.sites)
  N_total_treatment <- N_total_treatment[selected]
  N_fixed_treatment <- N_fixed_treatment[selected]

  Likelihood <- function(lambda2) {
    # Target function
    aa <- 1 - delta/3
    bb <- lambda2*aa + (1-lambda2)*delta
    #nt <- 0
    #for(i in 1:num.sites){
    #  nt <- nt + (N_total_treatment[i]- N_fixed_treatment[i])*log(bb) + N_fixed_treatment[i]*log(1 - bb)
    #}
    nt <- sum(apply(cbind(N_fixed_treatment,N_total_treatment),1,function(x) {(x[2]-x[1])*log(bb)+x[1]*log(1-bb)}))
    -1*nt
  }
  res <- optim(initial_lambda2, Likelihood, method = "Brent", lower = 0, upper = 1)
  return(res$par)

}
