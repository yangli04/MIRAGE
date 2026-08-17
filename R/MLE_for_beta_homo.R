MLE_for_beta_homo <- function(N_fixed_treatment, N_total_treatment, lambda1.est, lambda2.est, delta = 0.001, initial_beta = 0.2){

  Likelihood <- function(beta) {  # Target function
    aa <- lambda1.est*beta + lambda2.est*(1-beta)
    bb <- aa*(1-delta/3) + (1-aa)*delta
    nt <- (N_total_treatment - N_fixed_treatment)*log(bb) + N_fixed_treatment*log(1-bb)
    - 1*nt
  }

  res <- optim(initial_beta, Likelihood, method = "Brent", lower = 10^-6, upper = 1)
  return(res$par)

}


Bayesian_inference_for_beta_homo <- function(N_fixed_treatment, N_total_treatment, lambda1.est, lambda2.est, motif.prior.prob, delta = 0.001, initial_beta = 0.2){

  options(scipen = 16)

  aa <- lambda2.est
  bb <- aa*(1-delta/3) + (1-aa)*delta

  unmethyl_likelihood <- dbinom(N_total_treatment-N_fixed_treatment, N_total_treatment, bb)

  homo_Likelihood <- function(N_fixed_treatment, N_total_treatment, lambda1.est, lambda2.est){
    prob1 <- pbinom(N_total_treatment+1-N_fixed_treatment, N_total_treatment + 1, (1-delta*4/3)*lambda2.est) - pbinom(N_total_treatment+1-N_fixed_treatment, N_total_treatment+1, (1-delta*4/3)*lambda1.est)
    normalizing.factor <- (1-delta*4/3)*(lambda1.est -lambda2.est )*(N_total_treatment + 1)
    prob1/normalizing.factor
  }

  methyl_likelihood <- homo_Likelihood(N_fixed_treatment, N_total_treatment, lambda1.est, lambda2.est)

  aa <- (1-motif.prior.prob)*unmethyl_likelihood
  bb <- motif.prior.prob*methyl_likelihood

  posterior <- aa/(aa+bb)

  est_beta <- MLE_for_beta_homo(N_fixed_treatment, N_total_treatment, lambda1.est, lambda2.est, delta, initial_beta)

  return(c(est_beta, posterior))
}



