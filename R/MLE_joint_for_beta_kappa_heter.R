MLE_joint_for_beta_kappa_heter <- function(N_fixed_treatment, N_total_treatment, N_fixed_control, N_total_control, lambda1.est, lambda2.est, delta = 0.001, initial_beta = 0.2, initial_kappa = 0.5){

  Likelihood <- function(param) {

    aa <- lambda1.est*param[1] + lambda2.est*(1-param[1])
    bb <- aa*(1-delta/3) + (1-aa)*delta
    cc <- (1-param[2])*(1-delta/3)
    # Target function
    nt <- (N_total_treatment - N_fixed_treatment)*log(bb*param[2] + cc) + N_fixed_treatment*log(1-bb*param[2]-cc)
    nt <- nt + (N_total_control - N_fixed_control)*log(delta*param[2]+ cc) + N_fixed_control*log(1-delta*param[2]-cc)
    -1*nt
  }

  res <- optim(par=c(initial_beta, initial_kappa), fn=Likelihood, method = "L-BFGS-B", lower = c(10^-6, 10^-6), upper = c(1-10^-6, 1-10^-6))
  res$par <- round(res$par, 5)

  return(res$par)

}


MLE_for_beta_heter <- function(N_fixed_treatment, N_total_treatment, kappa, lambda1.est, lambda2.est,  delta = 0.001, initial_beta = 0.5){

  Likelihood <- function(beta) {  # Target function
    aa <- lambda1.est*beta + lambda2.est*(1-beta)
    bb <- aa*(1-delta/3) + (1-aa)*delta
    cc <- (1-kappa)*(1-delta/3)
    dd <- bb*kappa + cc
    nt <- (N_total_treatment - N_fixed_treatment)*log(dd) + N_fixed_treatment*log(1-dd)
    - 1*nt
  }

  res <- optim(initial_beta, Likelihood, method = "Brent", lower = 10^-6, upper = 1)
  return(res$par)


}


MLE_for_kappa_heter <- function(N_fixed_control, N_total_control, delta = 0.001, initial_kappa = 0.5){

  Likelihood <- function(kappa) {  # Target function
    bb <- (1-kappa)*(1-delta/3) + kappa*delta
    nt <- (N_total_control - N_fixed_control)*log(bb) + N_fixed_control*log(1-bb)
    - 1*nt
  }

  res <- optim(initial_kappa, Likelihood, method = "Brent", lower = 10^-6, upper = 1-10^-6)
  return(res$par)

}


Bayesian_inference_for_beta_kappa_heter <- function(N_fixed_treatment, N_total_treatment, N_fixed_control, N_total_control, lambda1.est, lambda2.est, motif.prior.prob, delta = 0.001){

  kappa_est <- MLE_for_kappa_heter(N_fixed_control, N_total_control)

  bb <- lambda2.est*(1-delta/3) + (1-lambda2.est)*delta
  cc <- (1-kappa_est)*(1-delta/3)
  dd <- bb*kappa_est + cc
  ee <- lambda1.est*(1-delta/3) + (1-lambda1.est)*delta

  unmethyl_likelihood <- dbinom(N_total_treatment-N_fixed_treatment, N_total_treatment, dd)

  prob2 <- pbinom(N_total_treatment+1-N_fixed_treatment, N_total_treatment+1, dd) - pbinom(N_total_treatment+1-N_fixed_treatment, N_total_treatment+1, cc + ee*kappa_est)
  normalizing.factor <- (1-delta*4/3)*(lambda1.est - lambda2.est)*(N_total_treatment+1)
  methyl_likelihood <- prob2/normalizing.factor

  aa <- (1-motif.prior.prob)*unmethyl_likelihood
  bb <- motif.prior.prob*methyl_likelihood
  posterior <- aa/(aa+bb)

  est_beta <- MLE_for_beta_heter(N_fixed_treatment, N_total_treatment, kappa_est, lambda1.est, lambda2.est)

  return(c(kappa_est, est_beta, posterior))
}




#
# Bayesian_inference_for_beta_kappa_heter <- function(N_fixed_treatment, N_total_treatment, N_fixed_control, N_total_control, lambda1.est, lambda2.est, motif.prior.prob, delta = 0.001, initial_beta = 0.2, initial_kappa = 0.5, num = 20){
#
#   para_beta <- seq(10^-3, 1, length.out = num)
#   para_kappa <- seq(10^-3, 1, length.out = num)
#
#   para_grid <- cbind(rep(para_beta, num), rep(para_kappa, each = num))
#
#   heter_Likelihood <- function(N_fixed_treatment, N_total_treatment, N_fixed_control, N_total_control, lambda2.est){
#     prob1 <- pbinom(N_total_control+1-N_fixed_control, N_total_control+1, 1-delta/3) - pbinom(N_total_control+1-N_fixed_control, N_total_control+1, delta)
#     prob2 <- pbinom(N_total_treatment+1-N_fixed_treatment, N_total_treatment+1, 1-delta/3) - pbinom(N_total_treatment+1-N_fixed_treatment, N_total_treatment+1, (1-delta/3)*lambda2.est + (1-lambda2.est)*delta)
#     normalizing.factor <- (1-delta*4/3)^2*(1 - lambda2.est)*(N_total_control+1)*(N_total_treatment+1)
#     prob1*prob2/normalizing.factor
#   }
#
#   unmethyl_likelihood <- heter_Likelihood(N_fixed_treatment, N_total_treatment, N_fixed_control, N_total_control, lambda2.est)
#   aa <- (1-motif.prior.prob)*unmethyl_likelihood
#
#   heter_methyal_approximate_Likelihood <- function(N_fixed_treatment, N_total_treatment, N_fixed_control, N_total_control, lambda1.est, lambda2.est){
#
#     single_likelihood <- function(beta, kappa){
#       aa <- lambda1.est*beta + lambda2.est*(1-beta)
#       bb <- aa*(1-delta/3) + (1-aa)*delta
#       cc <- (1-kappa)*(1-delta/3)
#       p_plus <- cc + bb*kappa
#       p_minus <- cc + delta*kappa
#       prob_plus <- dbinom(N_total_treatment-N_fixed_treatment, N_total_treatment, p_plus)
#       prob_minus <- dbinom(N_total_control-N_fixed_control, N_total_control, p_minus)
#       return(prob_plus*prob_minus)
#     }
#
#     all_param_likelihood <- apply(para_grid, 1, function(z){
#       single_likelihood(z[1], z[2])
#     })
#
#     approximated_likelihood <- mean(all_param_likelihood)
#     return(approximated_likelihood)
#    }
#
#   methyl_likelihood <- heter_methyal_approximate_Likelihood(N_fixed_treatment, N_total_treatment, N_fixed_control, N_total_control, lambda1.est, lambda2.est)
#   bb <- motif.prior.prob*methyl_likelihood
#   posterior <- aa/(aa+bb)
#
#   est_param <- MLE_joint_for_beta_kappa_heter(N_fixed_treatment, N_total_treatment, N_fixed_control, N_total_control, lambda1.est, lambda2.est, delta, initial_beta, initial_kappa)
#
#   return(c(est_param, posterior))
# }
#
#
#
