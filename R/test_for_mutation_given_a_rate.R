# log_lik_ratio_test is referred to methylSig (https://github.com/sartorlab/methylSig/blob/devel/R/diff_binomial.R)

log_lik_ratio_test <- function(treatment_total_count,treatment_fixed_count,control_total_count,control_fixed_count) {
  control_X_count <- control_total_count-control_fixed_count
  treatment_X_count <- treatment_total_count-treatment_fixed_count
  ratio <- 2 * (control_X_count * log(control_X_count / control_total_count + 1e-100)
    + control_fixed_count * log(control_fixed_count / control_total_count + 1e-100)
    + treatment_X_count * log(treatment_X_count / treatment_total_count + 1e-100)
    + treatment_fixed_count * log(treatment_fixed_count / treatment_total_count + 1e-100)
    - (control_X_count+treatment_X_count) * log((control_X_count+treatment_X_count) / (control_total_count+treatment_total_count) + 1e-100)
    - (control_fixed_count+treatment_fixed_count) * log((control_fixed_count+treatment_fixed_count) / (control_total_count+treatment_total_count) + 1e-100)
  )
  p <- pchisq(ratio, 1, lower.tail=FALSE)
  return(c(ratio,p))
}

test_for_mutation_given_a_rate <- function(treatment_fixed_count, treatment_total_count, control_fixed_count, control_total_count, mutation_rate){

  binom_test_stats <- apply(cbind(treatment_total_count-treatment_fixed_count, treatment_total_count), 1, function(x){
    aa <- binom.test(x[1], x[2], p = mutation_rate, alternative = "greater")
    return(c(aa$p.value, aa$estimate))
  })
  binom_fdr <- p.adjust(binom_test_stats[1,], method="BH")
  fisher_test_pvalue <- apply(cbind(treatment_total_count-treatment_fixed_count,treatment_fixed_count,control_total_count-control_fixed_count,control_fixed_count),1,function(x) {fisher.test(matrix(x,nrow=2),alternative="greater")$p.value})
  fisher_fdr <- p.adjust(fisher_test_pvalue, method="BH")
  binom_lrt <- apply(cbind(treatment_total_count,treatment_fixed_count,control_total_count,control_fixed_count),1,function(x) {log_lik_ratio_test(x[1],x[2],x[3],x[4])})
  lrt_fdr <- p.adjust(binom_lrt[2,], method="BH")
  res <- data.frame(binom_p = binom_test_stats[1, ], est_X_rate = binom_test_stats[2, ], binom_fdr = binom_fdr, fisher_p = fisher_test_pvalue, fisher_fdr = fisher_fdr, lrt = binom_lrt[1,], lrt_p = binom_lrt[2,], lrt_fdr = lrt_fdr)
  return(res)

}


