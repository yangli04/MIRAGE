#' MIRAGE: Mutation-Encoded Inference of RNA Activity
#'
#' MIRAGE estimates site-level RNA activity from matched treatment and control
#' read-count tables produced by mutation-encoded RNA assays (for example
#' BACS, SHAPE-MaP, and PAR-CLIP). It jointly models the on-target
#' conversion rate, the background conversion rate, sequencing error, and
#' allele status, and returns a latent activity level \code{beta_est} for
#' every site.
#'
#' The main entry points are \code{\link{estimate_inference_with_empirical}}
#' (motif-free inference), \code{\link{estimate_inference_with_prior}} together
#' with \code{\link{compute_prior}} (motif-aware inference), and
#' \code{\link{plot_signal_scatter}} (diagnostics).
#'
#' @docType package
#' @name MIRAGE-package
#' @aliases MIRAGE
#' @keywords internal
#' @import doParallel
#' @importFrom foreach foreach %dopar% getDoParRegistered getDoParWorkers
#' @importFrom parallel makeCluster stopCluster
#' @importFrom dplyr left_join
#' @importFrom utils globalVariables
#' @importFrom stats binom.test dbinom fisher.test median optim p.adjust pbinom pchisq quantile
"_PACKAGE"

globalVariables(c("i", ".data"))
