#' Compute motif prior probabilities
#'
#' \code{compute_prior} computes prior probabilities for 5-mer motif sequences.
#'
#' @param motif_freq_exp A \code{data.frame} with columns "motif" and "freq". "motif" stores 5-mer motif sequences and "freq" stores the frequency of each motif among independently defined signal sites.
#' @param motif_freq_bg A \code{data.frame} with columns "motif", "tx_count", and "genome_count". "tx_count" and "genome_count" are motif counts in the transcriptome and genome, respectively.
#' @param target_motif A \code{string} representing the target motif pattern. MIRAGE expands supported degenerate bases such as N, R, Y, W, and H. The default is "DRACH"; set this explicitly for the assay being analyzed.
#' @param target_base_index A positive integer numeric value indicating the position of the target base within \code{target_motif}. Default is 3.
#' @param prob_methylation A positive numeric value indicating the genome- or transcriptome-wide prior probability that the target base carries true signal. This argument name is retained for backward compatibility. Default is 0.007.
#'
#' @return A \code{data.frame} including the information from \code{motif_freq_exp} and \code{motif_freq_bg}, plus the prior probability for each 5-mer motif sequence.
#'
#' @details \code{compute_prior} calculates the prior probability that a site is a true signal site given its local 5-mer motif sequence. Using the property of conditional likelihood, we compute the prior probability as follows.
#' Let \eqn{P(S=1|D_i)} denote the prior probability that site \eqn{i} carries true signal given local motif \eqn{D_i}, we have
#' \deqn{P(S=1|D_i=j) = \frac{P(D_i=j| S=1) P (S=1)}{ P(D_i=j)}}
#' where \eqn{P(D_i=j|S=1)} is the likelihood for each motif among true signal sites and can be calculated using motif frequencies from an independent signal set.
#' \eqn{P(S=1)} is the overall prior probability that the target base carries true signal.
#' \eqn{P(D_i=j)} is the likelihood of observing that motif in the transcriptome.
#'
#' @examples
#' # Build motif frequency tables from a PAR-CLIP count table.
#' count_table <- read.delim(system.file("extdata", "pos.read.count.example.parclip.txt",
#'                                       package = "MIRAGE"))
#' motif_counts <- table(count_table$motif)
#' motif_freq_exp <- data.frame(motif = names(motif_counts),
#'                              freq = as.numeric(motif_counts) / sum(motif_counts))
#'
#' motif_freq_bg <- data.frame(motif = names(motif_counts),
#'                             tx_count = as.numeric(motif_counts),
#'                             genome_count = as.numeric(motif_counts))
#'
#' # Compute the prior probability for each 5-mer motif sequence.
#' motif_prior <- compute_prior(motif_freq_exp, motif_freq_bg, target_motif = "NNUNN",
#'                              target_base_index = 3, prob_methylation = 0.01)
#'
#' @export

compute_prior <- function(motif_freq_exp, motif_freq_bg, target_motif = "DRACH", target_base_index = 3, prob_methylation = 0.007){

  motif_freq_exp[, 1] <- gsub("T", "U", motif_freq_exp[, 1])
  motif_freq_bg[, 1] <- gsub("T", "U", motif_freq_bg[, 1])
  target_motif <- gsub("T", "U", target_motif)

  target_base <- strsplit(target_motif,"")[[1]][target_base_index]
  if(!(target_base %in% c("A","C","U","G"))) stop("The target base must be one of A, C, U or G!")
  reg_exprs <- paste0(paste0(rep("[AUCG]",target_base_index-1),collapse=""),target_base,paste0(rep("[AUCG]",nchar(target_motif)-target_base_index),collapse=""))
  motif_freq_exp <- motif_freq_exp[grepl(reg_exprs,motif_freq_exp$motif),]
  motif_freq_bg <- motif_freq_bg[grepl(reg_exprs,motif_freq_bg[,1]),]

  motif_freq_all <- data.frame(motif = motif_freq_bg[, 1], tx_freq = motif_freq_bg[, 2]/sum(motif_freq_bg[, 2]), genome_freq = motif_freq_bg[, 3]/sum(motif_freq_bg[, 3]))

  motif_seqs <- apply(generate_motif_sequences(target_motif), 1, function(x){paste(x, collapse = "")})

  motif_freq_exp$motif_type <- ifelse(motif_freq_exp$motif %in% motif_seqs,target_motif,paste0("non",target_motif))
  freq.tab <- merge(motif_freq_all, motif_freq_exp, by.x = "motif", by.y = "motif", all = T)

  freq.tab$prior_methylated <- freq.tab$freq*prob_methylation/freq.tab$tx_freq
  return(freq.tab)

}
