degenerate_bases <- list(
  A = c("A"),
  C = c("C"),
  G = c("G"),
  U = c("U"),
  R = c("A", "G"),
  Y = c("C", "U"),
  S = c("G", "C"),
  W = c("A", "U"),
  K = c("G", "U"),
  M = c("A", "C"),
  B = c("C", "G", "U"),
  D = c("A", "G", "U"),
  H = c("A", "C", "U"),
  V = c("A", "C", "G"),
  N = c("A", "C", "G", "U")
)

generate_motif_sequences <- function(mer) {
  bases <- strsplit(mer, "")[[1]]
  sequences <- expand.grid(lapply(bases, function(base) degenerate_bases[[base]]))
  colnames(sequences) <- 1:5
  return(sequences)
}
