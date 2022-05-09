# Rejet CC

rejet_cc <- function(x) {
  # fonction rejet
  p <- exp(-beta * c(1, 0, 1))
  p <- p / sum(p)
  while (verif_e(x) == FALSE) {
    x <- matrix(sample(-1:1, length(x), replace = TRUE, prob = p), nrow = nrow(x))
  }
  return(x)
}