# Metropolis Hastings compose chimique

metropolis_hastings_cc <- function(n_iter, beta, uniform = FALSE, x_0) {
  # fonction metropolis hastings
  if (uniform) {
    p <- rep(1 / 3, 3)
  } else {
    p <- exp(-beta * c(1, 0, 1))
    p <- p / sum(p)
  }
  j <- sample(seq_len(length(x_0)), n_iter, replace = TRUE)
  v <- sample(-1:1, n_iter, prob = p, replace = TRUE)
  u <- runif(n_iter)
  return(mapping(x_0, n_iter, j, v, u, beta, uniform))
}
