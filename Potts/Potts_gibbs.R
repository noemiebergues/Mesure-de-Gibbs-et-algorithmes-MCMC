# Gibbs Potts

gibbs_potts <- function(n_iter, beta, k, random, x_0) {
  # fonction gibbs
  taille <- length(x_0)
  u <- 1:taille
  for (j in 1:n_iter) {
    if (random) { # random scan
      u <- sample(1:taille, taille)
    }
    for (i in u) {
      x_0[i] <- sample(1:k, 1, prob = prob_cond(x_0, i, beta, k))
    }
  }
  return(x_0)
}
