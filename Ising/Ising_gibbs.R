# Gibbs Ising

gibbs_ising <- function(n_iter, beta, random, x_0) {
  # fonction gibbs
  taille <- length(x_0)
  u <- 1:taille
  for (j in 1:n_iter) {
    if (random) {
      # random scan
      u <- sample(1:taille, taille)
    }
    for (i in u) {
      x_0[i] <- sample(c(-1, 1), 1, prob = prob_cond(x_0, i, beta))
    }
  }
  return(x_0)
}

gibbs_ising_h <- function(n_iter, beta, random, x_0) {
  # fonction gibbs qui renvoie le hamiltonien de chaque iteration
  taille <- length(x_0)
  u <- 1:taille
  hamil <- numeric(n_iter)
  for (j in 1:n_iter) {
    if (random) { # random scan
      u <- sample(1:taille, taille)
    }
    for (i in u) {
      x_0[i] <- sample(c(-1, 1), 1, prob = prob_cond(x_0, i, beta))
    }
    hamil[j] <- hamiltonien(x_0)
  }
  return(hamil)
}
