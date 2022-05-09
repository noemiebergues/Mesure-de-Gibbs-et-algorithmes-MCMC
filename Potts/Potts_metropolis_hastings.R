# Metropolis-Hastings Potts

metropolis_hastings_potts <- function(n_iter, beta, k, x_0) {
  # fonction metropolis hastings
  j <- sample(seq_len(length(x_0)), n_iter, replace = TRUE) # choix du site
  v <- sample(1:k, n_iter, replace = TRUE) # choix de la nouvelle valeur du site
  u <- runif(n_iter)
  return(mapping(x_0, n_iter, j, v, u, beta))
}
