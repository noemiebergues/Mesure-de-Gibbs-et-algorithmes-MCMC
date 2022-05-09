# Metropolis Hastings Ising

metropolis_hastings_ising <- function(n_iter, beta, x_0) {
  # fonction metropolis hastings
  j <- sample(seq_len(length(x_0)), n_iter, replace = TRUE) # choix du site
  v <- sample(c(-1, 1), n_iter, replace = TRUE) # choix de la nouvelle valeur du site
  u <- runif(n_iter)
  return(mapping(x_0, n_iter, j, v, u, beta))
}

metropolis_hastings_ising_h <- function(n_iter, beta, x_0) {
  # fonction metropolis hastings qui renvoie le hamiltonien de chaque iteration
  j <- sample(seq_len(length(x_0)), n_iter, replace = TRUE) # choix du site
  v <- sample(c(-1, 1), n_iter, replace = TRUE) # choix de la nouvelle valeur du site
  u <- runif(n_iter)
  return(mapping_h(x_0, n_iter, j, v, u, beta))
}
