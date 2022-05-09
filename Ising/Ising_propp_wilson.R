# Propp Wilson Ising

propp_wilson_ising <- function(pas, beta, h, w) {
  # fonction propp wilson
  taille <- h * w
  m1 <- matrix(rep(-1, taille), h)
  m2 <- matrix(rep(1, taille), h)
  n_iter <- 0
  j <- c()
  v <- c()
  u <- c()
  while (identical(m1, m2) == FALSE) {
    n_iter <- n_iter + pas
    m1 <- matrix(rep(-1, taille), h)
    m2 <- matrix(rep(1, taille), h)
    j <- c(sample(seq_len(taille), pas, replace = TRUE), j) # choix du site
    v <- c(sample(c(-1, 1), pas, replace = TRUE), v) # choix de la nouvelle valeur du site
    u <- c(runif(pas), u)
    m1 <- mapping(m1, n_iter, j, v, u, beta)
    m2 <- mapping(m2, n_iter, j, v, u, beta)
  }
  return(m1)
}

propp_wilson_ising_h <- function(pas, n_iter, beta, h, w) {
  # fonction propp wilson qui renvoie le hamiltonien de chaque iteration
  taille <- h * w
  m1 <- matrix(rep(-1, taille), h)
  m2 <- matrix(rep(1, taille), h)
  n_iter <- 0
  j <- c()
  v <- c()
  u <- c()
  while (identical(m1, m2) == FALSE) {
    n_iter <- n_iter + pas
    m1 <- matrix(rep(-1, taille), h)
    m2 <- matrix(rep(1, taille), h)
    j <- c(sample(seq_len(taille), pas, replace = TRUE), j) # choix du site
    v <- c(sample(c(-1, 1), pas, replace = TRUE), v) # choix de la nouvelle valeur du site
    u <- c(runif(pas), u)
    m1 <- mapping(m1, n_iter, j, v, u, beta)
    m2 <- mapping(m2, n_iter, j, v, u, beta)
  }
  j <- sample(seq_len(taille), n_iter, replace = TRUE)
  v <- sample(c(-1, 1), n_iter, replace = TRUE)
  u <- runif(n_iter)
  return(mapping_h(m1, n_iter, j, v, u, beta))
}
