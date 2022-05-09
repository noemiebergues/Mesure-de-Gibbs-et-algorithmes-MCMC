# Propp-Wilson compose chimique

propp_wilson_cc <- function(pas, beta, h, w) {
  # fonction propp wilson
  taille <- h * w
  m1 <- matrix(rep(-1, taille), h)
  m2 <- matrix(rep(1, taille), h)
  nb_iter <- 0
  j <- c()
  v <- c()
  u <- c()
  p <- exp(-beta * c(1, 0, 1))
  p <- p / sum(p)
  while (identical(m1, m2) == FALSE) {
    nb_iter <- nb_iter + pas
    m1 <- matrix(rep(-1, taille), h)
    m2 <- matrix(rep(1, taille), h)
    j <- c(sample(seq_len(taille), pas, replace = TRUE), j) # choix du site
    v <- c(sample(-1:1, pas, prob = p, replace = TRUE), v) # choix de la nouvelle valeur du site
    u <- c(runif(pas), u)
    m1 <- mapping(m1, nb_iter, j, v, u, beta, FALSE)
    m2 <- mapping(m2, nb_iter, j, v, u, beta, FALSE)
  }
  return(m1)
}
