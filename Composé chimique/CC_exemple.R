# Exemple compose chimique

rm(list = ls())
set.seed(23)

verif_e <- function(x) {
  # fonction qui verifie si la matrice x appartient à E
  return(sum(x[-nrow(x), ] * x[-1, ] == -1) + sum(x[, -ncol(x)] * x[, -1] == -1) == 0)
}

verif_e_loc <- function(x, i, j) {
  # fonction qui vérifie si la contrainte de E est vérifiée pour le site i de x qui prend la valeur j et ses voisins
  return(sum(j * voisins(x, i) == -1) == 0)
}

voisins <- function(x, i) {
  # fonction qui renvoie les voisins du site i de la matrice x
  h <- nrow(x)
  v <- rep(2, 4)
  if (i %% h != 1) {
    v[1] <- x[i - 1]
  }
  if (i %% h != 0) {
    v[2] <- x[i + 1]
  }
  if (i > h) {
    v[3] <- x[i - h]
  }
  if (i < length(x) - h + 1) {
    v[4] <- x[i + h]
  }
  return(v[v != 2])
}

hamiltonien <- function(x) {
  # hamiltonien
  return(sum(abs(x)))
}

prob_cond <- function(x, i, beta) {
  # calcul des probas conditionnelles
  p <- exp(-beta * c(1, 0, 1))
  p <- p * (colSums(voisins(x, i) %*% t(-1:1) == -1) == 0)
  return(p / sum(p))
}

mapping <- function(x, n_iter, j, v, u, beta, uniform = TRUE) {
  test <- TRUE
  for (i in 1:n_iter) {
    if (uniform) {
      test <- u[i] <= exp(-beta * (abs(v[i]) - abs(x[j[i]])))
    }
    if (test & verif_e_loc(x, j[i], v[i])) {
      x[j[i]] <- v[i]
    }
  }
  return(x)
}

rinit_cc <- function(h, w) {
  x <- matrix(rep(0, h * w), h)
  i <- 0
  while (i < h) {
    x[i, ] <- rep(sample(-1:1, 1, replace = TRUE, prob = c(2 / 5, 1 / 5, 2 / 5)), w)
    i <- i + 2
  }
  return(x)
}

rejet_cc <- function(x) {
  # fonction rejet
  p <- exp(-beta * c(1, 0, 1))
  p <- p / sum(p)
  while (verif_e(x) == FALSE) {
    x <- matrix(sample(-1:1, length(x), replace = TRUE, prob = p), nrow = nrow(x))
  }
  return(x)
}

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

gibbs_cc <- function(n_iter, beta, random, x_0) {
  # fonction gibbs
  taille <- length(x_0)
  u <- 1:taille
  for (j in 1:n_iter) {
    if (random) {
      # random scan
      u <- sample(1:taille, taille)
    }
    for (i in u) {
      x_0[i] <- sample(-1:1, 1, prob = prob_cond(x_0, i, beta))
    }
  }
  return(x_0)
}

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

# initialisation des paramètres
w <- 10
h <- 20
beta <- 0.5
init <- rinit_cc(h, w)

# Rejet
x_1 <- rejet_cc(matrix(sample(-1:1, 10 * 5, replace = TRUE), 10))

# Metropolis Hastings
x_2 <- metropolis_hastings_cc(n_iter = h * w * 1000, beta, uniform = FALSE, x_0 = init)

# Gibbs
x_3 <- gibbs_cc(n_iter = 1000, beta, random = FALSE, x_0 = init)

# Propp Wilson
x_4 <- propp_wilson_cc(pas = 1000, beta, h, w)

# Affichage des realisations du rejet, de MH, Gibbs et PW
par(mfrow = c(2, 2))
image(x_1, axes = FALSE)
image(x_2, axes = FALSE)
image(x_3, axes = FALSE)
image(x_4, axes = FALSE)
