# Exemple Potts

rm(list = ls())
set.seed(23)

voisins <- function(x, i) {
  # fonction qui renvoie les voisins du site i de la matrice x
  h <- nrow(x)
  v <- rep(-1, 4)
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
  return(v[v != -1])
}

hamiltonien <- function(x) {
  # hamiltonien
  return(-sum(x[-1, ] == x[-nrow(x), ]) - sum(x[, -1] == x[, -ncol(x)]))
}

prob_cond <- function(x, i, beta, k) {
  # calcul des probas conditionnelles
  p <- numeric(k)
  vois <- voisins(x, i)
  for (j in 1:k) {
    p[j] <- exp(beta * sum(vois == j))
  }
  return(p / sum(p))
}

mapping <- function(x, n_iter, j, v, u, beta) {
  for (i in 1:n_iter) {
    vois <- voisins(x, j[i])
    if (u[i] <= exp(beta * sum((vois == v[i]) - (vois == x[j[i]])))) {
      x[j[i]] <- v[i]
    }
  }
  return(x)
}

rinit_potts <- function(h, w, k) {
  # initialisation d'une matrice dans E
  return(matrix(sample(1:k, h * w, replace = TRUE), h))
}

metropolis_hastings_potts <- function(n_iter, beta, k, x_0) {
  # fonction metropolis hastings
  j <- sample(seq_len(length(x_0)), n_iter, replace = TRUE) # choix du site
  v <- sample(1:k, n_iter, replace = TRUE) # choix de la nouvelle valeur du site
  u <- runif(n_iter)
  return(mapping(x_0, n_iter, j, v, u, beta))
}

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

# initialisation des parametres
w <- 30
h <- 20
beta <- 0.8
k <- 3
init <- rinit_potts(h, w, k)

# Metropolis-Hastings
x_1 <- metropolis_hastings_potts(n_iter = h * w * 1000, beta, k, x_0 = init)

# Gibbs
x_2 <- gibbs_potts(n_iter = 1000, beta, k, random = FALSE, x_0 = init)

# Affichage des realisations de MH et Gibbs
par(mfrow = c(1, 2))
image(x_1, axes = FALSE)
image(x_2, axes = FALSE)
