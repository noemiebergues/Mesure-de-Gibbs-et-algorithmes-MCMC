# Modele compose chimique

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
  # initialisation d'un matrice dans E
  x <- matrix(rep(0, h * w), h)
  i <- 0
  while (i < h) {
    x[i, ] <- rep(sample(-1:1, 1, replace = TRUE, prob = c(2 / 5, 1 / 5, 2 / 5)), w)
    i <- i + 2
  }
  return(x)
}
