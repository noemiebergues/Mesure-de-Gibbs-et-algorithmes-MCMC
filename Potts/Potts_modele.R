# Modele Potts

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
