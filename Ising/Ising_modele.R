# Modele Ising

voisins <- function(x, i) {
  # fonction qui renvoie les voisins du site i de la matrice x
  h <- nrow(x)
  v <- rep(0, 4)
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
  return(v[v != 0])
}

hamiltonien <- function(x) {
  # hamiltonien
  return(-sum(x[-1, ] * x[-nrow(x), ]) - sum(x[, -1] * x[, -ncol(x)]))
}

prob_cond <- function(x, i, beta) {
  # calcul des probas conditionnelles
  p <- exp(beta * c(-1, 1) * sum(voisins(x, i)))
  return(p / sum(p))
}

mapping <- function(x, n_iter, j, v, u, beta) {
  for (i in 1:n_iter) {
    if (u[i] <= exp(beta * sum(voisins(x, j[i]) * (v[i] - x[j[i]])))) {
      x[j[i]] <- v[i]
    }
  }
  return(x)
}

mapping_h <- function(x, n_iter, j, v, u, beta) {
  # fonction mapping qui renvoie le hamiltonien de chaque iteration
  hamil <- numeric(n_iter)
  hamil[1] <- hamiltonien(x)
  for (i in 2:n_iter) {
    temp <- sum(voisins(x, j[i]) * (v[i] - x[j[i]]))
    if (u[i] <= exp(beta * temp)) {
      hamil[i] <- hamil[i - 1] - temp
      x[j[i]] <- v[i]
    } else {
      hamil[i] <- hamil[i - 1]
    }
  }
  return(hamil)
}

rinit_ising <- function(h, w) {
  # initialisation d'un matrice dans E
  return(matrix(sample(c(-1, 1), h * w, replace = TRUE), h))
}
