# EXEMPLE SUDOKU
rm(list = ls())
set.seed(27)

empty_cells <- function(x) {
  # indice des cases vides par carre
  empty_cells <- which(x == 0)
  col <- (empty_cells - 1) %/% 9
  row <- (empty_cells - 1) %% 9
  return(split(empty_cells, col %/% 3 * 3 + row %/% 3))
}

missing_values <- function(x) {
  # valeurs manquantes par carre
  given_cells <- which(x > 0)
  col <- (given_cells - 1) %/% 9
  row <- (given_cells - 1) %% 9
  given_cells <- split(given_cells, col %/% 3 * 3 + row %/% 3)
  return(sapply(1:9, function(i) {
    setdiff(1:9, x[given_cells[[i]]])
  }))
}

hamiltonien <- function(x) {
  # hamiltonien : nombre de replicats par lignes et colonnes
  s <- 0
  for (i in 1:9) {
    s <- s + sum(duplicated(x[i, ]) + duplicated(x[, i]))
  }
  return(s)
}

recuit_simule_sudoku <- function(s, n_iter, T_) {
  # fonction recuit simule sudoku
  empty <- empty_cells(s)
  missing <- missing_values(s)
  for (i in 1:9) {
    # remplissage de la grille
    s[empty[[i]]] <- sample(missing[[i]])
  }
  h_old <- hamiltonien(s)
  j <- sample(1:9, n_iter, replace = TRUE) # on choisit un carre uniformement
  u <- runif(n_iter)
  for (i in 1:n_iter) {
    T_ <- T_ * 0.99999
    if (h_old == 2) {
      # on remonte la temperature car min local
      T_ <- 0.5
    }
    s_new <- s
    z <- sample(empty[[j[i]]], 2) # on choisit 2 cases a echanger
    s_new[z] <- s_new[z] %*% matrix(c(0, 1, 1, 0), 2) # on echange
    h_new <- hamiltonien(s_new)
    if (u[i] <= exp((1 / T_) * (h_old - h_new))) {
      s <- s_new
      h_old <- h_new
    }
    if (h_old == 0) {
      # on a trouve la solution
      break
    }
  }
  return(s)
}

# Remplissage de la grille de sudoku
d <- matrix(0, 9, 9)
d[1, c(4, 9)] <- c(9, 6)
d[2, c(3, 6, 7, 8)] <- c(3, 8, 4, 5)
d[3, c(5, 7)] <- c(5, 9)
d[4, c(3, 5, 6, 7, 9)] <- c(1, 4, 6, 3, 2)
d[5, c(1, 9)] <- c(6, 4)
d[6, c(1, 3, 4, 5, 7)] <- c(9, 4, 8, 3, 1)
d[7, c(3, 5)] <- c(2, 8)
d[8, c(2, 3, 4, 7)] <- c(5, 6, 2, 7)
d[9, c(1, 6)] <- c(8, 3)

# initialisation des parametres
n_iter <- 5000000
T_ <- 0.5
d <- recuit_simule_sudoku(d, n_iter, T_)

# Affichage de la grille apres recuit simule
print(d)
