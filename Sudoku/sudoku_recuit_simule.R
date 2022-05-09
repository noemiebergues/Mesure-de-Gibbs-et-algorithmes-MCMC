# Recuit simule sudoku

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
