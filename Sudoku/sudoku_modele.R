# MODELE SUDOKU

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
  # hamiltonien : nombre de réplicats par lignes et colonnes
  s <- 0
  for (i in 1:9) {
    s <- s + sum(duplicated(x[i, ]) + duplicated(x[, i]))
  }
  return(s)
}
