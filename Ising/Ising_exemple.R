# Exemple Ising

rm(list = ls())
set.seed(23)

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


metropolis_hastings_ising <- function(n_iter, beta, x_0) {
  # fonction metropolis hastings
  j <- sample(seq_len(length(x_0)), n_iter, replace = TRUE) # choix du site
  v <- sample(c(-1, 1), n_iter, replace = TRUE) # choix de la nouvelle valeur du site
  u <- runif(n_iter)
  return(mapping(x_0, n_iter, j, v, u, beta))
}


gibbs_ising <- function(n_iter, beta, random, x_0) {
  # fonction gibbs
  taille <- length(x_0)
  u <- 1:taille
  for (j in 1:n_iter) {
    if (random) {
      # random scan
      u <- sample(1:taille, taille)
    }
    for (i in u) {
      x_0[i] <- sample(c(-1, 1), 1, prob = prob_cond(x_0, i, beta))
    }
  }
  return(x_0)
}

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

metropolis_hastings_ising_h <- function(n_iter, beta, x_0) {
  # fonction metropolis hastings qui renvoie le hamiltonien de chaque iteration
  j <- sample(seq_len(length(x_0)), n_iter, replace = TRUE) # choix du site
  v <- sample(c(-1, 1), n_iter, replace = TRUE) # choix de la nouvelle valeur du site
  u <- runif(n_iter)
  return(mapping_h(x_0, n_iter, j, v, u, beta))
}

gibbs_ising_h <- function(n_iter, beta, random, x_0) {
  # fonction gibbs qui renvoie le hamiltonien de chaque iteration
  taille <- length(x_0)
  u <- 1:taille
  hamil <- numeric(n_iter)
  for (j in 1:n_iter) {
    if (random) { # random scan
      u <- sample(1:taille, taille)
    }
    for (i in u) {
      x_0[i] <- sample(c(-1, 1), 1, prob = prob_cond(x_0, i, beta))
    }
    hamil[j] <- hamiltonien(x_0)
  }
  return(hamil)
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

h_sample <- function(h, w, beta, n, n_burn, lag, sampler = "gibbs") {
  # retourne un echantillon du hamiltonien de realisations du modele
  if (sampler == "gibbs") {
    # en utilisant gibbs
    keep <- n_burn + seq(lag, (n * lag), lag)
    init <- rinit_ising(h, w)
    hamil <- gibbs_ising_h(n_burn + n * lag, beta, FALSE, init)
  } else if (sampler == "mh") {
    # en utilisant metropolis hastings
    keep <- n_burn + seq(lag, (n * lag), lag)
    init <- rinit_ising(h, w)
    hamil <- metropolis_hastings_ising_h(n_burn + n * lag, beta, init)
  } else {
    # en utilisant propp wilson
    keep <- seq(lag, (n * lag), lag)
    pas <- 1000
    hamil <- propp_wilson_ising_h(pas, n * lag, beta, h, w)
  }
  return(hamil[keep])
}

descente_gradient <- function(epsilon, b_0, x_obs, n_MC, n_burn, lag, sampler) {
  # descente de gradient
  h <- nrow(x_obs)
  w <- ncol(x_obs)
  h_obs <- hamiltonien(x_obs)
  hamil <- h_sample(h, w, b_0, n_MC, n_burn, lag, sampler)
  grad <- mean(hamil) - h_obs
  b <- c(b_0, b_0 + grad / abs(max(hamil)))
  i <- 2
  while (abs(grad) > epsilon) {
    hamil <- h_sample(h, w, b[i], n_MC, n_burn, lag, sampler)
    grad <- mean(hamil) - h_obs
    b <- c(b, b[i] + grad / (i * abs(max(hamil))))
    i <- i + 1
  }
  return(b)
}

algo_echange <- function(x_obs, n_iter, n_burn, n_MC, d, b_0 = runif(1, max = d), lag = 100, sampler = "gibbs") {
  # fonction algo echange
  h <- nrow(x_obs)
  w <- ncol(x_obs)
  h_obs <- hamiltonien(x_obs)
  b <- c(b_0, runif(n_iter - 1, 0, d))
  u <- runif(n_iter)
  for (i in 2:n_iter) {
    hamil <- h_sample(h, w, beta, n_MC, n_burn, lag, sampler)
    est_mc <- mean(exp(hamil * (b[i] - b[i - 1])))
    if (u[i] > est_mc * exp(h_obs * (b[i - 1] - b[i]))) {
      b[i] <- b[i - 1]
    }
  }
  return(b)
}

# initialisation des parametres
w <- 30
h <- 20
beta <- 0.3
d <- 1
init <- rinit_ising(h, w)

# Metropolis-Hastings
x_1 <- metropolis_hastings_ising(n_iter = h * w * 1000, beta, x_0 = init)

# Gibbs
x_2 <- gibbs_ising(n_iter = 1000, beta, random = FALSE, x_0 = init)

# Propp Wilson
x_3 <- propp_wilson_ising(pas = 1000, beta, h, w)

# Affichage des realisations de MH, Gibbs et Propp Wilson
par(mfrow = c(2, 2))
image(x_1, axes = FALSE)
image(x_2, axes = FALSE)
image(x_3, axes = FALSE)

# initialisation d'une observation
x_obs <- x_3

# Descente de gradient
b_1 <- descente_gradient(epsilon = 1, b_0 = 0.8, x_obs, n_MC = 1000, n_burn = 1000 * h * w, lag = 100, "mh")
cat("l'estimateur de Robbins Monro est ", b_1[length(b_1)], "\n")

# Affichage du chemin pris par beta lors de la descente de gradient:
par(mfrow = c(1, 1))
t_1 <- 1:length(b_1)
plot(t_1, b_1, type = "l", xlab = "nombre d'itérations", ylab = "beta")

# Algorithme d'echange
n_iter <- 5000
b_2 <- algo_echange(x_obs, n_iter, n_burn = 1000 * h * w, n_MC = 1, d, b_0 = runif(1, max = d), lag = 100, sampler = "mh")

# Trajectoire de la chaine
par(mfrow = c(1, 2))
t_2 <- 1:length(b_2)
plot(t_2, b_2, type = "l", xlab = "nombre d'itérations", ylab = "beta")

# Affichage de la loi a posteriori
b_3 <- b_2[1000:n_iter]
hist(b_3, xlab = "beta", freq = FALSE, xlim = c(0.25, 0.4), breaks = 15)
cat("l'estimateur bayésien est ", mean(b_3))
