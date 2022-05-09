# Descente de gradient Ising

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
