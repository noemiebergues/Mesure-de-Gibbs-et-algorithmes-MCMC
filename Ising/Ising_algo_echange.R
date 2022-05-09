# Algorithme echange Ising

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
