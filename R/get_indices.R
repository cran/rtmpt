
get_indices <- function(log_lik, Nchains, Nsample) {
  
  if (is.character(log_lik)) {
    log_lik <- as.matrix(fread(file=log_lik,skip=0))
    my_waic <- waic(log_lik)
    my_loo <- loo(x = log_lik, r_eff = relative_eff(exp(log_lik), chain_id = rep(1:Nchains, each = Nsample), cores = Nchains), cores = Nchains)
  } else if (is.matrix(log_lik)) {
    my_waic <- waic(x = log_lik)
    my_loo <- loo(x = log_lik, r_eff = relative_eff(exp(log_lik), chain_id = rep(1:Nchains, each = Nsample), cores = Nchains), cores = Nchains)
  }
  
  indices <- list()
  indices$WAIC <- my_waic
  indices$LOO <- my_loo
  indices$LogLik <- log_lik
  
  return(indices)
  
}