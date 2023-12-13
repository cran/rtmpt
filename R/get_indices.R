
#' @importFrom data.table ":=" ".SD" copy fread
get_indices_x <- function(log_lik, Nchains, Nsample, df) {
  
  indices_list <- lapply(unique(df$subj), function(x) which(df$subj == x))
  ll_dt <- NULL
  if (is.character(log_lik)) {
    ll_dt <- fread(file = log_lik, skip = 0, header = FALSE)
    ll_mat <- as.matrix(copy(ll_dt))
    # ll_dt[, (paste0("s", unique(df$subj))) := lapply(indices_list, FUN = function(i) {rowMeans(ll_dt[, .SD[, i, with = FALSE]])})]
    # ll_mat2 <- as.matrix(copy(ll_dt[, .SD, .SDcols = paste0("s", unique(df$subj))]))
    my_waic <- suppressWarnings(waic(ll_mat))
    my_loo <- loo(x = ll_mat, save_psis = TRUE, cores = Nchains,
                  r_eff = relative_eff(exp(ll_mat), chain_id = rep(1:Nchains, each = Nsample), cores = Nchains))
    # my_loo2 <- loo(x = ll_mat2, save_psis = TRUE, cores = Nchains,
    #                r_eff = relative_eff(exp(ll_mat2), chain_id = rep(1:Nchains, each = Nsample), cores = Nchains))
  } else if (is.matrix(log_lik)) {
    # log_lik2 <- matrix(unlist(lapply(indices_list, FUN = function(i) {rowMeans(log_lik[, i])})), ncol = length(indices_list))
    my_waic <- suppressWarnings(waic(x = log_lik))
    my_loo <- loo(x = log_lik, save_psis = TRUE, cores = Nchains,
                  r_eff = relative_eff(exp(log_lik), chain_id = rep(1:Nchains, each = Nsample), cores = Nchains))
    # my_loo2 <- loo(x = log_lik2, save_psis = TRUE, cores = Nchains,
    #               r_eff = relative_eff(exp(log_lik2), chain_id = rep(1:Nchains, each = Nsample), cores = Nchains))
  }
  
  indices <- list()
  indices$WAIC <- my_waic
  indices$LOO <- my_loo
  # indices$LOO2 <- my_loo2
  indices$LogLik <- log_lik
  if (is.character(log_lik)) indices$LogLik <- ll_dt
  
  
  return(indices)
  
}
