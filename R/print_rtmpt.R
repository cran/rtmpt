

#' @importFrom data.table as.data.table copy
StddevCorr2Cov <- function(mat, dim2, corr, sds) {
  
  dt <- as.data.table(mat)
  ind1 <- NULL
  for (i in 1:(dim2-1)) ind1 <- c(ind1, rep(i, dim2-i))
  ind2 <- NULL
  for (i in 2:dim2) ind2 <- c(ind2, i:dim2)
  cnams_corr <- colnames(dt)[corr]
  cnams_sds <- colnames(dt)[sds]
  for (i in seq_along(corr)) dt[, (cnams_corr[i]) := (get(cnams_corr[i])*(get(cnams_sds[ind1[i]]) * get(cnams_sds[ind2[i]])))]
  for (i in seq_along(sds)) dt[, (cnams_sds[i]) := get(cnams_sds[i])^2]
  return(as.matrix(copy(dt)))
  
  
  # if (!is.matrix(mat)) mat <- as.matrix(mat)
  # if ((dim2*(dim2-1)/2+dim2) != dim(mat)[2]) stop("wrong dim2 or 2nd dim of mat!")
  # 
  # matcov <- matrix(NA, ncol = dim(mat)[2], nrow = dim(mat)[1])
  # SDCorr <- matrix(0, ncol = dim2, nrow = dim2)
  # for (i in 1:dim(mat)[1]) {
  #   SDCorr[lower.tri(diag(dim2, nrow = dim2), diag = TRUE)] <- mat[i,]
  #   SDs <- diag(SDCorr)
  #   Corr <- SDCorr - diag(SDs, nrow = dim2) + diag(dim2)
  #   upperCorr <- t(Corr)-diag(diag(Corr))
  #   Corr <- Corr + upperCorr
  #   Cov <- (SDs %*% t(SDs)) * Corr
  #   matcov[i,] <- Cov[lower.tri(Cov, diag = TRUE)]
  # }
  # return(matcov)
  
}

#' @importFrom data.table as.data.table copy
Cov2StddevCorr <- function(mat, dim2, covs, vars) {
  
  dt <- as.data.table(mat)
  ind1 <- NULL
  for (i in 1:(dim2-1)) ind1 <- c(ind1, rep(i, dim2-i))
  ind2 <- NULL
  for (i in 2:dim2) ind2 <- c(ind2, i:dim2)
  cnams_covs <- colnames(dt)[covs]
  cnams_vars <- colnames(dt)[vars]
  for (i in seq_along(covs)) dt[, (cnams_covs[i]) := (get(cnams_covs[i])/sqrt(get(cnams_vars[ind1[i]]) * get(cnams_vars[ind2[i]])))]
  for (i in seq_along(vars)) dt[, (cnams_vars[i]) := sqrt(get(cnams_vars[i]))]
  return(as.matrix(copy(dt)))
  
  # if (!is.matrix(mat)) mat <- as.matrix(mat)
  # if ((dim2*(dim2-1)/2+dim2) != dim(mat)[2]) stop("wrong dim2 or 2nd dim of mat!")
  # 
  # matcorr <- matrix(NA, ncol = dim(mat)[2], nrow = dim(mat)[1])
  # CoVar <- matrix(0, ncol = dim2, nrow = dim2)
  # for (i in 1:dim(mat)[1]) {
  #   CoVar[lower.tri(diag(dim2, nrow = dim2), diag = TRUE)] <- mat[i,]
  #   tCoVar <- t(CoVar)
  #   diag(tCoVar) <- rep(0, dim2)
  #   CoVar <- CoVar + tCoVar
  #   SDs <- sqrt(diag(CoVar))
  #   tmp <- diag(dim2)*1/sqrt(diag(CoVar))
  #   Corr <- tmp %*% CoVar %*% tmp
  #   diag(Corr) <- SDs
  #   matcorr[i,] <- Corr[lower.tri(Corr, diag = TRUE)]
  # }
  # return(matcorr)
  
}

#' @importFrom stats dt
int_dt <- function(x, mu, sd, defr) {
  x * dt((x-mu)/sd, df = defr) / sd
}
#' @importFrom stats integrate
integrate_wrapper <- function(mu, sd, defr) {
  integrate(int_dt, 0, Inf, mu=mu, sd=sd, defr=defr)$value
}
vectorized_integrate <- Vectorize(integrate_wrapper)

#' @importFrom stats pt
#' @importFrom data.table as.data.table
expect_t = function(location, scale, df){
  
  dt_loc <- as.data.table(location)
  cnams_loc <- colnames(copy(dt_loc))
  dt_loc[, ("scale") := exp(scale)]
  dt_loc[, ("df") := df]
  
  for (cn in cnams_loc) {
    
    dt_loc[, (cn) := (vectorized_integrate(mu = get(cn), sd = scale, defr = df) / pt(-get(cn)/scale, df=df, lower.tail=FALSE))]
    
    # dt_loc[, (cn) := (stats::integrate(int_dt, 0, Inf, mu=get(cn), sd=scale, defr=df)$value / stats::pt(-get(cn)/scale, df=df, lower.tail=F))]
    # stats::integrate(int_dt, 0, Inf, mu=shift, sd=scale, defr=df)$value / stats::pt(-shift/scale, df=df, lower.tail=F)
  }
  
  if (length(cnams_loc) == 1) {
    ret <- dt_loc[, get(cnams_loc)]
  } else ret <- as.matrix(copy(dt_loc[, mget(cnams_loc)]))
  
  return(ret)
}


#' @importFrom coda effectiveSize gelman.diag
#' @importFrom stats dnorm pnorm
writeSummaryERTMPT <- function(x, keep, ...) {

  prob_names <- colnames(x$specs$model$params$probs[1, which(is.na(x$specs$model$params$probs[1,]))])
  tauM_names <- colnames(x$specs$model$params$taus[1, which(is.na(x$specs$model$params$taus[1,]))])
  tauP_names <- colnames(x$specs$model$params$taus[2, which(is.na(x$specs$model$params$taus[2,]))])
  tau_names <- c(tauM_names, tauP_names)
  Nprobs <- sum(is.na(x$specs$model$params$probs))
  Ntau_m <- sum(is.na(x$specs$model$params$taus[1,]))
  Ntau_p <- sum(is.na(x$specs$model$params$taus[2,]))
  Ntaus <- Ntau_m + Ntau_p
  Nparams <- Nprobs+Ntaus
  Nresps <- length(unique(x$specs$model$responses$MAP))
  Ngroups <- x$specs$n.groups
  Nsubj <- x$specs$n.subj
  subj_exist <- exists("subj", x$specs$transformation)
  subj_labels <- ifelse((keep & subj_exist), x$specs$transformation$group$new, 0:(Nsubj-1))
  group_exist <- exists("group", x$specs$transformation)
  group_labels <- ifelse((keep & group_exist), x$specs$transformation$subj$new, 0:(Ngroups-1))
  Names <- c("Mean", "SD", "2.5%", "50%", "97.5%", "Naive SE", "Time-series SE", "n.eff", "Rhat", "R_95%")
  
  samp <- x$samples
  
  # omega squared
  ind_omega2 <- Nparams*Ngroups + Nparams*(Nparams+1)/2 + Nparams*Nsubj + Nresps*Ngroups + 1
  
  # probabilities
  ind_probs <- 1:(Nprobs*Ngroups)
  
  # process times
  ind_taus <- (Nprobs*Ngroups+1):(Nprobs*Ngroups+Ntaus*Ngroups)
  
  # SIGMA
  ind_SIGMA <- (Nparams*Ngroups+1):(Nparams*Ngroups+Nparams*(Nparams+1)/2)
  mmm <- matrix(NA, ncol = 2, nrow = Nparams*(Nparams+1)/2)
  count=1
  for (i in 1:Nparams) for(j in i:Nparams) {mmm[count,] <- c(i,j); count = count + 1}
  corrs <- which(mmm[,1]!=mmm[,2])
  sds <- which(mmm[,1]==mmm[,2])
  
  # motor times
  aa <- Nparams*Ngroups + Nparams*(Nparams+1)/2 + Nparams*Nsubj + 1
  ind_resps <- aa:(aa+Nresps*Ngroups-1)
  
  # GAMMA
  ind_GAMMA <- (aa+Nresps*Ngroups+1):(aa+Nresps*Ngroups+1+Nresps*(Nresps+1)/2-1)
  mmmR <- matrix(NA, ncol = 2, nrow = Nresps*(Nresps+1)/2)
  count=1
  for (i in 1:Nresps) for(j in i:Nresps) {mmmR[count,] <- c(i,j); count = count + 1}
  corrsR <- which(mmmR[,1]!=mmmR[,2])
  sdsR <- which(mmmR[,1]==mmmR[,2])
  
  # sigma squared
  bb <- aa+Nresps*Ngroups+1+Nresps*(Nresps+1)/2+Nresps*Nsubj
  ind_sigma2 <- bb:(bb+Nsubj-1)
  
  # transform samp
  delta <- samp[, ind_resps]
  temp <- list()
  for (n in 1:x$specs$n.chains) {
    for (s in 1:Nsubj) {
      temp[[s]] <- samp[[n]][, ind_resps] + dnorm(-samp[[n]][, ind_resps]/sqrt(samp[[n]][, ind_sigma2[s]])) / 
        (1-pnorm(-samp[[n]][, ind_resps]/sqrt(samp[[n]][, ind_sigma2[s]]))) * sqrt(samp[[n]][, ind_sigma2[s]])
      if (s == 1) {delta[[n]] <- temp[[s]]} else delta[[n]] <- delta[[n]] + temp[[s]]
    }
    delta[[n]] <- delta[[n]] / Nsubj * 1000
    samp[[n]][, ind_probs] <- pnorm(samp[[n]][, ind_probs])
    samp[[n]][, ind_taus] <- 1000/(samp[[n]][, ind_taus])
    samp[[n]][, ind_SIGMA] <- StddevCorr2Cov(samp[[n]][, ind_SIGMA], Nparams, corrs, sds)
    samp[[n]][, ind_resps] <- 1000*(samp[[n]][, ind_resps])
    samp[[n]][, ind_GAMMA] <- StddevCorr2Cov(samp[[n]][, ind_GAMMA], Nresps, corrsR, sdsR)
  }
  
  # summary stats
  su <- summary(samp)
  su_orig <- summary(x$samples)
  sudelta <- summary(delta)
  n_eff <- effectiveSize(samp)
  n_eff_orig <- effectiveSize(x$samples)
  n_eff_delta <- effectiveSize(delta)
  R_hat <- gelman.diag(samp, multivariate = FALSE)
  R_hat_delta <- gelman.diag(delta, multivariate = FALSE)
  
  
  # output
  
  omega2 <- as.numeric(c(su$statistics[ind_omega2, c(1,2)], su$quantiles[ind_omega2, c(1,3,5)], 
                         su$statistics[ind_omega2, c(3,4)], n_eff[ind_omega2], R_hat$psrf[ind_omega2,]))
  omega2_mat <- matrix(omega2, ncol = length(Names))
  colnames(omega2_mat) <- Names
  rownames(omega2_mat) <- "omega_squared"
  
  ind_main <- c(ind_probs, ind_taus)
  main_mat <- cbind(su$statistics[ind_main, c(1,2)], su$quantiles[ind_main, c(1,3,5)], 
                    su$statistics[ind_main, c(3,4)], n_eff[ind_main], R_hat$psrf[ind_main,])
  if (Nresps > 1 | Ngroups > 1) {
	main_mat <- rbind(main_mat, cbind(sudelta$statistics[, c(1,2)], sudelta$quantiles[, c(1,3,5)], 
                                    sudelta$statistics[, c(3,4)], n_eff_delta, R_hat_delta$psrf))
  } else {
    main_mat <- rbind(main_mat, matrix(c(sudelta$statistics[c(1,2)], sudelta$quantiles[c(1,3,5)], 
	                                     sudelta$statistics[c(3,4)], n_eff_delta, R_hat_delta$psrf), nrow = 1))
  }
  colnames(main_mat) <- Names
  rnams <- c("theta_", "E(tau_minus_", "E(tau_plus_", "E(delta_")
  if (Ngroups == 1) {
    theta_names <- paste0(rnams[1], prob_names)
    tau_M_names <- paste0(rnams[2], tauM_names, ")")
    tau_P_names <- paste0(rnams[3], tauP_names, ")")
    delta_names <- paste0(rnams[4], rep(paste0("R", 0:(Nresps-1))), ")")
  } else {
    theta_names <- paste0(rnams[1], prob_names)
    tau_M_names <- paste0(rnams[2], tauM_names)
    tau_P_names <- paste0(rnams[3], tauP_names)
    delta_names <- paste0(rnams[4], rep(paste0("R", 0:(Nresps-1))))
    theta_names <- paste0(rep(theta_names, Ngroups), rep(paste0("[", group_labels, "]"), each=Nprobs))
    tau_M_names <- paste0(rep(tau_M_names, Ngroups), rep(paste0("[", group_labels, "]"), each=Ntau_m), ")")
    tau_P_names <- paste0(rep(tau_P_names, Ngroups), rep(paste0("[", group_labels, "]"), each=Ntau_p), ")")
    delta_names <- paste0(rep(delta_names, Ngroups), rep(paste0("[", group_labels, "]"), each=Nresps), ")")
  }
  rownames(main_mat) <- c(theta_names, tau_M_names, tau_P_names, delta_names)

  
  ind_main <- c(ind_main, ind_resps)
  orig_mat <- cbind(su_orig$statistics[ind_main, c(1,2)], su_orig$quantiles[ind_main, c(1,3,5)], 
                    su_orig$statistics[ind_main, c(3,4)], n_eff_orig[ind_main], 
                    x$diags$R_hat$psrf[ind_main,])
  colnames(orig_mat) <- Names

  
  CorrP_mat <- cbind(su_orig$statistics[ind_SIGMA[corrs], c(1,2)], su_orig$quantiles[ind_SIGMA[corrs], c(1,3,5)], 
                     su_orig$statistics[ind_SIGMA[corrs], c(3,4)], n_eff_orig[ind_SIGMA[corrs]], 
                     x$diags$R_hat$psrf[ind_SIGMA[corrs],])
  colnames(CorrP_mat) <- Names
  
  if (length(corrsR)==1) {
    
    CorrR <- c(su_orig$statistics[ind_GAMMA[corrsR], c(1,2)], su_orig$quantiles[ind_GAMMA[corrsR], c(1,3,5)], 
               su_orig$statistics[ind_GAMMA[corrsR], c(3,4)], n_eff_orig[ind_GAMMA[corrsR]], 
               x$diags$R_hat$psrf[ind_GAMMA[corrsR],])
    CorrR_mat <- matrix(CorrR, ncol = length(Names))
    colnames(CorrR_mat) <- Names
    rownames(CorrR_mat) <- names(n_eff_orig[ind_GAMMA[corrsR]])
  } else if (length(corrsR>1)) {
    CorrR_mat <- cbind(su_orig$statistics[ind_GAMMA[corrsR], c(1,2)], su_orig$quantiles[ind_GAMMA[corrsR], c(1,3,5)], 
                       su_orig$statistics[ind_GAMMA[corrsR], c(3,4)], n_eff_orig[ind_GAMMA[corrsR]], 
                       x$diags$R_hat$psrf[ind_GAMMA[corrsR],])
    colnames(CorrR_mat) <- Names
  }
  if (!exists("CorrR_mat")) CorrR_mat <- NULL
  
  sdsP_mat <- cbind(su_orig$statistics[ind_SIGMA[sds], c(1,2)], su_orig$quantiles[ind_SIGMA[sds], c(1,3,5)], 
                    su_orig$statistics[ind_SIGMA[sds], c(3,4)], n_eff_orig[ind_SIGMA[sds]], 
                    x$diags$R_hat$psrf[ind_SIGMA[sds],])
  colnames(sdsP_mat) <- Names
  
  if (length(sdsR)==1) {
    
    sds_R <- c(su_orig$statistics[ind_GAMMA[sdsR], c(1,2)], su_orig$quantiles[ind_GAMMA[sdsR], c(1,3,5)], 
               su_orig$statistics[ind_GAMMA[sdsR], c(3,4)], n_eff_orig[ind_GAMMA[sdsR]], 
               x$diags$R_hat$psrf[ind_GAMMA[sdsR],])
    sdsR_mat <- matrix(sds_R, ncol = length(Names))
    colnames(sdsR_mat) <- Names
    rownames(sdsR_mat) <- names(n_eff_orig[ind_GAMMA[sdsR]])
  } else if (length(sdsR>1)) {
    sdsR_mat <- cbind(su_orig$statistics[ind_GAMMA[sdsR], c(1,2)], su_orig$quantiles[ind_GAMMA[sdsR], c(1,3,5)], 
                      su_orig$statistics[ind_GAMMA[sdsR], c(3,4)], n_eff_orig[ind_GAMMA[sdsR]], 
                      x$diags$R_hat$psrf[ind_GAMMA[sdsR],])
    colnames(sdsR_mat) <- Names
  }

  
  CovP_mat <- cbind(su$statistics[ind_SIGMA[corrs], c(1,2)], su$quantiles[ind_SIGMA[corrs], c(1,3,5)], 
                    su$statistics[ind_SIGMA[corrs], c(3,4)], n_eff[ind_SIGMA[corrs]], 
                    R_hat$psrf[ind_SIGMA[corrs],])
  colnames(CovP_mat) <- Names
  
  if (length(corrsR)==1) {
    
    CovR <- c(su$statistics[ind_GAMMA[corrsR], c(1,2)], su$quantiles[ind_GAMMA[corrsR], c(1,3,5)], 
              su$statistics[ind_GAMMA[corrsR], c(3,4)], n_eff[ind_GAMMA[corrsR]], 
              R_hat$psrf[ind_GAMMA[corrsR],])
    CovR_mat <- matrix(CovR, ncol = length(Names))
    colnames(CovR_mat) <- Names
    rownames(CovR_mat) <- names(n_eff_orig[ind_GAMMA[corrsR]])
  } else if (length(corrsR>1)) {
    CovR_mat <- cbind(su$statistics[ind_GAMMA[corrsR], c(1,2)], su$quantiles[ind_GAMMA[corrsR], c(1,3,5)], 
                      su$statistics[ind_GAMMA[corrsR], c(3,4)], n_eff[ind_GAMMA[corrsR]], 
                      R_hat$psrf[ind_GAMMA[corrsR],])
    colnames(CovR_mat) <- Names
  }
  if (!exists("CovR_mat")) CovR_mat <- NULL
  
  VarP_mat <- cbind(su$statistics[ind_SIGMA[sds], c(1,2)], su$quantiles[ind_SIGMA[sds], c(1,3,5)], 
                    su$statistics[ind_SIGMA[sds], c(3,4)], n_eff[ind_SIGMA[sds]], 
                    R_hat$psrf[ind_SIGMA[sds],])
  colnames(VarP_mat) <- Names
  
  if (length(sdsR)==1) {
    
    VarR <- c(su$statistics[ind_GAMMA[sdsR], c(1,2)], su$quantiles[ind_GAMMA[sdsR], c(1,3,5)], 
              su$statistics[ind_GAMMA[sdsR], c(3,4)], n_eff[ind_GAMMA[sdsR]], 
              R_hat$psrf[ind_GAMMA[sdsR],])
    VarR_mat <- matrix(VarR, ncol = length(Names))
    colnames(VarR_mat) <- Names
    rownames(VarR_mat) <- names(n_eff_orig[ind_GAMMA[sdsR]])
  } else if (length(sdsR>1)) {
    VarR_mat <- cbind(su$statistics[ind_GAMMA[sdsR], c(1,2)], su$quantiles[ind_GAMMA[sdsR], c(1,3,5)], 
                      su$statistics[ind_GAMMA[sdsR], c(3,4)], n_eff[ind_GAMMA[sdsR]], 
                      R_hat$psrf[ind_GAMMA[sdsR],])
    colnames(VarR_mat) <- Names
  }
  
  elementsP <- su_orig$statistics[ind_SIGMA, 1]
  CORR_SD_P <- matrix(0, ncol = Nparams, nrow = Nparams)
  CORR_SD_P[lower.tri(CORR_SD_P, diag = TRUE)] <- elementsP
  CORR_SD_P <- CORR_SD_P + t(CORR_SD_P) - diag(elementsP[sds])

  elementsR <- su_orig$statistics[ind_GAMMA, 1]
  if (length(elementsR)==1) {
    CORR_SD_R <- elementsR
  } else {
    CORR_SD_R <- matrix(0, ncol = Nresps, nrow = Nresps)
    CORR_SD_R[lower.tri(CORR_SD_R, diag = TRUE)] <- elementsR
    CORR_SD_R <- CORR_SD_R + t(CORR_SD_R) - diag(elementsR[sdsR])
  }
  
  elements_covP <- su$statistics[ind_SIGMA, 1]
  COV_VAR_P <- matrix(0, ncol = Nparams, nrow = Nparams)
  COV_VAR_P[lower.tri(COV_VAR_P, diag = TRUE)] <- elements_covP
  COV_VAR_P <- COV_VAR_P + t(COV_VAR_P) - diag(elements_covP[sds])
  
  
  elements_covR <- su$statistics[ind_GAMMA, 1]
  if (length(elements_covR)==1) {
    COV_VAR_R <- elements_covR
  } else {
    COV_VAR_R <- matrix(0, ncol = Nresps, nrow = Nresps)
    COV_VAR_R[lower.tri(COV_VAR_R, diag = TRUE)] <- elements_covR
    COV_VAR_R <- COV_VAR_R + t(COV_VAR_R) - diag(elements_covR[sdsR])
  }
  
  summary_list <- list(resid_var = omega2_mat, 
                       transformed_pars = main_mat,
                       orig_pars = orig_mat, 
                       Corrs = list(process=CorrP_mat, motor=CorrR_mat), 
                       SDs = list(process=sdsP_mat, motor=sdsR_mat),
                       Covs = list(process=CovP_mat, motor=CovR_mat), 
                       Vars = list(process=VarP_mat, motor=VarR_mat), 
                       CorrSD = list(process=CORR_SD_P, motor=CORR_SD_R), 
                       CovVar = list(process=COV_VAR_P, motor=COV_VAR_R))
  return(summary_list)
}







get_mcmc_process <- function(list_param, ip, ipar, nam, sampm, mcmc_names, n.iter) {
  if (is.na(list_param[[ipar]][1, ip])) {
    ind_prob <- which(mcmc_names == paste0(nam, names(list_param[[ipar]][ip])))
    return(sampm[, ind_prob])
  } else {
    if (is.numeric(list_param[[ipar]][1, ip])) return(rep(list_param[[ipar]][1, ip], n.iter))
    if (is.character(list_param[[ipar]][1, ip])) {
      ind <- which(names(list_param[[ipar]])==list_param[[ipar]][1, ip])
      if (is.na(list_param[[ipar]][1, ind])) {
        ind_prob <- which(mcmc_names == paste0(nam, names(list_param[[ipar]][ind])))
        return(sampm[, ind_prob])
      } else {
        return(rep(list_param[[ipar]][1, ind], n.iter))
      }
    }
  }
}

#' @importFrom data.table fifelse
log1pem1 <- function(z) {
  fifelse(abs(z) < 1.0e-2, log(-expm1(z)), log1p(-exp(z)))
}

#' @importFrom data.table fifelse
logdiff <- function(xa, xb) {
  fifelse(xb >= xa, -Inf, fifelse(xb <= -Inf, xa, xa + log1pem1(xb - xa)))
}

prob_lower <- function(a, v, w) {
  len <- length(a)
  vm <- fifelse(v >= 0, v, -v)
  wm <- fifelse(v >= 0, w, 1-w)
  tmp <- exp( logdiff(rep(0, len), -2*vm*a*(1-wm)) - logdiff(2*vm*a*wm, -2*vm*a*(1-wm)) )
  if (any(tmp < 0)) stop("bla")
  return(fifelse(v >= rep(0, len), tmp, 1-tmp))
}

#' @importFrom data.table fifelse
expect_time_lower <- function(a, v, w) {
  temp <- NULL
  amw <- a * (1-w)
  return(fifelse(abs(v) < rep(1e-5, length(v)), (a^2 - amw^2) / 3, (a / tanh(a*v) - amw / tanh(v*amw)) / v))
}

prep_transform <- function(a, b, loc_o, scale) {
  range <- b - a
  loc <- (loc_o-a) / range
  loc <- log(loc / (1-loc))
  scale <- (loc_o + scale - a) / range
  scale <- log(scale / (1-scale))
  scale <- scale - loc
  return(list(a = a, b = b, location = loc, scale = scale, range = range))
}

#' @importFrom data.table fifelse
inv_scale_logit <- function(ttheta, a, b, location, scale, range) {
  dims <- dim(ttheta)
  z <- scale * ttheta + location
  return(matrix(fifelse(z < -700, a, a + range / (1+exp(-z))), nrow = dims[1], ncol = dims[2]))
}

scale_logit <- function(theta, a, b, location, scale, range) {
  tmp <- (theta - a) / range
  tmp <- log(tmp / (1-tmp))
  return((tmp-location)/scale)
}

#' @importFrom coda effectiveSize gelman.diag
#' @importFrom stats dnorm pnorm
writeSummaryDRTMPT <- function(x, keep, ...) {
  
  summary_list <- list()

  
  # define some names
  thresh_names <- colnames(x$specs$model$params$threshold[1, which(is.na(x$specs$model$params$threshold[1,]))])
  drift_names <- colnames(x$specs$model$params$driftrate[1, which(is.na(x$specs$model$params$driftrate[1,]))])
  start_names <- colnames(x$specs$model$params$startpoint[1, which(is.na(x$specs$model$params$startpoint[1,]))])

  
  # define numbers
  Nthresh <- sum(is.na(x$specs$model$params$threshold[1,]))
  Ndrift <- sum(is.na(x$specs$model$params$driftrate[1,]))
  Nstart <- sum(is.na(x$specs$model$params$startpoint[1,]))
  Nparams <- Nthresh + Ndrift + Nstart
  Nresps <- length(unique(x$specs$model$responses$MAP))
  Ngroups <- x$specs$n.groups
  Nsubj <- x$specs$n.subj
  
  
  # define flags and labels
  subj_exist <- exists("subj", x$specs$transformation)
  subj_labels <- ifelse((keep & subj_exist), x$specs$transformation$group$new, 0:(Nsubj-1))
  group_exist <- exists("group", x$specs$transformation)
  group_labels <- ifelse((keep & group_exist), x$specs$transformation$subj$new, 0:(Ngroups-1))
  Names <- c("Mean", "SD", "2.5%", "50%", "97.5%", "Naive SE", "Time-series SE", "n.eff", "Rhat", "R_95%")
  
  
  # copy samples for transformations
  samp <- x$samples
  
  
  # mcmc colnames
  mcmc_names <- colnames(samp[[1]])
  
  
  # GET INDICES
  ## omega squared
  ind_Omega2 <- which(mcmc_names == "Omega2") #Nparams*Ngroups + Nparams*Nsubj + Nresps*Ngroups + Nresps*Nsubj + Nsubj + Nparams*(Nparams+1)/2 + Nresps*(Nresps+1)/2 + 1
  
  ind_thresh <- ind_drift <- ind_start <- NULL
  ## thresholds
  if (Nthresh) ind_thresh <- which(grepl(pattern = "mu_a", x = mcmc_names))
  ## drift rate
  if (Ndrift) ind_drift <- which(grepl(pattern = "mu_nu", x = mcmc_names))
  ## drift rate
  if (Nstart) ind_start <- which(grepl(pattern = "mu_omega", x = mcmc_names))
  
  ## SIGMA
  ind_SIGMA <- which(grepl(pattern = "s_", x = mcmc_names)) #(Nparams*Ngroups+1):(Nparams*Ngroups+Nparams*(Nparams+1)/2)
  mmm <- matrix(NA, ncol = 2, nrow = Nparams*(Nparams+1)/2)
  count=1
  for (i in 1:Nparams) for(j in i:Nparams) {mmm[count,] <- c(i,j); count = count + 1}
  covs <- which(mmm[,1]!=mmm[,2])
  vars <- which(mmm[,1]==mmm[,2])
  
  ## motor times
  ind_resps <- which(grepl(pattern = "mu_gamma_R", x = mcmc_names)) # aa:(aa+Nresps*Ngroups-1)
  
  ## GAMMA
  ind_GAMMA <- which(grepl(pattern = "g_gamma", x = mcmc_names)) # (aa+Nresps*Ngroups+1):(aa+Nresps*Ngroups+1+Nresps*(Nresps+1)/2-1)
  mmmR <- matrix(NA, ncol = 2, nrow = Nresps*(Nresps+1)/2)
  count=1
  for (i in 1:Nresps) for(j in i:Nresps) {mmmR[count,] <- c(i,j); count = count + 1}
  covsR <- which(mmmR[,1]!=mmmR[,2])
  varsR <- which(mmmR[,1]==mmmR[,2])
  
  ## sigma squared
  ind_sigma <- which(grepl(pattern = "sigma", x = mcmc_names)) # bb:(bb+Nsubj-1)
  
  # transform samp
  ## prepare transformation of a v w
  logit_args <- prep_transform(a = c(.01, -100, .001), b = c(100, 100, .999), loc_o = c(.8, 0, .5), scale = c(.2, 1, .1))
  list_param <- x$specs$model$params
  Nprobs_all <- length(list_param$threshold)
  tmp_thresh <- tmp_drift <- tmp_start <- samp[, 1:(Nprobs_all*Ngroups)]
  NisNA <- numeric(Nprobs_all)
  for (ip in 1:Nprobs_all) {
    NisNA[ip] <- sum(c(is.na(list_param$threshold[ip]), is.na(list_param$driftrate[ip]), is.na(list_param$startpoint[ip])))
  }
  Nprobs_est <- sum(NisNA>0)
  theta <- tau_m <- tau_p <- samp[, 1:(Nprobs_est*Ngroups)]
  delta <- samp[, ind_resps]
  temp <- list()
  for (n in 1:x$specs$n.chains) {
    for (s in 1:Nsubj) {
      temp[[s]] <- expect_t(location = x$samples[[n]][, ind_resps], 
                            scale = samp[[n]][, ind_sigma[s]], 
                            df = x$specs$prior_params$delta_df)
      if (s == 1) {delta[[n]] <- temp[[s]]} else delta[[n]] <- delta[[n]] + temp[[s]]
    }
    delta[[n]] <- delta[[n]] / Nsubj * 1000
    
    if (Nthresh) samp[[n]][, ind_thresh] <- inv_scale_logit(ttheta = as.matrix(x$samples[[n]][, ind_thresh]), a = logit_args$a[1], b = logit_args$b[1], location = logit_args$location[1], 
                                                             scale = logit_args$scale[1], range = logit_args$range[1])
    if (Ndrift) samp[[n]][, ind_drift] <- inv_scale_logit(ttheta = as.matrix(x$samples[[n]][, ind_drift]), a = logit_args$a[2], b = logit_args$b[2], location = logit_args$location[2], 
                                                           scale = logit_args$scale[2], range = logit_args$range[2])
    if (Nstart) samp[[n]][, ind_start] <- inv_scale_logit(ttheta = as.matrix(x$samples[[n]][, ind_start]), a = logit_args$a[3], b = logit_args$b[3], location = logit_args$location[3], 
                                                           scale = logit_args$scale[3], range = logit_args$range[3])
    jp <- 0
    for (ip in 1:Nprobs_all) {
      if (NisNA[ip] > 0) {
        tmp_thresh[[n]][, ip] <- get_mcmc_process(list_param, ip, 1, "mu_a_", as.matrix(samp[[n]][, 1:Nparams]), mcmc_names, x$specs$n.iter)
        tmp_drift[[n]][, ip] <- get_mcmc_process(list_param, ip, 2, "mu_nu_", as.matrix(samp[[n]][, 1:Nparams]), mcmc_names, x$specs$n.iter)
        tmp_start[[n]][, ip] <- get_mcmc_process(list_param, ip, 3, "mu_omega_", as.matrix(samp[[n]][, 1:Nparams]), mcmc_names, x$specs$n.iter)
        jp <- jp + 1
        theta[[n]][, jp] <- prob_lower(as.numeric(tmp_thresh[[n]][, ip]), 
                                       -as.numeric(tmp_drift[[n]][, ip]), 
                                       1-as.numeric(tmp_start[[n]][, ip]))
        tau_m[[n]][, jp] <- expect_time_lower(as.numeric(tmp_thresh[[n]][, ip]), 
                                              as.numeric(tmp_drift[[n]][, ip]), 
                                              as.numeric(tmp_start[[n]][, ip]))
        tau_p[[n]][, jp] <- expect_time_lower(as.numeric(tmp_thresh[[n]][, ip]), 
                                              -as.numeric(tmp_drift[[n]][, ip]), 
                                              1-as.numeric(tmp_start[[n]][, ip]))
      }
    }
    samp[[n]][, ind_SIGMA] <- Cov2StddevCorr(x$samples[[n]][, ind_SIGMA], Nparams, covs = covs, vars = vars)
    samp[[n]][, ind_resps] <- 1000*(x$samples[[n]][, ind_resps])
    samp[[n]][, ind_GAMMA] <- Cov2StddevCorr(x$samples[[n]][, ind_GAMMA], Nresps, covs = covsR, vars = varsR)
  }
  
  
  # summary stats
  su <- summary(samp)
  su_orig <- summary(x$samples)
  sudelta <- summary(delta)
  sutheta <- summary(theta)
  sutau_m <- summary(tau_m)
  sutau_p <- summary(tau_p)
  n_eff <- effectiveSize(samp)
  n_eff_orig <- effectiveSize(x$samples)
  n_eff_delta <- effectiveSize(delta)
  n_eff_theta <- effectiveSize(theta)
  n_eff_tau_m <- effectiveSize(tau_m)
  n_eff_tau_p <- effectiveSize(tau_p)
  R_hat <- gelman.diag(samp, multivariate = FALSE)
  R_hat_delta <- gelman.diag(delta, multivariate = FALSE)
  R_hat_theta <- gelman.diag(theta, multivariate = FALSE)
  R_hat_tau_m <- gelman.diag(tau_m, multivariate = FALSE)
  R_hat_tau_p <- gelman.diag(tau_p, multivariate = FALSE)
  
  
  # output
  ## residual variance
  Omega2 <- as.numeric(c(su$statistics[ind_Omega2, c(1,2)], su$quantiles[ind_Omega2, c(1,3,5)], 
                         su$statistics[ind_Omega2, c(3,4)], n_eff[ind_Omega2], R_hat$psrf[ind_Omega2,]))
  Omega2_mat <- matrix(Omega2, ncol = length(Names))
  colnames(Omega2_mat) <- Names
  rownames(Omega2_mat) <- "Omega_squared"
  summary_list$resid_var <- Omega2_mat
  
  
  ## main transformed params
  ind_main <- c(ind_thresh, ind_drift, ind_start)
  main_mat <- cbind(su$statistics[ind_main, c(1,2)], su$quantiles[ind_main, c(1,3,5)],
                    su$statistics[ind_main, c(3,4)], n_eff[ind_main], R_hat$psrf[ind_main,])
  main_mat <- rbind(main_mat, cbind(sutheta$statistics[, c(1,2)], sutheta$quantiles[, c(1,3,5)], 
                                    sutheta$statistics[, c(3,4)], n_eff_theta, R_hat_theta$psrf))
  main_mat <- rbind(main_mat, cbind(sutau_m$statistics[, c(1,2)], sutau_m$quantiles[, c(1,3,5)], 
                                    sutau_m$statistics[, c(3,4)], n_eff_tau_m, R_hat_tau_m$psrf))
  main_mat <- rbind(main_mat, cbind(sutau_p$statistics[, c(1,2)], sutau_p$quantiles[, c(1,3,5)], 
                                    sutau_p$statistics[, c(3,4)], n_eff_tau_p, R_hat_tau_p$psrf))
  if (Nresps > 1 | Ngroups > 1) {
    main_mat <- rbind(main_mat, cbind(sudelta$statistics[, c(1,2)], sudelta$quantiles[, c(1,3,5)],
                                      sudelta$statistics[, c(3,4)], n_eff_delta, R_hat_delta$psrf))
  } else {
    main_mat <- rbind(main_mat, matrix(c(sudelta$statistics[c(1,2)], sudelta$quantiles[c(1,3,5)],
                                         sudelta$statistics[c(3,4)], n_eff_delta, R_hat_delta$psrf), nrow = 1))
  }
  colnames(main_mat) <- Names
  rnams <- c("mean_a_", "mean_nu_", "mean_omega_", "E(theta_", "E(tau_minus_", "E(tau_plus", "E(delta_")
  a_names <- nu_names <- omega_names <- NULL
  prob_names <- names(x$specs$model$params$threshold)[which(NisNA>0)]
  if (Ngroups == 1) {
    if (Nthresh) a_names <- paste0(rnams[1], thresh_names)
    if (Ndrift) nu_names <- paste0(rnams[2], drift_names)
    if (Nstart) omega_names <- paste0(rnams[3], start_names)
    theta_names <- paste0(rnams[4], prob_names, ")")
    tau_M_names <- paste0(rnams[5], prob_names, ")")
    tau_P_names <- paste0(rnams[6], prob_names, ")")
    delta_names <- paste0(rnams[7], rep(paste0("R", 0:(Nresps-1))), ")")
  } else {
    if (Nthresh) a_names <- paste0(rnams[1], thresh_names)
    if (Ndrift) nu_names <- paste0(rnams[2], drift_names)
    if (Nstart) omega_names <- paste0(rnams[3], start_names)
    theta_names <- paste0(rnams[4], prob_names)
    tau_M_names <- paste0(rnams[5], prob_names)
    tau_P_names <- paste0(rnams[6], prob_names)
    delta_names <- paste0(rnams[7], rep(paste0("R", 0:(Nresps-1))))
    if (Nthresh) a_names <- paste0(rep(a_names, Ngroups), rep(paste0("[", group_labels, "]"), each=Nthresh))
    if (Ndrift) nu_names <- paste0(rep(nu_names, Ngroups), rep(paste0("[", group_labels, "]"), each=Ndrift))
    if (Nstart) omega_names <- paste0(rep(omega_names, Ngroups), rep(paste0("[", group_labels, "]"), each=Nstart))
    theta_names <- paste0(rep(theta_names, Ngroups), rep(paste0("[", group_labels, "]"), each=Nprobs_est), ")")
    tau_M_names <- paste0(rep(tau_M_names, Ngroups), rep(paste0("[", group_labels, "]"), each=Nprobs_est), ")")
    tau_P_names <- paste0(rep(tau_P_names, Ngroups), rep(paste0("[", group_labels, "]"), each=Nprobs_est), ")")
    delta_names <- paste0(rep(delta_names, Ngroups), rep(paste0("[", group_labels, "]"), each=Nresps), ")")
  }
  rownames(main_mat) <- c(a_names, nu_names, omega_names, theta_names, tau_M_names, tau_P_names, delta_names)
  summary_list$transformed_pars <- main_mat
  
  
  ## main original params
  ind_main <- c(ind_main, ind_resps)
  orig_mat <- cbind(su_orig$statistics[ind_main, c(1,2)], su_orig$quantiles[ind_main, c(1,3,5)], 
                    su_orig$statistics[ind_main, c(3,4)], n_eff_orig[ind_main], 
                    x$diags$R_hat$psrf[ind_main,])
  colnames(orig_mat) <- Names
  summary_list$orig_pars <- orig_mat
  

  ## covariance
  CovP_mat <- cbind(su_orig$statistics[ind_SIGMA[covs], c(1,2)], su_orig$quantiles[ind_SIGMA[covs], c(1,3,5)], 
                    su_orig$statistics[ind_SIGMA[covs], c(3,4)], n_eff_orig[ind_SIGMA[covs]], 
                    x$diags$R_hat$psrf[ind_SIGMA[covs],])
  colnames(CovP_mat) <- Names
  if (length(covsR)==1) {
    CovR <- c(su_orig$statistics[ind_GAMMA[covsR], c(1,2)], su_orig$quantiles[ind_GAMMA[covsR], c(1,3,5)], 
              su_orig$statistics[ind_GAMMA[covsR], c(3,4)], n_eff_orig[ind_GAMMA[covsR]], 
              x$diags$R_hat$psrf[ind_GAMMA[covsR],])
    CovR_mat <- matrix(CovR, ncol = length(Names))
    colnames(CovR_mat) <- Names
    rownames(CovR_mat) <- names(n_eff_orig[ind_GAMMA[covsR]])
  } else if (length(covsR>1)) {
    CovR_mat <- cbind(su_orig$statistics[ind_GAMMA[covsR], c(1,2)], su_orig$quantiles[ind_GAMMA[covsR], c(1,3,5)], 
                      su_orig$statistics[ind_GAMMA[covsR], c(3,4)], n_eff_orig[ind_GAMMA[covsR]], 
                      x$diags$R_hat$psrf[ind_GAMMA[covsR],])
    colnames(CovR_mat) <- Names
  }
  if (!exists("CovR_mat")) CovR_mat <- NULL
  summary_list$Covs = list(process=CovP_mat, motor=CovR_mat)
  
  
  ## variances
  varsP_mat <- cbind(su_orig$statistics[ind_SIGMA[vars], c(1,2)], su_orig$quantiles[ind_SIGMA[vars], c(1,3,5)], 
                     su_orig$statistics[ind_SIGMA[vars], c(3,4)], n_eff_orig[ind_SIGMA[vars]], 
                     x$diags$R_hat$psrf[ind_SIGMA[vars],])
  colnames(varsP_mat) <- Names
  if (length(varsR)==1) {
    vars_R <- c(su_orig$statistics[ind_GAMMA[varsR], c(1,2)], su_orig$quantiles[ind_GAMMA[varsR], c(1,3,5)], 
                su_orig$statistics[ind_GAMMA[varsR], c(3,4)], n_eff_orig[ind_GAMMA[varsR]], 
                x$diags$R_hat$psrf[ind_GAMMA[varsR],])
    varsR_mat <- matrix(vars_R, ncol = length(Names))
    colnames(varsR_mat) <- Names
    rownames(varsR_mat) <- names(n_eff_orig[ind_GAMMA[varsR]])
  } else if (length(varsR>1)) {
    varsR_mat <- cbind(su_orig$statistics[ind_GAMMA[varsR], c(1,2)], su_orig$quantiles[ind_GAMMA[varsR], c(1,3,5)], 
                       su_orig$statistics[ind_GAMMA[varsR], c(3,4)], n_eff_orig[ind_GAMMA[varsR]], 
                       x$diags$R_hat$psrf[ind_GAMMA[varsR],])
    colnames(varsR_mat) <- Names
  }
  summary_list$Vars = list(process=varsP_mat, motor=varsR_mat)
  
  
  ## correlations
  CorrP_mat <- cbind(su$statistics[ind_SIGMA[covs], c(1,2)], su$quantiles[ind_SIGMA[covs], c(1,3,5)], 
                     su$statistics[ind_SIGMA[covs], c(3,4)], n_eff[ind_SIGMA[covs]], 
                     R_hat$psrf[ind_SIGMA[covs],])
  colnames(CorrP_mat) <- Names
  if (length(covsR)==1) {
    CorrR <- c(su$statistics[ind_GAMMA[covsR], c(1,2)], su$quantiles[ind_GAMMA[covsR], c(1,3,5)], 
               su$statistics[ind_GAMMA[covsR], c(3,4)], n_eff[ind_GAMMA[covsR]], 
               R_hat$psrf[ind_GAMMA[covsR],])
    CorrR_mat <- matrix(CorrR, ncol = length(Names))
    colnames(CorrR_mat) <- Names
    rownames(CorrR_mat) <- names(n_eff_orig[ind_GAMMA[covsR]])
  } else if (length(covsR>1)) {
    CorrR_mat <- cbind(su$statistics[ind_GAMMA[covsR], c(1,2)], su$quantiles[ind_GAMMA[covsR], c(1,3,5)], 
                       su$statistics[ind_GAMMA[covsR], c(3,4)], n_eff[ind_GAMMA[covsR]], 
                       R_hat$psrf[ind_GAMMA[covsR],])
    colnames(CorrR_mat) <- Names
  }
  if (!exists("CorrR_mat")) CorrR_mat <- NULL
  summary_list$Corrs = list(process=CorrP_mat, motor=CorrR_mat)
  
  
  ## standard deviation
  sdsP_mat <- cbind(su$statistics[ind_SIGMA[vars], c(1,2)], su$quantiles[ind_SIGMA[vars], c(1,3,5)], 
                    su$statistics[ind_SIGMA[vars], c(3,4)], n_eff[ind_SIGMA[vars]], 
                    R_hat$psrf[ind_SIGMA[vars],])
  colnames(sdsP_mat) <- Names
  if (length(varsR)==1) {
    sds_R <- c(su$statistics[ind_GAMMA[varsR], c(1,2)], su$quantiles[ind_GAMMA[varsR], c(1,3,5)], 
               su$statistics[ind_GAMMA[varsR], c(3,4)], n_eff_orig[ind_GAMMA[varsR]], 
               x$diags$R_hat$psrf[ind_GAMMA[varsR],])
    sdsR_mat <- matrix(sds_R, ncol = length(Names))
    colnames(sdsR_mat) <- Names
    rownames(sdsR_mat) <- names(n_eff_orig[ind_GAMMA[varsR]])
  } else if (length(varsR>1)) {
    sdsR_mat <- cbind(su$statistics[ind_GAMMA[varsR], c(1,2)], su$quantiles[ind_GAMMA[varsR], c(1,3,5)], 
                      su$statistics[ind_GAMMA[varsR], c(3,4)], n_eff_orig[ind_GAMMA[varsR]], 
                      x$diags$R_hat$psrf[ind_GAMMA[varsR],])
    colnames(sdsR_mat) <- Names
  }
  summary_list$SDs = list(process=sdsP_mat, motor=sdsR_mat)
  
  
  ## standard-deviation-correlation matrix
  elementsP <- su$statistics[ind_SIGMA, 1]
  CORR_SD_P <- matrix(0, ncol = Nparams, nrow = Nparams)
  CORR_SD_P[lower.tri(CORR_SD_P, diag = TRUE)] <- elementsP
  CORR_SD_P <- CORR_SD_P + t(CORR_SD_P) - diag(elementsP[vars])
  
  elementsR <- su$statistics[ind_GAMMA, 1]
  if (length(elementsR)==1) {
    CORR_SD_R <- elementsR
  } else {
    CORR_SD_R <- matrix(0, ncol = Nresps, nrow = Nresps)
    CORR_SD_R[lower.tri(CORR_SD_R, diag = TRUE)] <- elementsR
    CORR_SD_R <- CORR_SD_R + t(CORR_SD_R) - diag(elementsR[varsR])
  }
  summary_list$CorrSD = list(process=CORR_SD_P, motor=CORR_SD_R)
  
  
  ## Variance-covariance matrix
  elements_covP <- su_orig$statistics[ind_SIGMA, 1]
  COV_VAR_P <- matrix(0, ncol = Nparams, nrow = Nparams)
  COV_VAR_P[lower.tri(COV_VAR_P, diag = TRUE)] <- elements_covP
  COV_VAR_P <- COV_VAR_P + t(COV_VAR_P) - diag(elements_covP[vars])
  
  elements_covR <- su_orig$statistics[ind_GAMMA, 1]
  if (length(elements_covR)==1) {
    COV_VAR_R <- elements_covR
  } else {
    COV_VAR_R <- matrix(0, ncol = Nresps, nrow = Nresps)
    COV_VAR_R[lower.tri(COV_VAR_R, diag = TRUE)] <- elements_covR
    COV_VAR_R <- COV_VAR_R + t(COV_VAR_R) - diag(elements_covR[varsR])
  }
  summary_list$CovVar = list(process=COV_VAR_P, motor=COV_VAR_R)
  
  
  return(summary_list)
  
}



printSummaryERTMPT <- function(x, ...) {
  cat("\nCall: \n")
  print(x$call)
  cat("\n\n# Residual variance:\n")
  print(round(x$resid_var , x$round))
  cat("\n\n# Transformed main parameters (probabilities, process times in ms, and motor times in ms):\n")
  print(round(x$transformed_pars, x$round))
  cat("\n\n# Group-level parameters on original scale:\n")
  print(round(x$orig_pars, x$round))
  cat("\n\n# CORRELATIONS\n")
  cat("## Process-related:\n")
  print(round(x$Corrs$process, x$round))
  if (!is.null(x$Corrs$motor)) {
    cat("## Motor-related:\n")
    print(round(x$Corrs$motor, x$round))
  }
  cat("\n\n# STANDARD DEVIATIONS\n")
  cat("## Process-related:\n")
  print(round(x$SDs$process, x$round))
  cat("## Motor-related:\n")
  print(round(x$SDs$motor, x$round))
  cat("\n# NOTE: The same is available for variances and covariances\n")
  cat("\n\n# CORRELATION MATRICES WITH STANDARD DEVIATIONS ON DIAGONAL\n")
  cat("## Process-related:\n")
  print(round(x$CorrSD$process, x$round))
  cat("## Motor-related:\n")
  print(round(x$CorrSD$motor, x$round))
  cat("\n# NOTE: Covariance matrices are also available\n")
}

printSummaryDRTMPT <- function(x, ...) {
  cat("\nCall: \n")
  print(x$call)
  cat("\n\n# Residual variance:\n")
  print(round(x$resid_var , x$round))
  cat("\n\n# Transformed main parameters (probabilities, process times in ms, and motor times in ms):\n")
  print(round(x$transformed_pars, x$round))
  cat("\n\n# Group-level parameters on original scale:\n")
  print(round(x$orig_pars, x$round))
  cat("\n\n# CORRELATIONS\n")
  cat("## Process-related:\n")
  print(round(x$Corrs$process, x$round))
  if (!is.null(x$Corrs$motor)) {
    cat("## Motor-related:\n")
    print(round(x$Corrs$motor, x$round))
  }
  cat("\n\n# STANDARD DEVIATIONS\n")
  cat("## Process-related:\n")
  print(round(x$SDs$process, x$round))
  cat("## Motor-related:\n")
  print(round(x$SDs$motor, x$round))
  cat("\n# NOTE: The same is available for variances and covariances\n")
  cat("\n\n# CORRELATION MATRICES WITH STANDARD DEVIATIONS ON DIAGONAL\n")
  cat("## Process-related:\n")
  print(round(x$CorrSD$process, x$round))
  cat("## Motor-related:\n")
  print(round(x$CorrSD$motor, x$round))
  cat("\n# NOTE: Covariance matrices are also available\n")
}


#' @export
print.summary.ertmpt_fit <- function(x,  ...){
  printSummaryERTMPT(x)
}

#' @export
print.summary.rtmpt_fit <- function(x,  ...){
  printSummaryERTMPT(x)
}

#' @export
print.summary.drtmpt_fit <- function(x,  ...){
  printSummaryDRTMPT(x)
}


#' @export
summary.ertmpt_fit <- function(object, round=3, ...){
  summa <- object$summary
  summa$call <- object$specs$call
  summa$round <- round
  class(summa) <- "summary.ertmpt_fit"
  return(summa)
}

#' @export
summary.rtmpt_fit <- function(object, round=3, ...){
  summa <- object$summary
  summa$call <- object$specs$call
  summa$round <- round
  class(summa) <- "summary.ertmpt_fit"
  return(summa)
}

#' @export
summary.drtmpt_fit <- function(object, round=3, ...){
  summa <- object$summary
  summa$call <- object$specs$call
  summa$round <- round
  class(summa) <- "summary.drtmpt_fit"
  return(summa)
}


#' @export
print.ertmpt_fit <- function(x, ...) {
  cat("\nFUNCTION CALL\n\n")
  print(x$specs$call)
  
  cat("\n\nMEDIAN OF THE GROUP-LEVEL PARAMETERS\n\n")
  rowNms <- names(x$specs$model$params$probs)
  whole <- as.data.frame(matrix(data = rep(0, 3*length(rowNms)*x$specs$n.groups), nrow = length(rowNms)*x$specs$n.groups))
  Group_labels <- 0:(x$specs$n.groups-1)
  if (!is.null(x$specs$call[["old_label"]])) if (x$specs$call$old_label) Group_labels <- as.character(x$specs$transformation$group$old)
  rownames(whole) <- paste0(rep(rowNms, x$specs$n.groups), "[", rep(Group_labels, each = length(rowNms)), "]")
  colnames(whole) <- c("mu_probs", "mu_tau_minus", "mu_tau_plus")
  index <- which(is.na(x$specs$model$params$probs))
  if (x$specs$n.groups>1) for (i in 1:(x$specs$n.groups-1)) {index <- c(index, index+length(rowNms))}
  whole[index,1] <- x$diags$mu_probs[,3]
  if (length(which(!is.na(x$specs$model$params$probs)))>0) {
    index <- which(!is.na(x$specs$model$params$probs))
    index_orig <- index
    if (x$specs$n.groups>1) for (i in 1:(x$specs$n.groups-1)) {index <- c(index, index+length(rowNms))}
    if (any(x$specs$model$params$probs[1,index_orig] %in% rowNms)) {
      tmp <- numeric(length = length(index_orig)*x$specs$n.groups)
      if (all(x$specs$model$params$probs[1,index_orig] %in% rowNms)) {
        
        nms <- as.character(x$specs$model$params$probs[1,index_orig])
        ind_nms <- sapply(X = nms, FUN = function(i) {which(names(x$specs$model$params$probs) %in% i)})
        ind_mu_probs <- sapply(X = ind_nms, FUN = function(i) sum(is.na(x$specs$model$params$probs[1, 1:i])))
        if (x$specs$n.groups == 1) {
          tmp <- x$diags$mu_probs[ind_mu_probs,3]
        } else {
          ind_mu_probs_g <- ind_mu_probs; for ( i in 1:(x$specs$n.groups-1) ) {ind_mu_probs_g <- c(ind_mu_probs_g, ind_mu_probs_g+length(x$diags$mu_probs[,3])/x$specs$n.groups)}
          tmp <- x$diags$mu_probs[ind_mu_probs_g,3]
        }
        
      } else {
        ind_cnst <- which(!x$specs$model$params$probs[1,index_orig] %in% rowNms)
        ind_eql <- which(x$specs$model$params$probs[1,index_orig] %in% rowNms)
        nms <- as.character(x$specs$model$params$probs[1,index_orig[ind_eql]])
        ind_nms <- sapply(X = nms, FUN = function(i) {which(names(x$specs$model$params$probs) %in% i)})
        ind_mu_probs <- sapply(X = ind_nms, FUN = function(i) sum(is.na(x$specs$model$params$probs[1, 1:i])))
        if (x$specs$n.groups == 1) {
          tmp[ind_cnst] <- x$specs$model$params$probs[1,index_orig[ind_cnst]]
          tmp[ind_eql] <- x$diags$mu_probs[ind_mu_probs,3]
        } else {
          ind_cnst_g <- ind_cnst; for ( i in 1:(x$specs$n.groups-1) ) {ind_cnst_g <- c(ind_cnst_g, ind_cnst_g+length(index_orig))}
          ind_eql_g <- ind_eql; for ( i in 1:(x$specs$n.groups-1) ) {ind_eql_g <- c(ind_eql_g, ind_eql_g+length(index_orig))}
          tmp[ind_cnst_g] <- x$specs$model$params$probs[1,index_orig[ind_cnst]]
          ind_mu_probs_g <- ind_mu_probs; for ( i in 1:(x$specs$n.groups-1) ) {ind_mu_probs_g <- c(ind_mu_probs_g, ind_mu_probs_g+length(x$diags$mu_probs[,3])/x$specs$n.groups)}
          tmp[ind_eql_g] <- x$diags$mu_probs[ind_mu_probs_g,3]
        }
        
      }
      whole[index,1] <- tmp
    } else whole[index,1] <- rep(as.numeric(x$specs$model$params$probs[1,index_orig]), x$specs$n.groups)
  }
  
  index <- which(is.na(x$specs$model$params$taus[1,]))
  if (x$specs$n.groups>1) for (i in 1:(x$specs$n.groups-1)) {index <- c(index, index+length(rowNms))}
  whole[index,2] <- x$diags$mu_tau_minus[,3]
  if (length(which(!is.na(x$specs$model$params$taus[1,])))>0) {
    index <- which(!is.na(x$specs$model$params$taus[1,]))
    index_orig <- index
    if (x$specs$n.groups>1) for (i in 1:(x$specs$n.groups-1)) {index <- c(index, index+length(rowNms))}
    if (any(x$specs$model$params$taus[1,index_orig] %in% rowNms)) {
      tmp <- numeric(length = length(index_orig)*x$specs$n.groups)
      if (all(x$specs$model$params$taus[1,index_orig] %in% rowNms)) {
        
        nms <- as.character(x$specs$model$params$taus[1,index_orig])
        ind_nms <- sapply(X = nms, FUN = function(i) {which(rowNms %in% i)})
        ind_mu_tau_minus <- sapply(X = ind_nms, FUN = function(i) sum(is.na(x$specs$model$params$taus[1, 1:i])))
        if (x$specs$n.groups == 1) {
          tmp <- x$diags$mu_tau_minus[ind_mu_tau_minus,3]
        } else {
          ind_mu_tau_minus_g <- ind_mu_tau_minus; for ( i in 1:(x$specs$n.groups-1) ) {ind_mu_tau_minus_g <- c(ind_mu_tau_minus_g, ind_mu_tau_minus_g+length(x$diags$mu_tau_minus[,3])/x$specs$n.groups)}
          tmp <- x$diags$mu_tau_minus[ind_mu_tau_minus_g,3]
        }
        
      } else {
        ind_zero <- which(!x$specs$model$params$taus[1,index_orig] %in% rowNms)
        ind_eql <- which(x$specs$model$params$taus[1,index_orig] %in% rowNms)
        nms <- as.character(x$specs$model$params$taus[1,index_orig[ind_eql]])
        ind_nms <- sapply(X = nms, FUN = function(i) {which(rowNms %in% i)})
        ind_mu_tau_minus <- sapply(X = ind_nms, FUN = function(i) sum(is.na(x$specs$model$params$taus[1, 1:i])))
        if (x$specs$n.groups == 1) {
          tmp[ind_zero] <- x$specs$model$params$taus[1,index_orig[ind_zero]]
          tmp[ind_eql] <- x$diags$mu_tau_minus[ind_mu_tau_minus,3]
        } else {
          ind_zero_g <- ind_zero; for ( i in 1:(x$specs$n.groups-1) ) {ind_zero_g <- c(ind_zero_g, ind_zero_g+length(index_orig))}
          ind_eql_g <- ind_eql; for ( i in 1:(x$specs$n.groups-1) ) {ind_eql_g <- c(ind_eql_g, ind_eql_g+length(index_orig))}
          tmp[ind_zero_g] <- x$specs$model$params$taus[1,index_orig[ind_zero]]
          ind_mu_tau_minus_g <- ind_mu_tau_minus; for ( i in 1:(x$specs$n.groups-1) ) {ind_mu_tau_minus_g <- c(ind_mu_tau_minus_g, ind_mu_tau_minus_g+length(x$diags$mu_tau_minus[,3])/x$specs$n.groups)}
          tmp[ind_eql_g] <- x$diags$mu_tau_minus[ind_mu_tau_minus_g,3]
        }
        
      }
      whole[index,2] <- tmp
    } else whole[index,2] <- rep(as.numeric(x$specs$model$params$taus[1,index_orig]), x$specs$n.groups)
  }
  
  index <- which(is.na(x$specs$model$params$taus[2,]))
  if (x$specs$n.groups>1) for (i in 1:(x$specs$n.groups-1)) {index <- c(index, index+length(rowNms))}
  whole[index,3] <- x$diags$mu_tau_plus[,3]
  if (length(which(!is.na(x$specs$model$params$taus[2,])))>0) {
    index <- which(!is.na(x$specs$model$params$taus[2,]))
    index_orig <- index
    if (x$specs$n.groups>1) for (i in 1:(x$specs$n.groups-1)) {index <- c(index, index+length(rowNms))}
    if (any(x$specs$model$params$taus[2,index_orig] %in% rowNms)) {
      tmp <- numeric(length = length(index_orig)*x$specs$n.groups)
      if (all(x$specs$model$params$taus[2,index_orig] %in% rowNms)) {
        
        nms <- as.character(x$specs$model$params$taus[2,index_orig])
        ind_nms <- sapply(X = nms, FUN = function(i) {which(rowNms %in% i)})
        ind_mu_tau_plus <- sapply(X = ind_nms, FUN = function(i) sum(is.na(x$specs$model$params$taus[2, 1:i])))
        if (x$specs$n.groups == 1) {
          tmp <- x$diags$mu_tau_plus[ind_mu_tau_plus,3]
        } else {
          ind_mu_tau_plus_g <- ind_mu_tau_plus; for ( i in 1:(x$specs$n.groups-1) ) {ind_mu_tau_plus_g <- c(ind_mu_tau_plus_g, ind_mu_tau_plus_g+length(x$diags$mu_tau_plus[,3])/x$specs$n.groups)}
          tmp <- x$diags$mu_tau_plus[ind_mu_tau_plus_g,3]
        }
        
      } else {
        ind_zero <- which(!x$specs$model$params$taus[2,index_orig] %in% rowNms)
        ind_eql <- which(x$specs$model$params$taus[2,index_orig] %in% rowNms)
        nms <- as.character(x$specs$model$params$taus[2,index_orig[ind_eql]])
        ind_nms <- sapply(X = nms, FUN = function(i) {which(rowNms %in% i)})
        ind_mu_tau_plus <- sapply(X = ind_nms, FUN = function(i) sum(is.na(x$specs$model$params$taus[2, 1:i])))
        if (x$specs$n.groups == 1) {
          tmp[ind_zero] <- x$specs$model$params$taus[2,index_orig[ind_zero]]
          tmp[ind_eql] <- x$diags$mu_tau_plus[ind_mu_tau_plus,3]
        } else {
          ind_zero_g <- ind_zero; for ( i in 1:(x$specs$n.groups-1) ) {ind_zero_g <- c(ind_zero_g, ind_zero_g+length(index_orig))}
          ind_eql_g <- ind_eql; for ( i in 1:(x$specs$n.groups-1) ) {ind_eql_g <- c(ind_eql_g, ind_eql_g+length(index_orig))}
          tmp[ind_zero_g] <- x$specs$model$params$taus[2,index_orig[ind_zero]]
          ind_mu_tau_plus_g <- ind_mu_tau_plus; for ( i in 1:(x$specs$n.groups-1) ) {ind_mu_tau_plus_g <- c(ind_mu_tau_plus_g, ind_mu_tau_plus_g+length(x$diags$mu_tau_plus[,3])/x$specs$n.groups)}
          tmp[ind_eql_g] <- x$diags$mu_tau_plus[ind_mu_tau_plus_g,3]
        }
        
      }
      whole[index,3] <- tmp
    } else whole[index,3] <- rep(as.numeric(x$specs$model$params$taus[2,index_orig]), x$specs$n.groups)
  }
  
  cat("Process-related parameters:\n")
  print(whole)
  cat("\n* NOTE 1: Process completion times in ms.")
  if (any(!is.na(x$specs$model$params$probs))) {
    pass <- paste0("\n* NOTE 2: Constants are also displayed in the process parameter table", 
                   if(any(!is.na(x$specs$model$params$taus))) {"\n\t  as well as zeros for suppressed process completion times"}, ".")
    cat(pass)
  } else if (any(!is.na(x$specs$model$params$taus))) {
    pass <- paste0("\n* NOTE 2: Zeros for suppressed process completion times are also displayed\n\t  in the process parameter table.")
  }
  cat("\n---------------------------\n")
  
  cat("\nEncoding plus motor execution parameter(s):\n")
  print(data.frame(mu_gamma=x$diags$mu_gamma[,3]))
  cat("\n* NOTE 1: Encoding plus motor execution time(s) in ms.")
  cat("\n-------------------------------------------\n\n")
}

#' @export
print.rtmpt_fit <- function(x, ...) {
  cat("\nFUNCTION CALL\n\n")
  print(x$specs$call)
  
  cat("\n\nMEDIAN OF THE GROUP-LEVEL PARAMETERS\n\n")
  rowNms <- names(x$specs$model$params$probs)
  whole <- as.data.frame(matrix(data = rep(0, 3*length(rowNms)*x$specs$n.groups), nrow = length(rowNms)*x$specs$n.groups))
  Group_labels <- 0:(x$specs$n.groups-1)
  if (!is.null(x$specs$call[["old_label"]])) if (x$specs$call$old_label) Group_labels <- as.character(x$specs$transformation$group$old)
  rownames(whole) <- paste0(rep(rowNms, x$specs$n.groups), "[", rep(Group_labels, each = length(rowNms)), "]")
  colnames(whole) <- c("mu_probs", "mu_tau_minus", "mu_tau_plus")
  index <- which(is.na(x$specs$model$params$probs))
  if (x$specs$n.groups>1) for (i in 1:(x$specs$n.groups-1)) {index <- c(index, index+length(rowNms))}
  whole[index,1] <- x$diags$mu_probs[,3]
  if (length(which(!is.na(x$specs$model$params$probs)))>0) {
    index <- which(!is.na(x$specs$model$params$probs))
    index_orig <- index
    if (x$specs$n.groups>1) for (i in 1:(x$specs$n.groups-1)) {index <- c(index, index+length(rowNms))}
    if (any(x$specs$model$params$probs[1,index_orig] %in% rowNms)) {
      tmp <- numeric(length = length(index_orig)*x$specs$n.groups)
      if (all(x$specs$model$params$probs[1,index_orig] %in% rowNms)) {
        
        nms <- as.character(x$specs$model$params$probs[1,index_orig])
        ind_nms <- sapply(X = nms, FUN = function(i) {which(names(x$specs$model$params$probs) %in% i)})
        ind_mu_probs <- sapply(X = ind_nms, FUN = function(i) sum(is.na(x$specs$model$params$probs[1, 1:i])))
        if (x$specs$n.groups == 1) {
          tmp <- x$diags$mu_probs[ind_mu_probs,3]
        } else {
          ind_mu_probs_g <- ind_mu_probs; for ( i in 1:(x$specs$n.groups-1) ) {ind_mu_probs_g <- c(ind_mu_probs_g, ind_mu_probs_g+length(x$diags$mu_probs[,3])/x$specs$n.groups)}
          tmp <- x$diags$mu_probs[ind_mu_probs_g,3]
        }
        
      } else {
        ind_cnst <- which(!x$specs$model$params$probs[1,index_orig] %in% rowNms)
        ind_eql <- which(x$specs$model$params$probs[1,index_orig] %in% rowNms)
        nms <- as.character(x$specs$model$params$probs[1,index_orig[ind_eql]])
        ind_nms <- sapply(X = nms, FUN = function(i) {which(names(x$specs$model$params$probs) %in% i)})
        ind_mu_probs <- sapply(X = ind_nms, FUN = function(i) sum(is.na(x$specs$model$params$probs[1, 1:i])))
        if (x$specs$n.groups == 1) {
          tmp[ind_cnst] <- x$specs$model$params$probs[1,index_orig[ind_cnst]]
          tmp[ind_eql] <- x$diags$mu_probs[ind_mu_probs,3]
        } else {
          ind_cnst_g <- ind_cnst; for ( i in 1:(x$specs$n.groups-1) ) {ind_cnst_g <- c(ind_cnst_g, ind_cnst_g+length(index_orig))}
          ind_eql_g <- ind_eql; for ( i in 1:(x$specs$n.groups-1) ) {ind_eql_g <- c(ind_eql_g, ind_eql_g+length(index_orig))}
          tmp[ind_cnst_g] <- x$specs$model$params$probs[1,index_orig[ind_cnst]]
          ind_mu_probs_g <- ind_mu_probs; for ( i in 1:(x$specs$n.groups-1) ) {ind_mu_probs_g <- c(ind_mu_probs_g, ind_mu_probs_g+length(x$diags$mu_probs[,3])/x$specs$n.groups)}
          tmp[ind_eql_g] <- x$diags$mu_probs[ind_mu_probs_g,3]
        }
        
      }
      whole[index,1] <- tmp
    } else whole[index,1] <- rep(as.numeric(x$specs$model$params$probs[1,index_orig]), x$specs$n.groups)
  }
  
  index <- which(is.na(x$specs$model$params$taus[1,]))
  if (x$specs$n.groups>1) for (i in 1:(x$specs$n.groups-1)) {index <- c(index, index+length(rowNms))}
  whole[index,2] <- x$diags$mu_tau_minus[,3]
  if (length(which(!is.na(x$specs$model$params$taus[1,])))>0) {
    index <- which(!is.na(x$specs$model$params$taus[1,]))
    index_orig <- index
    if (x$specs$n.groups>1) for (i in 1:(x$specs$n.groups-1)) {index <- c(index, index+length(rowNms))}
    if (any(x$specs$model$params$taus[1,index_orig] %in% rowNms)) {
      tmp <- numeric(length = length(index_orig)*x$specs$n.groups)
      if (all(x$specs$model$params$taus[1,index_orig] %in% rowNms)) {
        
        nms <- as.character(x$specs$model$params$taus[1,index_orig])
        ind_nms <- sapply(X = nms, FUN = function(i) {which(rowNms %in% i)})
        ind_mu_tau_minus <- sapply(X = ind_nms, FUN = function(i) sum(is.na(x$specs$model$params$taus[1, 1:i])))
        if (x$specs$n.groups == 1) {
          tmp <- x$diags$mu_tau_minus[ind_mu_tau_minus,3]
        } else {
          ind_mu_tau_minus_g <- ind_mu_tau_minus; for ( i in 1:(x$specs$n.groups-1) ) {ind_mu_tau_minus_g <- c(ind_mu_tau_minus_g, ind_mu_tau_minus_g+length(x$diags$mu_tau_minus[,3])/x$specs$n.groups)}
          tmp <- x$diags$mu_tau_minus[ind_mu_tau_minus_g,3]
        }
        
      } else {
        ind_zero <- which(!x$specs$model$params$taus[1,index_orig] %in% rowNms)
        ind_eql <- which(x$specs$model$params$taus[1,index_orig] %in% rowNms)
        nms <- as.character(x$specs$model$params$taus[1,index_orig[ind_eql]])
        ind_nms <- sapply(X = nms, FUN = function(i) {which(rowNms %in% i)})
        ind_mu_tau_minus <- sapply(X = ind_nms, FUN = function(i) sum(is.na(x$specs$model$params$taus[1, 1:i])))
        if (x$specs$n.groups == 1) {
          tmp[ind_zero] <- x$specs$model$params$taus[1,index_orig[ind_zero]]
          tmp[ind_eql] <- x$diags$mu_tau_minus[ind_mu_tau_minus,3]
        } else {
          ind_zero_g <- ind_zero; for ( i in 1:(x$specs$n.groups-1) ) {ind_zero_g <- c(ind_zero_g, ind_zero_g+length(index_orig))}
          ind_eql_g <- ind_eql; for ( i in 1:(x$specs$n.groups-1) ) {ind_eql_g <- c(ind_eql_g, ind_eql_g+length(index_orig))}
          tmp[ind_zero_g] <- x$specs$model$params$taus[1,index_orig[ind_zero]]
          ind_mu_tau_minus_g <- ind_mu_tau_minus; for ( i in 1:(x$specs$n.groups-1) ) {ind_mu_tau_minus_g <- c(ind_mu_tau_minus_g, ind_mu_tau_minus_g+length(x$diags$mu_tau_minus[,3])/x$specs$n.groups)}
          tmp[ind_eql_g] <- x$diags$mu_tau_minus[ind_mu_tau_minus_g,3]
        }
        
      }
      whole[index,2] <- tmp
    } else whole[index,2] <- rep(as.numeric(x$specs$model$params$taus[1,index_orig]), x$specs$n.groups)
  }
  
  index <- which(is.na(x$specs$model$params$taus[2,]))
  if (x$specs$n.groups>1) for (i in 1:(x$specs$n.groups-1)) {index <- c(index, index+length(rowNms))}
  whole[index,3] <- x$diags$mu_tau_plus[,3]
  if (length(which(!is.na(x$specs$model$params$taus[2,])))>0) {
    index <- which(!is.na(x$specs$model$params$taus[2,]))
    index_orig <- index
    if (x$specs$n.groups>1) for (i in 1:(x$specs$n.groups-1)) {index <- c(index, index+length(rowNms))}
    if (any(x$specs$model$params$taus[2,index_orig] %in% rowNms)) {
      tmp <- numeric(length = length(index_orig)*x$specs$n.groups)
      if (all(x$specs$model$params$taus[2,index_orig] %in% rowNms)) {
        
        nms <- as.character(x$specs$model$params$taus[2,index_orig])
        ind_nms <- sapply(X = nms, FUN = function(i) {which(rowNms %in% i)})
        ind_mu_tau_plus <- sapply(X = ind_nms, FUN = function(i) sum(is.na(x$specs$model$params$taus[2, 1:i])))
        if (x$specs$n.groups == 1) {
          tmp <- x$diags$mu_tau_plus[ind_mu_tau_plus,3]
        } else {
          ind_mu_tau_plus_g <- ind_mu_tau_plus; for ( i in 1:(x$specs$n.groups-1) ) {ind_mu_tau_plus_g <- c(ind_mu_tau_plus_g, ind_mu_tau_plus_g+length(x$diags$mu_tau_plus[,3])/x$specs$n.groups)}
          tmp <- x$diags$mu_tau_plus[ind_mu_tau_plus_g,3]
        }
        
      } else {
        ind_zero <- which(!x$specs$model$params$taus[2,index_orig] %in% rowNms)
        ind_eql <- which(x$specs$model$params$taus[2,index_orig] %in% rowNms)
        nms <- as.character(x$specs$model$params$taus[2,index_orig[ind_eql]])
        ind_nms <- sapply(X = nms, FUN = function(i) {which(rowNms %in% i)})
        ind_mu_tau_plus <- sapply(X = ind_nms, FUN = function(i) sum(is.na(x$specs$model$params$taus[2, 1:i])))
        if (x$specs$n.groups == 1) {
          tmp[ind_zero] <- x$specs$model$params$taus[2,index_orig[ind_zero]]
          tmp[ind_eql] <- x$diags$mu_tau_plus[ind_mu_tau_plus,3]
        } else {
          ind_zero_g <- ind_zero; for ( i in 1:(x$specs$n.groups-1) ) {ind_zero_g <- c(ind_zero_g, ind_zero_g+length(index_orig))}
          ind_eql_g <- ind_eql; for ( i in 1:(x$specs$n.groups-1) ) {ind_eql_g <- c(ind_eql_g, ind_eql_g+length(index_orig))}
          tmp[ind_zero_g] <- x$specs$model$params$taus[2,index_orig[ind_zero]]
          ind_mu_tau_plus_g <- ind_mu_tau_plus; for ( i in 1:(x$specs$n.groups-1) ) {ind_mu_tau_plus_g <- c(ind_mu_tau_plus_g, ind_mu_tau_plus_g+length(x$diags$mu_tau_plus[,3])/x$specs$n.groups)}
          tmp[ind_eql_g] <- x$diags$mu_tau_plus[ind_mu_tau_plus_g,3]
        }
        
      }
      whole[index,3] <- tmp
    } else whole[index,3] <- rep(as.numeric(x$specs$model$params$taus[2,index_orig]), x$specs$n.groups)
  }
  
  cat("Process-related parameters:\n")
  print(whole)
  cat("\n* NOTE 1: Process completion times in ms.")
  if (any(!is.na(x$specs$model$params$probs))) {
    pass <- paste0("\n* NOTE 2: Constants are also displayed in the process parameter table", 
                   if(any(!is.na(x$specs$model$params$taus))) {"\n\t  as well as zeros for suppressed process completion times"}, ".")
    cat(pass)
  } else if (any(!is.na(x$specs$model$params$taus))) {
    pass <- paste0("\n* NOTE 2: Zeros for suppressed process completion times are also displayed\n\t  in the process parameter table.")
  }
  cat("\n---------------------------\n")
  
  cat("\nEncoding plus motor execution parameter(s):\n")
  print(data.frame(mu_gamma=x$diags$mu_gamma[,3]))
  cat("\n* NOTE 1: Encoding plus motor execution time(s) in ms.")
  cat("\n-------------------------------------------\n\n")
}


#' @export
print.drtmpt_fit <- function(x, ...) {
  cat("\nFUNCTION CALL\n\n")
  print(x$specs$call)
  
  cat("\n\nMEDIAN OF THE GROUP-LEVEL PARAMETERS\n\n")
  rowNms <- names(x$specs$model$params$threshold)
  whole <- as.data.frame(matrix(data = rep(0, 3*length(rowNms)*x$specs$n.groups), nrow = length(rowNms)*x$specs$n.groups))
  Group_labels <- 0:(x$specs$n.groups-1)
  if (!is.null(x$specs$call[["old_label"]])) if (x$specs$call$old_label) Group_labels <- as.character(x$specs$transformation$group$old)
  rownames(whole) <- paste0(rep(rowNms, x$specs$n.groups), "[", rep(Group_labels, each = length(rowNms)), "]")
  colnames(whole) <- c("mu_a", "mu_nu", "mu_omega")
  index <- which(is.na(x$specs$model$params$threshold[1,]))
  if (x$specs$n.groups>1) for (i in 1:(x$specs$n.groups-1)) {index <- c(index, index+length(rowNms))}
  whole[index,1] <- x$diags$mu_a[,3]
  if (length(which(!is.na(x$specs$model$params$threshold[1,])))>0) {
    index <- which(!is.na(x$specs$model$params$threshold[1,]))
    index_orig <- index
    if (x$specs$n.groups>1) for (i in 1:(x$specs$n.groups-1)) {index <- c(index, index+length(rowNms))}
    if (any(x$specs$model$params$threshold[1,index_orig] %in% rowNms)) {
      tmp <- numeric(length = length(index_orig)*x$specs$n.groups)
      if (all(x$specs$model$params$threshold[1,index_orig] %in% rowNms)) {
        
        nms <- as.character(x$specs$model$params$threshold[1,index_orig])
        ind_nms <- sapply(X = nms, FUN = function(i) {which(names(x$specs$model$params$threshold[1,]) %in% i)})
        ind_mu_a <- sapply(X = ind_nms, FUN = function(i) sum(is.na(x$specs$model$params$threshold[1, 1:i])))
        if (x$specs$n.groups == 1) {
          tmp <- x$diags$mu_a[ind_mu_a,3]
        } else {
          ind_mu_a_g <- ind_mu_a; for ( i in 1:(x$specs$n.groups-1) ) {ind_mu_a_g <- c(ind_mu_a_g, ind_mu_a_g+length(x$diags$mu_a[,3])/x$specs$n.groups)}
          tmp <- x$diags$mu_a[ind_mu_a_g,3]
        }
        
      } else {
        ind_cnst <- which(!x$specs$model$params$threshold[1,index_orig] %in% rowNms)
        ind_eql <- which(x$specs$model$params$threshold[1,index_orig] %in% rowNms)
        nms <- as.character(x$specs$model$params$threshold[1,index_orig[ind_eql]])
        ind_nms <- sapply(X = nms, FUN = function(i) {which(names(x$specs$model$params$threshold[1,]) %in% i)})
        ind_mu_a <- sapply(X = ind_nms, FUN = function(i) sum(is.na(x$specs$model$params$threshold[1, 1:i])))
        if (x$specs$n.groups == 1) {
          tmp[ind_cnst] <- x$specs$model$params$threshold[1,index_orig[ind_cnst]]
          tmp[ind_eql] <- x$diags$mu_a[ind_mu_a,3]
        } else {
          ind_cnst_g <- ind_cnst; for ( i in 1:(x$specs$n.groups-1) ) {ind_cnst_g <- c(ind_cnst_g, ind_cnst_g+length(index_orig))}
          ind_eql_g <- ind_eql; for ( i in 1:(x$specs$n.groups-1) ) {ind_eql_g <- c(ind_eql_g, ind_eql_g+length(index_orig))}
          tmp[ind_cnst_g] <- x$specs$model$params$threshold[1,index_orig[ind_cnst]]
          ind_mu_a_g <- ind_mu_a; for ( i in 1:(x$specs$n.groups-1) ) {ind_mu_a_g <- c(ind_mu_a_g, ind_mu_a_g+length(x$diags$mu_a[,3])/x$specs$n.groups)}
          tmp[ind_eql_g] <- x$diags$mu_a[ind_mu_a_g,3]
        }
        
      }
      whole[index,1] <- tmp
    } else whole[index,1] <- rep(as.numeric(x$specs$model$params$threshold[1,index_orig]), x$specs$n.groups)
  }
  
  index <- which(is.na(x$specs$model$params$driftrate[1,]))
  if (x$specs$n.groups>1) for (i in 1:(x$specs$n.groups-1)) {index <- c(index, index+length(rowNms))}
  whole[index,2] <- x$diags$mu_nu[,3]
  if (length(which(!is.na(x$specs$model$params$driftrate[1,])))>0) {
    index <- which(!is.na(x$specs$model$params$driftrate[1,]))
    index_orig <- index
    if (x$specs$n.groups>1) for (i in 1:(x$specs$n.groups-1)) {index <- c(index, index+length(rowNms))}
    if (any(x$specs$model$params$driftrate[1,index_orig] %in% rowNms)) {
      tmp <- numeric(length = length(index_orig)*x$specs$n.groups)
      if (all(x$specs$model$params$driftrate[1,index_orig] %in% rowNms)) {
        
        nms <- as.character(x$specs$model$params$driftrate[1,index_orig])
        ind_nms <- sapply(X = nms, FUN = function(i) {which(rowNms %in% i)})
        ind_mu_nu <- sapply(X = ind_nms, FUN = function(i) sum(is.na(x$specs$model$params$driftrate[1, 1:i])))
        if (x$specs$n.groups == 1) {
          tmp <- x$diags$mu_nu[ind_mu_nu,3]
        } else {
          ind_mu_nu_g <- ind_mu_nu; for ( i in 1:(x$specs$n.groups-1) ) {ind_mu_nu_g <- c(ind_mu_nu_g, ind_mu_nu_g+length(x$diags$mu_nu[,3])/x$specs$n.groups)}
          tmp <- x$diags$mu_nu[ind_mu_nu_g,3]
        }
        
      } else {
        ind_zero <- which(!x$specs$model$params$driftrate[1,index_orig] %in% rowNms)
        ind_eql <- which(x$specs$model$params$driftrate[1,index_orig] %in% rowNms)
        nms <- as.character(x$specs$model$params$driftrate[1,index_orig[ind_eql]])
        ind_nms <- sapply(X = nms, FUN = function(i) {which(rowNms %in% i)})
        ind_mu_nu <- sapply(X = ind_nms, FUN = function(i) sum(is.na(x$specs$model$params$driftrate[1, 1:i])))
        if (x$specs$n.groups == 1) {
          tmp[ind_zero] <- x$specs$model$params$driftrate[1,index_orig[ind_zero]]
          tmp[ind_eql] <- x$diags$mu_nu[ind_mu_nu,3]
        } else {
          ind_zero_g <- ind_zero; for ( i in 1:(x$specs$n.groups-1) ) {ind_zero_g <- c(ind_zero_g, ind_zero_g+length(index_orig))}
          ind_eql_g <- ind_eql; for ( i in 1:(x$specs$n.groups-1) ) {ind_eql_g <- c(ind_eql_g, ind_eql_g+length(index_orig))}
          tmp[ind_zero_g] <- x$specs$model$params$driftrate[1,index_orig[ind_zero]]
          ind_mu_nu_g <- ind_mu_nu; for ( i in 1:(x$specs$n.groups-1) ) {ind_mu_nu_g <- c(ind_mu_nu_g, ind_mu_nu_g+length(x$diags$mu_nu[,3])/x$specs$n.groups)}
          tmp[ind_eql_g] <- x$diags$mu_nu[ind_mu_nu_g,3]
        }
        
      }
      whole[index,2] <- tmp
    } else whole[index,2] <- rep(as.numeric(x$specs$model$params$driftrate[1,index_orig]), x$specs$n.groups)
  }
  
  index <- which(is.na(x$specs$model$params$startpoint[1,]))
  if (x$specs$n.groups>1) for (i in 1:(x$specs$n.groups-1)) {index <- c(index, index+length(rowNms))}
  whole[index,3] <- x$diags$mu_omega[,3]
  if (length(which(!is.na(x$specs$model$params$startpoint[1,])))>0) {
    index <- which(!is.na(x$specs$model$params$startpoint[1,]))
    index_orig <- index
    if (x$specs$n.groups>1) for (i in 1:(x$specs$n.groups-1)) {index <- c(index, index+length(rowNms))}
    if (any(x$specs$model$params$startpoint[1,index_orig] %in% rowNms)) {
      tmp <- numeric(length = length(index_orig)*x$specs$n.groups)
      if (all(x$specs$model$params$startpoint[1,index_orig] %in% rowNms)) {
        
        nms <- as.character(x$specs$model$params$startpoint[1,index_orig])
        ind_nms <- sapply(X = nms, FUN = function(i) {which(rowNms %in% i)})
        ind_mu_omega <- sapply(X = ind_nms, FUN = function(i) sum(is.na(x$specs$model$params$startpoint[1, 1:i])))
        if (x$specs$n.groups == 1) {
          tmp <- x$diags$mu_omega[ind_mu_omega,3]
        } else {
          ind_mu_omega_g <- ind_mu_omega; for ( i in 1:(x$specs$n.groups-1) ) {ind_mu_omega_g <- c(ind_mu_omega_g, ind_mu_omega_g+length(x$diags$mu_omega[,3])/x$specs$n.groups)}
          tmp <- x$diags$mu_omega[ind_mu_omega_g,3]
        }
        
      } else {
        ind_zero <- which(!x$specs$model$params$startpoint[1,index_orig] %in% rowNms)
        ind_eql <- which(x$specs$model$params$startpoint[1,index_orig] %in% rowNms)
        nms <- as.character(x$specs$model$params$startpoint[1,index_orig[ind_eql]])
        ind_nms <- sapply(X = nms, FUN = function(i) {which(rowNms %in% i)})
        ind_mu_omega <- sapply(X = ind_nms, FUN = function(i) sum(is.na(x$specs$model$params$startpoint[1, 1:i])))
        if (x$specs$n.groups == 1) {
          tmp[ind_zero] <- x$specs$model$params$startpoint[1,index_orig[ind_zero]]
          tmp[ind_eql] <- x$diags$mu_omega[ind_mu_omega,3]
        } else {
          ind_zero_g <- ind_zero; for ( i in 1:(x$specs$n.groups-1) ) {ind_zero_g <- c(ind_zero_g, ind_zero_g+length(index_orig))}
          ind_eql_g <- ind_eql; for ( i in 1:(x$specs$n.groups-1) ) {ind_eql_g <- c(ind_eql_g, ind_eql_g+length(index_orig))}
          tmp[ind_zero_g] <- x$specs$model$params$startpoint[1,index_orig[ind_zero]]
          ind_mu_omega_g <- ind_mu_omega; for ( i in 1:(x$specs$n.groups-1) ) {ind_mu_omega_g <- c(ind_mu_omega_g, ind_mu_omega_g+length(x$diags$mu_omega[,3])/x$specs$n.groups)}
          tmp[ind_eql_g] <- x$diags$mu_omega[ind_mu_omega_g,3]
        }
        
      }
      whole[index,3] <- tmp
    } else whole[index,3] <- rep(as.numeric(x$specs$model$params$startpoint[1,index_orig]), x$specs$n.groups)
  }
  
  cat("Process-related parameters:\n")
  print(whole)
  cat("\n* NOTE 1: Process parameters are displayed on the natural scale of diffusion parameters.")
  if (any(!is.na(x$specs$model$params$threshold)) | any(!is.na(x$specs$model$params$driftrate)) | any(!is.na(x$specs$model$params$startpoint))) {
    pass <- paste0("\n* NOTE 2: Constants are also displayed in the process parameter table", ".")
    cat(pass)
  }
  cat("\n---------------------------\n")
  
  cat("\nEncoding plus motor execution parameter(s):\n")
  print(data.frame(mu_gamma=x$diags$mu_gamma[,3]))
  cat("\n* NOTE 1: Encoding plus motor execution time(s) in ms.")
  cat("\n-------------------------------------------\n\n")
}
