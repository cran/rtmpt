StddevCorr2Cov <- function(mat, dim2) {
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  if ((dim2*(dim2-1)/2+dim2) != dim(mat)[2]) stop("wrong dim2 or 2nd dim of mat!")
  
  matcov <- matrix(NA, ncol = dim(mat)[2], nrow = dim(mat)[1])
  SDCorr <- matrix(0, ncol = dim2, nrow = dim2)
  for (i in 1:dim(mat)[1]) {
    SDCorr[lower.tri(diag(dim2, nrow = dim2), diag = TRUE)] <- mat[i,]
    SDs <- diag(SDCorr)
    Corr <- SDCorr - diag(SDs, nrow = dim2) + diag(dim2)
    upperCorr <- t(Corr)-diag(diag(Corr))
    Corr <- Corr + upperCorr
    Cov <- (SDs %*% t(SDs)) * Corr
    matcov[i,] <- Cov[lower.tri(Cov, diag = TRUE)]
  }
  return(matcov)
  
}

#' @importFrom coda effectiveSize gelman.diag
#' @importFrom stats dnorm pnorm
writeSummaryRTMPT <- function(x, keep, ...) {

  prob_names <- colnames(x$specs$model$params$probs[1, which(is.na(x$specs$model$params$probs[1,]))])
  tauM_names <- colnames(x$specs$model$params$taus[1, which(is.na(x$specs$model$params$taus[1,]))])
  tauP_names <- colnames(x$specs$model$params$taus[2, which(is.na(x$specs$model$params$taus[2,]))])
  tau_names <- c(tauM_names, tauP_names)
  Nprobs <- sum(is.na(x$specs$model$params$probs))
  Ntaus <- sum(is.na(x$specs$model$params$taus))
  Nparams <- Nprobs+Ntaus
  Nresps <- length(unique(x$specs$model$responses$RESP))
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
    samp[[n]][, ind_SIGMA] <- StddevCorr2Cov(samp[[n]][, ind_SIGMA], Nparams)
    samp[[n]][, ind_resps] <- 1000*(samp[[n]][, ind_resps])
    samp[[n]][, ind_GAMMA] <- StddevCorr2Cov(samp[[n]][, ind_GAMMA], Nresps)
  }
  
  # summary stats
  su <- summary(samp)
  su_orig <- summary(x$samples)
  sudelta <- summary(delta)
  n_eff <- effectiveSize(samp)
  n_eff_orig <- effectiveSize(x$samples)
  n_eff_delta <- effectiveSize(delta)
  R_hat <- gelman.diag(samp)
  R_hat_delta <- gelman.diag(delta)
  
  
  # output
  
  omega2 <- as.numeric(c(su$statistics[ind_omega2, c(1,2)], su$quantiles[ind_omega2, c(1,3,5)], 
                         su$statistics[ind_omega2, c(3,4)], n_eff[ind_omega2], R_hat$psrf[ind_omega2,]))
  omega2_mat <- matrix(omega2, ncol = length(Names))
  colnames(omega2_mat) <- Names
  rownames(omega2_mat) <- "omega_squared"
  
  ind_main <- c(ind_probs, ind_taus)
  main_mat <- cbind(su$statistics[ind_main, c(1,2)], su$quantiles[ind_main, c(1,3,5)], 
                    su$statistics[ind_main, c(3,4)], n_eff[ind_main], R_hat$psrf[ind_main,])
  if (Nresps > 1) {
	main_mat <- rbind(main_mat, cbind(sudelta$statistics[, c(1,2)], sudelta$quantiles[, c(1,3,5)], 
                                    sudelta$statistics[, c(3,4)], n_eff_delta, R_hat_delta$psrf))
  } else {
    main_mat <- rbind(main_mat, matrix(c(sudelta$statistics[c(1,2)], sudelta$quantiles[c(1,3,5)], 
	                                     sudelta$statistics[c(3,4)], n_eff_delta, R_hat_delta$psrf), nrow = 1))
  }
  
  colnames(main_mat) <- Names
  rnams <- c("theta_", "tau_minus_", "tau_plus_", "delta_")
  theta_names <- paste0(rnams[1], prob_names)
  tau_M_names <- paste0(rnams[2], tauM_names)
  tau_P_names <- paste0(rnams[3], tauP_names)
  delta_names <- paste0(rnams[4], rep(paste0("R", 0:(Nresps-1))))
  if(Ngroups>1) {
    theta_names <- paste0(rep(theta_names, Ngroups), rep(paste0("[", group_labels, "]"), each=Nprobs))
    tau_M_names <- paste0(rep(tau_M_names, Ngroups), rep(paste0("[", group_labels, "]"), each=Nprobs))
    tau_P_names <- paste0(rep(tau_P_names, Ngroups), rep(paste0("[", group_labels, "]"), each=Nprobs))
    delta_names <- paste0(rep(delta_names, Ngroups), rep(paste0("[", group_labels, "]"), each=Nprobs))
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
  
  
  elements_covR <- su_orig$statistics[ind_GAMMA, 1]
  if (length(elements_covR)==1) {
    COV_VAR_R <- elements_covR
  } else {
    COV_VAR_R <- matrix(0, ncol = Nresps, nrow = Nresps)
    COV_VAR_R[lower.tri(COV_VAR_R, diag = TRUE)] <- elements_covR
    COV_VAR_R <- COV_VAR_R + t(COV_VAR_R) - diag(elements_covR[sdsR])
  }
  
  summary_list <- list(resid_var = omega2_mat, 
                       main_pars = main_mat,
                       orig_pars = orig_mat, 
                       Corrs = list(process=CorrP_mat, motor=CorrR_mat), 
                       SDs = list(process=sdsP_mat, motor=sdsR_mat),
                       Covs = list(process=CovP_mat, motor=CovR_mat), 
                       Vars = list(process=VarP_mat, motor=VarR_mat), 
                       CorrSD = list(process=CORR_SD_P, motor=CORR_SD_R), 
                       CovVar = list(process=COV_VAR_P, motor=COV_VAR_R))
  return(summary_list)
}



printSummaryRTMPT <- function(x, ...) {
  cat("\nCall: \n")
  print(x$call)
  cat("\n\n# Residual variance:\n")
  print(round(x$resid_var , x$round))
  cat("\n\n# Transformed main parameters (probabilities, process times in ms, and motor times in ms):\n")
  print(round(x$main_pars, x$round))
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
  cat("\n\n# CORRELATION MATRICES WITH STANDARD DEVIATIO ON DIAGONAL\n")
  cat("## Process-related:\n")
  print(round(x$CorrSD$process, x$round))
  cat("## Motor-related:\n")
  print(round(x$CorrSD$motor, x$round))
  cat("\n# NOTE: Covariance matrices are also available\n")
}


#' @export
print.summary.rtmpt_fit <- function(x,  ...){
  printSummaryRTMPT(x)
}


#' @export
summary.rtmpt_fit <- function(object, round=3, ...){
  summa <- object$summary
  summa$call <- object$specs$call
  summa$round <- round
  class(summa) <- "summary.rtmpt_fit"
  return(summa)
}


#' @export
print.rtmpt_fit <- function(x, ...) {
  cat("\nFUNCTION CALL\n\n")
  print(x$specs$call)
  
  cat("\n\nMEDIAN OF THE GROUP-LEVEL PARAMETERS\n\n")
  rowNms <- names(x$specs$model$params$probs)
  whole <- as.data.frame(matrix(data = rep(0, 3*length(rowNms)*x$specs$n.groups), nrow = length(rowNms)*x$specs$n.groups))
  rownames(whole) <- paste0(rep(rowNms, x$specs$n.groups), "[", rep(0:(x$specs$n.groups-1), each = length(rowNms)), "]")
  colnames(whole) <- c("mu_probs", "mu_tau_minus", "mu_tau_plus")
  index <- which(is.na(x$specs$model$params$probs))
  if (x$specs$n.groups>1) for (i in 1:(x$specs$n.groups-1)) {index <- c(index, index+length(rowNms))}
  whole[index,1] <- x$diags$mu_probs[,3]
  if (length(which(!is.na(x$specs$model$params$probs)))>0) {
    index <- which(!is.na(x$specs$model$params$probs))
    index_orig <- index
    if (x$specs$n.groups>1) for (i in 1:(x$specs$n.groups-1)) {index <- c(index, index+length(rowNms))}
    whole[index,1] <- rep(x$specs$model$params$probs[1,index_orig], x$specs$n.groups)
  }
  index <- which(is.na(x$specs$model$params$taus[1,]))
  if (x$specs$n.groups>1) for (i in 1:(x$specs$n.groups-1)) {index <- c(index, index+length(rowNms))}
  whole[index,2] <- x$diags$mu_tau_minus[,3]
  index <- which(is.na(x$specs$model$params$taus[2,]))
  if (x$specs$n.groups>1) for (i in 1:(x$specs$n.groups-1)) {index <- c(index, index+length(rowNms))}
  whole[index,3] <- x$diags$mu_tau_plus[,3]
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

