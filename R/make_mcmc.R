
# function that gives variablenames/columnnames
labelnames <- function(data_info) {
  
  Gr <- data_info$Ngroups
  
  probs_string <- data_info$probs_string
  minus_string <- data_info$minus_string
  plus_string <- data_info$plus_string
  
  P <- length(probs_string)
  Pm <- length(minus_string)
  Pp <- length(plus_string)

  Ptot <- P+Pm+Pp
  Plam <- Ptot-P

  S <- data_info$Nsubj
  R <- data_info$Nresps
  
  
  # lable variables
  sig <- c("", "minus", "plus")
  Prime <- c("alpha_prime", "beta_prime")
  all_string <- c(probs_string, minus_string, plus_string)
  
  
  label <- c()
  index <- 0
  ## group mean probs
  for (gr in 1:Gr) {
    for (p in 1:P) label[index+p] <- paste0("mu_alpha_", probs_string[p], if (Gr>1) paste0("[", gr-1, "]"))
    index <- index + P
  }
  ## group mean rates
  for (gr in 1:Gr) {
    for (p in 1:Pm) label[index+p] <- paste0("exp_mu_beta_", minus_string[p], "_", sig[2], if (Gr>1) paste0("[", gr-1, "]"))
    index <- index + Pm
    for (p in 1:Pp) label[index+p] <- paste0("exp_mu_beta_", plus_string[p], "_", sig[3], if (Gr>1) paste0("[", gr-1, "]"))
    index <- index + Pp
  }
  # upper triangular of alpha' and beta' (SD and correlation)
  for (i in 1:( Ptot )) {
    for (j in i:( Ptot )) {
      label[index+j] <- paste0(ifelse(i<=P, Prime[1], Prime[2]), "_", all_string[i], 
                               ifelse(i<=P, sig[1], paste0("_", ifelse(i<=P+Pm, sig[2], sig[3]))), "__",
                               ifelse(j<=P, Prime[1], Prime[2]), "_", all_string[j],
                               ifelse(j<=P, sig[1], paste0("_", ifelse(j<=P+Pm, sig[2], sig[3]))))
    }
    index <- index + (Ptot-i)
  }
  index <- index + Ptot
  ## alpha'
  for (s in 1:S) label <- c(label, paste0(Prime[1], "_", probs_string, "[", s-1, "]"))
  index <- index + (S*P)
  ## beta'
  beta_ind <- (index+1):(index+S*Plam)
  for (s in 1:S) label <- c(label, paste0(Prime[2], "_", all_string[(P+1):Ptot], "_", 
                                          c(rep(sig[2], Pm), rep(sig[3], Pp)), "[", s-1, "]"))
  index <- index + (S*Plam)
  ##  group mean encoding & response execution times
  for (gr in 1:Gr) {
    for (r in 1:R) label[index+r] <- paste0("mu_gamma", if(R>1) paste0("_R", r-1), if (Gr>1) paste0("[", gr-1, "]"))
    index <- index + R
  }
  ## omega^2
  label[index+1] <- paste("omega_square")
  index <- index + 1
  ## upper triangular of gamma' (SD and correlation)
  for (i in 1:R) {
    for (j in i:R) label[index+j] <- paste0("gamma_prime_", "R", i-1, "__", "gamma_prime_", "R", j-1)
    index <- index + (R-i)
  }
  index <- index + R
  ## gamma'
  for (s in 1:S) label <- c(label, paste0("gamma_prime", if(R>1) paste0("_R", (1:2)-1), "[", s-1, "]"))
  index <- index + (S*R)
  ## sigma square per participant
  for (s in 1:S) label[index+s] <- paste0("sigma_square[", s-1, "]")
  index <- index + S
  
  label <- c(label, "log_lik")
  index = index + 1
  
  return(list(label=label, beta_ind=beta_ind))

}

labelnames_d <- function(data_info) {
  
  Gr <- data_info$Ngroups
  
  thresh_string <- data_info$thresh_string
  drift_string <- data_info$drift_string
  start_string <- data_info$start_string
  
  Pt <- length(thresh_string)
  Pd <- length(drift_string)
  Ps <- length(start_string)
  
  Ptot <- Pt+Pd+Ps
  
  S <- data_info$Nsubj
  R <- data_info$Nresps
  
  
  # lable variables
  sig <- c("a", "nu", "omega")
  Prime <- c("a_prime", "nu_prime", "omega_prime")
  all_string <- c(thresh_string, drift_string, start_string)
  
  
  label <- c()
  index <- 0
  
  ## prepare group labels
  if (Gr > 1) {
    gr_labels_t <- paste0("[", rep(0:(Gr-1), each = Pt), "]")
    gr_labels_d <- paste0("[", rep(0:(Gr-1), each = Pd), "]")
    gr_labels_s <- paste0("[", rep(0:(Gr-1), each = Ps), "]")
    gr_labels_r <- paste0("[", rep(0:(Gr-1), each = R), "]")
  } else {
    gr_labels_t <- rep("", Pt)
    gr_labels_d <- rep("", Pd)
    gr_labels_s <- rep("", Ps)
    gr_labels_r <- rep("", R)
  }
  
  ## prepare subject labels
  subj_labels <- paste0("[", 0:(S-1), "]")
  
  ## group mean thresholds
  if (Pt) label <- paste0("mu_a_", rep(thresh_string, Gr), gr_labels_t)
  
  ## group mean drift rates
  if (Pd) label <- c(label, paste0("mu_nu_", rep(drift_string, Gr), gr_labels_d))
  
  ## group mean starting points
  if (Ps) label <- c(label, paste0("mu_omega_", rep(start_string, Gr), gr_labels_s))
  
  ## person-level process params
  label <- c(label, as.vector(sapply(1:S, FUN = function(x) {
    c(if (Pt) {paste0(Prime[1], "_", thresh_string, subj_labels[x])},
      if (Pd) {paste0(Prime[2], "_", drift_string, subj_labels[x])},
      if (Ps) {paste0(Prime[3], "_", start_string, subj_labels[x])})
  }, simplify = TRUE)))
  
  ## group mean motor times
  label <- c(label, paste0("mu_gamma_R", if(R>1) rep(paste0("_R", 0:(R-1)), Gr), gr_labels_r))
  
  ## person-level motor params
  label <- c(label, as.vector(sapply(1:S, FUN = function(x) {
    paste0("gamma_prime_R", if (R > 1) {0:(R-1)}, subj_labels[x])
  }, simplify = TRUE)))
  
  ## scale params
  label <- c(label, paste0("log(sigma", subj_labels, ")"))
  
  ## SIGMA
  s1 <- rep(unlist(sapply(1:3, function(i) rep(Prime[i], c(Pt, Pd, Ps)[i]))), Ptot)
  s2 <- rep(all_string, Ptot)
  s3 <- unlist(sapply(1:3, function(i) rep(Prime[i], (Ptot*c(Pt, Pd, Ps))[i])))
  s4 <- rep(all_string, each = Ptot)
  label <- c(label,
             matrix(paste0("s_", s1, "_", s2, "__", s3, "_", s4), ncol = Ptot)[lower.tri(x = matrix(NA, ncol = Ptot, nrow = Ptot), diag = TRUE)])
  
  ## GAMMA
  g2 <- rep(0:(R-1), R)
  g4 <- rep(0:(R-1), each = R)
  label <- c(label,
             matrix(paste0("g_gamma_prime_", g2, "__", "gamma_prime_", g4), ncol = R)[lower.tri(x = matrix(NA, ncol = R, nrow = R), diag = TRUE)])
  
  ## OMEGA2
  label <- c(label, "Omega2")
  
  # label <- c(label, "log_lik")
  
  return(label)
  
}

labelnames_keep <- function(data_info) {
  
  Gr <- data_info$Ngroups
  
  probs_string <- data_info$probs_string
  minus_string <- data_info$minus_string
  plus_string <- data_info$plus_string
  
  P <- length(probs_string)
  Pm <- length(minus_string)
  Pp <- length(plus_string)
  
  Ptot <- P+Pm+Pp
  Plam <- Ptot-P
  
  S <- data_info$Nsubj
  R <- data_info$Nresps
  
  
  # lable variables
  sig <- c("", "minus", "plus")
  Prime <- c("alpha_prime", "beta_prime")
  all_string <- c(probs_string, minus_string, plus_string)
  
  group_flag <- FALSE
  if (exists("group", data_info$transformation)) {
    group_flag <- TRUE
  }
  subj_flag <- FALSE
  if (exists("subj", data_info$transformation)) {
    subj_flag <- TRUE
  }
  
  
  label <- c()
  index <- 0
  ## group mean probs
  for (gr in 1:Gr) {
    if (group_flag) gr_label <- data_info$transformation$group$old[gr] else gr_label <- gr-1
    for (p in 1:P) label[index+p] <- paste0("mu_alpha_", probs_string[p], if (Gr>1) paste0("[", gr_label, "]"))
    index <- index + P
  }
  ## group mean rates
  for (gr in 1:Gr) {
    if (group_flag) gr_label <- data_info$transformation$group$old[gr] else gr_label <- gr-1
    for (p in 1:Pm) label[index+p] <- paste0("exp_mu_beta_", minus_string[p], "_", sig[2], if (Gr>1) paste0("[", gr_label, "]"))
    index <- index + Pm
    for (p in 1:Pp) label[index+p] <- paste0("exp_mu_beta_", plus_string[p], "_", sig[3], if (Gr>1) paste0("[", gr_label, "]"))
    index <- index + Pp
  }
  # upper triangular of alpha' and beta' (SD and correlation)
  for (i in 1:( Ptot )) {
    for (j in i:( Ptot )) {
      label[index+j] <- paste0(ifelse(i<=P, Prime[1], Prime[2]), "_", all_string[i], 
                               ifelse(i<=P, sig[1], paste0("_", ifelse(i<=P+Pm, sig[2], sig[3]))), "__",
                               ifelse(j<=P, Prime[1], Prime[2]), "_", all_string[j],
                               ifelse(j<=P, sig[1], paste0("_", ifelse(j<=P+Pm, sig[2], sig[3]))))
    }
    index <- index + (Ptot-i)
  }
  index <- index + Ptot
  ## alpha'
  for (s in 1:S) {
    if (subj_flag) s_label <- data_info$transformation$subj$old[s] else s_label <- s-1
    label <- c(label, paste0(Prime[1], "_", probs_string, "[", s_label, "]"))
  }
  index <- index + (S*P)
  ## beta'
  beta_ind <- (index+1):(index+S*Plam)
  for (s in 1:S) {
    if (subj_flag) s_label <- data_info$transformation$subj$old[s] else s_label <- s-1
    label <- c(label, paste0(Prime[2], "_", all_string[(P+1):Ptot], "_", 
                             c(rep(sig[2], Pm), rep(sig[3], Pp)), "[", s_label, "]"))
  }
  index <- index + (S*Plam)
  ##  group mean encoding & response execution times
  for (gr in 1:Gr) {
    if (group_flag) gr_label <- data_info$transformation$group$old[gr] else gr_label <- gr-1
    for (r in 1:R) label[index+r] <- paste0("mu_gamma", if(R>1) paste0("_R", r-1), if (Gr>1) paste0("[", gr_label, "]"))
    index <- index + R
  }
  ## omega^2
  label[index+1] <- paste("omega_square")
  index <- index + 1
  ## upper triangular of gamma' (SD and correlation)
  for (i in 1:R) {
    for (j in i:R) label[index+j] <- paste0("gamma_prime_", "R", i-1, "__", "gamma_prime_", "R", j-1)
    index <- index + (R-i)
  }
  index <- index + R
  ## gamma'
  for (s in 1:S) {
    if (subj_flag) s_label <- data_info$transformation$subj$old[s] else s_label <- s-1
    label <- c(label, paste0("gamma_prime", if(R>1) paste0("_R", (1:2)-1), "[", s_label, "]"))
  }
  index <- index + (S*R)
  ## sigma square per participant
  for (s in 1:S) {
    if (subj_flag) s_label <- data_info$transformation$subj$old[s] else s_label <- s-1
    label[index+s] <- paste0("sigma_square[", s_label, "]")
  }
  index <- index + S
  
  label <- c(label, "log_lik")
  index = index + 1
  
  return(list(label=label, beta_ind=beta_ind))
  
}


labelnames_keep_d <- function(data_info) {
  
  Gr <- data_info$Ngroups
  
  thresh_string <- data_info$thresh_string
  drift_string <- data_info$drift_string
  start_string <- data_info$start_string
  
  Pt <- length(thresh_string)
  Pd <- length(drift_string)
  Ps <- length(start_string)
  
  Ptot <- Pt+Pd+Ps
  
  S <- data_info$Nsubj
  R <- data_info$Nresps
  
  
  # lable variables
  sig <- c("a", "nu", "omega")
  Prime <- c("a_prime", "nu_prime", "omega_prime")
  all_string <- c(thresh_string, drift_string, start_string)
  
  group_flag <- FALSE
  if (exists("group", data_info$transformation)) {
    group_flag <- TRUE
  }
  subj_flag <- FALSE
  if (exists("subj", data_info$transformation)) {
    subj_flag <- TRUE
  }
  
  
  label <- c()
  index <- 0
  
  ## prepare group labels
  if (Gr > 1 && group_flag) {
    gr_labels_t <- paste0("[", rep(data_info$transformation$group$old, each = Pt), "]")
    gr_labels_d <- paste0("[", rep(data_info$transformation$group$old, each = Pd), "]")
    gr_labels_s <- paste0("[", rep(data_info$transformation$group$old, each = Ps), "]")
    gr_labels_r <- paste0("[", rep(data_info$transformation$group$old, each = R), "]")
  } else if (Gr > 1 && !group_flag) {
    gr_labels_t <- paste0("[", rep(0:(Gr-1), each = Pt), "]")
    gr_labels_d <- paste0("[", rep(0:(Gr-1), each = Pd), "]")
    gr_labels_s <- paste0("[", rep(0:(Gr-1), each = Ps), "]")
    gr_labels_r <- paste0("[", rep(0:(Gr-1), each = R), "]")
  } else if (Gr <= 1) {
    gr_labels_t <- rep("", Pt)
    gr_labels_d <- rep("", Pd)
    gr_labels_s <- rep("", Ps)
    gr_labels_r <- rep("", R)
  }
  
  ## prepare subject labels
  if (group_flag) {
    subj_labels <- paste0("[", data_info$transformation$subj$old, "]")
  } else {
    subj_labels <- paste0("[", 0:(S-1), "]")
  }
  
  ## group mean thresholds
  if (Pt) label <- paste0("mu_a_", rep(thresh_string, Gr), gr_labels_t)
  
  ## group mean drift rates
  if (Pd) label <- c(label, paste0("mu_nu_", rep(drift_string, Gr), gr_labels_d))
  
  ## group mean starting points
  if (Ps) label <- c(label, paste0("mu_omega_", rep(start_string, Gr), gr_labels_s))
  
  ## person-level process params
  label <- c(label, as.vector(sapply(1:S, FUN = function(x) {
    c(if (Pt) {paste0(Prime[1], "_", thresh_string, subj_labels[x])},
      if (Pd) {paste0(Prime[2], "_", drift_string, subj_labels[x])},
      if (Ps) {paste0(Prime[3], "_", start_string, subj_labels[x])})
  }, simplify = TRUE)))
  
  ## group mean motor times
  label <- c(label, paste0("mu_gamma_R", if(R>1) rep(paste0("_R", 0:(R-1)), Gr), gr_labels_r))
  
  ## person-level motor params
  label <- c(label, as.vector(sapply(1:S, FUN = function(x) {
    paste0("gamma_prime_R", if (R > 1) {0:(R-1)}, subj_labels[x])
  }, simplify = TRUE)))
  
  ## scale params
  label <- c(label, paste0("log(sigma", subj_labels, ")"))
  
  ## SIGMA
  s1 <- rep(unlist(sapply(1:3, function(i) rep(Prime[i], c(Pt, Pd, Ps)[i]))), Ptot)
  s2 <- rep(all_string, Ptot)
  s3 <- unlist(sapply(1:3, function(i) rep(Prime[i], (Ptot*c(Pt, Pd, Ps))[i])))
  s4 <- rep(all_string, each = Ptot)
  label <- c(label,
             matrix(paste0("s_", s1, "_", s2, "__", s3, "_", s4), ncol = Ptot)[lower.tri(x = matrix(NA, ncol = Ptot, nrow = Ptot), diag = TRUE)])
  
  ## GAMMA
  g2 <- rep(0:(R-1), R)
  g4 <- rep(0:(R-1), each = R)
  label <- c(label,
             matrix(paste0("g_gamma_prime_", g2, "__", "gamma_prime_", g4), ncol = R)[lower.tri(x = matrix(NA, ncol = R, nrow = R), diag = TRUE)])
  
  ## OMEGA2
  label <- c(label, "Omega2")
  
  # label <- c(label, "log_lik")
  
  return(label)
  
}




#' @importFrom data.table as.data.table
#' @importFrom coda as.mcmc
make_mcmc_list <- function(file, infofile, Nchains, Nsamples, data_info, keep) {
  
  # read text file with chains
  temp <- c()
  if (is.character(file)) {
    temp <- as.vector(read.table(file = file, header = FALSE, nrows = 1))
    dt <- fread(file=file,skip=1)
  } else if (is.data.frame(file) || is.matrix(file)) {
    temp <- dim(file)
    dt <- as.data.table(file)
  }
  rm(file); gc()

  
  # specify parameters used in MCMC
  start <- 1
  end <- Nsamples
  npar <- temp[2]
  
  
  # generate MCMC-list for coda
  vec <- vector("list", Nchains)
  for (i in 1:Nchains) {
    vec[[i]] <- dt[(i-1)*end + start:end]
    vec[[i]] <- as.mcmc(vec[[i]], start = start, end = end, thin = 1)
  }
  rm(dt); gc()
  samples <- as.mcmc.list(vec, start = start, end = end, thin = 1)
  rm(vec); gc()
  
  
  # name chain columns
  # if(npar != length(varnames(samples))) stop("Number of columns in chains is wrong.")       ###---###
  # varstring <- system(paste("gawk 'END {print}' ", infofile), intern = TRUE)
  # var_char <- strsplit(varstring, "\\s+")[[1]]
  
  
  # label the samples
  if (keep) {
    label_list <- labelnames_keep(data_info = data_info)
    varnames(samples) <- label_list$label
  } else {
    label_list <- labelnames(data_info = data_info)
    varnames(samples) <- label_list$label
  }  
  
  # change beta_prime to log-scale
  for (nc in 1:Nchains) {
    for (ind in label_list$beta_ind) {
	    samples[[nc]][,ind] <- log(samples[[nc]][,ind])
  	}
  }
  gc()
  return(samples)
  
}



#' @importFrom data.table as.data.table fread setnames
#' @importFrom coda as.mcmc
make_mcmc_list_d <- function(file, infofile, Nchains, Nsamples, data_info, keep) {
  
  # label the samples
  dt_colnames <- NULL
  if (keep) {
    dt_colnames <- labelnames_keep_d(data_info = data_info)
  } else {
    dt_colnames <- labelnames_d(data_info = data_info)
  }
  gc()
  
  
  # read text file with chains
  temp <- c()
  if (is.character(file)) {
    temp <- as.vector(read.table(file = file, header = FALSE, nrows = 1))
    dt <- fread(file=file, skip=1, header = FALSE, col.names = dt_colnames)
  } else if (is.data.frame(file) || is.matrix(file)) {
    temp <- dim(file)
    dt <- as.data.table(file)
    setnames(dt, dt_colnames)
  }
  rm(file); gc()
  
  
  # read text file with chains
  # temp <- as.vector(read.table(file = file, header = FALSE, nrows = 1))
  # dt <- fread(file = file, skip = 1, header = FALSE, col.names = dt_colnames, blank.lines.skip = TRUE)
  # gc()
  
  
  # specify parameters used in MCMC
  factor <- nrow(dt)/Nsamples/Nchains
  start <- 1 + (factor-1)*Nsamples
  end <- Nsamples + (factor-1)*Nsamples
  npar <- temp[2]
  
  
  # generate MCMC-list for coda
  vec <- vector("list", Nchains)
  for (i in 1:Nchains) {
    vec[[i]] <- dt[(i-1)*end + start:end]
    vec[[i]] <- as.mcmc(vec[[i]], start = start, end = end, thin = 1)
  }
  rm(dt); gc()
  samples <- as.mcmc.list(vec, start = start, end = end, thin = 1)
  rm(vec); gc()
  
  
  return(samples)
  
}

