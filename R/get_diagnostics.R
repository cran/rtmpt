
get_diags <- function(diag_file, data_info, keep) {
  
  # needed variables
  diag_names <- c("MUs per group", "rho minus per group", "rho plus per group",
                  "SIG", "Restpars", "RSIG", "Daten", "DIC1", "frequencies",
                  "latencies", "group-tests mu", "group-tests pho",
                  "group-tests residual", "within-group residuals")
  diag_lines <- readLines(con = diag_file)

  Nall <- data_info$Nprobs+data_info$Nminus+data_info$Nplus
  
  # flags
  group_flag <- FALSE
  subj_flag <- FALSE
  if (keep) {
    if (exists("group", data_info$transformation)) {
      group_flag <- TRUE
    }
    if (exists("subj", data_info$transformation)) {
      subj_flag <- TRUE
    }
  }
  if (group_flag) gr_labels <- data_info$transformation$group$old[1:(data_info$Ngroups)] else gr_labels <- (1:(data_info$Ngroups))-1 
  
  
  # list to return
  diag_list <- list()
  
  
  # some functions
  get_indices <- function(pattern, N, add = 0) {
    start <- grep(pattern = pattern, x = diag_lines)[1]
    indices <- (start+1+add):(start+N+add)
    return(indices)
  }
  save_HDI2list <- function(indices, names, skip) {
    dfx <- as.data.frame(matrix(data = NA, ncol = length(indices), nrow = 5))
    colnames(dfx) <- names
    rownames(dfx) <- c("low_HDI_99", "low_HDI_95", "median", "upp_HDI_95", "upp_HDI_99")
    for (i in 1:length(indices)) {
      value_char <- strsplit(diag_lines[indices[i]], "\\s+")[[1]]
      dfx[,i] <- as.numeric(value_char[value_char != ""][(1:5)+skip])
    }
    return(dfx)
  }
  save_test2list <- function(indices, names) {
    dfx <- as.data.frame(matrix(data = NA, ncol = length(indices), nrow = 5))
    colnames(dfx) <- names
    rownames(dfx) <- c("mean.1", "mean.2", "p", "low_HDI_95", "upp_HDI_95")
    value_char <- c(strsplit(diag_lines[indices[1]+1], "\\s+")[[1]], strsplit(diag_lines[indices[1]+3], "\\s+")[[1]])
    for (i in 1:length(indices)) dfx[,i] <- as.numeric(value_char[value_char != ""])
    return(dfx)
  }
  save_test3list <- function(indices, names) {
    dfx <- as.data.frame(matrix(data = NA, ncol = length(indices), nrow = 5))
    colnames(dfx) <- names
    rownames(dfx) <- c("observed", "expected", "ppp", "low_HDI_95", "upp_HDI_95")
    value_char <- c(strsplit(diag_lines[indices[1]+1], "\\s+")[[1]], strsplit(diag_lines[indices[1]+3], "\\s+")[[1]])
    for (i in 1:length(indices)) dfx[,i] <- as.numeric(value_char[value_char != ""])
    return(dfx)
  }

  
  # prepare indices
  probs_ind <- get_indices(pattern = diag_names[1], N = data_info$Ngroups*data_info$Nprobs)
  lambda_minus_ind <- get_indices(pattern = diag_names[2], N = data_info$Ngroups*data_info$Nminus)
  lambda_plus_ind <- get_indices(pattern = diag_names[3], N = data_info$Ngroups*data_info$Nplus)
  SIGMA_ind <- get_indices(pattern = diag_names[4], N = Nall*(Nall+1)/2)
  resps_ind <- get_indices(pattern = diag_names[5], N = data_info$Ngroups*data_info$Nresps)
  omega2_ind <- get_indices(pattern = diag_names[5], N = 1, add = data_info$Ngroups*data_info$Nresps)
  TAU_ind <- get_indices(pattern = diag_names[6], N = data_info$Nresps*(data_info$Nresps+1)/2)
  dat_ind <- get_indices(pattern = diag_names[7], N = 1)
  DIC_ind <- get_indices(pattern = diag_names[8], N = 2)
  freq_ind <- get_indices(pattern = diag_names[9], N = 1, add = -1)
  time_ind <- get_indices(pattern = diag_names[10], N = 1, add = -1)

  # read and save
  diag_list$mu_probs <- t(save_HDI2list(indices = probs_ind, names = paste0(data_info$probs_string, "[",rep(gr_labels, each=data_info$Nprobs), "]"), skip = 1))
  diag_list$mu_tau_minus <- t(save_HDI2list(indices = lambda_minus_ind, names = paste0(data_info$minus_string, "[",rep(gr_labels, each=data_info$Nminus), "]"), skip = 1))
  diag_list$mu_tau_plus <- t(save_HDI2list(indices = lambda_plus_ind, names = paste0(data_info$plus_string, "[",rep(gr_labels, each=data_info$Nplus), "]"), skip = 1))

  process_chars <- c(paste0("alpha_prm_", data_info$probs_string), 
                     paste0("beta_prm_", data_info$minus_string, "_minus"), 
                     paste0("beta_prm_", data_info$plus_string, "_plus"))
  SIGMA_names <- rep(NA, Nall*(Nall+1)/2)
  ind = 0
  for (i in 1:Nall) {
    for (j in i:Nall) {
      ind = ind + 1
      SIGMA_names[ind] <- paste0(process_chars[i], "..", process_chars[j])
    }
  }
  diag_list$SD_CORR_Proc <- t(save_HDI2list(indices = SIGMA_ind, names = SIGMA_names, skip = 2))
  diag_list$mu_gamma <- t(save_HDI2list(indices = resps_ind, names = paste0(paste0("R", 0:(data_info$Nresps-1)), "[",rep(gr_labels, each=data_info$Nresps), "]"), skip = 0))
  diag_list$omega2 <- t(save_HDI2list(indices = omega2_ind, names = "omega2", skip = 0))

  resps_chars <- paste0("gamma_prm_", 0:(data_info$Nresps-1))
  DELTA_names <- rep(NA, data_info$Nresps*(data_info$Nresps+1)/2)
  ind = 0
  for (i in 1:data_info$Nresps) {
    for (j in i:data_info$Nresps) {
      ind = ind + 1
      DELTA_names[ind] <- paste0(resps_chars[i], "..", resps_chars[j])
    }
  }
  diag_list$SD_CORR_Motor <- t(save_HDI2list(indices = TAU_ind, names = DELTA_names, skip = 2))
  diag_list$RT <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 3))
  colnames(diag_list$RT) <- c("mean", "variance", "residual_var")
  val_char <- strsplit(diag_lines[dat_ind[1]], "\\s+")[[1]]
  diag_list$RT[1,] <- as.numeric(val_char[val_char != ""])
  # diag_list$CORR <- save_HDI2list(indices = dat_ind+(1:2), names = c("beta_prm", "alphabeta_prm"), skip = 1)

  val_char <- strsplit(diag_lines[DIC_ind[1]], "\\s+")[[1]]
  DIC <- as.numeric(val_char[val_char != ""])
  diag_list$DIC1 <- data.frame(DIC=DIC[1],pd = DIC[3])
  val_char <- strsplit(diag_lines[DIC_ind[2]], "\\s+")[[1]]
  diag_list$DIC2 <- data.frame(DIC=DIC[2],pv = as.numeric(val_char[val_char != ""]))
  diag_list$PostPredCheck_Frequencies <- t(save_test3list(indices = freq_ind, names = "freq"))
  diag_list$PostPredCheck_Latencies <- t(save_test3list(indices = time_ind, names = "time"))

    
  # group-wise
  if (data_info$Ngroups > 1) {
    # prepare indices
    mu_ind <- which(diag_lines == diag_names[11])
    rho_ind <- which(diag_lines == diag_names[12])
    residual_ind <- which(diag_lines == diag_names[13])
    wi_residuals_ind <- which(diag_lines == diag_names[14])
    
    
    # read and save
    try(diag_list$group.test.probs <- save_test2list(indices = mu_ind, names = data_info$probs_string), silent = TRUE)
    try(diag_list$group.test.tau <- save_test2list(indices = rho_ind, names = c(paste0(data_info$minus_string, "_minus"), paste0(data_info$plus_string, "_plus"))), silent = TRUE)
    try(diag_list$group.test.gamma <- save_test2list(indices = residual_ind, names = paste0("R", 0:(data_info$Nresps-1))), silent = TRUE)
    try(group.test.within.residuals <- save_test2list(indices = wi_residuals_ind, names = paste0("residual[", 0:(data_info$Ngroups-1), "]")), silent = TRUE)
  }
  

  return(diag_list)
  
}





get_diags_d <- function(diag_file, data_info, keep, DIC) {
  
  # needed variables
  diag_names <- c("Thresholds per group", "Drift rates per group", "Start point per group",
                  "SIG", "Mean motor/encoding", "Sig motor/encoding", "Residual variance", "RT: mean, variance, residual variance", 
                  "frequencies", "latencies", "DIC2"
                  #, "group-tests mu", "group-tests pho", "group-tests residual", "within-group residuals"
                  )
  diag_lines <- readLines(con = diag_file)
  
  Nall <- data_info$Nthreshold + data_info$Ndriftrate + data_info$Nstartpoint
  
  # flags
  group_flag <- FALSE
  subj_flag <- FALSE
  if (keep) {
    if (exists("group", data_info$transformation)) {
      group_flag <- TRUE
    }
    if (exists("subj", data_info$transformation)) {
      subj_flag <- TRUE
    }
  }
  if (group_flag) gr_labels <- data_info$transformation$group$old[1:(data_info$Ngroups)] else gr_labels <- (1:(data_info$Ngroups))-1 
  
  
  # list to return
  diag_list <- list()
  
  
  # some functions
  get_indices <- function(pattern, N, add = 0) {
    start <- grep(pattern = pattern, x = diag_lines)[1]
    indices <- (start+1+add):(start+N+add)
    return(indices)
  }
  save_HDI2list <- function(indices, names, skip) {
    dfx <- as.data.frame(matrix(data = NA, ncol = length(indices), nrow = 5))
    colnames(dfx) <- names
    rownames(dfx) <- c("low_HDI_99", "low_HDI_95", "median", "upp_HDI_95", "upp_HDI_99")
    for (i in 1:length(indices)) {
      value_char <- strsplit(diag_lines[indices[i]], "\\s+")[[1]]
      dfx[,i] <- as.numeric(value_char[value_char != ""][(1:5)+skip])
    }
    return(dfx)
  }
  save_test2list <- function(indices, names) {
    dfx <- as.data.frame(matrix(data = NA, ncol = length(indices), nrow = 5))
    colnames(dfx) <- names
    rownames(dfx) <- c("mean.1", "mean.2", "p", "low_HDI_95", "upp_HDI_95")
    value_char <- c(strsplit(diag_lines[indices[1]+1], "\\s+")[[1]], strsplit(diag_lines[indices[1]+3], "\\s+")[[1]])
    for (i in 1:length(indices)) dfx[,i] <- as.numeric(value_char[value_char != ""])
    return(dfx)
  }
  save_test3list <- function(indices, names) {
    dfx <- as.data.frame(matrix(data = NA, ncol = length(indices), nrow = 5))
    colnames(dfx) <- names
    rownames(dfx) <- c("observed", "expected", "ppp", "low_HDI_95", "upp_HDI_95")
    value_char <- c(strsplit(diag_lines[indices[1]+1], "\\s+")[[1]], strsplit(diag_lines[indices[1]+3], "\\s+")[[1]])
    for (i in 1:length(indices)) dfx[,i] <- as.numeric(value_char[value_char != ""])
    return(dfx)
  }
  
  
  # prepare indices
  thresh_ind <- get_indices(pattern = diag_names[1], N = data_info$Ngroups*data_info$Nthreshold)
  drift_ind <- get_indices(pattern = diag_names[2], N = data_info$Ngroups*data_info$Ndriftrate)
  start_ind <- get_indices(pattern = diag_names[3], N = data_info$Ngroups*data_info$Nstartpoint)
  SIGMA_ind <- get_indices(pattern = diag_names[4], N = Nall*(Nall+1)/2)
  resps_ind <- get_indices(pattern = diag_names[5], N = data_info$Ngroups*data_info$Nresps)
  GAMMA_ind <- get_indices(pattern = diag_names[6], N = data_info$Nresps*(data_info$Nresps+1)/2)
  Omega2_ind <- get_indices(pattern = diag_names[7], N = 1)
  dat_ind <- get_indices(pattern = diag_names[8], N = 1)
  freq_ind <- get_indices(pattern = diag_names[9], N = 1, add = -1)
  time_ind <- get_indices(pattern = diag_names[10], N = 1, add = -1)
  if (DIC) DIC_ind <- get_indices(pattern = diag_names[11], N = 2)
  
  # read and save
  if (data_info$Nthreshold) diag_list$mu_a <- t(save_HDI2list(indices = thresh_ind, names = paste0(data_info$thresh_string, "[",rep(gr_labels, each=data_info$Nthreshold), "]"), skip = 2))
  if (data_info$Ndriftrate) diag_list$mu_nu <- t(save_HDI2list(indices = drift_ind, names = paste0(data_info$drift_string, "[",rep(gr_labels, each=data_info$Ndriftrate), "]"), skip = 2))
  if (data_info$Nstartpoint) diag_list$mu_omega <- t(save_HDI2list(indices = start_ind, names = paste0(data_info$start_string, "[",rep(gr_labels, each=data_info$Nstartpoint), "]"), skip = 2))
  
  process_chars <- c(paste0("a_prm_", data_info$thresh_string), 
                     paste0("nu_prm_", data_info$drift_string), 
                     paste0("omega_prm_", data_info$start_string))
  SIGMA_names <- rep(NA, Nall*(Nall+1)/2)
  ind = 0
  for (i in 1:Nall) {
    for (j in i:Nall) {
      ind = ind + 1
      SIGMA_names[ind] <- paste0(process_chars[i], "..", process_chars[j])
    }
  }
  
  #### NOT SURE WHETHER IT REALLY IS SD AND CORR !!!!!!!!!!!!
  
  diag_list$VAR_COVAR_Proc <- t(save_HDI2list(indices = SIGMA_ind, names = SIGMA_names, skip = 2))
  diag_list$mu_gamma <- t(save_HDI2list(indices = resps_ind, names = paste0(paste0("R", 0:(data_info$Nresps-1)), "[",rep(gr_labels, each=data_info$Nresps), "]"), skip = 0))
  diag_list$Omega2 <- t(save_HDI2list(indices = Omega2_ind, names = "Omega2", skip = 0))
  
  resps_chars <- paste0("gamma_prm_", 0:(data_info$Nresps-1))
  DELTA_names <- rep(NA, data_info$Nresps*(data_info$Nresps+1)/2)
  ind = 0
  for (i in 1:data_info$Nresps) {
    for (j in i:data_info$Nresps) {
      ind = ind + 1
      DELTA_names[ind] <- paste0(resps_chars[i], "..", resps_chars[j])
    }
  }
  diag_list$VAR_COVAR_Motor <- t(save_HDI2list(indices = GAMMA_ind, names = DELTA_names, skip = 2))
  diag_list$RT <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 3))
  colnames(diag_list$RT) <- c("mean", "variance", "residual_var")
  val_char <- strsplit(diag_lines[dat_ind[1]], "\\s+")[[1]]
  diag_list$RT[1,] <- as.numeric(val_char[val_char != ""])
  # diag_list$CORR <- save_HDI2list(indices = dat_ind+(1:2), names = c("beta_prm", "alphabeta_prm"), skip = 1)
  
  if (DIC) {
    val_char <- strsplit(diag_lines[DIC_ind[1]], "\\s+")[[1]]
    DIC_vals <- as.numeric(val_char[val_char != ""])
    # diag_list$DIC1 <- data.frame(DIC=DIC[1],pd = DIC[3])
    # val_char <- strsplit(diag_lines[DIC_ind[2]], "\\s+")[[1]]
    diag_list$DIC2 <- data.frame(DIC=DIC_vals[1],pd = DIC_vals[2])
  }
  diag_list$PostPredCheck_Frequencies <- t(save_test3list(indices = freq_ind, names = "freq"))
  diag_list$PostPredCheck_Latencies <- t(save_test3list(indices = time_ind, names = "time"))
  
  
  # group-wise
  # if (data_info$Ngroups > 1) {
  #   # prepare indices
  #   mu_ind <- which(diag_lines == diag_names[11])
  #   rho_ind <- which(diag_lines == diag_names[12])
  #   residual_ind <- which(diag_lines == diag_names[13])
  #   wi_residuals_ind <- which(diag_lines == diag_names[14])
  #   
  #   
  #   # read and save
  #   try(diag_list$group.test.probs <- save_test2list(indices = mu_ind, names = data_info$probs_string), silent = TRUE)
  #   try(diag_list$group.test.tau <- save_test2list(indices = rho_ind, names = c(paste0(data_info$minus_string, "_minus"), paste0(data_info$plus_string, "_plus"))), silent = TRUE)
  #   try(diag_list$group.test.gamma <- save_test2list(indices = residual_ind, names = paste0("R", 0:(data_info$Nresps-1))), silent = TRUE)
  #   try(group.test.within.residuals <- save_test2list(indices = wi_residuals_ind, names = paste0("residual[", 0:(data_info$Ngroups-1), "]")), silent = TRUE)
  # }
  
  
  return(diag_list)
  
}
