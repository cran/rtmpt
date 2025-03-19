
#' Simulate data from an RT-MPT model
#' 
#' Simulate data from RT-MPT models using \code{ertmpt_model} objects. The difference to \code{\link{sim_ertmpt_data}} is that here only scalars are allowed. This makes it usable for 
#'   simulation-based calibration (SBC; Talts et al., 2018). You can specify the random seed, number of subjects, number of trials, and some
#'   parameters (same as \code{prior_params} from \code{\link{fit_ertmpt}}).
#'
#' @param model A list of the class \code{ertmpt_model}.
#' @param seed Random seed number.
#' @param n.subj <- Number of subjects.
#' @param n.trials <- Number of trials per tree.
#' @param params Named list of parameters from which the data will be generated. This must be the same named list as \code{prior_params} from 
#'   \code{\link{fit_ertmpt}} and has the same defaults. It is not recommended to use the defaults since they lead to many probabilities close or 
#'   equal to \code{0} and/or \code{1} and to RTs close or equal to \code{0}. Allowed parameters are:
#'   \itemize{
#'     \item \code{mean_of_exp_mu_beta}: This is the expected exponential rate (\code{E(exp(beta)) = E(lambda)}) and 
#'           \code{1/mean_of_exp_mu_beta} is the expected process time (\code{1/E(exp(beta)) = E(tau)}). The default
#'           mean is set to \code{10}, such that the expected process time is \code{0.1} seconds.
#'     \item \code{var_of_exp_mu_beta}: The group-specific variance of the exponential rates. Since
#'           \code{exp(mu_beta)} is Gamma distributed, the rate of the distribution is just mean divided by variance and
#'           the shape is the mean times the rate. The default is set to \code{100}.
#'     \item \code{mean_of_mu_gamma}: This is the expected \emph{mean parameter} of the encoding and response execution times,
#'           which follow a normal distribution truncated from below at zero, so \code{E(mu_gamma) < E(gamma)}. The default is \code{0}.
#'     \item \code{var_of_mu_gamma}: The group-specific variance of the \emph{mean parameter}. Its default is \code{10}.
#'     \item \code{mean_of_omega_sqr}: This is the expected residual variance (\code{E(omega^2)}). The default is \code{0.005}.
#'     \item \code{var_of_omega_sqr}: The variance of the residual variance (\code{Var(omega^2)}). The default is
#'           \code{0.01}. The default of the mean and variance is equivalent to a shape and rate of \code{0.0025} and 
#'           \code{0.5}, respectivly.
#'     \item \code{df_of_sigma_sqr}: degrees of freedom for the individual variance of the response executions. The
#'           individual variance follows a scaled inverse chi-squared distribution with \code{df_of_sigma_sqr} degrees of freedom and
#'           \code{omega^2} as scale. \code{2} is the default and it should be an integer.
#'     \item \code{sf_of_scale_matrix_SIGMA}: The original scaling matrix (S) of the (scaled) inverse Wishart distribution for the process 
#'           related parameters is an identity matrix \code{S=I}. \code{sf_of_scale_matrix_SIGMA} is a scaling factor, that scales this 
#'           matrix (\code{S=sf_of_scale_matrix_SIGMA*I}). Its default is \code{1}.
#'     \item \code{sf_of_scale_matrix_GAMMA}: The original scaling matrix (S) of the (scaled) inverse Wishart distribution for the encoding and
#'           motor execution parameters is an identity matrix \code{S=I}. \code{sf_of_scale_matrix_GAMMA} is a scaling factor that scales 
#'           this matrix (\code{S=sf_of_scale_matrix_GAMMA*I}). Its default is \code{1}.
#'     \item \code{prec_epsilon}: This is epsilon in the paper. It is the precision of mu_alpha and all xi (scaling parameter
#'           in the scaled inverse Wishart distribution). Its default is also \code{1}.
#'     \item \code{add_df_to_invWish}: If \code{P} is the number of parameters or rather the size of the scale matrix used in the (scaled)
#'           inverse Wishart distribution then \code{add_df_to_invWish} is the number of degrees of freedom that can be added to it. So
#'           \code{DF = P + add_df_to_invWish}. The default for \code{add_df_to_invWish} is \code{1}, such that the correlations are uniformly 
#'           distributed within \code{[-1, 1]}.
#'   }
#' @return A list of the class \code{ertmpt_sim} containing 
#'   \itemize{
#'     \item \code{data}: the data.frame with the simulated data,
#'     \item \code{gen_list}: a list containing lists of the group-level and subject-specific parameters for the process-related parameters and the motor-related
#'                            parameters, and the trial-specific probabilities, process-times, and motor-times,
#'     \item \code{specs}: some specifications like the model, seed number, etc.,
#'   }
#' @references
#'   Talts, S., Betancourt, M., Simpson, D., Vehtari, A., & Gelman, A. (2018). Validating Bayesian inference algorithms with simulation-based calibration. \emph{arXiv preprint arXiv:1804.06788}.
#' @examples 
#' ########################################################################################
#' # Detect-Guess variant of the Two-High Threshold model.
#' # The encoding and motor execution times are assumed to be different for each response.
#' ########################################################################################
#' 
#' mdl_2HTM <- "
#' # targets
#' do+(1-do)*g     ; 0
#' (1-do)*(1-g)    ; 1
#'
#' # lures
#' (1-dn)*g        ; 0
#' dn+(1-dn)*(1-g) ; 1
#' 
#' # do: detect old; dn: detect new; g: guess
#' "
#' 
#' model <- to_ertmpt_model(mdl_file = mdl_2HTM)
#' 
#' params <- list(mean_of_exp_mu_beta = 10, 
#'                var_of_exp_mu_beta = 10, 
#'                mean_of_mu_gamma = 0.5, 
#'                var_of_mu_gamma = 0.0025, 
#'                mean_of_omega_sqr = 0.005, 
#'                var_of_omega_sqr = 0.000025,
#'                df_of_sigma_sqr = 10, 
#'                sf_of_scale_matrix_SIGMA = 0.1, 
#'                sf_of_scale_matrix_GAMMA = 0.01, 
#'                prec_epsilon = 10,
#'                add_df_to_invWish = 5)
#' 
#' sim_dat <- rtmpt:::sim_ertmpt_data_SBC(model, seed = 123, n.subj = 40, 
#'                                       n.trials = 30, params = params)
#' 
#' @author Raphael Hartmann
#' @importFrom truncnorm rtruncnorm
#' @importFrom stats rgamma rnorm pnorm rexp runif
sim_ertmpt_data_SBC<- function(model, 
                              seed,
                              n.subj,
                              n.trials,
                              params = NULL) {
  
  
  # some controls ----
  model_elmnts <- c("lines", "params", "responses")
  if (!is.list(model)) stop("\"model\" must be a list.")
  if (!all(model_elmnts  %in% names(model))) stop("\"model\" must contain \"", model_elmnts[which(!(model_elmnts %in% names(model)))[1]], "\".")
  
  if (!is.numeric(seed)) stop("\"seed\" must be numeric.")
  
  if (n.subj < 2) stop("\"n.subj\" must be larger than or equal to two.")
  if (n.trials < 2) stop("\"n.trials\" must be larger than or equal to two.")
  if (n.trials < 30) warning("\"n.trials\" is recommended to be larger than 30.")
  
  if (!is.null(params) && !is.list(params)) stop("\"params\" must be a list.")
  
  
  
  # argument preparation ----
  Nproc <- length(model$params$probs[1,])
  Nprob <- sum(is.na(model$params$probs[1,]))
  Nminus <- sum(is.na(model$params$taus[1,]))
  Nplus <- sum(is.na(model$params$taus[2,]))
  Nsubj <- n.subj
  Nresp <- length(unique(model$responses$MAP))
  NTrials <- n.trials
  Executions <- model$responses$MAP+1
  some_const <- any(!is.na(model$params$probs[1,]))
  PositionProb <- ifelse(some_const, which(!is.na(model$params$probs[1,])), 0)
  ConstProb <- ifelse(some_const, model$params$probs[1, PositionProb], 0)
  SupprMinus <- ifelse(any(!is.na(model$params$taus[1,])), which(!is.na(model$params$taus[1,])), 0)
  SupprPlus <- ifelse(any(!is.na(model$params$taus[2,])), which(!is.na(model$params$taus[2,])), 0)
  
  
  
  # produce infofile ----
  mdl_txt <- gsub("\\\\", "/", tempfile(pattern = "model", tmpdir = tempdir(), fileext = ".txt"))
  #mdl_txt <- paste0(directory, "/model.txt")
  mdl_info <- gsub("\\\\", "/", tempfile(pattern = "model", tmpdir = tempdir(), fileext = ".info"))
  #mdl_info <- paste0(directory, "/model.info")
  infofile <- try(get_infofile(model, mdl_txt = mdl_txt, mdl_info = mdl_info))
  if(is(infofile, "try-error")) stop("problem with S4 routines in makin infofile.\n")
  
  
  # parameters ----
  param_names <- c("mean_of_exp_mu_beta", "var_of_exp_mu_beta", "mean_of_mu_gamma", "var_of_mu_gamma",
                   "mean_of_omega_sqr", "var_of_omega_sqr", "sf_of_scale_matrix_SIGMA", "sf_of_scale_matrix_GAMMA",
                   "df_of_sigma_sqr", "prec_epsilon", "add_df_to_invWish")
  if (is.null(params)) params <- list()
  if (!("mean_of_exp_mu_beta" %in% names(params))) params$mean_of_exp_mu_beta <- 10
  if (!("var_of_exp_mu_beta" %in% names(params))) params$var_of_exp_mu_beta <- 100
  if (!("mean_of_mu_gamma" %in% names(params))) params$mean_of_mu_gamma <- 0
  if (!("var_of_mu_gamma" %in% names(params))) params$var_of_mu_gamma <- 10
  if (!("mean_of_omega_sqr" %in% names(params))) params$mean_of_omega_sqr <- 0.005
  if (!("var_of_omega_sqr" %in% names(params))) params$var_of_omega_sqr <- 0.01
  if (!("df_of_sigma_sqr" %in% names(params))) params$df_of_sigma_sqr <- 2
  if (!("sf_of_scale_matrix_SIGMA" %in% names(params))) params$sf_of_scale_matrix_SIGMA <- 1
  if (!("sf_of_scale_matrix_GAMMA" %in% names(params))) params$sf_of_scale_matrix_GAMMA <- 1
  if (!("prec_epsilon" %in% names(params))) params$prec_epsilon <- 1
  if (!("add_df_to_invWish" %in% names(params))) params$add_df_to_invWish <- 1
  if (params$add_df_to_invWish %% 1 != 0) stop("\"add_df_to_invWish\" must be an integer!")
  if (!all(names(params) %in% param_names)) {
    ind <- which(!(names(params) %in% param_names))
    stop("The following \"params\" are not valid parameters: ", paste0(names(params)[ind], sep = " "))
  }
  
  
  if (length(params$mean_of_exp_mu_beta) == 1) {
    params$mean_of_exp_mu_beta <- rep(params$mean_of_exp_mu_beta, 2*Nproc)
  } else stop("The length of \"mean_of_exp_mu_beta\" must be 1")
  if (length(params$var_of_exp_mu_beta) == 1) {
    params$var_of_exp_mu_beta <- rep(params$var_of_exp_mu_beta, 2*Nproc)
  } else stop("The length of \"var_of_exp_mu_beta\" must be 1")
  if (length(params$mean_of_mu_gamma) == 1) {
    params$mean_of_mu_gamma <- rep(params$mean_of_mu_gamma, Nresp)
  } else stop("The length of \"mean_of_mu_gamma\" must be 1")
  if (length(params$var_of_mu_gamma) == 1) {
    params$var_of_mu_gamma <- rep(params$var_of_mu_gamma, Nresp)
  } else stop("The length of \"var_of_mu_gamma\" must be 1")
  if (length(params$mean_of_omega_sqr) != 1 | length(params$var_of_omega_sqr) != 1) stop("mean and variance of omega_sqr must have length 1.")
  if (length(params$df_of_sigma_sqr) != 1) stop("\"df_of_sigma_sqr\" must have length 1.")
  if (length(params$sf_of_scale_matrix_SIGMA) != 1) stop("\"sf_of_scale_matrix_SIGMA\" must have length 1.")
  if (length(params$sf_of_scale_matrix_GAMMA) != 1) stop("\"sf_of_scale_matrix_GAMMA\" must have length 1.")
  if (length(params$prec_epsilon) != 1) stop("\"prec_epsilon\" must have length 1.")
  if (length(params$add_df_to_invWish) != 1) stop("\"add_df_to_invWish\" must have length 1.")
  
  
  e_mu_beta_rate = params$mean_of_exp_mu_beta / params$var_of_exp_mu_beta
  e_mu_beta_shape = params$mean_of_exp_mu_beta * e_mu_beta_rate
  mu_gamma_mean = params$mean_of_mu_gamma
  mu_gamma_var = params$var_of_mu_gamma
  omega_sqr_rate = params$mean_of_omega_sqr / params$var_of_omega_sqr
  omega_sqr_shape = params$mean_of_omega_sqr * omega_sqr_rate
  sigma_sqr_df = params$df_of_sigma_sqr
  SF_cov_SIG = params$sf_of_scale_matrix_SIGMA
  SF_cov_LAM = params$sf_of_scale_matrix_GAMMA
  epsilon_prec = params$prec_epsilon
  add_COV_df = params$add_df_to_invWish - 1
   
  
  
  # ASSIGN PROCESS PARAMS ----
  ProcessAssign <- function(Nproc, Nprob, Nminus, Nplus, Nsubj, epsilon = 1, SF_P = 1, mu_beta_shape = 1, mu_beta_rate = .1, k = 0) {
    if (!is.numeric(Nproc)) stop("Nproc needs to be numeric")
    if (Nproc < 2) stop("Nproc must be larger than 1")
    Nproc <- as.integer(Nproc)
    
    if (!is.numeric(Nprob)) stop("Nprob needs to be numeric")
    if (Nprob < 2) stop("Nprob must be larger than 1")
    Nprob <- as.integer(Nprob)
    
    if (!is.numeric(Nminus)) stop("Nminus needs to be numeric")
    if (Nminus < 1) stop("Nminus must be larger than 0")
    Nmi <- as.integer(Nminus)
    
    if (!is.numeric(Nplus)) stop("Nplus needs to be numeric")
    if (Nplus < 1) stop("Nplus must be larger than 0")
    Npl <- as.integer(Nplus)
    
    if (max(Nprob, Nmi, Npl) > Nproc) stop("Nproc must be number of potential probability \n
                                           parameters in the model and therefore cannot \n
                                           be smaller than Nprob, Nminus, or Nplus")
    
    if (!is.numeric(Nsubj)) stop("Nsubj needs to be numeric")
    if (Nsubj < 1) stop("Nsubj must be larger than 0")
    Nsu <- as.integer(Nsubj)
    
    ## number of total parameter
    Npars <- as.integer(Nprob+Nmi+Npl)
    
    ## zetas
    zeta <- rnorm(n = Npars, mean = 1.0, sd = epsilon^(-1/2))
    
    ## Covariance Matrix for doubleprimes
    S <- rinvwishart(nu = Npars+1+k, S = SF_P*diag(Npars))
    
    ## doubleprimes
    doubleprimes <- rmvn(n = Nsu, mu = rep(0, Npars), Sigma = S)
    
    ## primes
    primes <- zeta*t(doubleprimes)
    
    ## mean alpha
    mu_alpha <- t(rmvn(n = 1, mu = rep(0, Nprob), 
                       Sigma = epsilon^(-1)*diag(Nprob)))
    
    ## mean beta
    mu_beta <- matrix(log(rgamma(n = (Nmi+Npl), shape = mu_beta_shape, 
                                 rate = mu_beta_rate)), ncol = 1)
    
    ## alpha_beta
    alpha_beta <- c(mu_alpha, mu_beta) + primes
    
    ## probabilities
    probs <- pnorm(q = alpha_beta[1:Nprob, ])
    
    ## rates
    rates_minus <- exp(alpha_beta[(Nprob+1):(Nprob+Nmi), ])
    rates_plus <- exp(alpha_beta[(Nprob+Nmi+1):(Nprob+Nmi+Npl), ])
    
    ## generate list
    processList <- list(Nproc = Nproc, Nprob = Nprob, Nminus = Nmi, Nplus = Npl, 
                        Nsubj = Nsu, zeta = zeta, S_doubleprime = S, 
                        doubleprimes = doubleprimes,
                        primes = primes, mu_alpha = mu_alpha,
                        mu_beta = mu_beta, alpha_beta = alpha_beta,
                        probs = probs, rates_minus = rates_minus,
                        rates_plus = rates_plus)
    return(processList)
  }
  
  
  
  # Process Suppression ----
  Suppress_n_Constant <- function(ProcessList, ConstProb = 0, PositionProb = 0, SupprMinus = 0, SupprPlus = 0) {
    if (!is.numeric(ConstProb)) stop("ConstProb needs to be numeric")
    if (any((ConstProb < 0) || (ConstProb >= 1))) stop("ConstProb must be larger than 0 and smaller than 1")
    
    if (!is.numeric(PositionProb)) stop("PositionProb needs to be numeric")
    if (any(PositionProb < 0)) stop("PositionProb must be larger or equal to 0")
    PositionProb <- as.integer(PositionProb)
    
    if (length(ConstProb) != length(PositionProb)) stop("Length of ConstProb and PositionProb do not match")
    
    if (!is.numeric(SupprMinus)) stop("SupprMinus needs to be numeric")
    if (any(SupprMinus < 0) || any(SupprMinus > ProcessList$Nproc)) stop("SupprMinus must be larger or equal to 0 and \n
                                                                         smaller than the potential number of processes")
    SupprMinus <- as.integer(SupprMinus)
    
    if (!is.numeric(SupprPlus)) stop("SupprPlus needs to be numeric")
    if (any(SupprPlus < 0) || any(SupprPlus > ProcessList$Nproc)) stop("SupprPlus must be larger or equal to 0 and \n
                                                                       smaller than the potential number of processes")
    SupprPlus <- as.integer(SupprPlus)
    
    DiffProb <- ProcessList$Nproc-ProcessList$Nprob
    DiffMinus <- ProcessList$Nproc-ProcessList$Nminus
    DiffPlus <- ProcessList$Nproc-ProcessList$Nplus
    
    if (DiffProb != length(ConstProb) && (DiffProb != 0)) stop("Length of ConstProb not right")
    if (DiffMinus != length(SupprMinus) && (DiffMinus != 0)) stop("Length of SupprMinus not right")
    if (DiffPlus != length(SupprPlus) && (DiffPlus != 0)) stop("Length of SupprPlus not right")
    
    ProcessListNew <- ProcessList
    
    if (DiffProb == 0) {
      # message("Nothing is done for the probability parameters")
    } else if (any(ConstProb == 0)) {
      stop("Give allowed real(s) for ConstProb and also integer(s) for PositionProb)")
    } else {
      probs <- matrix(NA, nrow = ProcessList$Nproc, ncol = ProcessList$Nsubj)
      subtr <- 0
      for (p in 1:ProcessList$Nproc) {
        if (any(p == PositionProb)) {
          index <- which(PositionProb == p)
          probs[p, ] <- ConstProb[index]
          subtr <- subtr + 1
        } else {
          probs[p, ] <- ProcessList$probs[p-subtr, ]
        }
      }
      ProcessListNew$probs <- probs
    }
    
    
    if (DiffMinus == 0) {
      # message("Nothing is done for the rate_minus parameters")
    } else if (any(SupprMinus == 0)) {
      stop("Give allowed integer(s) for SupprMinus")
    } else {
      rates_minus <- matrix(NA, nrow = ProcessList$Nproc, ncol = ProcessList$Nsubj)
      subtr <- 0
      for (p in 1:ProcessList$Nproc) {
        if (any(p == SupprMinus)) {
          index <- which(SupprMinus == p)
          rates_minus[p, ] <- 0
          subtr <- subtr + 1
        } else {
          rates_minus[p, ] <- ProcessList$rates_minus[p-subtr, ]
        }
      }
      ProcessListNew$rates_minus <- rates_minus
    }
    
    if (DiffPlus == 0) {
      # message("Nothing is done for the rate_plus parameters")
    } else if (any(SupprPlus == 0)) {
      stop("Give allowed integer(s) for SupprPlus")
    } else {
      rates_plus <- matrix(NA, nrow = ProcessList$Nproc, ncol = ProcessList$Nsubj)
      subtr <- 0
      for (p in 1:ProcessList$Nproc) {
        if (any(p == SupprPlus)) {
          index <- which(SupprPlus == p)
          rates_plus[p, ] <- 0
          subtr <- subtr + 1
        } else {
          rates_plus[p, ] <- ProcessList$rates_plus[p-subtr, ]
        }
      }
      ProcessListNew$rates_plus <- rates_plus
    }
    
    return(ProcessListNew)
  }
  
  
  
  # ASSIGN RESPONSE PARAMS ----
  ResponseAssign <- function(Nresp, Nsubj, epsilon = 1, SF_R = 1, w2_shape = 1, w2_rate = 200, mu_gamma_mean = 0, mu_gamma_var = 1, df = 2, k = 0) {
    if (!is.numeric(Nresp)) stop("Nresp needs to be numeric")
    if (Nresp < 1) stop("Nresp must be larger than 0")
    Nresp <- as.integer(Nresp)
    
    if (!is.numeric(Nsubj)) stop("Nsubj needs to be numeric")
    if (Nsubj < 1) stop("Nsubj must be larger than 0")
    Nsu <- as.integer(Nsubj)
    
    ## zetas
    zeta <- rnorm(n = Nresp, mean = 1.0, sd = epsilon^(-1/2))
    
    ## Covariance Matrix for doubleprimes
    S <- rinvwishart(nu = Nresp+1+k, S = SF_R*diag(Nresp))
    
    ## doubleprimes
    doubleprimes <- rmvn(n = Nsu, mu = rep(0, Nresp), Sigma = S)
    
    ## primes
    primes <- zeta*t(doubleprimes)
    
    ## mean gamma
    mu_gamma <- t(rmvn(n = 1, mu = mu_gamma_mean, 
                       Sigma = mu_gamma_var*diag(Nresp)))
    
    ## gammas
    gamma <- c(mu_gamma) + primes
    
    ## omega_square
    w2 <- rgamma(n = 1, shape = w2_shape, rate = w2_rate)
    #w2 <- rgamma(n = 1, shape = 1, rate = 200)
    #w2 <- rgamma(n = 1, shape = .0025, rate = .5)
    
    ## sigma
    sigma2 <- rinvchisq(n = Nsu, df = df, scale = w2)
    #sigma2 <- rinvgamma(n = Nsu, shape = 1, scale = w2)
    
    ## generate list
    executionlist <- list(Nresp = Nresp, zeta = zeta, 
                          Nsubj = Nsu, S_doubleprime = S, 
                          doubleprimes = doubleprimes,
                          primes = primes, mu_gamma = mu_gamma,
                          gamma = gamma, omega_square = w2, df = df,
                          sigma_square = sigma2)
    return(executionlist)
  }
  
  
  
  # CALCULATE BRANCH PROBABILITIES ----
  CalcBranchProbs <- function(infolist, ProcessList) {
    NiB <- matrix(infolist$NodeInBranch, ncol = infolist$MaxNode*infolist$MaxBranch)
    NTree <- infolist$NTrees
    MNode <- infolist$MaxNode
    MBran <- infolist$MaxBranch
    NiB_list <- list()
    for (i in 1:MNode) {
      name <- paste0("node_", i)
      NiB_list[[name]] <- NiB[ , ((i-1)*MBran+1):(i*MBran)]
    }
    probabilitylist <- list()
    
    for (s in 1:ProcessList$Nsubj) {
      
      BranchProbs <- matrix(rep(1, MBran*2*NTree), ncol = MBran)
      unused <- matrix(rep(1, MBran*2*NTree), ncol = MBran)
      for (t in 1:NTree) {
        corr <- (t-1)*2+1
        fals <- 2*t
        for (b in 1:MBran) {
          for (n in 1:MNode) {
            
            # correct branches
            if (NiB_list[[n]][corr,b] == 1) {
              index <- infolist$tree2node[[t]]$nodes[n]
              BranchProbs[corr,b] <- BranchProbs[corr,b] * ProcessList$probs[index,s]
              unused[corr,b] <- 0
            } else if (NiB_list[[n]][corr,b] == -1) {
              index <- infolist$tree2node[[t]]$nodes[n]
              BranchProbs[corr,b] <- BranchProbs[corr,b] * (1-ProcessList$probs[index,s])
              unused[corr,b] <- 0
            }
            
            # false branches
            if (NiB_list[[n]][fals,b] == 1) {
              index <- infolist$tree2node[[t]]$nodes[n]
              BranchProbs[fals,b] <- BranchProbs[fals,b] * ProcessList$probs[index,s]
              unused[fals,b] <- 0
            } else if (NiB_list[[n]][fals,b] == -1) {
              index <- infolist$tree2node[[t]]$nodes[n]
              BranchProbs[fals,b] <- BranchProbs[fals,b] * (1-ProcessList$probs[index,s])
              unused[fals,b] <- 0
            }
            
          }
          
          # cleaning BranchProbs
          if (unused[corr,b] == 1) BranchProbs[corr,b] <- 0
          if (unused[fals,b] == 1) BranchProbs[fals,b] <- 0
        }
      }
      if (all(BranchProbs==0)) cat("Zero matrix produced in s=", s, "\n")
      probabilitylist[[s]] <- BranchProbs
      
    }
    probslist <- list()
    probslist$probs <- probabilitylist
    probslist$NiB_list <- NiB_list
    
    return(probslist)
  }
  
  
  
  # ASSIGN PROCESS TIMES AND CATEGORIES ----
  Data_List <- function(NTrials, infolist, ProbsList, ProcessList, ExecutionList, Executions) {
    if (!is.numeric(NTrials)) stop("NTrials needs to be numeric")
    if (any(NTrials < 1)) stop("NTrials must be larger than 0")
    if (length(NTrials)==1) NTr <- as.integer(NTrials) else NTr <- as.integer(max(NTrials))
    
    if (!is.numeric(Executions)) stop("Executions needs to be numeric")
    if (any((Executions < 1)) || any((Executions > ExecutionList$Nresp))) stop("Executions must be larger than 0 but not \n
                                                                               be larger than the number of responses")
    Executions <- as.integer(Executions)
    if (length(Executions) < infolist$NCateg) {
      message("The length of Executions does no equal the number of categories.\n
              The vector Executions will be replicated until its length is equal.")
      ExLen <- length(Executions)
      if (infolist$NCateg %% ExLen > 0) warning("Length of Executions should be equal to the number \n
                                                of categories or at least be a factor of it. This will\n
                                                caus weird results now.")
      factor <- ceiling(infolist$NCateg / ExLen)
      Executions <- rep(Executions, factor)[infolist$NCateg]
    }
    
    NTree <- dim(ProbsList$probs[[1]])[1]/2
    MBran <- dim(ProbsList$probs[[1]])[2]
    
    BranchList <- list()
    Data_List <- list()
    for (s in 1:ProcessList$Nsubj) {
      Subj_List <- list(Times = NULL, CAT = NULL, RT = NULL)
      sigma <- sqrt(ExecutionList$sigma_square[s])
      Branch <- matrix(NA, ncol = NTr, nrow = NTree)
      CumulProbs <- t(apply(X = matrix(as.vector(t(ProbsList$probs[[s]])), nrow = NTree, byrow = TRUE), MARGIN = 1, FUN = cumsum))
      U <- matrix(runif(n = NTree*NTr), nrow = NTree)
      for (t in 1:NTree) {
        for (x in 1:NTr) {
          for (b in 1:(2*MBran)) {
            if (U[t, x] <= CumulProbs[t, b]) {
              Branch[t, x] <- b
              break()
            }
          }
          if (is.na(Branch[t, x])) {
            cat("U=",U[t,x], " s=", s, " t=", t, " x=", x, " CumProb=", CumulProbs[t, 2*MBran],"\n")
            ProbsList$probs[[s]]
          }
        }
        
      }
      raw_CAT <- ifelse(Branch > MBran, 2, 1)
      if (any(is.na(Branch))) stop("here is a NA\n")
      CAT <- matrix(NA, nrow = NTree, ncol = NTr)
      for(t in 1:dim(raw_CAT)[1]) CAT[t, ] <- raw_CAT[t, ] + (t-1)*2
      # PT_cat <- list()
      # PT_cat$cat <- CAT
      Subj_List$CAT <- CAT
      BranchList[[s]] <- Branch
      
      PT_Tree <- list()
      for (t in 1:NTree) {
        Categ <- (t-1)*2+1
        PT_TrialParam <- matrix(0, ncol = (2*infolist$NParam+ExecutionList$Nresp), nrow = NTr)
        colnames(PT_TrialParam) <- c(paste0(infolist$ProbNames, "-"), paste0(infolist$ProbNames, "+"), paste0("R", as.character(1:ExecutionList$Nresp)))
        for (x in 1:NTr) {
          add <- 0
          bran <- BranchList[[s]][t, x]
          if (bran > MBran) {
            add <- add+1
            bran <- bran - MBran
          }
          ## Process Times
          for (n in 1:infolist$MaxNode) {
            
            index <- infolist$tree2node[[t]]$nodes[n]
            if (ProbsList$NiB_list[[n]][Categ+add, bran] == 1) {
              rate <- ProcessList$rates_plus[index, s]
              #if(rate==0) warning("rate 0 produced")
              PT_TrialParam[x, index+infolist$NParam] <- PT_TrialParam[x, index+infolist$NParam] + ifelse(rate==0, 0, rexp(n = 1, rate = rate))
            } 
            if (ProbsList$NiB_list[[n]][Categ+add, bran] == -1) {
              rate <- ProcessList$rates_minus[index, s]
              #if(rate==0) warning("rate 0 produced")
              PT_TrialParam[x, index] <- PT_TrialParam[x, index] + ifelse(rate==0, 0, rexp(n = 1, rate = rate))
            }
          }
          ## Execution Times
          resp <- Executions[Categ+add]
          gamma <- ExecutionList$gamma[resp, s]
          PT_TrialParam[x, 2*ProcessList$Nproc+resp] <- rtruncnorm(n = 1, a = 0, b = Inf, mean = gamma, sd = sigma)
        }
        PT_Tree$tree[[t]] <- PT_TrialParam
      }
      Subj_List$RT$tree <- lapply(X = PT_Tree$tree, FUN = function(x) {apply(X = x, MARGIN = 1, FUN = function(y) {sum(y)})})
      Subj_List$Times <- PT_Tree
      Data_List$Subj[[s]] <- Subj_List
    }
    
    return(Data_List)
    
    
  }
  
  
  
  # MAKE DATA FRAME ----
  make_DF <- function(Data_List) {
    NTree <- dim(Data_List$Subj[[1]]$CAT)[1]
    NTr <- dim(Data_List$Subj[[1]]$CAT)[2]
    NSubj <- length(Data_List$Subj)
    DF <- data.frame(subj = rep(NA, NSubj*NTr*NTree), group = NA, tree = NA, cat = NA, rt = NA)
    DF$subj <- rep(1:NSubj, each = (NTr*NTree)) - 1
    DF$group <- 0
    DF$tree <- rep(rep(1:NTree, each = (NTr)), NSubj) - 1
    DF$cat <- c(sapply(X = Data_List$Subj, FUN = function(x) {c(t(x$CAT))})) - 1
    DF$rt <- round(1000*c(sapply(X = Data_List$Subj, FUN = function(x) {unlist(x$RT$tree)})))
    return(DF)
  }
  
  
  set.seed(seed)
  
  infolist <- readinfofile(infofile = infofile)
  if (infolist$NCateg/infolist$NTrees > 2) stop("Currently this function can only be used with models containing two response categories per tree")
  ProcessList <- ProcessAssign(Nproc = Nproc, Nprob = Nprob, Nminus = Nminus, Nplus = Nplus, 
                               Nsubj = Nsubj, epsilon = epsilon_prec, SF_P = SF_cov_SIG, 
                               mu_beta_shape = e_mu_beta_shape, mu_beta_rate = e_mu_beta_rate, 
                               k = add_COV_df)
  ProcessList <- Suppress_n_Constant(ProcessList, ConstProb, PositionProb, SupprMinus, SupprPlus)
  ExecutionList <- ResponseAssign(Nresp = Nresp, Nsubj = Nsubj, epsilon = epsilon_prec, 
                                  SF_R = SF_cov_LAM, w2_shape = omega_sqr_shape, 
                                  w2_rate = omega_sqr_rate, mu_gamma_mean = mu_gamma_mean, 
                                  mu_gamma_var = mu_gamma_var, df = sigma_sqr_df, k = add_COV_df)
  ProbsList <- CalcBranchProbs(infolist, ProcessList)
  PT <- Data_List(NTrials = NTrials, infolist = infolist, 
                  ProbsList = ProbsList, 
                  ProcessList = ProcessList, 
                  ExecutionList = ExecutionList, Executions = Executions)
  DF <- make_DF(PT)
  
  gen_list <- list(model_info = infolist, process_list = ProcessList, motor_list = ExecutionList, 
                   probs_list = ProbsList, RT_list = PT)
  
  specs <- list(model = model, seed = seed, params = params, n.subj = n.subj, n.trials = n.trials, call = match.call())
  
  
  sim_list <- list(data_frame = DF, gen_list = gen_list, specs = specs)
  
  class(sim_list) <- "ertmpt_sim"
  
  return(sim_list)
  
}


