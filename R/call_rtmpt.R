
#' Fit Exponential-RT-MPT Models
#' 
#' Given model and data, this function calls an altered version of the C++ program by Klauer and Kellen (2018) to sample from
#'   the posterior distribution via a Metropolis-Gibbs sampler and storing it in an mcmc.list called \code{samples}. 
#'   Posterior predictive checks developed by Klauer (2010), deviance information criterion (DIC; Spiegelhalter et al., 2002),
#'   99\% and 95\% highest density intervals (HDI) together with the median will be provided for the main parameters in a list 
#'   called \code{diags}. Optionally, the \code{indices} widely applicable information criterion (WAIC; Watanabe, 2010; Vehtari et al., 2017) and 
#'   leave-one-out cross-validation (LOO; Vehtari et al., 2017) can be saved. Additionally the log-likelihood (\code{LogLik}) can also be stored. 
#'   Some specifications of the function call are also saved in \code{specs}.
#'
#' @param model A list of the class \code{ertmpt_model}.
#' @param data Optimally, a list of class \code{ertmpt_data}. Also possible is a \code{data.frame} or a 
#'   path to the text file. Both, \code{data.frame} and the text file must contain the column names "subj", 
#'   "group", "tree", "cat", and "rt" preferably but not necessarily in this order. The values of the latter must 
#'   be in milliseconds. It is always advised to use \code{\link{to_ertmpt_data}} first, which gives back an \code{ertmpt_data} list
#'   with informations about the changes in the data, that were needed.
#' @param n.chains Number of chains to use. Default is 4. Must be larger than 1 and smaller or equal to 16.
#' @param n.iter Number of samples per chain. Default is 5000.
#' @param n.burnin Number of warm-up samples. Default is 200.
#' @param n.thin Thinning factor. Default is 1.
#' @param Rhat_max Maximal Potential scale reduction factor: A lower threshold that needs to be reached before the actual sampling starts. Default is 1.05
#' @param Irep Every \code{Irep} samples an interim state with the current maximal potential scale reduction
#'   factor is shown. Default is 1000. The following statements must hold true for \code{Irep}:
#'   \itemize{
#'     \item \code{n.burnin} is smaller than or equal to \code{Irep},
#'     \item \code{Irep} is a multiple of \code{n.thin} and
#'     \item \code{n.iter} is a multiple of \code{Irep / n.thin}.
#'   }
#' @param prior_params Named list with prior parameters. All parameters have default values, that lead to uninformative priors. Vectors are not allowed.
#'   Allowed parameters are:
#'   \itemize{
#'     \item \code{mean_of_exp_mu_beta}: This is the a priori expected exponential rate (\code{E(exp(beta)) = E(lambda)}) and 
#'           \code{1/mean_of_exp_mu_beta} is the a priori expected process time (\code{1/E(exp(beta)) = E(tau)}). The default
#'           mean is set to \code{10}, such that the expected a priori process time is \code{0.1} seconds.
#'     \item \code{var_of_exp_mu_beta}: The a priori group-specific variance of the exponential rates. Since
#'           \code{exp(mu_beta)} is Gamma distributed, the rate of the distribution is just mean divided by variance and
#'           the shape is the mean times the rate. The default is set to \code{100}.
#'     \item \code{mean_of_mu_gamma}: This is the a priori expected \emph{mean parameter} of the encoding and response execution times,
#'           which follow a normal distribution truncated from below at zero, so \code{E(mu_gamma) < E(gamma)}. The default is \code{0}.
#'     \item \code{var_of_mu_gamma}: The a priori group-specific variance of the \emph{mean parameter}. Its default is \code{10}.
#'     \item \code{mean_of_omega_sqr}: This is the a priori expected residual variance (\code{E(omega^2)}). Its distribution
#'           differs from the one used in the paper. Here it is a Gamma distribution instead of an improper one. The default
#'           is \code{0.005}.
#'     \item \code{var_of_omega_sqr}: The a priori variance of the residual variance (\code{Var(omega^2)}). The default is
#'           \code{0.01}. The default of the mean and variance is equivalent to a shape and rate of \code{0.0025} and 
#'           \code{0.5}, respectivly.
#'     \item \code{df_of_sigma_sqr}: A priori degrees of freedom for the individual variance of the response executions. The
#'           individual variance has a scaled inverse chi-squared prior with \code{df_of_sigma_sqr} degrees of freedom and
#'           \code{omega^2} as scale. \code{2} is the default and it should be an integer.
#'     \item \code{sf_of_scale_matrix_SIGMA}: The original scaling matrix (S) of the (scaled) inverse Wishart distribution for the process 
#'           related parameters is an identity matrix \code{S=I}. \code{sf_of_scale_matrix_SIGMA} is a scaling factor, that scales this 
#'           matrix (\code{S=sf_of_scale_matrix_SIGMA*I}). Its default is \code{1}.
#'     \item \code{sf_of_scale_matrix_GAMMA}: The original scaling matrix (S) of the (scaled) inverse Wishart distribution for the encoding and
#'           motor execution parameters is an identity matrix \code{S=I}. \code{sf_of_scale_matrix_GAMMA} is a scaling factor, that scales 
#'           this matrix (\code{S=sf_of_scale_matrix_GAMMA*I}). Its default is \code{1}.
#'     \item \code{prec_epsilon}: This is epsilon in the paper. It is the precision of mu_alpha and all xi (scaling parameter
#'           in the scaled inverse Wishart distribution). Its default is also \code{1}.
#'     \item \code{add_df_to_invWish}: If \code{P} is the number of parameters or rather the size of the scale matrix used in the (scaled)
#'           inverse Wishart distribution then \code{add_df_to_invWish} is the number of degrees of freedom that can be added to it. So
#'           \code{DF = P + add_df_to_invWish}. The default for \code{add_df_to_invWish} is \code{1}, such that the correlations are uniformly 
#'           distributed within \code{[-1, 1]}.
#'   }
#' @param indices Model selection indices. If set to \code{TRUE} the log-likelihood for each iteration and trial will be stored temporarily
#'   and with that the WAIC and LOO will be calculated via the \code{loo} package. If you want to have this log-likelihood matrix stored in the
#'   output of this function, you can set \code{save_log_lik} to \code{TRUE}. The default for \code{indices} is \code{FALSE}.
#' @param save_log_lik If set to \code{TRUE} and \code{indices = TRUE} the log-likelihood matrix for each iteration and trial will
#'   be saved in the output as a matrix. Its default is \code{FALSE}.
#' @param old_label If set to \code{TRUE} the old labels of "subj" and "group" of the data will be used in the elements of the output list. Default is \code{FALSE}.
#' @return A list of the class \code{ertmpt_fit} containing 
#'   \itemize{
#'     \item \code{samples}: the posterior samples as an \code{mcmc.list} object,
#'     \item \code{diags}: some diagnostics like deviance information criterion, posterior predictive checks for the frequencies and latencies, 
#'                         potential scale reduction factors, and also the 99\% and 95\% HDIs and medians for the group-level parameters,
#'     \item \code{specs}: some model specifications like the model, arguments of the model call, and information about the data transformation,
#'     \item \code{indices} (optional): if enabled, WAIC and LOO,
#'     \item \code{LogLik} (optional): if enabled, the log-likelihood matrix used for WAIC and LOO. 
#'     \item \code{summary} includes posterior mean and median of the main parameters.
#'   }
#' @references
#' Hartmann, R., Johannsen, L., & Klauer, K. C. (2020). rtmpt: An R package for fitting response-time extended multinomial processing tree models. 
#'   \emph{Behavior Research Methods, 52}(3), 1313–1338. 
#' 
#' Hartmann, R., & Klauer, K. C. (2020). Extending RT-MPTs to enable equal process times. \emph{Journal of Mathematical Psychology, 96}, 102340.
#' 
#' Klauer, K. C. (2010). Hierarchical multinomial processing tree models: A latent-trait approach. \emph{Psychometrika, 75(1)}, 70-98.
#' 
#' Klauer, K. C., & Kellen, D. (2018). RT-MPTs: Process models for response-time distributions based on multinomial processing trees with 
#'   applications to recognition memory. \emph{Journal of Mathematical Psychology, 82}, 111-130.
#' 
#' Spiegelhalter, D. J., Best, N. G., Carlin, B. P., & Van Der Linde, A. (2002). Bayesian measures of model complexity and fit. 
#'   \emph{Journal of the royal statistical society: Series b (statistical methodology), 64(4)}, 583-639.
#'   
#' Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. 
#'   \emph{Statistics and Computing, 27(5)}, 1413-1432.
#' 
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory. 
#'   \emph{Journal of Machine Learning Research, 11(Dec)}, 3571-3594.
#' @examples 
#' ####################################################################################
#' # Detect-Guess variant of the Two-High Threshold model.
#' # The encoding and motor execution times are assumed to be equal for each response.
#' ####################################################################################
#' 
#' mdl_2HTM <- "
#' # targets
#' do+(1-do)*g
#' (1-do)*(1-g)
#'
#' # lures
#' (1-dn)*g
#' dn+(1-dn)*(1-g)
#' 
#' # do: detect old; dn: detect new; g: guess
#' "
#' 
#' model <- to_ertmpt_model(mdl_file = mdl_2HTM)
#' 
#' data_file <- system.file("extdata/data.txt", package="rtmpt")
#' data <- read.table(file = data_file, header = TRUE)
#' data_list <- to_ertmpt_data(raw_data = data, model = model)
#' \donttest{
#' # This might take some time
#' ertmpt_out <- fit_ertmpt(model = model, data = data_list, Rhat_max = 1.1)
#' ertmpt_out
#' }
#' # Type ?SimData for another working example.
#' @author Raphael Hartmann
#' @useDynLib "rtmpt", .registration=TRUE
#' @export
#' @importFrom coda as.mcmc as.mcmc.list "varnames<-" gelman.diag
#' @importFrom data.table as.data.table fread
#' @importFrom loo loo waic relative_eff
#' @importFrom methods as callGeneric new
#' @importFrom stats runif
#' @importFrom utils read.table write.table
#' @importFrom stringr str_split str_detect str_c str_replace
#' @importFrom methods is
fit_ertmpt <- function(model, 
                       data,
                       n.chains = 4, 
                       n.iter = 5000, 
                       n.burnin = 200, 
                       n.thin = 1,
                       Rhat_max = 1.05, 
                       Irep = 1000, 
                       prior_params = NULL,
                       indices = FALSE, 
                       save_log_lik = FALSE,
                       old_label = FALSE) {
  
  # internal argument name changes
  Nchains <- n.chains
  Nsamples <- n.iter
  Nwarmup <- n.burnin
  thin <- n.thin

  # some controls
  model_elmnts <- c("lines", "params", "responses")
  if (!is.list(model)) stop("\"model\" must be a list.")
  if (!all(model_elmnts  %in% names(model))) stop("\"model\" must contain \"", model_elmnts[which(!(model_elmnts %in% names(model)))[1]], "\".")
  
  if (Nchains < 2 || Nchains > 16) stop("\"n.chains\" must be larger than 1 and lower or equal to 16.")
  
  if (Irep < 1) stop("\"Irep\" must be larger than or equal to one.")
  if (thin < 1) stop("\"n.thin\" must be larger than or equal to one.")
  
  if (Irep %% thin != 0) stop("\"Irep\" must be a multiple of \"n.thin\".")
  if (Nsamples %% (Irep/thin) != 0) stop("\"n.iter\" must be a multiple of \"Irep\" / \"n.thin\".")
  if (Nsamples < Irep/thin) stop("\"n.iter\" must be greater or equal to \"Irep\" / \"n.thin\" = ", Irep*thin, ".")
  if (Nwarmup > Irep) stop("\"n.burnin\" must be smaller than or equal to \"Irep\".")
  
  if (Rhat_max < 1) stop("\"Rhat_max\" must be larger than or equal to one.")
  
  if (!is.logical(indices)) stop("\"indices\" must either be TRUE or FALSE.")
  # if (!is.logical(bridge)) stop("\"bridge\" must either be TRUE or FALSE.")
  bridge = FALSE
  if (!is.logical(save_log_lik)) stop("\"save_log_lik\" must either be TRUE or FALSE.")
  if (!is.logical(old_label)) stop("\"old_label\" must either be TRUE or FALSE.")
  if (!is.null(prior_params) && !is.list(prior_params)) stop("\"prior_params\" must be a list.")

  # prepare data
  keep_data_path <- FALSE
  if (is.data.frame(data)) {
    temp_data <- to_ertmpt_data(data, model)
    data_frame <- temp_data$data
  	if("transformation" %in% names(temp_data)) {
  	  transformation <- temp_data$transformation
	  } else transformation <- list()
  } else if (is.character(data)) {
    temp_data <- to_ertmpt_data(read.table(file = data, header = TRUE), model)
    data_frame <- temp_data$data
  	if("transformation" %in% names(temp_data)) {
  	  transformation <- temp_data$transformation
  	} else transformation <- list()
    keep_data_path <- TRUE
  } else if (inherits(data, c("ertmpt_data", "rtmpt_data"))) {
    data_frame <- data$data
  	if("transformation" %in% names(data)) {
  	  transformation <- data$transformation
  	} else transformation <- list()
  }
  rm(data); gc()
  
  data_elmnts <- c("subj", "group", "tree", "cat", "rt")
  if (!all(data_elmnts %in% names(data_frame))) stop("\"data\" must contain \"", data_elmnts[which(!(data_elmnts %in% names(data_frame)))[1]], "\".")
  if (!all(data_elmnts == names(data_frame))) {
    df <- data_frame
    data_frame <- df[,match(data_elmnts, names(data_frame))]
  }
  if (any(is.na(data_frame)) || min(data_frame) < 0) stop("All values in \"data\" need to be larger than or equal to zero and must not be NA.")
  if (max(data_frame[,1]+1) != length(unique(data_frame[,1]))) stop("\"max(data$subj+1)\" must equal \"length(unique(data$subj))\". You might want to use to_ertmpt_data().")
  if (max(data_frame[,2]+1) != length(unique(data_frame[,2]))) stop("\"max(data$group+1)\" must equal \"length(unique(data$group))\". You might want to use to_ertmpt_data().")
  if (max(data_frame[,3]+1) != length(unique(data_frame[,3]))) stop("\"max(data$tree+1)\" must equal \"length(unique(data$tree))\". You might want to use to_ertmpt_data().")
  if (max(data_frame[,4]+1) != length(unique(data_frame[,4]))) stop("\"max(data$cat+1)\" must equal \"length(unique(data$cat))\". You might want to use to_ertmpt_data().")
  
  temp_dir <- gsub("\\\\", "/", tempdir())
  data_path <- gsub("\\\\", "/", tempfile(pattern = "ertmpt_data", tmpdir = tempdir(), fileext = ".txt"))
  write.table(x = data_frame, file = data_path, sep = " ", 
              row.names = FALSE, col.names = TRUE)
  diag_path <- gsub("\\\\", "/", tempfile(pattern = "diagnosis", tmpdir = tempdir(), fileext = ".out"))
  chains_path <- gsub("\\\\", "/", tempfile(pattern = "chains_raw", tmpdir = tempdir(), fileext = ".out"))
  loglik_path <- gsub("\\\\", "/", tempfile(pattern = "log_likelihood", tmpdir = tempdir(), fileext = ".out"))
  mdl_txt <- gsub("\\\\", "/", tempfile(pattern = "model", tmpdir = tempdir(), fileext = ".txt"))
  mdl_info <- gsub("\\\\", "/", tempfile(pattern = "model", tmpdir = tempdir(), fileext = ".info"))

  
  # produce infofile
  infofile <- try(get_infofile(model, mdl_txt = mdl_txt, mdl_info = mdl_info))
  if(is(infofile, "try-error")) stop("problem with S4 routines in makin infofile.\n")
  
  # # functions to use
  # varnames <- function(infofile) {
  #   varstring <- system(paste("gawk 'END {print}' ", infofile), intern = TRUE)
  #   var_char <- strsplit(varstring, "\\s+")[[1]]
  #   return(var_char)
  # }
  # 
  # var_char <- varnames(infofile = infofile)
  # if (!all(var_char == names(model$params$probs))) stop("Fatal error: Wrong order of parameters in model file.")
  
  
  # prior parameters
  prior_names <- c("mean_of_exp_mu_beta", "var_of_exp_mu_beta", "mean_of_mu_gamma", "var_of_mu_gamma",
                   "mean_of_omega_sqr", "var_of_omega_sqr", "sf_of_scale_matrix_SIGMA", "sf_of_scale_matrix_GAMMA",
                   "df_of_sigma_sqr", "prec_epsilon", "add_df_to_invWish")
  if (is.null(prior_params)) prior_params <- list()
  if (!("mean_of_exp_mu_beta" %in% names(prior_params))) prior_params$mean_of_exp_mu_beta <- 10
  if (!("var_of_exp_mu_beta" %in% names(prior_params))) prior_params$var_of_exp_mu_beta <- 100
  if (!("mean_of_mu_gamma" %in% names(prior_params))) prior_params$mean_of_mu_gamma <- 0
  if (!("var_of_mu_gamma" %in% names(prior_params))) prior_params$var_of_mu_gamma <- 10
  if (!("mean_of_omega_sqr" %in% names(prior_params))) prior_params$mean_of_omega_sqr <- 0.005
  if (!("var_of_omega_sqr" %in% names(prior_params))) prior_params$var_of_omega_sqr <- 0.01
  if (!("df_of_sigma_sqr" %in% names(prior_params))) prior_params$df_of_sigma_sqr <- 2
  if (!("sf_of_scale_matrix_SIGMA" %in% names(prior_params))) prior_params$sf_of_scale_matrix_SIGMA <- 1
  if (!("sf_of_scale_matrix_GAMMA" %in% names(prior_params))) prior_params$sf_of_scale_matrix_GAMMA <- 1
  if (!("prec_epsilon" %in% names(prior_params))) prior_params$prec_epsilon <- 1
  if (!("add_df_to_invWish" %in% names(prior_params))) prior_params$add_df_to_invWish <- 1
  if (prior_params$add_df_to_invWish %% 1 != 0) stop("\"add_df_to_invWish\" must be an integer!")
  if (!all(names(prior_params) %in% prior_names)) {
    ind <- which(!(names(prior_params) %in% prior_names))
    stop("The following \"prior_params\" are not valid parameters: ", paste0(names(prior_params)[ind], sep = " "))
  }
  
  
  # prior transformation
  rate_exp_mu_beta <- prior_params$mean_of_exp_mu_beta / prior_params$var_of_exp_mu_beta
  shape_exp_mu_beta <- prior_params$mean_of_exp_mu_beta * rate_exp_mu_beta
  rate_omega_sqr <- prior_params$mean_of_omega_sqr / prior_params$var_of_omega_sqr
  shape_omega_sqr <- prior_params$mean_of_omega_sqr * rate_omega_sqr
  add_df_invWish <- prior_params$add_df_to_invWish - 1

  
  # process names and number of precesses in total
  proc_names <- names(model$params$probs)
  nprocs <- length(proc_names)
  
  # prepare arguments for RTMPT
  CHAR <- c(Data = data_path, 
            Model = infofile, 
            Raus = chains_path, 
            Diagnostics = diag_path,
            loglike = loglik_path)
  
  REAL <- c(Rhat_max)
  
  REAL2 <- sapply(X = model$params$probs, FUN = function(x) {ifelse(is.na(x) | x %in% proc_names, -1, x)})  # constants in probs

  INTEGER <- c(NoThr = Nchains, 
               BurnIn = Nwarmup, 
               Thin = thin, 
               SampSize = Nchains*Nsamples,
               Irep = Irep,
               nKERN = dim(model$responses)[1], 
               nRESP = length(unique(model$responses$MAP)))
  INTEGER2 <- model$responses$MAP  # cat2resp
  
  BOOL1 <- sapply(X = model$params[["taus"]]["minus",], FUN = function(x) {ifelse(is.na(x) | x %in% proc_names, 1, 0)})  # tau minus
  BOOL2 <- sapply(X = model$params[["taus"]]["plus",], FUN = function(x) {ifelse(is.na(x) | x %in% proc_names, 1, 0)})  # tau plus
  
  BOOL3 <- c(loglik = indices, bridge = FALSE)#bridge)
  
  REAL3 <- c(prior_params$df_of_sigma_sqr, 
             shape_omega_sqr, 
             rate_omega_sqr, 
             prior_params$mean_of_mu_gamma, 
             prior_params$var_of_mu_gamma, 
             prior_params$prec_epsilon, 
             shape_exp_mu_beta, 
             rate_exp_mu_beta, 
             prior_params$sf_of_scale_matrix_SIGMA,
             prior_params$sf_of_scale_matrix_GAMMA)
  
  INTEGER3 <- add_df_invWish
  

  F2K <- integer(3*nprocs)
  for(i in 1:nprocs) {
    theta <- model$params$probs[i]
    F2K[i] <- ifelse(theta %in% proc_names, which(proc_names==as.character(theta)), 
                     ifelse(is.na(theta), i, -1))
    tauminus <- model$params[["taus"]]["minus", i]
    F2K[i+nprocs] <- ifelse(tauminus %in% proc_names, nprocs + which(proc_names==as.character(tauminus)), 
                            ifelse(BOOL1[i], nprocs + i, -1))
    tauplus <- model$params[["taus"]]["plus", i]
    F2K[i+2*nprocs] <- ifelse(tauplus %in% proc_names, 2*nprocs + which(proc_names==as.character(tauplus)), 
                              ifelse(BOOL2[i], 2*nprocs + i, -1))
  }
  INTEGER5 <- unique(F2K[F2K != -1])-1

  
  K2F <- integer(3*nprocs)
  cntk2f <- 1
  for(i in 1:(3*nprocs)) {
    tmp <- NULL
    if (i <= nprocs) {
      tmp <- model$params$probs[i]
    } else if (i <= 2*nprocs) {
      tmp <- model$params[["taus"]]["minus", i-nprocs]
    } else {tmp <- model$params[["taus"]]["plus", i-2*nprocs]}
    
    if(is.na(tmp)) {
      K2F[i] <- cntk2f
      cntk2f <- cntk2f + 1
    } else if(tmp == 0 | (tmp < 1 & tmp > 0)) {
      K2F[i] <- -1
    } else if(tmp %in% proc_names) {
      K2F[i] <- K2F[which(proc_names == as.character(tmp)) + (i>nprocs)*nprocs + (i>2*nprocs)*nprocs]
    }
  }
  INTEGER4 <- K2F-1*(K2F!=-1)
  
  
  
  # call C++ function RTMPT
  out <- .Call("ertmpt_fit", 
               as.numeric(REAL), 
               as.numeric(REAL2),
               as.numeric(REAL3),
               as.character(CHAR), 
               as.integer(INTEGER), 
               as.integer(INTEGER2),
               as.integer(INTEGER3),
               as.integer(INTEGER4),
               as.integer(INTEGER5),
               as.logical(BOOL1), 
               as.logical(BOOL2), 
               as.logical(BOOL3))

  
  # get data information
  data_info <- list(Nsubj = length(unique(data_frame[,1])), 
                    Ngroups = length(unique(data_frame[,2])), 
                    Nprobs = sum(is.na(model$params$probs)), 
                    Nminus = sum(is.na(model$params$taus[1,])), 
                    Nplus = sum(is.na(model$params$taus[2,])), 
                    Nresps = length(unique(model$responses$MAP)), 
                    # epsilon = 1,
                    probs_string = names(model$params$probs[which(is.na(model$params$probs))]), 
                    minus_string = names(model$params$probs[which(is.na(model$params$taus[1,]))]),
                    plus_string = names(model$params$probs[which(is.na(model$params$taus[2,]))]),
                    transformation = transformation)
  

  
  # prepare output list
  ertmpt <- list()
  ertmpt$samples <- make_mcmc_list(file = out$pars_samples, infofile = infofile, 
                                  Nchains = Nchains, Nsamples = Nsamples, 
                                  data_info = data_info, keep = old_label)
  out$pars_samples <- NULL; gc()
  
  
  # diagnostics
  ertmpt$diags <- get_diags(diag_file = diag_path, data_info = data_info, keep = old_label)
  if (file.exists(diag_path)) file.remove(diag_path)
  ertmpt$diags$R_hat <- gelman.diag(ertmpt$samples, multivariate = FALSE)

  
  # specs
  infos <- readinfofile(infofile)
  ertmpt$specs <- list(model = model, n.chains = Nchains, n.iter = Nsamples, n.burnin = Nwarmup, n.thin = thin, 
                      n.groups = data_info$Ngroups, n.subj = data_info$Nsubj, Irep = Irep, Rhat_max = Rhat_max, 
                      prior_params = prior_params, infolist = infos, call = match.call(), 
					  n.groups = max(data_frame[,2]+1))
  if(exists("transformation")) {
    ertmpt$specs$transformation <- transformation
  }

  
  # WAIC & LOO
  if (file.exists(loglik_path)) {
    if (indices) {
      suppressWarnings( temp <- get_indices_x(loglik_path, Nchains, Nsamples, data_frame) )
      ertmpt$indices <- list(WAIC = temp$WAIC, LOO = temp$LOO)
      if (save_log_lik) ertmpt$LogLik <- temp$LogLik
      rm(temp)
    }
  }
  file.remove(loglik_path)
  
  
  # remove data file
  if (!keep_data_path && file.exists(data_path)) file.remove(data_path)
  
  
  # remove info file and chains
  if (file.exists(infofile)) file.remove(infofile)
  if (file.exists(chains_path)) file.remove(chains_path)
  
  
  # summary
  ertmpt$summary <- writeSummaryERTMPT(ertmpt, keep = old_label)
  
  
  # output
  class(ertmpt) <- "ertmpt_fit"
  gc()
  return(ertmpt)
  
}






#' Fit Diffusion-RT-MPT Models
#'
#' Given model and data, this function a Hamiltonian MCMC sampler and stores the samples in an mcmc.list called \code{samples}.
#'   Posterior predictive checks developed by Klauer (2010), deviance information criterion (DIC; Spiegelhalter et al., 2002),
#'   99\% and 95\% highest density intervals (HDI) together with the median will be provided for the main parameters in a list
#'   called \code{diags}. Optionally, the \code{indices} widely applicable information criterion (WAIC; Watanabe, 2010; Vehtari et al., 2017) and
#'   leave-one-out cross-validation (LOO; Vehtari et al., 2017) can be saved. Additionally the log-likelihood (\code{LogLik}) can also be stored.
#'   Some specifications of the function call are also saved in \code{specs}.
#'
#' @param model A list of the class \code{drtmpt_model}.
#' @param data Optimally, a list of class \code{drtmpt_data}. Also possible is a \code{data.frame} or a
#'   path to the text file. Both, \code{data.frame} and the text file must contain the column names "subj",
#'   "group", "tree", "cat", and "rt" preferably but not necessarily in this order. The values of the latter must
#'   be in milliseconds. It is always advised to use \code{\link{to_drtmpt_data}} first, which gives back a
#'   \code{drtmpt_data} list with information about the changes in the data, that were needed.
#' @param n.chains Number of chains to use. Default is 4. Must be larger than 1 and smaller or equal to 16.
#' @param n.iter Number of samples per chain. Default is 1000.
#' @param n.phase1 Number of samples for phase 1 (adaptation phase). Default is 1000.
#' @param n.phase2 Number of samples for phase 2. Default is 2000.
#' @param n.thin Thinning factor. Default is 1.
#' @param Rhat_max Maximal Potential scale reduction factor: A lower threshold that needs to be reached before the actual sampling starts. Default is 1.05
#' @param Irep Every \code{Irep} samples an interim state with the current maximal potential scale reduction
#'   factor is shown. Default is 1000. The following statements must hold true for \code{Irep}:
#'   \itemize{
#'     \item \code{n.phase1} is a multiple of \code{Irep}
#'     \item \code{n.phase2} is a multiple of \code{Irep}
#'     \item \code{n.phase1} is smaller than or equal to \code{n.phase2},
#'     \item \code{Irep} is a multiple of \code{n.thin} and
#'     \item \code{n.iter} is a multiple of \code{Irep / n.thin}.
#'   }
#' @param prior_params Named list with prior parameters. All parameters have default values, that lead to uninformative priors. Vectors are not allowed.
#'   Allowed parameters are:
#'   \itemize{
#'     \item \code{prec_epsilon}: prior precision for population means of process-related parameters.
#'       Default is \code{1.0}.
#'     \item \code{delta_df}: degrees of freedom of t-distribution for motor times. Default is \code{10}.
#'     \item \code{delta_mu}: mean of t-distribution for motor times. Default is \code{0.5}.
#'     \item \code{delta_scale}: scale of t-distribution for motor times. Default is \code{1.0}.
#'     \item \code{SIGMA_Corr_eta}: shape parameter for LKJ distribution for process-related parameters. Default is \code{4.0}.
#'     \item \code{SIGMA_SD_rho}: scale parameter of half-Cauchy distribution for process-related parameters. Default is \code{2.5}.
#'     \item \code{GAMMA_Corr_eta}: shape parameter for LKJ distribution for motor-related parameters. Default is \code{4.0}.
#'     \item \code{GAMMA_SD_rho}: scale parameter of half-Cauchy distribution for motor-related parameters. Default is \code{0.5}.
#'     \item \code{Omega2_alpha}: shape parameter of gamma distribution for residual variance. Default is \code{0.0025}.
#'     \item \code{Omega2_beta}: rate parameter of gamma distribution for residual variance. Default is \code{0.5}.
#'   }
#' @param flags Either NULL or a list of
#'   \itemize{
#'     \item \code{old_label} If set to \code{TRUE} the old labels of "subj" and "group" of the data will be used in the elements of the output list.
#'       Default is \code{FALSE}.
#'     \item \code{indices} Model selection indices. If set to \code{TRUE} the log-likelihood for each iteration and trial will be
#'       stored temporarily and with that the WAIC and LOO will be calculated via the \code{loo} package. If you want to have this
#'       log-likelihood matrix stored in the output of this function, you can set \code{loglik} to \code{TRUE}. Default for
#'       \code{indices} is \code{FALSE}.
#'     \item \code{loglik} If set to \code{TRUE} and \code{indices = TRUE} the log-likelihood matrix for each iteration and trial will
#'       be saved in the output as a matrix. Default is \code{FALSE}.
#'     \item \code{random_init} If set to \code{TRUE} the initial values are randomly drawn. If \code{FLASE} maximum likelihood is used
#'       for initial values.
#'   }
#' @param control Either NULL or a list of
#'   \itemize{
#'     \item \code{maxthreads} for the ML estimation of the initial values and the calculation of the DIC values one can use more than
#'       \code{n.chains} threads for parallelization. Default is 4 like \code{n.chians}. \code{maxthreads} must be larger or equal to
#'       \code{n.chains}.
#'     \item \code{maxtreedepth1_3} maxtree-depth of the no-U-turn algorithm in Phases 1 to 3
#'     \item \code{maxtreedepth4} maxtree-depth of the no-U-turn algorithm in Phases 4
#'   }
#' @return A list of the class \code{drtmpt_fit} containing
#'   \itemize{
#'     \item \code{samples}: the posterior samples as an \code{mcmc.list} object,
#'     \item \code{diags}: some diagnostics like deviance information criterion, posterior predictive checks for the frequencies and latencies,
#'                         potential scale reduction factors, and also the 99\% and 95\% HDIs and medians for the group-level parameters,
#'     \item \code{specs}: some model specifications like the model, arguments of the model call, and information about the data transformation,
#'     \item \code{indices} (optional): if enabled, WAIC and LOO,
#'     \item \code{LogLik} (optional): if enabled, the log-likelihood matrix used for WAIC and LOO.
#'     \item \code{summary} includes posterior mean and median of the main parameters.
#'   }
#' @references
#' Klauer, K. C. (2010). Hierarchical multinomial processing tree models: A latent-trait approach. \emph{Psychometrika, 75(1)}, 70-98.
#'
#' Spiegelhalter, D. J., Best, N. G., Carlin, B. P., & Van Der Linde, A. (2002). Bayesian measures of model complexity and fit.
#'   \emph{Journal of the royal statistical society: Series b (statistical methodology), 64(4)}, 583-639.
#'
#' Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC.
#'   \emph{Statistics and Computing, 27(5)}, 1413-1432.
#'
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory.
#'   \emph{Journal of Machine Learning Research, 11(Dec)}, 3571-3594.
#'
#' @examples
#' ####################################################################################
#' # Detect-Guess variant of the Two-High Threshold model.
#' # The encoding and motor execution times are assumed to be equal for each response.
#' ####################################################################################
#'
#' mdl_2HTM <- "
#' # targets
#' do+(1-do)*g
#' (1-do)*(1-g)
#'
#' # lures
#' (1-dn)*g
#' dn+(1-dn)*(1-g)
#'
#' # do: detect old; dn: detect new; g: guess
#' "
#'
#' model <- to_drtmpt_model(mdl_file = mdl_2HTM)
#'
#' data_file <- system.file("extdata/data.txt", package="rtmpt")
#' data <- read.table(file = data_file, header = TRUE)
#' data_list <- to_drtmpt_data(raw_data = data, model = model)
#' \donttest{
#' # This might take some time
#' drtmpt_out <- fit_drtmpt(model = model, data = data_list, Rhat_max = 1.1)
#' drtmpt_out
#' }
#' @author Raphael Hartmann
#' @useDynLib "rtmpt", .registration=TRUE
#' @export
#' @importFrom coda as.mcmc as.mcmc.list "varnames<-" gelman.diag
#' @importFrom data.table as.data.table fread
#' @importFrom loo loo waic relative_eff
#' @importFrom methods as callGeneric new
#' @importFrom stats runif
#' @importFrom utils read.table write.table
#' @importFrom stringr str_split str_detect str_c str_replace
#' @importFrom methods is
fit_drtmpt <- function(model,
                       data,
                       n.chains = 4,
                       n.iter = 1000,
                       n.phase1 = 1000,
                       n.phase2 = 2000,
                       n.thin = 1,
                       Rhat_max = 1.1,
                       Irep = 1000,
                       prior_params = NULL,
                       flags = NULL,
                       control = NULL) {
  
  
  # CHECKS
  model_elmnts <- c("lines", "params", "responses")
  if (!is.list(model)) stop("\"model\" must be a list.")
  if (!all(model_elmnts  %in% names(model))) stop("\"model\" must contain \"", model_elmnts[which(!(model_elmnts %in% names(model)))[1]], "\".")
  
  if (n.chains < 2 || n.chains > 16) stop("\"n.chains\" must be larger than 1 and lower or equal to 16.")
  
  if (Irep < 1) stop("\"Irep\" must be larger than or equal to one.")
  if (n.thin < 1) stop("\"n.thin\" must be larger than or equal to one.")
  
  if (n.phase1 %% Irep) stop("\"n.phase1\" must be a multiple of \"Irep\".")
  if (n.phase2 %% Irep) stop("\"n.phase2\" must be a multiple of \"Irep\".")
  if (n.phase1 > n.phase2) stop("\"n.phase2\" must be larger or equal to \"n.phase1\".")
  
  if (Irep %% n.thin != 0) stop("\"Irep\" must be a multiple of \"n.thin\".")
  if (n.iter %% (Irep/n.thin) != 0) stop("\"n.iter\" must be a multiple of \"Irep\" / \"n.thin\".")
  if (n.iter < Irep/n.thin) stop("\"n.iter\" must be greater or equal to \"Irep\" / \"n.thin\" = ", Irep*n.thin, ".")
  
  if (Rhat_max < 1) stop("\"Rhat_max\" must be larger than or equal to one.")
  
  if (!is.null(prior_params) && !is.list(prior_params)) stop("\"prior_params\" must be a list or NULL.")
  if (!is.null(flags) && !is.list(flags)) stop("\"flags\" must be a list or NULL.")
  if (!is.null(control) && !is.list(control)) stop("\"control\" must be a list or NULL.")
  
  if (!is.null(flags)) {
    if ("old_label" %in% names(flags)) if (!is.logical(flags$old_label)) stop("\"flags$old_label\" must either be TRUE or FALSE.")
    if ("indices" %in% names(flags)) if (!is.logical(flags$indices)) stop("\"flags$indices\" must either be TRUE or FALSE.")
    if ("loglik" %in% names(flags)) if (!is.logical(flags$loglik)) stop("\"flags$loglik\" must either be TRUE or FALSE.")
    if ("random_init" %in% names(flags)) if (!is.logical(flags$random_init)) stop("\"flags$random_init\" must either be TRUE or FALSE.")
  }
  
  if (!is.null(control)) {
    if ("maxthreads" %in% control) if (!is.numeric(control$maxthreads)) stop("\"control$maxthreads\" must be numeric.")
    if ("maxtreedepth1_3" %in% control) if (!is.numeric(control$maxtreedepth1_3)) stop("\"control$maxtreedepth1_3\" must be numeric.")
    if ("maxtreedepth4" %in% control) if (!is.numeric(control$maxtreedepth4)) stop("\"control$maxtreedepth4\" must be numeric.")
  }
  
  
  # SET FLAGS
  if (!is.list(flags)) flags <- list()
  if (!"old_label" %in% names(flags)) {
    old_label <- FALSE
  } else {
    old_label <- flags$old_label
  }
  if (!"indices" %in% names(flags)) flags$indices <- FALSE
  if (!"loglik" %in% names(flags)) flags$loglik <- FALSE
  if (!"random_init" %in% names(flags)) flags$random_init <- FALSE
  
  
  # SET CONTROLS
  if (!is.list(control)) control <- list()
  if (!"maxthreads" %in% names(flags)) control$maxthreads <- 4
  if (!"maxtreedepth1_3" %in% names(flags)) control$maxtreedepth1_3 <- 5
  if (!"maxtreedepth4" %in% names(flags)) control$maxtreedepth4 <- 9
  
  
  # PREPARE DATA
  keep_data_path <- FALSE
  data_frame <- NULL
  if (is.data.frame(data)) {
    temp_data <- to_drtmpt_data(data, model)
    data_frame <- temp_data$data
    if("transformation" %in% names(temp_data)) {
      transformation <- temp_data$transformation
    } else transformation <- list()
  } else if (is.character(data)) {
    temp_data <- to_drtmpt_data(read.table(file = data, header = TRUE), model)
    data_frame <- temp_data$data
    if("transformation" %in% names(temp_data)) {
      transformation <- temp_data$transformation
    } else transformation <- list()
    keep_data_path <- TRUE
  } else if (inherits(data, "drtmpt_data")) {
    data_frame <- data$data
    if("transformation" %in% names(data)) {
      transformation <- data$transformation
    } else transformation <- list()
  }
  rm(data); gc()
  
  data_elmnts <- c("subj", "group", "tree", "cat", "rt")
  if (!all(data_elmnts %in% names(data_frame))) stop("\"data\" must contain \"", data_elmnts[which(!(data_elmnts %in% names(data_frame)))[1]], "\".")
  if (!all(data_elmnts == names(data_frame))) {
    df <- data_frame
    data_frame <- df[,match(data_elmnts, names(data_frame))]
  }
  if (any(is.na(data_frame)) || min(data_frame) < 0) stop("All values in \"data\" need to be larger than or equal to zero and must not be NA.")
  if (max(data_frame[,1]+1) != length(unique(data_frame[,1]))) stop("\"max(data$subj+1)\" must equal \"length(unique(data$subj))\". You might want to use to_rtmpt_data().")
  if (max(data_frame[,2]+1) != length(unique(data_frame[,2]))) stop("\"max(data$group+1)\" must equal \"length(unique(data$group))\". You might want to use to_rtmpt_data().")
  if (max(data_frame[,3]+1) != length(unique(data_frame[,3]))) stop("\"max(data$tree+1)\" must equal \"length(unique(data$tree))\". You might want to use to_rtmpt_data().")
  if (max(data_frame[,4]+1) != length(unique(data_frame[,4]))) stop("\"max(data$cat+1)\" must equal \"length(unique(data$cat))\". You might want to use to_rtmpt_data().")
  
  
  # PREPARE PATHS
  # Inputs
  model_path <- gsub("\\\\", "/", tempfile(pattern = "model", tmpdir = tempdir(), fileext = ".info"))
  mdl_path <- gsub("\\\\", "/", tempfile(pattern = "model", tmpdir = tempdir(), fileext = ".txt"))
  data_path <- gsub("\\\\", "/", tempfile(pattern = "data", tmpdir = tempdir(), fileext = ".txt"))
  write.table(x = data_frame, file = data_path, sep = " ", row.names = FALSE, col.names = TRUE)
  # Outputs
  raus_path <- gsub("\\\\", "/", tempfile(pattern = "raus", tmpdir = tempdir(), fileext = ".out"))
  ll_path <- gsub("\\\\", "/", tempfile(pattern = "loglik", tmpdir = tempdir(), fileext = ".out"))
  cont_path <- gsub("\\\\", "/", tempfile(pattern = "continue", tmpdir = tempdir(), fileext = ".out"))
  means_path <- gsub("\\\\", "/", tempfile(pattern = "means", tmpdir = tempdir(), fileext = ".out"))
  tests_path <- gsub("\\\\", "/", tempfile(pattern = "tests", tmpdir = tempdir(), fileext = ".out"))
  rand_path <- gsub("\\\\", "/", tempfile(pattern = "random", tmpdir = tempdir(), fileext = ".out"))
  
  
  # PREPARE INFOFILE
  infofile <- try(get_infofile(model, mdl_txt = mdl_path, mdl_info = model_path))
  if(is(infofile, "try-error")) stop("problem with S4 routines in makin infofile.\n")
  
  
  # Prior Parameter
  prior_names <- c("prec_epsilon", "delta_df", "delta_mu", "delta_scale",
                   "SIGMA_Corr_eta", "SIGMA_SD_rho", "GAMMA_Corr_eta", "GAMMA_SD_rho",
                   "Omega2_alpha", "Omega2_beta")
  if (is.null(prior_params)) prior_params <- list()
  if (!("prec_epsilon" %in% names(prior_params))) prior_params$prec_epsilon <- 1.0
  if (!("delta_df" %in% names(prior_params))) prior_params$delta_df <- 10
  if (!("delta_mu" %in% names(prior_params))) prior_params$delta_mu <- 0.5
  if (!("delta_scale" %in% names(prior_params))) prior_params$delta_scale <- 1.0
  if (!("SIGMA_Corr_eta" %in% names(prior_params))) prior_params$SIGMA_Corr_eta <- 4.0
  if (!("SIGMA_SD_rho" %in% names(prior_params))) prior_params$SIGMA_SD_rho <- 2.5
  if (!("GAMMA_Corr_eta" %in% names(prior_params))) prior_params$GAMMA_Corr_eta <- 4.0
  if (!("GAMMA_SD_rho" %in% names(prior_params))) prior_params$GAMMA_SD_rho <- 0.5
  if (!("Omega2_alpha" %in% names(prior_params))) prior_params$Omega2_alpha <- 0.0025
  if (!("Omega2_beta" %in% names(prior_params))) prior_params$Omega2_beta <- 0.5
  if (!is.numeric(unlist(prior_params))) stop("prior_params must be a list of numerics and integers")
  if (prior_params$delta_df %% 1 != 0) stop("delta_df must be an integer")
  
  
  # PROCESS NAMES AND NUMBER
  proc_names <- names(model$params$threshold)
  nprocs <- length(proc_names)
  
  
  # KERN2FREE
  K2F <- integer(3*nprocs)
  cntk2f <- integer(3)
  consts <- rep(NaN, 3*nprocs)
  comp <- rep(FALSE, 3*nprocs)
  cntcomp <- integer(3)
  for (i in 1:3) {
    cntk2f[i] <- 0
    for (j in 1:nprocs) {
      tmp <- NULL
      if (i == 1) {
        tmp <- model$params$threshold[1, j]
      } else if (i == 2) {
        tmp <- model$params$driftrate[1, j]
      } else if (i == 3) {
        tmp <- model$params$startpoint[1, j]
      }
      
      ind <- (i-1)*nprocs+j
      if(is.na(tmp) | is.numeric(tmp)) {
        K2F[ind] <- cntk2f[i]
        cntk2f[i] <- cntk2f[i] + 1
        if (is.na(tmp)) {
          comp[(cntk2f[i]-1)*3 + i] <- TRUE
        }
        if (is.numeric(tmp)) {
          consts[(cntk2f[i]-1)*3 + i] <- tmp
        }
      } else if(tmp %in% proc_names) {
        K2F[ind] <- K2F[which(proc_names == as.character(tmp)) + (i>=2)*nprocs + (i==3)*nprocs]
      }
    }
  }
  
  
  # PATHS ARGUMENT
  CHAR1 <- as.character(c(Data = data_path,
                          Model = infofile,
                          Raus = raus_path,
                          LogLik = ll_path,
                          Continue = cont_path,
                          Means = means_path,
                          Tests = tests_path,
                          Random = rand_path,
                          tmpdir = paste0(tempdir(), "/")))
  
  
  # MCMC SETTINGS ARGUMENT
  INTEGER1 <- as.integer(c(IREP = Irep,
                           PHASE1 = n.phase1,
                           PHASE2 = n.phase2,
                           THIN = n.thin,
                           NOTHREADS = n.chains,
                           SAMPLE_SIZE = n.chains*n.iter,
                           MAXTHREADS = control$maxthreads,
                           kernpar = nprocs,
                           respno = length(unique(model$responses$MAP)),
                           cat2resp = model$responses$MAP))
  REAL1 <- c(Rhat_max = Rhat_max)
  
  
  # FLAGS ARGUMENT
  BOOL1 <- as.integer(c(DIC_FL = flags$indices,
                        LOG_LIK_FL = flags$loglik,
                        ML_INIT_FL = !flags$random_init))
  
  
  # PRIORS ARGUMENT
  INTEGER2 <- as.integer(c(degf = prior_params$delta_df))
  REAL2 <- c(PRIOR = prior_params$prec_epsilon,
             etat = prior_params$SIGMA_Corr_eta,
             taut = prior_params$SIGMA_SD_rho,
             etar = prior_params$GAMMA_Corr_eta,
             taur = prior_params$GAMMA_SD_rho,
             mu_prior = prior_params$delta_mu,
             rsd = prior_params$delta_scale,
             prioralpha = prior_params$Omega2_alpha,
             priorbeta = prior_params$Omega2_beta)
  
  
  # HMC OPTIONS ARGUMENT
  INTEGER3 <- as.integer(c(maxtreedepth1_3 = control$maxtreedepth1_3,
                           maxtreedepth4 = control$maxtreedepth4))
  
  
  # CONTINUATION ARGUMENT
  INTEGER4 <- as.integer(c(goon = FALSE,
                           ADDITION = 0))
  
  
  # CONSTANTS AND EQUALIZATION
  INTEGER5 <- as.integer(c(kern2free = K2F,
                           comp = comp,
                           ifree = cntk2f))
  REAL3 <- c(consts = consts)
  
  
  # C++ FUNCTION DRTMPT CALL
  out <- .Call("drtmpt_fit", CHAR1, INTEGER1, REAL1, BOOL1,
               INTEGER2, REAL2,
               INTEGER3, INTEGER4,
               REAL3, INTEGER5)
  file.remove(data_path)
  
  
  # DATA INFORMATION
  data_info <- list(Nsubj = length(unique(data_frame[,1])),
                    Ngroups = length(unique(data_frame[,2])),
                    Nthreshold = sum(is.na(model$params$threshold[1,])),
                    Ndriftrate = sum(is.na(model$params$driftrate[1,])),
                    Nstartpoint = sum(is.na(model$params$startpoint[1,])),
                    Nresps = length(unique(model$responses$MAP)),
                    thresh_string = names(model$params$threshold[which(is.na(model$params$threshold[1,]))]),
                    drift_string = names(model$params$driftrate[which(is.na(model$params$driftrate[1,]))]),
                    start_string = names(model$params$startpoint[which(is.na(model$params$startpoint[1,]))]),
                    transformation = transformation)
  
  
  # PREPARE OUTPUT LIST
  drtmpt <- list()
  drtmpt$samples <- make_mcmc_list_d(file = out$pars_samples, infofile = infofile,
                                     Nchains = n.chains, Nsamples = n.iter,
                                     data_info = data_info, keep = old_label)
  file.remove(raus_path)
  
  
  # DIAGNOSTICS
  drtmpt$diags <- get_diags_d(diag_file = tests_path, data_info = data_info, keep = old_label, DIC = flags$indices)
  drtmpt$diags$R_hat <- gelman.diag(drtmpt$samples, multivariate = FALSE)
  
  
  # SPECS
  infos <- readinfofile(infofile)
  drtmpt$specs <- list(model = model, n.chains = n.chains, n.iter = n.iter, n.phase1 = n.phase1, n.phase2 = n.phase2,
                       n.thin = n.thin, n.groups = data_info$Ngroups, n.subj = data_info$Nsubj, Irep = Irep,
                       Rhat_max = Rhat_max, prior_params = prior_params, infolist = infos, call = match.call())
  if(exists("transformation")) {
    drtmpt$specs$transformation <- transformation
  }
  
  
  # COMPUTE WAIC AND LOO-IC
  if (flags$indices) {
    if (file.exists(ll_path)) {
      if (flags$indices) {
        temp <- NULL
        suppressWarnings( temp <- get_indices_x(log_lik = ll_path, Nchains = n.chains, Nsample = n.iter, df = data_frame) )
        drtmpt$indices <- list(WAIC = temp$WAIC, LOO = temp$LOO)
        if (flags$loglik) drtmpt$LogLik <- temp$LogLik
        rm(temp)
      }
      file.remove(ll_path)
    }
  }
  
  
  # CONTINUATION
  if (FALSE) {
    drtmpt$continue <- readLines(con = cont_path)
  }
  file.remove(cont_path)
  
  
  # SUMMARY
  drtmpt$summary <- writeSummaryDRTMPT(x = drtmpt, keep = old_label)
  
  
  # CLEAN-UP
  file.remove(means_path)
  file.remove(tests_path)
  file.remove(rand_path)
  file.remove(infofile)
  file.remove(mdl_path)
  
  
  # OUTPUT
  class(drtmpt) <- "drtmpt_fit"
  gc()
  return(drtmpt)
  
}

