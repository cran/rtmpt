
#' Simulation-based calibration for RT-MPT models
#'
#' Simulate data from RT-MPT models using \code{ertmpt_model} objects. The difference to \code{\link{sim_ertmpt_data}} is that here only scalars are allowed. This makes it usable for
#'   simulation-based calibration (SBC; Talts et al., 2018). You can specify the random seed, number of subjects, number of trials, and some
#'   parameters (same as \code{prior_params} from \code{\link{fit_ertmpt}}).
#'
#' @param model A list of the class \code{ertmpt_model}.
#' @param seed Random seed number.
#' @param n.eff_samples Number of effective samples. Default is 99, leading to 100 possible ranks (from 0 to 99).
#' @param n.chains Number of chains to use. Default is 4. Must be larger than 1 and smaller or equal to 16.
#' @param n.iter Number of samples per chain. Default is 5000. Must be larger or equal to \code{n.eff_samples}.
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
#' @param n.subj Number of subjects. Default is 40.
#' @param n.trials Number of trials per tree. Default is 30.
#' @param prior_params Named list of parameters from which the data will be generated. This must be the same named list as \code{prior_params} from
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
#' @param sim_list Object of class \code{ertmpt_sim}. This is also an output object. Can be used to re-fit the model if \code{n.eff_samples} was not achieved in a previous fitting attempt.
#'   It will then use the data stored in this object. Default is NULL and this object will be created anew.
#' @return A list of the class \code{ertmpt_sbc} containing
#'   \itemize{
#'     \item \code{ranks}: the rank statistic for all parameters,
#'     \item \code{sim_list}: an object of the class \code{ertmpt_sim},
#'     \item \code{fit_list}: an object of the class \code{ertmpt_fit},
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
#' d+(1-d)*g     ; 0
#' (1-d)*(1-g)   ; 1
#'
#' # lures
#' (1-d)*g       ; 0
#' d+(1-d)*(1-g) ; 1
#'
#' # d: detect; g: guess
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
#' \donttest{
#' R = 2 # typically 2000 with n.eff_samples = 99, but this will run many days
#' rank_mat <- matrix(NA, ncol = 393, nrow = 2)
#' for (r in 1:R) {
#'   SBC_out <- fit_ertmpt_SBC(model, seed = r*123, prior_params = params,
#'                            n.eff_samples = 99, n.thin = 5,
#'                            n.iter = 5000, n.burnin = 2000, Irep = 5000)
#'   rank_mat[r, ] <- SBC_out$ranks
#' }
#' }
#' 
#' @author Raphael Hartmann
#' @export
#' @importFrom coda effectiveSize varnames
fit_ertmpt_SBC <- function(model,
                          seed,
                          n.eff_samples = 99,
                          n.chains = 4,
                          n.iter = 5000,
                          n.burnin = 200,
                          n.thin = 1,
                          Rhat_max = 1.05,
                          Irep = 1000,
                          n.subj = 40,
                          n.trials = 30,
                          prior_params = NULL,
                          sim_list = NULL) {


  model_elmnts <- c("lines", "params", "responses")

  if (is.null(sim_list)) {
    # some controls ----

    if (!is.list(model)) stop("\"model\" must be a list.")
    if (!all(model_elmnts  %in% names(model))) stop("\"model\" must contain \"", model_elmnts[which(!(model_elmnts %in% names(model)))[1]], "\".")

    if (!is.numeric(seed)) stop("\"seed\" must be numeric.")

    if (n.subj < 2) stop("\"n.subj\" must be larger than or equal to two.")
    if (n.trials < 2) stop("\"n.trials\" must be larger than or equal to two.")
    if (n.trials < 30) warning("\"n.trials\" is recommended to be larger than 30.")

    if (!is.null(prior_params) && !is.list(prior_params)) stop("\"prior_params\" must be a list.")

    if (!is.numeric(n.eff_samples) | n.eff_samples%%1!=0 | n.eff_samples<1) stop("\"n.eff_samples\" must be a natural number larger than one.")

    # generate data ----
    sim_list <- sim_ertmpt_data_SBC(model = model, seed = seed, n.subj = n.subj, n.trials = n.trials, params = prior_params)

  } else {
    model <- sim_list$specs$model

    seed <- sim_list$specs$seed

    n.subj <- sim_list$specs$n.subj

    n.trials <- sim_list$specs$n.trials

    if (!is.null(prior_params) && !is.list(prior_params)) stop("\"prior_params\" must be a list.")

    if (!is.numeric(n.eff_samples) | n.eff_samples%%1!=0 | n.eff_samples<1) stop("\"n.eff_samples\" must be a natural number larger than one.")

    if (!inherits(sim_list, c("ertmpt_sim", "rtmpt_sim"))) stop("\"sim_list\" is not of class \"ertmpt_sim\".")
  }


  # fitting model to data ----
  fit_list <- fit_ertmpt(model = sim_list$specs$model, data = sim_list$data_frame, n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin,
                        n.thin = n.thin, Rhat_max = Rhat_max, Irep = Irep, prior_params = prior_params,
                        indices = FALSE, save_log_lik = FALSE)


  # preparations
  cat("\n\ncalculating minimal effective sample size: ")
  min_eff <- min(coda::effectiveSize(fit_list$samples))
  cat(min_eff, "\n\n")

  nprobs <- sum(is.na(model$params$probs[1,]))
  nminus <- sum(is.na(model$params$taus[1,]))
  nplus <- sum(is.na(model$params$taus[2,]))
  npars <- nprobs + nminus + nplus
  nresp <- length(unique(model$responses$MAP))

  nparams <- length(coda::varnames(fit_list$samples)) - 1
  rankmat <- matrix(NA, ncol = nparams, nrow = 1)

  sample_L <- function() {

    chains <- fit_list$samples[,-(nparams+1)]
    varnam <- coda::varnames(chains)
    samp_all <- list()

    nn <- floor(n.chains*n.iter / n.eff_samples)

    samp_ind <- seq(nn, n.chains*n.iter, by = nn)
    if (length(samp_ind) > n.eff_samples) samp_ind <- samp_ind[1:n.eff_samples]
    # samp <- matrix(NA, nrow = n.eff_samples, ncol = nparams)
    chains_mat <- as.matrix(chains[[1]])
    for (i in 2:n.chains) {
      chains_mat <- rbind(chains_mat, as.matrix(chains[[i]]))
    }
    samp <- chains_mat[samp_ind,]

    # transform std dev to var
    if (npars == 1) mat_corrAB <- as.matrix(samp[, ind_ncorrAB]) else mat_corrAB <- samp[, ind_ncorrAB]
    samp[, ind_ncorrAB] <- StddevCorr2Cov(mat_corrAB, npars)
    if (nresp == 1) mat_corrC <- as.matrix(samp[, ind_ncorrC]) else mat_corrC <- samp[, ind_ncorrC]
    samp[, ind_ncorrC] <- StddevCorr2Cov(mat_corrC, nresp)

    colnames(samp) <- varnam
    samp_all$samp <- samp
    return(samp_all)
  }


  # calculate ranks ----
  n_run_low <- 1
  n_run_upp <- sum(is.na(model$params$probs[1,]))
  ind_nprobs <- n_run_low:n_run_upp

  n_run_low <- n_run_upp + 1
  n_run_upp <- n_run_upp + sum(is.na(model$params$taus[1,]))
  ind_nminus <- n_run_low:n_run_upp

  n_run_low <- n_run_upp + 1
  n_run_upp <- n_run_upp + sum(is.na(model$params$taus[2,]))
  ind_nplus <- n_run_low:n_run_upp

  n_run_low <- n_run_upp + 1
  n_run_upp <- n_run_upp + npars*(npars+1)/2
  ind_ncorrAB <- n_run_low:n_run_upp

  n_run_low <- n_run_upp + 1
  n_run_upp <- n_run_upp + n.subj*nprobs
  ind_nalpha_prim <- n_run_low:n_run_upp

  n_run_low <- n_run_upp + 1
  n_run_upp <- n_run_upp + n.subj*(nminus+nplus)
  ind_nbeta_prim <- n_run_low:n_run_upp

  n_run_low <- n_run_upp + 1
  n_run_upp <- n_run_upp + nresp
  ind_nresp <- n_run_low:n_run_upp

  n_run_low <- n_run_upp + 1
  n_run_upp <- n_run_upp + 1
  ind_nomega2 <- n_run_low:n_run_upp

  n_run_low <- n_run_upp + 1
  n_run_upp <- n_run_upp + nresp*(nresp+1)/2
  ind_ncorrC <- n_run_low:n_run_upp

  n_run_low <- n_run_upp + 1
  n_run_upp <- n_run_upp + n.subj*nresp
  ind_ngamma_prim <- n_run_low:n_run_upp

  n_run_low <- n_run_upp + 1
  n_run_upp <- n_run_upp + n.subj
  ind_nsigma2 <- n_run_low:n_run_upp


  StddevCorr2Cov <- function(mat, dim2) {
    if ((dim2*(dim2-1)/2+dim2) != dim(mat)[2]) stop("wrong dim2 or 2nd dim of mat!")

    matcov <- matrix(NA, ncol = dim(mat)[2], nrow = dim(mat)[1])
    if (dim2 == 1) {
      matcov[,1] <- mat[,1]*mat[,1]
    } else {
      SDCorr <- matrix(0, ncol = dim2, nrow = dim2)
      for (i in 1:dim(mat)[1]) {
        SDCorr[lower.tri(diag(dim2), diag = TRUE)] <- mat[i,]
        SDs <- diag(SDCorr)
        Corr <- SDCorr - diag(SDs) + diag(dim2)
        upperCorr <- t(Corr)-diag(diag(Corr))
        Corr <- Corr + upperCorr
        Cov <- (SDs %*% t(SDs)) * Corr
        matcov[i,] <- Cov[lower.tri(Cov, diag = TRUE)]
      }
    }
    return(matcov)

  }

  sample_L <- function() {

    chains <- fit_list$samples[,-(nparams+1)]
    varnam <- coda::varnames(chains)
    samp_all <- list()

    nn <- floor(n.chains*n.iter / n.eff_samples)

    samp_ind <- seq(nn, n.chains*n.iter, by = nn)
    if (length(samp_ind) > n.eff_samples) samp_ind <- samp_ind[1:n.eff_samples]
    # samp <- matrix(NA, nrow = n.eff_samples, ncol = nparams)
    chains_mat <- as.matrix(chains[[1]])
    for (i in 2:n.chains) {
      chains_mat <- rbind(chains_mat, as.matrix(chains[[i]]))
    }
    samp <- chains_mat[samp_ind,]

    # transform std dev to var
    if (npars == 1) mat_corrAB <- as.matrix(samp[, ind_ncorrAB]) else mat_corrAB <- samp[, ind_ncorrAB]
    samp[, ind_ncorrAB] <- StddevCorr2Cov(mat_corrAB, npars)
    if (nresp == 1) mat_corrC <- as.matrix(samp[, ind_ncorrC]) else mat_corrC <- samp[, ind_ncorrC]
    samp[, ind_ncorrC] <- StddevCorr2Cov(mat_corrC, nresp)

    colnames(samp) <- varnam
    samp_all$samp <- samp
    return(samp_all)
  }

  groundtruth <- function() {
    GT <- sim_list$gen_list
    grnd_trth <- GT$process_list$mu_alpha
    grnd_trth <- c(grnd_trth, exp(GT$process_list$mu_beta))
    l_zetaAB <- length(GT$process_list$zeta)
    grnd_trth <- c(grnd_trth, (diag(GT$process_list$zeta, nrow = l_zetaAB) %*% GT$process_list$S_doubleprime %*% diag(GT$process_list$zeta, nrow = l_zetaAB))[FALSE==upper.tri(diag(npars))])
    grnd_trth <- c(grnd_trth, as.vector(GT$process_list$primes[1:nprobs,]))
    grnd_trth <- c(grnd_trth, as.vector(GT$process_list$primes[(nprobs+1):npars,]))
    grnd_trth <- c(grnd_trth, GT$motor_list$mu_gamma)
    grnd_trth <- c(grnd_trth, GT$motor_list$omega_square)
    l_zetaC <- length(GT$motor_list$zeta)
    grnd_trth <- c(grnd_trth, (diag(GT$motor_list$zeta, nrow = l_zetaC) %*% GT$motor_list$S_doubleprime %*% diag(GT$motor_list$zeta, nrow = l_zetaC))[FALSE==upper.tri(diag(nresp))])
    grnd_trth <- c(grnd_trth, as.vector(GT$motor_list$primes))
    grnd_trth <- c(grnd_trth, GT$motor_list$sigma_square)
    return(grnd_trth)
  }

  # rankmat <- matrix(NA, ncol = nparams, nrow = 1)
  sample <- sample_L()$samp
  parnames <- colnames(sample)
  grnd_trth <- groundtruth()
  colnames(rankmat) <- parnames
  if (min_eff < n.eff_samples) {
    warning("minimal effective sample size is smaller than \"n.eff_samples\"! Ranks will not be calculated! Try with increased \"n.thin\".")
    rankmat[1, ] <- rep(NA, nparams)
  } else {
    for (p in 1:nparams) {
      rankmat[1, p] <- sum(unlist(apply(X = sample, MARGIN = 1, FUN = function(x) {x[p]})) < grnd_trth[p])
    }
  }



  # specs ----
  specs <- list(n.eff_samples = n.eff_samples, call = match.call())



  # output ----
  sbc_list <- list(ranks = rankmat, sim_list = sim_list, fit_list = fit_list, specs = specs)

  class(sbc_list) <- "ertmpt_sbc"

  return(sbc_list)

}
