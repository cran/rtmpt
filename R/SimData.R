#' Data simulated from the restricted 2HTM
#'
#' Data set generated from a restricted Two-High Threshold model.
#' @details
#' Fourty subjects with thirty trials per condition (Studied items, new Items) were simulated.
#' @format A data frame with five variables:
#' \describe{
#'   \item{\code{subj}}{subjects number}
#'   \item{\code{group}}{group label of the subjects}
#'   \item{\code{tree}}{condition of the current trial}
#'   \item{\code{cat}}{observed response category}
#'   \item{\code{rt}}{observed response time in ms}
#' }
#' @examples
#' ###################################################################
#' # Detect-Guess variant of the restricted Two-High Threshold model.
#' ###################################################################
#'
#' head(SimData)
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
#' model <- to_rtmpt_model(mdl_file = mdl_2HTM)
#'
#' data <- to_rtmpt_data(raw_data = SimData, model = model)
#' \donttest{
#' # this might take some time to run
#' rtmpt_out <- fit_rtmpt(model = model, data = data)
#'
#' # convergence
#' ## traceplot and summary of the first six parameters
#' plot(rtmpt_out$samples[,1:6])
#' summary(rtmpt_out)
#' }
"SimData"
