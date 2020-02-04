
#' Set responses in a \code{rtmpt_model}
#' 
#' Change the responses for a tree and the categories within that tree.
#'
#' @param model A list of the class \code{rtmpt_model}.
#' @param tree Character or numerical value of the tree for which the responses 
#'   should be changed.
#' @param categories Character or numerical vector identifying category/ies within 
#'   the specified \code{tree} for which the responses should be changed.
#' @param values Numerical vector of length \code{length(categories)} providing the responses. Default is 0.
#' @return A list of the class \code{rtmpt_model}.
#' @examples
#' #########################################################################
#' # Detect-Guess variant of the Two-High Threshold model.
#' # The encoding and motor execution times will be set to different values 
#' #   for each response.
#' #########################################################################
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
#' model <- to_rtmpt_model(mdl_file = mdl_2HTM)
#' 
#' ## changing the model to have two different encoding and response execution 
#' ## times for "old" and "new" responses.
#' for(i in c(0,1)) model <- set_resps(model = model, tree = i, 
#'                                     categories = i*2+1, values = 1)
#'                                  
#' @seealso \code{\link{set_params}}
#' @author Raphael Hartmann
#' @export
set_resps <- function(model, tree, categories, values = 0) {
  
  
  if (!is.list(model)) stop("model must be a list.")
  if (!("lines" %in% names(model)) || !("params" %in% names(model)) || !("responses" %in% names(model))) stop("No valid model file.")
  resps_df <- model$responses
  if ( !("RESP" %in% names(resps_df)) ||
       !("TREE" %in% names(resps_df)) ||
       !("CAT" %in% names(resps_df)) ) stop("\"responses\" must contain columns \"TREE\", \"CAT\" and \"RESP\".")
  # if (!is.character(tree)) stop("\"tree\" must be a character.")
  if (length(tree) != 1) stop("\"tree\" must be of length one.")
  if (!is.vector(categories)) stop("\"categories\" must be a vector")
  if (!is.vector(values)) stop("\"values\" must be a vector")
  if (!(tree %in% resps_df$TREE)) stop("\"tree\" name must match names used in \"model$responses$TREE\".")
  # if (!is.character(categories)) stop("\"categories\" must be characters.")
  if (!is.numeric(values)) stop("\"values\" must be numerical.")
  if (length(values)!=1 & length(categories) != length(values)) stop("Length of \"categories\" and \"values\" must match.")
  if (length(values)==1 & length(categories) != length(values)) values <- rep(values, length(categories))
  if (any(!(categories %in% resps_df$CAT))) stop("\"categories\" do not match with the names used in \"model$responses$CAT\".")
  if (!all(categories == unique(categories))) {
    uniq_cats <- unique(categories)
    index <- sapply(X = uniq_cats, FUN = function(x) {which(categories == x)[1]})
    categories <- categories[index]
    values <- values[index]
    warning("\"categories\" is not a unique vector. Transformed \"categories\" into a unique vector and changed \"values\" accordingly.")
  } 
  
  tree_ind <- which(resps_df$TREE == tree)
  cat_ind <- unique(sapply(X = categories, FUN = function(x) {which(resps_df$CAT == x)}))
  ind <- intersect(cat_ind, tree_ind)
  resps_df$RESP[ind] <- values
  
  if(length(unique(resps_df$RESP)) != (max(resps_df$RESP)+1)) {stop("The elements of the responses in the model must range from 0 to max(\"model$responses$RESP\")")}
  
  
  model$responses <- resps_df
  return(model)
  
}
