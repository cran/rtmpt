
#' Set mapping between response categories and encoding plus motor execution times
#' 
#' Mapping response categories with encoding and motor execution times (deltas). Unlike the processes there are no names for 
#'   the different deltas and therefore a mapping from response categories to different deltas must be specified.
#'
#' @param model A list of the class \code{ertmpt_model} or \code{drtmpt_model}.
#' @param trees Character or numerical vector giving the trees
#' @param categories Character or numerical vector identifying category/ies within 
#'   the specified \code{trees} for which the deltas should be changed.
#' @param mappings Numerical vector of length \code{length(categories)} providing the mappings. Default is 0.
#' @return A list of the class \code{ertmpt_model}.
#' @examples
#' ###########################################################################
#' # Detect-Guess variant of the Two-High Threshold model.
#' # The encoding and motor execution times will be set to different responses 
#' ###########################################################################
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
#' ## changing the model to have two different encoding and motor execution 
#' ## times for "old" and "new" responses.
#' new_model <- delta2delta(model = model, trees = c(0, 1), 
#'                          categories = c(1,3), mappings = c(1,1))
#' new_model
#' 
#' 
#' model <- to_drtmpt_model(mdl_file = mdl_2HTM)
#' 
#' ## changing the model to have two different encoding and motor execution 
#' ## times for "old" and "new" responses.
#' new_model <- delta2delta(model = model, trees = c(0, 1), 
#'                          categories = c(1,3), mappings = c(1,1))
#' new_model
#'                                  
#' @seealso \code{\link{theta2const}}, \code{\link{tau2zero}}, \code{\link{theta2theta}}, and \code{\link{tau2tau}}, 
#' @author Raphael Hartmann
#' @export
delta2delta <- function(model, trees, categories, mappings = 0) {
  
  
  if (!inherits(model, c("ertmpt_model", "rtmpt_model", "drtmpt_model"))) stop("model must be of class \"ertmpt_model\", or \"drtmpt_model\".")
  if (!("lines" %in% names(model)) || !("params" %in% names(model)) || !("responses" %in% names(model))) stop("No valid model file.")

  resps_df <- model$responses
  if ( !("MAP" %in% names(resps_df)) ||
       !("TREE" %in% names(resps_df)) ||
       !("CAT" %in% names(resps_df)) ) stop("\"model$responses\" must contain columns \"TREE\", \"CAT\" and \"MAP\".")
  if (!is.vector(trees)) stop("\"trees\" must be a vector")
  if (!is.vector(categories)) stop("\"categories\" must be a vector")
  if (!is.vector(mappings)) stop("\"mappings\" must be a vector")
 
  if (any(!(trees %in% resps_df$TREE))) stop("\"trees\" names must match names used in \"model$responses$TREE\".")
  if (!is.numeric(mappings)) stop("\"mappings\" must be numerical.")

  if (length(mappings)!=1 & length(trees) != length(mappings)) stop("Length of \"trees\" and \"mappings\" must match or \"mappings\" must be of length one.")
  if (length(mappings)!=1 & length(categories) != length(mappings)) stop("Length of \"categories\" and \"mappings\" must match or \"mappings\" must be of length one.")
  if (length(trees)!=length(categories)) stop("Length of \"trees\" and \"categories\" must match")
  if (length(mappings)==1 & length(categories) != length(mappings)) mappings <- rep(mappings, length(categories))
  
  if (any(!(categories %in% resps_df$CAT))) stop("\"categories\" do not match with the names used in \"model$responses$CAT\".")
  if (!all(trees == unique(trees))) {
    uniq_trees <- unique(trees)
    index <- sapply(X = uniq_trees, FUN = function(x) {which(trees == x)[1]})
    trees <- trees[index]
    categories <- categories[index]
    mappings <- mappings[index]
    warning("\"trees\" is not a unique vector. Transformed \"trees\" into a unique vector and changed \"mappings\" accordingly.")
  } 

  combo <- data.frame(trees = trees, categories = categories, mappings = mappings)
  
  for (i in 1:length(trees)) {
    if (any(sapply(1:dim(resps_df)[1], function(j) all(resps_df[j, c(1,2)]==combo[i, 1:2])))) {
      ind <- which(sapply(1:dim(resps_df)[1], function(j) all(resps_df[j, c(1,2)]==combo[i, 1:2])))
      resps_df[ind, "MAP"] <- mappings[i]
    } else {
      warning("at least one combination of \"trees\" and \"categories\" does not exist and will be skipped.")
    }
    
  }
  
  if(length(unique(resps_df$MAP)) != (max(resps_df$MAP)+1)) {stop("The elements of the \"mappings\" in the model must start at 0 and always be increased by only one")}
  
  model$responses <- resps_df
  return(model)
  
}


#' @rdname delta2delta
#' @examples
#' 
#' ## changing the model to have two different encoding and response execution 
#' ## times for "old" and "new" responses.
#' new_model <- set_deltas_equal(model = model, trees = c(0, 1), 
#'                               categories = c(1,3), mappings = c(1,1))
#' new_model
#' @export
set_deltas_equal <- delta2delta














#' Set responses in an \code{ertmpt_model} or a \code{drtmpt_model}
#' 
#' Change the responses for a tree and the categories within that tree.
#'
#' @param model A list of the class \code{ertmpt_model} or \code{drtmpt_model}.
#' @param tree Character or numerical value of the tree for which the responses 
#'   should be changed.
#' @param categories Character or numerical vector identifying category/ies within 
#'   the specified \code{tree} for which the responses should be changed.
#' @param values Numerical vector of length \code{length(categories)} providing the responses. Default is 0.
#' @return A list of the class \code{ertmpt_model} or \code{drtmpt_model}.
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
#' model <- to_ertmpt_model(mdl_file = mdl_2HTM)
#' 
#' ## changing the model to have two different encoding and response execution 
#' ## times for "old" and "new" responses.
#' for(i in c(0,1)) model <- set_resps(model = model, tree = i, 
#'                                     categories = i*2+1, values = 1)
#' 
#' 
#' #' model <- to_drtmpt_model(mdl_file = mdl_2HTM)
#' 
#' ## changing the model to have two different encoding and response execution 
#' ## times for "old" and "new" responses.
#' for(i in c(0,1)) model <- set_resps(model = model, tree = i, 
#'                                     categories = i*2+1, values = 1)
#'                                  
#' @author Raphael Hartmann
#' @export
set_resps <- function(model, tree, categories, values = 0) {
  
  warning("this function is deprecated and will be removed in a future version. Please use delta2delta instead.")
  
  if (!is.list(model)) stop("model must be a list.")
  if (!("lines" %in% names(model)) || !("params" %in% names(model)) || !("responses" %in% names(model))) stop("No valid model file.")
  resps_df <- model$responses
  if ( !("MAP" %in% names(resps_df)) ||
       !("TREE" %in% names(resps_df)) ||
       !("CAT" %in% names(resps_df)) ) stop("\"model$responses\" must contain columns \"TREE\", \"CAT\" and \"MAP\".")
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
  resps_df$MAP[ind] <- values
  
  if(length(unique(resps_df$MAP)) != (max(resps_df$MAP)+1)) {stop("The elements of the responses in the model must range from 0 to max(\"model$responses$MAP\")")}
  
  
  model$responses <- resps_df
  return(model)
  
}
