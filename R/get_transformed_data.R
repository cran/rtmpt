
#' Transform data to be used in RT-MPT model fitting
#'
#' Transform data, such that it can be used in \code{\link{fit_rtmpt}}. This implies changing each value/label in
#'   "subj", "group", "tree", and "cat" to numbers such that it starts from zero (e.g. data$tree = c(1,1,3,3,2,2,...)
#'   will be changed to data$tree = c(0,0,2,2,1,1,...)) and the columns will be ordered in the right way.
#'   "rt" must be provided in milliseconds. If it has decimal places it will be rounded to a whole number.
#'   \code{\link{fit_rtmpt}} will automatically call this function if its input is not already an \code{rtmpt_data} list, 
#'   but it is advised to use it anyway because it provides information about the transformations of the data.
#'
#' @param raw_data \code{data.frame} or path to data containing columns "subj", "group", "tree", "cat", and "rt". 
#'   If not provided in this order it will be reordered and unused variables will be moved to the end of the
#'   new data frame.
#' @param model A list of the class \code{rtmpt_model}.
#' @return A list of the class \code{rtmpt_data} containing transformed data and information about 
#'   the transformation that has been done.
#' @examples 
#' ####################################################################################
#' # Detect-Guess variant of the Two-High Threshold model.
#' # The encoding and motor execution times are assumed to be equal for each response.
#' ####################################################################################
#' 
#' eqn_2HTM <- "
#' # CORE MPT EQN
#' # tree ; cat  ; mpt
#' target ; hit  ; do
#' target ; hit  ; (1-do)*g
#' target ; miss ; (1-do)*(1-g)
#'        
#'   lure ;  f_a ; (1-dn)*g
#'   lure ;  c_r ; dn
#'   lure ;  c_r ; (1-dn)*(1-g)
#' "
#' 
#' model <- to_rtmpt_model(eqn_file = eqn_2HTM)
#' 
#' data_file <- system.file("extdata/labeled_data.txt", package="rtmpt")
#' data <- read.table(file = data_file, header = TRUE)
#' 
#' data_list <- to_rtmpt_data(raw_data = data, model = model)
#' data_list 
#' 
#' @author Raphael Hartmann
#' @export
to_rtmpt_data <- function(raw_data, model) {
  
  if (is.data.frame(raw_data)) {raw_data <- raw_data
  } else if (is.character(raw_data)) raw_data <- read.table(file = raw_data, header = TRUE)
  
  data_elmnts <- c("subj", "group", "tree", "cat", "rt")
  if (!all(data_elmnts %in% names(raw_data))) stop("\"raw_data\" must contain \"", data_elmnts[which(!(data_elmnts %in% names(raw_data)))[1]], "\".")
  if (!all(data_elmnts == names(raw_data)[1:5])) {
    df <- raw_data
    unused_ind <- which(!(names(df) %in% data_elmnts))
    raw_data <- df[,c(match(data_elmnts, names(raw_data)), unused_ind)]
  }
  if (any(is.na(raw_data)) || min(raw_data$rt) < 0) stop("All values in \"raw_data\" must not be NA and \"rt\" must be larger than zero.")
  
  if (is.factor(raw_data$subj)) raw_data$subj <- as.character(raw_data$subj)
  if (is.factor(raw_data$group)) raw_data$group <- as.character(raw_data$group)
  if (is.factor(raw_data$tree)) raw_data$tree <- as.character(raw_data$tree)
  if (is.factor(raw_data$cat)) raw_data$cat <- as.character(raw_data$cat)
  
  df <- raw_data
  # colnames(df) <- data_elmnts
  
  df_ord <- df[with(df, order(df[,1], df[,2], df[,3], df[,4])), ]
  
  old_new <- list()
  
  
  # SUBJECTS
  if ( (min(df_ord$subj) != 0) || (length(unique(df_ord$subj)) != max(df_ord$subj+1)) ) {
    uniq_subj <- unique(df_ord$subj)
    corr_subj <- seq(0, (length(uniq_subj)-1))
    df_ord$subj <- sapply(X = df_ord$subj, FUN = function(x) {ind <- which(uniq_subj == x); return(corr_subj[ind])})
    old_new$subj <- data.frame(old = uniq_subj, new = corr_subj)
  }
  
  # GROUPS
  if ( (min(df_ord$group) != 0) || (length(unique(df_ord$group)) != max(df_ord$group+1)) ) {
    uniq_group <- unique(df_ord$group)
    corr_group <- seq(0, (length(uniq_group)-1))
    df_ord$group <- sapply(X = df_ord$group, FUN = function(x) {ind <- which(uniq_group == x); return(corr_group[ind])})
    old_new$group <- data.frame(old = uniq_group, new = corr_group)
  }
  
  # TREES
  if (length(unique(model$responses$TREE)) != length(unique(df_ord$tree))) stop("Number of trees in \"model\" and \"raw_data\" do not match.")
  
  if (is.numeric(model$responses$TREE) && is.numeric(df_ord$tree)) {
    
    if ( (min(df_ord$tree) != 0) || (length(unique(df_ord$tree)) != max(df_ord$tree+1)) ) {
      uniq_tree <- unique(df_ord$tree)
      corr_tree <- seq(0, (length(uniq_tree)-1))
      df_ord$tree <- sapply(X = df_ord$tree, FUN = function(x) {ind <- which(uniq_tree == x); return(corr_tree[ind])})
      old_new$tree <- data.frame(old = uniq_tree, new = corr_tree)
    }
    
  } else if (is.numeric(model$responses$TREE) && is.character(df_ord$tree)) {
    
    stop("Please change the tree-column in your \"raw_data\" frame to numerics, such that it matches the tree-column in the \"model\".")
  
  } else if (is.character(model$responses$TREE) && is.character(df_ord$tree)) {
    
    if (!all(unique(model$responses$TREE) %in% unique(df_ord$tree))) {stop("\"raw_data\" and \"model\" cannot be matched since names of trees in \"raw_data\" and \"model\" differ.")
    } else {
      uniq_tree <- unique(model$responses$TREE)
      corr_tree <- seq(0, (length(uniq_tree)-1))
      df_ord$tree <- sapply(X = df_ord$tree, FUN = function(x) {ind <- which(uniq_tree == x); return(corr_tree[ind])})
      old_new$tree <- data.frame(old = uniq_tree, new = corr_tree)
    }
    
  } else if (is.character(model$responses$TREE) && is.numeric(df_ord$tree)) {
    
    if(suppressWarnings(all(!is.na(as.numeric(model$responses$TREE))))) {
      model$responses$TREE <- as.numeric(model$responses$TREE)
    } else warning("Found tree-column in \"model\" to be character and tree-column of \"raw_data\" to be numeric. It will be assumed that sort(x = unique(raw_data$tree), decreasing = FALSE) correspond to unique(model$responses$TREE). 
                   In other words it will be assumed, that your model specification starts with the tree, that has the lowest number in \"raw_data\" and so on and ends with the tree, that has the highest number.")
    
    if ( (min(df_ord$tree) != 0) || (length(unique(df_ord$tree)) != max(df_ord$tree+1)) ) {
      uniq_tree <- unique(df_ord$tree)
      corr_tree <- seq(0, (length(uniq_tree)-1))
      df_ord$tree <- sapply(X = df_ord$tree, FUN = function(x) {ind <- which(uniq_tree == x); return(corr_tree[ind])})
      old_new$tree <- data.frame(old = uniq_tree, new = corr_tree)
    }
    
  }
  
  # CATEGORIES
  if (length(unique(model$responses$CAT)) != length(unique(df_ord$cat))) stop("Number of categories in \"model\" and \"raw_data\" do not match.")
  
  if (is.numeric(model$responses$CAT) && is.numeric(df_ord$cat)) {
    
    if ( (min(df_ord$cat) != 0) || (length(unique(df_ord$cat)) != max(df_ord$cat+1)) ) {
      uniq_cat <- unique(df_ord$cat)
      corr_cat <- seq(0, (length(uniq_cat)-1))
      df_ord$cat <- sapply(X = df_ord$cat, FUN = function(x) {ind <- which(uniq_cat == x); return(corr_cat[ind])})
      old_new$cat <- data.frame(old = uniq_cat, new = corr_cat)
    }
    
  } else if (is.numeric(model$responses$CAT) && is.character(df_ord$cat)) {
    
    stop("Please change the cat-column in your data frame to numerics, such that it matches the cat-column in the \"model\".")
    
  } else if (is.character(model$responses$CAT) && is.character(df_ord$cat)) {
    
    if (!all(unique(model$responses$CAT) %in% unique(df_ord$cat))) {stop("\"raw_data\" and \"model\" cannot be matched since names of categories in \"raw_data\" and \"model\" differ.")
    } else {
      uniq_cat <- unique(model$responses$CAT)
      corr_cat <- seq(0, (length(uniq_cat)-1))
      df_ord$cat <- sapply(X = df_ord$cat, FUN = function(x) {ind <- which(uniq_cat == x); return(corr_cat[ind])})
      old_new$cat <- data.frame(old = uniq_cat, new = corr_cat)
    }
    
  } else if (is.character(model$responses$CAT) && is.numeric(df_ord$cat)) {
    
    if(suppressWarnings(all(!is.na(as.numeric(model$responses$CAT))))) {
      model$responses$CAT <- as.numeric(model$responses$CAT)
    } else warning("Found cat-column in \"model\" to be character and cat-column of \"raw_data\" to be numeric. It will be assumed that sort(x = unique(raw_data$cat), decreasing = FALSE) correspond to unique(model$responses$CAT). 
                   In other words it will be assumed, that your model specification starts with the category, that has the lowest number in \"raw_data\" and so on and ends with the category, that has the highest number.")
    if ( (min(df_ord$cat) != 0) || (length(unique(df_ord$cat)) != max(df_ord$cat+1)) ) {
      uniq_cat <- unique(df_ord$cat)
      corr_cat <- seq(0, (length(uniq_cat)-1))
      df_ord$cat <- sapply(X = df_ord$cat, FUN = function(x) {ind <- which(uniq_cat == x); return(corr_cat[ind])})
      old_new$cat <- data.frame(old = uniq_cat, new = corr_cat)
    }
    
  }
  
  # RESPONSE TIME
  if (!all((df_ord$rt %% 1) != 0)) {
    df_ord$rt <- round(df_ord$rt)
  }
  # if (max(df_ord$rt) < 100) {
  #   df_ord$rt <- df_ord$rt*1000
  # }
  
  data_list <- list(data=df_ord)
  if (length(old_new) > 0) data_list$transformation <- old_new
  
  
  class(data_list) <- "rtmpt_data"
  
  return(data_list)
}

#' @export
print.rtmpt_data <- function(x, ...) {
  cat("\nDATA TRANSFORMATION OVERVIEW\n\n")
  
  cat("\nReordered variables:\nsubj, group, tree, cat, rt\n")
  cat("* NOTE1: Additional variables are attached next to these five.\n")
  cat("* NOTE2: To see your data frame use <object name>$data.\n")
  cat("--------------------\n")
  
  cat("\nTransformed variable(s):")
  if("transformation" %in% names(x)) {
    if("subj" %in% names(x$transformation)) {
      cat("\n\"subj\"\n")
      print(x$transformation$subj)
    }
    if("group" %in% names(x$transformation)) {
      cat("\n\"group\"\n")
      print(x$transformation$group)
    }
    if("tree" %in% names(x$transformation)) {
      cat("\n\"tree\"\n")
      print(x$transformation$tree)
    }
    if("cat" %in% names(x$transformation)) {
      cat("\n\"cat\"\n")
      print(x$transformation$cat)
    }
    cat("\n* NOTE: \"old\" refers to the used labels and \"new\" to the ones that will be used.\n")
  } else cat("* NOTE: No transformations needed.\n")
  cat("------------------------\n\n")
}

