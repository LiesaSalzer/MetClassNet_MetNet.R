library(tidyverse)
library(MetNet)
library(igraph)
library(reshape2)
library(Hmisc)
setwd("/Users/Liesa4/Library/Mobile Documents/com~apple~CloudDocs/Promotion/R/MetClassNet/MetNet")


#' Changes to MetNet:
#' structural() has additional list entry of matrix containing mass values of respective matches
structural <- function(x, transformation, ppm = 5, directed = FALSE) {
  
  if (!is.data.frame(transformation))
    stop("transformation is not a data.frame")
  if (!"group" %in% colnames(transformation))
    stop("transformation does not contain the column group")
  if (!"mass" %in% colnames(transformation))
    stop("transformation does not contain the column mass")
  if (!"mz" %in% colnames(x)) stop("x does not contain the column mz")
  
  if (!is.numeric(ppm)) stop("ppm is not numeric")
  
  mass <- x[, "mz"]
  mat <- matrix(0, nrow = length(mass), ncol = length(mass))
  rownames(mat) <- colnames(mat) <- mass
  
  ## create matrix which has rowmames per row
  mat <- apply(mat, 1, function(x) as.numeric(mass))
  
  ## calculate ppm deviation
  mat_1 <- mat / abs(ppm / 10 ^ 6 + 1)
  mat_2 <- mat / abs(ppm / 10 ^ 6 - 1)
  
  ## calculate difference between rownames and colnames
  ## (difference between features)
  
  mat_1 <- mat - t(mat_1) ## max
  mat_2 <- mat - t(mat_2) ## min
  
  if (!directed) {
    mat_1_abs <- abs(mat_1)
    mat_2_abs <- abs(mat_2)
    mat_1 <- ifelse(mat_1_abs <= mat_2_abs, mat_2_abs, mat_1_abs) ## max
    mat_2 <- ifelse(mat_1_abs > mat_2_abs, mat_2_abs, mat_1_abs) ## min
  }
  
  ## create three matrices to store result (additional to MetNet: mat_mass)
  mat <- matrix(0, nrow = length(mass), ncol = length(mass))
  mat_type <- matrix("", nrow = length(mass), ncol = length(mass))
  mat_mass <- matrix("", nrow = length(mass), ncol = length(mass))
  
  ## iterate through each column and check if the "mass" is in the interval
  ## defined by the m/z value and ppm
  for (i in seq_along(transformation[, "mass"])) {
    
    transformation_i <- transformation[i, ]
    ind_mat_1 <- which(mat_1 >= transformation_i[["mass"]])
    ind_mat_2 <- which(mat_2 <= transformation_i[["mass"]])
    
    ## get intersect from the two (indices where "mass" is in the interval)
    ind_hit <- intersect(ind_mat_1, ind_mat_2)
    
    ## write to these indices 1 and the "group"
    mat[ind_hit] <- 1
    mat_type[ind_hit] <- ifelse(nchar(mat_type[ind_hit]) != 0,
                                yes = paste(mat_type[ind_hit], transformation_i[["group"]],
                                            sep = "/"),
                                no = as.character(transformation_i[["group"]]))
    ## additional to MetNet:
    ## wirte to these indices 1 and the "mass"
    mat_mass[ind_hit] <- ifelse(nchar(mat_mass[ind_hit]) != 0,
                                yes = paste(mat_mass[ind_hit], transformation_i[["mass"]],
                                            sep = "/"),
                                no = as.numeric(transformation_i[["mass"]]))
    
  }
  
  rownames(mat) <- colnames(mat) <- rownames(x)
  rownames(mat_type) <- colnames(mat_type) <- rownames(x)
  rownames(mat_mass) <- colnames(mat_mass) <- rownames(x) #additional to MetNet
  
  return(list(mat, mat_type, mat_mass))
  
}

#'Changes to MetNet:
#' additional attribute in funtion: p; default is FALSE (so works like MetNet), if p is TRUE and model is spearman/pearson
#' than the output will contain lists of pearson/spearman containing the corresponding correlation values (INCLUDING positive 
#' and negative values) and p-values
statistical <- function(x, model, p = FALSE, ...) {
  
  ## check if model complies with the implemented model and return error
  ## if not so
  if (!(all(model %in% c("lasso", "randomForest", "clr", "aracne",
                         "pearson", "pearson_partial", "pearson_semipartial",
                         "spearman", "spearman_partial", "spearman_semipartial", "bayes"))))
    stop("'model' not implemented in statistical")
  
  ## check if x is numeric matrix and return error if not so
  if (mode(x) != "numeric") stop("x is not a numerical matrix")
  
  ## z-scale x and transpose
  x_z <- apply(x, 1, function(y) {
    (y - mean(y, na.rm = TRUE)) / sd(y, na.rm = TRUE)
  })
  x_z <- t(x_z)
  
  l <- list()
  
  ## add entry for lasso if "lasso" is in model
  if ("lasso" %in% model) {
    lasso <- lasso(x = x_z, ...)
    diag(lasso) <- NaN
    l <- addToList(l, "lasso", lasso)
    print("lasso finished")
  }
  
  ## add entry for randomForest if "randomForest" is in model
  if ("randomForest" %in% model) {
    randomForest <- randomForest(x = x, ...)
    diag(randomForest) <- NaN
    l <- addToList(l, "randomForest", randomForest)
    print("randomForest finished.")
  }
  
  ## calculate mutual information if "clr" or "aracne" is in model
  if (any(c("clr", "aracne") %in% model)) {
    mi_x_z <- mpmi::cmi(t(x_z))$bcmi
    rownames(mi_x_z) <- colnames(mi_x_z) <- rownames(x)
  }
  
  ## add entry for clr if "clr" is in model
  if ("clr" %in% model) {
    clr <- threeDotsCall("clr", mi = mi_x_z, ...)
    diag(clr) <- NaN
    l <- addToList(l, "clr", clr)
    print("clr finished.")
  }
  
  ## add entry for aracne if "aracne" is in model
  if ("aracne" %in% model) {
    aracne <- threeDotsCall("aracne", mi = mi_x_z, ...)
    diag(aracne) <- NaN
    l <- addToList(l, "aracne", aracne)
    print("aracne finished.")
  }
  
  ## add entry for pearson if "pearson" is in mode
  ## FALSE %in% p is default and corresponds to original MetNet function
  if ("pearson" %in% model & FALSE %in% p) {
    pearson <- threeDotsCall("correlation", x = x, type = "pearson", ...)
    diag(pearson) <- NaN
    l <- addToList(l, "pearson", pearson)
    print("pearson finished.")
  }
  
  
  ## add entry for pearson_partial if "pearson_partial" is in model
  if ("pearson_partial" %in% model) {
    pearson_partial <- threeDotsCall("correlation", x = x,
                                     type = "pearson_partial", ...)
    diag(pearson_partial) <- NaN
    l <- addToList(l, "pearson_partial", pearson_partial)
    print("pearson_partial finished.")
  }
  
  ## add entry for pearson_semipartial if "pearson_semipartial" is in model
  if ("pearson_semipartial" %in% model) {
    pearson_sp <- threeDotsCall("correlation", x = x,
                                type = "pearson_semipartial", ...)
    diag(pearson_sp) <- NaN
    l <- addToList(l, "pearson_semipartial", pearson_sp)
    print("pearson_semipartial finished.")
  }
  
  ## add entry for spearman if "spearman" is in model
  ## FALSE %in% p is default and corresponds to original MetNet function
  if ("spearman" %in% model & FALSE %in% p) {
    spearman <- threeDotsCall("correlation", x = x, type = "spearman", ...)
    diag(spearman) <- NaN
    l <- addToList(l, "spearman", spearman)
    print("spearman finished.")
  }
  
  ## add entry for spearman_partial if "spearman_partial" is in model
  if ("spearman_partial" %in% model) {
    spearman_partial <- threeDotsCall("correlation", x = x,
                                      type = "spearman_partial", ...)
    diag(spearman_partial) <- NaN
    l <- addToList(l, "spearman_partial", spearman_partial)
    print("spearman_partial finished.")
  }
  
  ## add entry for spearman_semipartial if "spearman_semipartial" is in model
  if ("spearman_semipartial" %in% model) {
    spearman_sp <- threeDotsCall("correlation", x = x,
                                 type = "spearman_semipartial", ...)
    diag(spearman_sp) <- NaN
    l <- addToList(l, "spearman_semipartial", spearman_sp)
    print("spearman_semipartial finished.")
  }
  
  ## add entry for pearson if "pearson" is in model and p=TRUE
  ## Changed to MetNet: if p=TRUE, negative and positive correlation values were calculated and corresponding p-values
  if ("pearson" %in% model & TRUE %in% p) {
    pearson <- threeDotsCall("correlation_p", x = x, type = "pearson", ...)
    #diag(pearson) <- NaN
    diag(pearson[[1]]) <- NaN
    diag(pearson[[2]]) <- NaN
    l <- addToList_p(l, "pearson", pearson)
    print("pearson finished.")
  }
  
  ## add entry for pearson if "pearson" is in model and P=TRUE
  ## Changed to MetNet: if p=TRUE, negative and positive correlation values were calculated and corresponding p-values
  if ("spearman" %in% model & TRUE %in% p) {
    spearman <- threeDotsCall("correlation_p", x = x, type = "spearman", ...)
    #diag(spearman) <- NaN
    diag(spearman[[1]]) <- NaN
    diag(spearman[[2]]) <- NaN
    l <- addToList_p(l, "spearman", spearman)
    print("spearman finished.")
  }
  
  
  ## add entry for bayes if "bayes" is in model
  if ("bayes" %in% model) {
    bayes <- threeDotsCall("bayes", x = x, ...)
    diag(bayes) <- NaN
    l <- addToList(l, "bayes", bayes)
    print("bayes finished.")
  }
  
  return(l)
}
#' correlation_p is an additional function (complementary to correlation()), and needed for positive/negative pearson/spearman
#' correlation value and p-value calculation
correlation_p <- function(x, type = "pearson", use = "pairwise.complete.obs") {
  
  ## for pearson/spearman correlation
  if (type %in% c("pearson", "spearman")) {
    cor_list <- rcorr(x = t(x), type = type)
  }
  
  
  names(cor_list) <- c("Correlation Value", "n", "p-Value")
  ## exclude "n" column
  cor_list <- cor_list[-2] 
  
  return(cor_list)
}
addToList <- function(l, name, object) {
  
  ## test validity of objects
  if (!is.list(l)) {
    stop("l is not a list")
  }
  
  if (!is.character(name)) {
    stop("name is not a character")
  }
  
  if (!is.matrix(object)) {
    stop("object is not a matrix")
  }
  
  ## add object to l
  new_index <- length(l) + 1
  l[[new_index]] <- object
  
  ## assign the name to the newly added entry
  names(l)[new_index] <- name
  
  return(l)
}#nothing changed
#'Not tested if object is a matrix (since object is a list) --> maybe same function can be used?
addToList_p <- function(l, name, object) {
  
  ## test validity of objects
  if (!is.list(l)) {
    stop("l is not a list")
  }
  
  if (!is.character(name)) {
    stop("name is not a character")
  }
  
  ## add object to l
  new_index <- length(l) + 1
  l[[new_index]] <- object
  
  ## assign the name to the newly added entry
  names(l)[new_index] <- name
  
  return(l)
}
threeDotsCall <- function(fun, ...) {
  
  formal_args <- formalArgs(fun)
  args <- list(...)
  if (any(duplicated(names(args)))) stop("duplicated args in ...")
  
  input <- args[names(args) %in% formal_args]
  
  ## call the function
  res <- do.call(fun, input)
  return(res)
}#nothing changed

#'Changes to MetNet: new attribute for type: "threshold_p"
#'A list is created instead of a single matrix as output containing 1/0 assigned values of 
#'model matrices (e.g. pearson and spearman) and consensus matrix.
#'
#'If "treshold_p" is selected in 'type' all values are assigned to 1 if their p-value 
#'is BELOW a defined threshold (defined in'args')
#'If "treshold" is selected in 'type' all values are assigned to 1 if their Correlation-value 
#'is ABOVE a defined threshold (defined in'args')
threshold <- function(statistical, type, args,
                      values = c("all", "min", "max"), ...) {
  
  l <- statistical
  ## args, either N for tops
  ## or a list of threshold
  if (any(duplicated(names(args)))) {
    stop("names(args) contain duplicated entries")
  }
  
  ##Changes to MetNet: new attribute for type: "threshold_p"
  if (!type %in% c("top1", "top2", "mean", "threshold", "threshold_p"))
    stop("type not in 'top1', 'top2', 'mean', 'threshold', 'threshold_p'")
  
  ## check args
  if (type %in% c("threshold")) {
    if (!(all(names(l) %in% names(args)))) {
      stop("'args' does not contain entries for all 'model's in ",
           "'statistical'")
    }
    
    if (!"threshold" %in% names(args) && length(args$threshold) != 1) {
      stop("'args' does not contain entry 'threshold' of length 1")
    }
  }
  
  ## check args
  if (type %in% c("threshold")) {
    if (!(all(names(l) %in% names(args)))) {
      stop("'args' does not contain entries for all 'model's in ",
           "'statistical'")
    }
    
    if (!"threshold" %in% names(args) && length(args$threshold) != 1) {
      stop("'args' does not contain entry 'threshold' of length 1")
    }
  }
  ## complementary to "threshold":
  if (type %in% c("threshold_p")) {
    if (!(all(names(l) %in% names(args)))) {
      stop("'args' does not contain entries for all 'model's in ",
           "'statistical'")
    }
    
    if (!"threshold_p" %in% names(args) && length(args$threshold) != 1) {
      stop("'args' does not contain entry 'threshold' of length 1")
    }
  }
  
  ## check match.arg for values
  values <- match.arg(values)
  
  if (type %in% c("top1", "top2", "mean")) {
    if (!("n"  %in% names(args) && length(args$n) == 1 &&
          is.numeric(args$n)))
      stop("args does not contain the numeric entry `n` of length 1")
  }
  
  
  
  
  if (type == "threshold" || type == "threshold_p") {
    ## iterate through the list and remove the links below or above the
    ## threshold and write to list
    l <- lapply(seq_along(l), function(x) {
      
      ## find corresponding model in l
      name_x <- names(l)[x]
      
      ## get corresponding threshold in args
      threshold_x <- args[[names(l)[x]]]
      
      ## Changed to MetNet
      if ("threshold" %in% type) {
        
        if("Correlation Value" %in% names(l[[name_x]][1])) {
          ## get corresponding adjacency matrix of Correlation Values in l
          ## is used when statistical was calculated with p=TRUE
          l_x <- l[[name_x]]$`Correlation Value`

          ## only assign 1 to values that are above the threshold
          ifelse(l_x > threshold_x, 1, 0)
        } 
        else{
          ## get corresponding adjacency matrix in l
          ## corresponds to MetNet function
          l_x <- l[[name_x]]
          ## for pearson/spearman correlation models (incl. partial and
          ## semi-partial), lasso, randomForest, clr, aracne and bayes higher
          ## values corresond to higher confidence
          ## only assign 1 to values that are above the threshold
          ifelse(l_x > threshold_x, 1, 0)
        }
      }
      else if("threshold_p" %in% type){
        ## get corresponding adjacency matrix of p-Values in l
        ## is used when statistical was calculated with p=TRUE
        l_x <- l[[name_x]]$`p-Value`
        ## only assign 1 to values that are below the threshold
        ifelse(l_x < threshold_x, 1, 0)
      }      
      
    })
    
    ## allow for compatibility of arguments
    ## calculate consenses from the binary matrices
    cons <- threeDotsCall(sna::consensus, dat = l, ...)
    
    ## threshold consensus that it is a binary matrix
    cons <- ifelse(cons >= args$threshold, 1, 0)
    
    rownames(cons) <- colnames(cons) <- colnames(l[[1]])
  } 
  else { ## if type is in "top1", "top2" or "mean"
    l_df <- lapply(seq_along(l), function(x) {
      
      ## find corresponding model in l
      name_x <- names(l)[x]
      
      ## get corresponding adjacency matrix in l
      l_x <- l[[name_x]]
      
      ## take the respective minimum or maximum depending on `values`,
      ## do not do anything if `values` is equal to `all`
      if (values %in% c("min", "max")) {
        
        ## get values from the lower triangle
        lower_tri <- l_x[lower.tri(l_x)]
        
        ## get values from the upper triangle (requires transposing)
        l_x_t <- t(l_x)
        upper_tri <- l_x_t[lower.tri(l_x_t)]
        
        ## get min of lower_tri and upper_tri
        if (values == "min") {
          values <- apply(rbind(lower_tri, upper_tri), 2, min)
        } else {
          values <- apply(rbind(lower_tri, upper_tri), 2, max)
        }
        
        ## write back to the matrix
        l_x[lower.tri(l_x)] <- values
        l_x <- t(l_x)
        l_x[lower.tri(l_x)] <- values
      }
      
      ## for pearson/spearman correlation (incl. partial and
      ## semi-partial), lasso, randomForest, clr, aracne and bayes
      ## higher values corresond to higher confidence
      if (grepl(name_x, pattern = "lasso|randomForest|bayes")) {
        ## set values that are equal to 0 to NaN (values that are 0)
        ## do not explain the variability
        res <- getLinks(l_x, exclude = "== 0")
      }
      if (grepl(name_x, pattern = "pearson|spearman|clr|aracne")) {
        res <- getLinks(l_x, exclude = NULL)
      }
      
      res
    })
    
    names(l_df) <- names(l)
    
    ## bind together the ranks of the models, stored in l_df
    ranks <- lapply(l_df, function(x) x$rank)
    ranks <- do.call("cbind", ranks)
    colnames(ranks) <- names(l_df)
    
    ## calculate the consensus information, i.e. either get the first or
    ## second top rank per row or calculate the average across rows
    ## depending on the type argument
    cons_val <- MetNet:::topKnet(ranks, type)
    
    ## bind row and col information with cons information
    row_col <- l_df[[1]][, c("row", "col")]
    ranks <- cbind(row_col, cons_val)
    
    ## get the top N features
    n <- args$n
    top_n <- sort(unique(cons_val))[1:n]
    ranks_top <- ranks[cons_val %in% top_n, ]
    
    ## write links in ranks_top to binary adjacency matrix cons
    cons <- matrix(0, nrow = ncol(l[[1]]), ncol = ncol(l[[1]]))
    rownames(cons) <- colnames(cons) <- colnames(l[[1]])
    cons[as.numeric(rownames(ranks_top))] <- 1
    
  }
  
  ## Changes to MetNet: A list is created as output containing 1/0 assigned values of 
  ## model matrices (e.g. pearson and spearman) and consensus matrix
  names(l) <- names(statistical)
  l[["Consensus"]] <- cons
  class(l[[3]]) <- "numeric"
  return(l)
}

#' Changes to MetNet: New attributes added, if model = "combined" (default) the result will be the same 
#' as in MetNet, except the output list items were named to "combined" and "Character" 
#' if model = "pearson" or "spearman" than also the corresponding weighted statistical 
#' adjacency matrix is required as attribute (weighted_statistical = XY)
#' The output in this case will be a list containing 4 listitems, where combination relies on the
#' unweighted adjacency matrix of either Pearson or Spearman. 
#' Moreover corresponding Correlation and p-values will be displayes as listitems
combine <- function(structural, statistical, threshold = 1, model = "combined", weighted_statistical) {
  
  ## Is changed since structural list is now lenght 3
  if (!is.list(structural) | length(structural) != 3)
    stop("structural is not a list of length 3")
  
  if (!is.matrix(structural[[1]]) | !is.numeric(structural[[1]]))
    stop("structural[[1]] is not a numeric matrix")
  
  if (!is.matrix(structural[[2]]) | !is.character(structural[[2]]))
    stop("structural[[2]] is not a character matrix")
  
  ## Additional to MetNet:
  if (!is.matrix(statistical[[3]]) | !is.numeric(statistical[[3]]))
    stop("statistical is not a numeric matrix")
  
  if (!all(rownames(structural[[1]]) == rownames(structural[[2]])))
    stop(c("rownames of structural[[1]] are not identical to rownames of ",
           "structural[[2]]"))
  
  if (!all(colnames(structural[[1]]) == colnames(structural[[2]])))
    stop(c("colnames of structural[[1]] are not identical to colnames of ",
           "structural[[2]]"))
  
  if (!all(rownames(structural[[1]]) == rownames(statistical)))
    stop("rownames are not identical")
  
  if (!all(colnames(structural[[1]]) == colnames(statistical)))
    stop("colnames are not identical")
  
  if (!is.numeric(threshold)) stop("threshold is not numeric")
  
  ## create list to store results
  res <- list()
  
  
  ## Changes to MetNet: Distinguish between default (model == "combined") and model = "pearson" or "spearman"
  if(model == "pearson"){
    
    ## create the first entry of the list
    ## sum the matrices structural and statistical, if the value is above
    ## threshold then assign 1, otherwise 0
    cons_num <- structural[[1]] + statistical[["pearson"]]
    cons_num <- ifelse(cons_num > threshold, 1, 0)
    
    ## if p-values have previously been calculated
    if("Correlation Value" %in% names(weighted_statistical[["pearson"]][1])){
    cons_corr <- ifelse(cons_num == 1, weighted_statistical[["pearson"]][["Correlation Value"]], "")
    cons_p <- ifelse(cons_num == 1, weighted_statistical[["pearson"]][["p-Value"]], "")}
    else {
    cons_corr <- ifelse(cons_num == 1, weighted_statistical[["pearson"]], "")
    cons_p <- NaN
     }}
    
  if(model == "spearman"){
    ## create the first entry of the list
    ## sum the matrices structural and statistical, if the value is above
    ## threshold then assign 1, otherwise 0
    cons_num <- structural[[1]] + statistical[["spearman"]]
    cons_num <- ifelse(cons_num > threshold, 1, 0)
    
    ## if p-values have previously been calculated
    if("Correlation Value" %in% names(weighted_statistical[["spearman"]][1])){
      cons_corr <- ifelse(cons_num == 1, weighted_statistical[["spearman"]][["Correlation Value"]], "")
      cons_p <- ifelse(cons_num == 1, weighted_statistical[["spearman"]][["p-Value"]], "")}
    else {
      cons_corr <- ifelse(cons_num == 1, weighted_statistical[["spearman"]], "")
      cons_p <- NaN
    }}
    
  
  if(model == "combined"){
    ## create the first entry of the list
    ## sum the matrices structural and statistical, if the value is above
    ## threshold then assign 1, otherwise 0
    cons_num <- structural[[1]] + statistical[[3]]
    cons_num <- ifelse(cons_num > threshold, 1, 0)
    cons_corr <- NaN
    cons_p <- NaN
  }
  
  ## create the second entry of the list
  ## if element in cons_num is equal to 1, take the element in structural[[2]]
  ## (the type of link), otherwise ""
  cons_char <- ifelse(cons_num == 1, structural[[2]], "")
  
  ## assign to list
  ## Compared to MetNet names were assigned
  res[[model]] <- cons_num
  res[["Character"]] <- cons_char
  
  ## assign Correlation and p-values to list if model is "pearson" or "spearman"
  if(model == "pearson" || model == "spearman"){
    res[["Correlation Value"]] <- cons_corr
    res[["p-Value"]] <- cons_p
  }

  
  return(res)
}



exportNet2gml <- function (x, from, ...) {
  if ("structural" %in% from) {
    
    mat <- x[[1]]
    mat_type <- x[[2]]
    mat_mass <- x[[3]]
    class(mat_mass) <- "numeric"
    net       <- 
      graph_from_adjacency_matrix(mat, mode = "undirected", weighted = T)
    net_type  <-
      graph_from_adjacency_matrix(mat_type, mode = "undirected", weighted = T)
    net_mass  <-
      graph_from_adjacency_matrix(mat_mass, mode = "undirected", weighted = T)
    net_comb  <- union(net, net_mass)
    names(edge_attr(net_comb))[1] <- "adj"
    names(edge_attr(net_comb))[2] <- "mass difference"
    
    #net_plot <- plot(net_type, edge.width = 5, vertex.label.cex = 0.5, edge.color = "grey")
    
    write_graph(net_comb, "structural_type.gml", format = c("gml"))
  }
  else if ("statistical+p" %in% from) {
    
    for (i in 1:length(x)) {
      cor_list <- x[[i]]
      ##Plot structural adjacency matrix and export to gml
      net_cor <-
        igraph::graph_from_adjacency_matrix(cor_list[[1]], mode = "undirected", weighted = T)
      net_p   <-
        igraph::graph_from_adjacency_matrix(cor_list[[2]], mode = "undirected", weighted = T)
      net_comb <- union(net_cor, net_p)
      names(edge_attr(net_comb))[1] <- "correlation"
      names(edge_attr(net_comb))[2] <- "p"
      # #net_plot <- plot(net_type, edge.width = 5, vertex.label.cex = 0.5, edge.color = "grey")
      q <- names(x[i])
      write_graph(net_comb, file = sprintf('statistical.%s.gml', q), format = c("gml"))
      
    }
    
  }
  else if ("combine" %in% from) {
    if ("pearson" %in% names(x[1]) | "spearman" %in% names(x[1]) ){
      class(x[[3]]) <- "numeric"
      class(x[[4]]) <- "numeric"
      net_cor <- graph_from_adjacency_matrix(x[[3]], mode = "undirected", weighted = T)
      net_p   <- graph_from_adjacency_matrix(x[[4]], mode = "undirected", weighted = T)
      net     <- union(net_cor, net_p)
      names(edge_attr(net))[1] <- "correlation"
      names(edge_attr(net))[2] <- "p"
    }
    else { #if "combined" or other model
      net       <- 
        graph_from_adjacency_matrix(x[[1]], mode = "undirected", weighted = T)
    }
    write_graph(net, "combined.gml", format = c("gml"))
    
  }
}

adjacency_list <- function(x, from){
  
  if (!(all(from %in% c("structural", "statistical", "threshold", "combine"))))
    stop("'from' not implemented in adjacency_list")
  
  if ("structural" %in% from) {
    
    
    x[[2]][upper.tri(x[[2]])] <- ''
    x[[3]][upper.tri(x[[3]])] <- ''
    
    list_type <- melt(x[[2]]) %>% filter(Var1 != Var2) %>% filter(value != '')
    list_mass <- melt(x[[3]]) %>% filter(Var1 != Var2) %>% filter(value != '')
    combine <- add_column(list_type,  `mass difference`= list_mass$value) %>% as.data.frame()
    return(combine)
  }
  else if ("statistical" %in% from) {
    
    for (i in seq_along(x)) {
      if (i == 1) {
        x[[i]][upper.tri(x[[i]])] <- ''
        list_corr <- melt(x[[i]]) %>% filter(Var1 != Var2) %>% filter(value != '') %>% 
          select(Var1, Var2, value) 
        colnames(list_corr) <- c("Feature1", "Feature2", names(x[i]))
        #return(list_corr)
      }
      if (i != 1){
        model = names(x[i])
        x[[i]][upper.tri(x[[i]])] <- ''
        list_corr2 <- melt(x[[i]]) %>% filter(Var1 != Var2) %>% filter(value != '')
        list_comb <- add_column(list_corr, list_corr2$value)
        list_comb <- as.data.frame(list_comb) 
        colnames(list_comb)[i+2] <- c(names(x[i]))
        list_corr <- list_comb
      } 
    } 
    return(list_corr)
  }
  else if ("combine" %in% from){
    x[[2]][upper.tri(x[[2]])] <- ''
    x[[3]][upper.tri(x[[3]])] <- ''
    x[[4]][upper.tri(x[[4]])] <- ''
    
    list_mass <- melt(x[[2]]) %>% filter(Var1 != Var2) %>% filter(value != '')
    list_corr <- melt(x[[3]]) %>% filter(Var1 != Var2) %>% filter(value != '')
    list_p    <- melt(x[[4]]) %>% filter(Var1 != Var2) %>% filter(value != '')
    listed <- add_column(list_mass, `Correlation Value` = list_corr$value)
    listed <- add_column(listed, `p-Value` = list_p$value)
return(listed)
    
  }
  
}

sum_mass <- function(adjacency_list){
  
  if("mass difference" %in% names(adjacency_list)){
  sum_mass <- adjacency_list %>% group_by(`mass difference`) %>% summarise(count=n()) %>%
    as.data.frame()
  sum_comb <- adjacency_list %>% group_by(`value`) %>% summarise(count=n()) %>%
    as.data.frame()  %>% add_column(sum_mass$`mass difference`)
  colnames(sum_comb) <- c("Type", "Counts", "Mass Difference")
  sum_comb <- sum_comb %>% select(Type, `Mass Difference`, Counts)}
  else{
    sum_comb <- adjacency_list %>% group_by(`value`) %>% summarise(count=n()) %>%
      as.data.frame()
    colnames(sum_comb) <- c("Type", "Counts")
    
  }
  
  plot_list <- ggplot(sum_comb, aes(x=Type, y=Counts, fill=Type)) + geom_bar(stat = "identity") + theme_minimal() + 
    labs(title = "Numbers of destinct type of a biochemical reaction")+ scale_fill_brewer(palette = "Blues") + theme(legend.position = "right")
  plot(plot_list)
  
  return(sum_comb)
}
