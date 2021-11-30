if (TRUE){
  
  debug <- 3
  
  true_bn <- bnrepository("asia")
  
  set.seed(1)
  x <- bnlearn::rbn(true_bn, n = 1e4)
  x <- as.data.frame(sapply(x, function(x) as.factor(as.integer(x) - 1L)),
                     stringsAsFactors = TRUE)
  
  eg <- expand.grid(restrict = c("ppc", "true", "cig"),
                    maximize = c(""),
                    test = bnlearn:::available.discrete.tests,
                    score = "bic",
                    alpha = c(1e-2),
                    undirected = c(FALSE),
                    sort_pval = c(TRUE),
                    path = c(10),
                    hgi = c(TRUE),
                    stringsAsFactors = FALSE)
  
  results <- vector("list", nrow(eg))
  
  for (i in seq_len(nrow(eg))){
    
    if (!is.null(results[[i]])) next
    
    list2env(eg[i,], envir = .GlobalEnv)
    if (restrict == "" && maximize == "") next
    
    cat(sprintf("%g) restrict = %s, maximize = %s, test = %s, score = %s, alpha = %g,\n    undirected = %s, sort_pval = %s, path = %g, hgi = %s\n",
                i, restrict, maximize, test, score, alpha, undirected, sort_pval, path, hgi))
    
    if (restrict %in% c("true", "cig"))
      temp_bn <- true_bn
    else
      temp_bn <- NULL
    
    result <- bnsl(x = x,
                   restrict = restrict, maximize = maximize,
                   restrict.args = list(test = test, alpha = alpha, 
                                        max.sx = 3, sort_pval = sort_pval),
                   maximize.args = list(score = score, maxp = 8),
                   undirected = undirected,
                   path = path, hgi = hgi, 
                   true_bn = temp_bn, debug = debug)
    
    hide <- capture.output(print(result))
    
    results[[i]] <- result
    
    ## TODO: write actual tests
  }
  all(sapply(results, function(x) if (!is.null(x)) sum(bnlearn::amat(x)) else 1) > 0)
}