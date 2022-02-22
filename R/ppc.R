# Partitioned PC (pPC) skeleton learning
# 
# Skeleton learning phase of the pPC algorithm.

ppc_skeleton <- function (x, cluster = NULL, whitelist, blacklist, test, alpha, B, complete,
                          max.sx = ncol(x)-2, max_wthn_sx = max.sx, max_btwn_sx = max.sx,
                          max_btwn_nbr = ncol(x)-2, sort_pval = TRUE, max_groups = 20, 
                          true_bn = NULL, debug = FALSE)
{
  times <- tests <- 
    c(partition = 0, within = 0, between = 0, complete = 0)
  tests[] <- 1 * bnlearn::test.counter()
  
  nodes <- names(x)
  nnodes <- length(nodes)
  mb <- structure(vector(length(nodes), mode = "list"), names = nodes)
  skeleton <- bnlearn:::subsets(nodes, 2)
  node.pairs <- apply(skeleton, 1, function(x) list(arc = x, 
                                                    max.adjacent = nnodes - 1))
  nbr.size = sapply(node.pairs, `[[`, "max.adjacent")
  max.dsep.size <- min(max.sx, length(nodes) - 2)
  max_wthn_sx <- min(max.dsep.size, max_wthn_sx)
  max_btwn_sx <- min(max.dsep.size, max_btwn_sx)
  max_btwn_nbr <- min(max_btwn_nbr, length(nodes) - 2)
  bnlearn:::check.logical(sort_pval)
  debug_cli(!is.numeric(max_groups) || (max_groups %% 1 != 0) || 
              is.infinite(max_groups), cli::cli_abort, 
            "max_groups must be a positive numeric integer")
  nlev <- sapply(x, function(y) length(unique(y)))
  key <- build_key(nodes = nodes)
  
  
  ## partition
  start_time_p <- Sys.time()
  test0 <- ifelse(max_groups <= 1, test,
                  ifelse(test %in% bnlearn:::available.discrete.tests,
                         "mi", "cor"))
  node.pairs <- bnlearn:::smartSapply(cl = cluster, node.pairs, ppc_heuristic,
                                      data = x, alpha = alpha, B = B, 
                                      whitelist = whitelist, blacklist = blacklist, 
                                      test = test0, dsep.size = 0, complete = complete, 
                                      true_bn = true_bn, debug = debug)
  
  if (max_groups > 1){
    distances <- get_distances(nodes = nodes, nnodes = nnodes, key = key,
                               skeleton = skeleton, node.pairs = node.pairs,
                               x = x, test = test, B = B,
                               alpha = alpha, complete = complete)
    group <- get_group(distances, nnodes = nnodes, max_groups = max_groups)
  } 
  else {
    group <- sapply(nodes, function(x) 1)
  }
  n_group <- length(unique(group))
  bool_btwn <- get_bool_btwn(group)
  
  ## attempt at memory efficiency by only keeping necessary information
  ## TODO: verify whether or not helpful
  node.pairs <- lapply(node.pairs, function(x) x[c("arc", "p.value", 
                                                   "dsep.set")])
  
  parentList <- node.pairs2parentList(node.pairs = node.pairs[!bool_btwn], 
                                      nodes = nodes, alpha = alpha)
  bool_pairs <- get_bool_pairs(min_size = 1, node.pairs, alpha, parentList)
  
  end_time_p <- Sys.time()
  times[1] <- as.numeric(end_time_p - start_time_p, unit = "secs")
  tests[1] <- bnlearn::test.counter() - tests[1]
  
  debug_cli(debug, cli::cli_alert,
            c("divided {ncol(x)} variables into {n_group} groups ",
              "with sizes {paste(table(group), collapse = ', ')} ",
              "in {prettyunits::pretty_sec(times[1])} with {tests[1]} calls"))
  
  
  ## within cluster skeleton estimation
  for (dsep.size in seq_len(max_wthn_sx)){
    
    ## consistently update bool_pairs to minimize the number of
    ## elements in node.pairs that we pass over
    bool_pairs <- bool_pairs & !bool_btwn
    
    node.pairs[bool_pairs] <- 
      bnlearn:::smartSapply(cluster, 
                            node.pairs[bool_pairs], ppc_heuristic, 
                            data = x, alpha = alpha, B = B, whitelist = whitelist, 
                            blacklist = blacklist, test = test, skeleton = NULL, 
                            dsep.size = dsep.size, complete = complete, 
                            sort_pval = sort_pval, group = NULL, 
                            parentList = parentList, node.pairs = node.pairs, 
                            key = key, true_bn = true_bn, debug = debug)
    
    parentList <- node.pairs2parentList(node.pairs = node.pairs[!bool_btwn], 
                                        nodes = nodes, alpha = alpha)
    bool_pairs <- get_bool_pairs(min_size = dsep.size+1, node.pairs, alpha, parentList)
    if (! any(bool_pairs))
      break
  }
  end_time_w <- Sys.time()
  times[2] <- as.numeric(end_time_w - end_time_p, unit = "secs")
  tests[2] <- bnlearn::test.counter() - sum(tests[1:2])
  
  debug_cli(debug, cli::cli_alert,
            c("estimated edges within {n_group} clusters ",
              "in {prettyunits::pretty_sec(times[2])} with {tests[2]} calls"))
  
  
  ## screen edges between clusters, if multiple clusters
  if (n_group > 1){
    
    ## first screen, of between cluster node pairs
    bool_pairs <- sapply(node.pairs, function(x) x$p.value <= alpha && 
                           any(sapply(parentList[x$arc], length) > 0))
    bool_pairs <- bool_pairs & bool_btwn
    node.pairs[bool_pairs] <- 
      bnlearn:::smartSapply(cluster, 
                            node.pairs[bool_pairs], btwn_heuristic, 
                            data = x, alpha = alpha, B = B, whitelist = whitelist, 
                            blacklist = blacklist, test = test, skeleton = NULL, 
                            dsep.size = max_btwn_sx, complete = complete, 
                            max_btwn_nbr = max_btwn_nbr, group = NULL, 
                            parentList = parentList, node.pairs = node.pairs, 
                            key = key, round = 1, true_bn = true_bn, debug = debug)
    
    parentList <- node.pairs2parentList(node.pairs = node.pairs,
                                        nodes = nodes, alpha = alpha)
    bool_pairs <- bool_btwn & get_bool_pairs(min_size = 1, 
                                             node.pairs, alpha, parentList)
    
    ## exclude node pairs that are not connected to 
    ## other between cluster node pairs, meaning they
    ## don't have updated neighbors to investigate
    if (any(bool_pairs)){
      btwn_skel <- t(sapply(node.pairs[bool_pairs], `[[`, "arc"))
      tab_btwn_skel <- table(btwn_skel)
      bool_connected <- apply(btwn_skel, 1, function(x) any(tab_btwn_skel[x] > 1))
      bool_pairs[bool_pairs] <- bool_pairs[bool_pairs] & bool_connected
    }
    
    
    ## second screen, now of between cluster edges
    if (any(bool_pairs)){
      node.pairs[bool_pairs] <- 
        bnlearn:::smartSapply(cluster, 
                              node.pairs[bool_pairs], btwn_heuristic, 
                              data = x, alpha = alpha, B = B, whitelist = whitelist, 
                              blacklist = blacklist, test = test, skeleton = NULL, 
                              dsep.size = max_btwn_sx, complete = complete, 
                              max_btwn_nbr = max_btwn_nbr, sort_pval = sort_pval, 
                              group = group, parentList = parentList, 
                              node.pairs = node.pairs, key = key, round = 2, 
                              true_bn = true_bn, debug = debug)
      
      parentList <- node.pairs2parentList(node.pairs = node.pairs,
                                          nodes = nodes, alpha = alpha)
    }
    end_time_b <- Sys.time()
    times[3] <- as.numeric(end_time_b - end_time_w, unit = "secs")
    tests[3] <- bnlearn::test.counter() - sum(tests[1:3])
    
    debug_cli(debug, cli::cli_alert,
              c("screened edges between {n_group} clusters ",
                "in {prettyunits::pretty_sec(times[3])} with {tests[3]} calls"))
    
    
    bool_pairs <- get_bool_pairs(min_size = 1, node.pairs, alpha, parentList)
    btwn_skel <- t(sapply(node.pairs[bool_pairs & bool_btwn], `[[`, "arc"))
    bool_pairs <- bool_pairs & get_affected(key[unique(as.character(btwn_skel))], nnodes)
    
    ## complete remaining tests
    for (dsep.size in seq(from = 1, to = max.dsep.size)){
      node.pairs[bool_pairs] <- 
        bnlearn:::smartSapply(cluster,
                              node.pairs[bool_pairs], ppc_heuristic,
                              data = x, alpha = alpha, B = B, whitelist = whitelist,
                              blacklist = blacklist, test = test, skeleton = NULL,
                              dsep.size = dsep.size, complete = complete, 
                              max_wthn_sx = max_wthn_sx, sort_pval = sort_pval, 
                              group = group, parentList = parentList, 
                              node.pairs = node.pairs, key = key, true_bn = true_bn, 
                              debug = debug)
      
      parentList <- node.pairs2parentList(node.pairs = node.pairs,
                                          nodes = nodes, alpha = alpha)
      bool_pairs <- bool_pairs & get_bool_pairs(min_size = 1, 
                                                node.pairs, alpha, parentList)
      if (! any(bool_pairs)){
        break
      }
    }
    end_time_c <- Sys.time()
    times[4] <- as.numeric(end_time_c - end_time_b, unit = 'secs')
    tests[4] <- bnlearn::test.counter() - sum(tests[1:4])
    
    debug_cli(debug, cli::cli_alert,
              c("completed conditional independence investigation ",
                "in {prettyunits::pretty_sec(times[4])} with {tests[4]} calls"))
  } else{
    
    tests[3:4] <- 0
  }
  skeleton = do.call(rbind, lapply(node.pairs, function(x) {
    if (x$p.value <= alpha) 
      return(x$arc)
    else return(NULL)
  }))
  skeleton <- bnlearn:::cache.structure(nodes, bnlearn:::arcs.rbind(skeleton, skeleton,
                                                                    reverse2 = TRUE))
  attr(skeleton, "dsep.set") <- node.pairs
  attr(skeleton, "learning") <- list(group = group,
                                     tests = tests,
                                     times = times)
  invisible(skeleton)
}



# Version of bnlearn:::pc.heuristic() for pPC.

ppc_heuristic <- function(pair, data, alpha, B, whitelist, blacklist, test, skeleton, 
                          dsep.size, complete, max_wthn_sx = dsep.size,
                          sort_pval = TRUE, group = NULL, parentList = NULL, 
                          node.pairs = NULL, key = NULL, true_bn = NULL, debug = FALSE) 
{
  arc = pair$arc
  
  ## partition step
  ## requires distance measure, even if test == "dsep"
  if (dsep.size == 0){
    
    a0 <- bnlearn:::indep.test(x = arc[1], y = arc[2], sx = character(0), 
                               data = data, test = ifelse(test == 'zf', 'cor', test), 
                               B = B, alpha = alpha, learning = FALSE, complete = complete)
    if (!is.null(true_bn))
      p.value <- bnlearn::dsep(bn = true_bn, x = arc[1], y = arc[2])
    else if (!is.null(whitelist) && bnlearn:::is.whitelisted(whitelist, arc, either = TRUE))
      p.value <- 0
    else if (!is.null(blacklist) && bnlearn:::is.blacklisted(blacklist, arc, both = TRUE))
      p.value <- 1
    else
      p.value <- a0[["p.value"]]
    
    return(list(arc = arc, p.value = p.value, 
                dsep.set = character(0), statistic = a0[["statistic"]]))
  }
  if (!is.null(whitelist) && bnlearn:::is.whitelisted(whitelist, arc, either = TRUE))
    return(list(arc = arc, p.value = 0, dsep.set = NULL))
  else if (!is.null(blacklist) && bnlearn:::is.blacklisted(blacklist, arc, both = TRUE))
    return(list(arc = arc, p.value = 1, dsep.set = NULL))
  if (!is.null(pair[["p.value"]]) && pair[["p.value"]] > alpha)
    return(pair)
  
  ## not partition step
  if (is.null(parentList)){
    
    nbr1 = union(skeleton[skeleton[, 2] == arc[1], 1], skeleton[skeleton[, 
                                                                         1] == arc[1], 2])
    nbr2 = union(skeleton[skeleton[, 2] == arc[2], 1], skeleton[skeleton[, 
                                                                         1] == arc[2], 2])
  } else{ 
    
    ## faster, or at least far more scalable
    nbr1 <- parentList[[arc[1]]]
    nbr2 <- parentList[[arc[2]]]
  }
  nbr1 = setdiff(nbr1, arc[2])
  nbr2 = setdiff(nbr2, arc[1])
  
  if ((length(nbr1) < dsep.size) && (length(nbr2) < dsep.size)) 
    return(pair)
  
  if (debug >= 3) {
    cat("----------------------------------------------------------------\n")
    cat("* investigating", arc[1], "-", arc[2], ", d-separating sets of size", 
        dsep.size, ".\n")
    cat("  > neighbours of", arc[1], ":", nbr1, "\n")
  }
  nnodes <- ncol(data)
  if (length(nbr1) >= dsep.size) {
    
    ## sort by p-values
    if (sort_pval && !is.null(node.pairs))
      nbr1 <- nbr1[order(
        sapply(nbr1, function(x) {
          node.pairs[[find_arc_index(key[c(arc[1], x)], 
                                     nnodes)]]$p.value
        })
      )]
    a1 = allsubs.test_(x = arc[1], y = arc[2], sx = nbr1, 
                       min = dsep.size, max = dsep.size, data = data, test = test, 
                       alpha = alpha, B = B, complete = complete, max_wthn_sx = max_wthn_sx,
                       group = group, true_bn = true_bn, debug = debug)
    if (a1["p.value"] > pair$p.value){
      pair$p.value <- a1["p.value"]
      pair$dsep.set <- attr(a1, "dsep.set")
    }
  }
  if (debug >= 3) 
    cat("  > neighbours of", arc[2], ":", nbr2, "\n")
  if (dsep.size == 1) 
    nbr2 = setdiff(nbr2, nbr1)
  if (pair$p.value <= alpha && (length(nbr2) >= dsep.size) && (dsep.size > 0) && !setequal(nbr1, 
                                                                                           nbr2)) {
    ## sort by p-values
    if (sort_pval && !is.null(node.pairs))
      nbr2 <- nbr2[order(
        sapply(nbr2, function(x) {
          node.pairs[[find_arc_index(key[c(arc[2], x)], nnodes)]]$p.value
        })
      )]
    a2 = allsubs.test_(x = arc[2], y = arc[1], sx = nbr2, 
                       min = dsep.size, max = dsep.size, data = data, test = test, 
                       alpha = alpha, B = B, complete = complete, max_wthn_sx = max_wthn_sx,
                       group = group, true_bn = true_bn, debug = debug)
    if (a2["p.value"] > pair$p.value){
      pair$p.value <- a2["p.value"]
      pair$dsep.set <- attr(a2, "dsep.set")
    }
  }
  return(pair)
}



# Version of bnlearn:::pc.heuristic() for screening edges between groups.

btwn_heuristic <- function (pair, data, alpha, B, whitelist, blacklist, 
                            test, skeleton, dsep.size, complete, max_btwn_nbr = ncol(data)-2,
                            sort_pval = TRUE, group = NULL, parentList = NULL, node.pairs = NULL, 
                            key = NULL, round = 1, true_bn = NULL, debug = FALSE)
{
  arc = pair$arc
  if (!is.null(whitelist) && bnlearn:::is.whitelisted(whitelist, arc, either = TRUE))
    return(list(arc = arc, p.value = 0, dsep.set = NULL))
  else if (!is.null(blacklist) && bnlearn:::is.blacklisted(blacklist, arc, both = TRUE))
    return(list(arc = arc, p.value = 1, dsep.set = NULL))
  if (!is.null(pair[["p.value"]]) && pair[["p.value"]] > alpha)
    return(pair)
  if (is.null(parentList)){
    nbr1 = union(skeleton[skeleton[, 2] == arc[1], 1], skeleton[skeleton[, 
                                                                         1] == arc[1], 2])
    nbr2 = union(skeleton[skeleton[, 2] == arc[2], 1], skeleton[skeleton[, 
                                                                         1] == arc[2], 2])
  } else{
    nbr1 <- parentList[[arc[1]]]
    nbr2 <- parentList[[arc[2]]]
  }
  nbr1 = setdiff(nbr1, arc[2])
  nbr2 = setdiff(nbr2, arc[1])
  if ((length(nbr1) == 0) && (length(nbr2) == 0))
    return(pair)
  if (debug >= 3) {
    cat("----------------------------------------------------------------\n")
    cat("* investigating", arc[1], "-", arc[2], ", d-separating sets of size", 
        dsep.size, ".\n")
    cat("  > neighbours of", arc[1], ":", nbr1, "\n")
  }
  max.dsep.size <- dsep.size
  nnodes <- ncol(data)
  
  ## first round of between-edge screening
  if (round == 1){
    nbr3 <- union(nbr1, nbr2)
    
    ## sort by p-value
    if (sort_pval && !is.null(node.pairs))
      nbr3 <- nbr3[order(
        sapply(nbr3, function(x) {
          node.pairs[[find_arc_index(key[c(arc[1], x)], nnodes)]]$p.value + 
            node.pairs[[find_arc_index(key[c(arc[2], x)], nnodes)]]$p.value
        })
      )]
    if (length(nbr3) > max_btwn_nbr)
      nbr3 <- nbr3[seq_len(max_btwn_nbr)]
    dsep.size <- min(length(nbr3), max.dsep.size)
    a3 = allsubs.test_(x = arc[1], y = arc[2], sx = nbr3,
                       min = dsep.size, max = dsep.size, data = data, test = test, 
                       alpha = alpha, B = B, complete = complete, 
                       group = NULL, true_bn = true_bn, debug = debug)
    if (a3["p.value"] > pair$p.value){
      pair$p.value <- a3["p.value"]
      pair$dsep.set <- attr(a3, "dsep.set")
    }
    
    ## second round of between-edge screening
  } else if (round == 2){
    if (length(nbr1) > 0){
      
      ## sort by p-value
      if (sort_pval && !is.null(node.pairs))
        nbr1 <- nbr1[order(
          sapply(nbr1, function(x) {
            node.pairs[[find_arc_index(key[c(arc[1], x)], nnodes)]]$p.value
          })
        )]
      if (length(nbr1) > max_btwn_nbr)
        nbr1 <- nbr1[seq_len(max_btwn_nbr)]
      dsep.size <- min(length(nbr1), max.dsep.size)
      a1 = allsubs.test_(x = arc[1], y = arc[2], sx = nbr1,
                         min = dsep.size, max = dsep.size, data = data, test = test, 
                         alpha = alpha, B = B, complete = complete, 
                         group = group, true_bn = true_bn, debug = debug)
      if (a1["p.value"] > pair$p.value){
        pair$p.value <- a1["p.value"]
        pair$dsep.set <- attr(a1, "dsep.set")
      }
    }
    if (pair$p.value <= alpha && length(nbr2) > 0 && !setequal(nbr1,
                                                               nbr2)){
      ## sort by p-value
      if (sort_pval && !is.null(node.pairs))
        nbr2 <- nbr2[order(
          sapply(nbr2, function(x) {
            node.pairs[[find_arc_index(key[c(arc[2], x)], nnodes)]]$p.value
          })
        )]
      if (length(nbr2) > max_btwn_nbr)
        nbr2 <- nbr2[seq_len(max_btwn_nbr)]
      dsep.size <- min(length(nbr2), max.dsep.size)
      a2 = allsubs.test_(x = arc[2], y = arc[1], sx = nbr2,
                         min = dsep.size, max = dsep.size, data = data, test = test, 
                         alpha = alpha, B = B, complete = complete, 
                         group = group, true_bn = true_bn, debug = debug)
      if (a2["p.value"] > pair$p.value){
        pair$p.value <- a2["p.value"]
        pair$dsep.set <- attr(a2, "dsep.set")
      }
    }
  }
  return(pair)
}



# Modified version of bnlearn:::allsubs.test() that switches
# to somesubs.test() for extended functionality.

allsubs.test_ <- function (x, y, sx, fixed = character(0), data, test, B = 0L, 
                           alpha = 1, min = 0, max = length(sx), complete, 
                           max_wthn_sx = max, group = NULL, true_bn = NULL, 
                           debug = FALSE)
{
  ## TODO: remove because doesn't store dsep.set for pc; perhaps is speed needed
  if (FALSE && is.null(true_bn) && is.null(group)){
    
    res <- .Call(bnlearn:::call_allsubs_test, x = x, y = y, sx = c(fixed, sx),
                 fixed = fixed, data = data, test = test, B = B, alpha = alpha,
                 min = as.integer(min), max = as.integer(max), complete = complete,
                 debug = debug >= 3)
  } else {
    
    res <- somesubs.test(x = x, y = y, sx = sx, fixed = fixed, 
                         data = data, test = test, B = B, alpha = alpha,
                         min = min, max = max, complete = complete,
                         max_wthn_sx = max_wthn_sx, group = group, 
                         true_bn = true_bn, debug = debug)
  }
  return(res)
}



# Modified version of bnlearn:::allsubs.test() that takes
# some inspiration from pcalg::skeleton(). Allows for:
# 1) investigation of some and not all subsets
# 2) investigation of d-separation using true_bn
# 3) stores dsep.set corresponding to the maximum p-value for PATH

somesubs.test <- function (x, y, sx, fixed = character(0), data, test, B = 0L, 
                           alpha = 1, min = 0, max = length(sx), complete,
                           max_wthn_sx = max, group = NULL, true_bn = NULL, debug = FALSE)
{
  ## TODO: rename dsep.set to just set, because it doesn't d-separate
  res <- c(p.value = 0, min.p.value = 0, max.p.value = 0)
  min <- min + length(fixed)
  max <- max + length(fixed)
  sx <- union(fixed, sx)
  p.value <- 0
  for (dsep.size in seq(min, max)){
    S <- seq_len(dsep.size)
    nextS <- list(wasLast = FALSE)
    repeat {
      if (
        dsep.size == 0 ||  # empty conditioning set
        (
          ## none fixed or includes all fixed
          (length(fixed) == 0 || all(fixed %in% sx[S])) &&
          (
            ## if no group, no restriction
            is.null(group) || 
            (
              ## TODO: within set size
              ## set larger than largest within set size, or
              dsep.size > max_wthn_sx ||
              ## different groups, or
              group[x] != group[y] ||  ## TODO: implement criteria (3.5) (b)
              ## same group but conditioning set contains other group(s)
              !all(group[sx[S]] %in% group[x])
            )
          )
        )
      ){
        if (is.null(true_bn)){
          
          ## conditional independence tests
          p.value <- .Call(bnlearn:::call_allsubs_test, x = x, y = y, sx = sx[S],
                           fixed = character(0), data = data, test = test, B = B, alpha = alpha,
                           min = dsep.size, max = dsep.size, complete = complete,
                           debug = debug >= 3)[1]
        } else {
          ## d-separation
          if (length(S)){
            p.value <- bnlearn::dsep(bn = true_bn, x = x, y = y, z = sx[S])
          } else{
            p.value <- bnlearn::dsep(bn = true_bn, x = x, y = y)
          }
          increment.test.counter_(1)
        }
        if (is.na(p.value)) 
          p.value <- 1L
        if (p.value > res["max.p.value"]){
          res["max.p.value"] <- p.value
          attr(res, "dsep.set") <- sx[S]
          if (res["min.p.value"] == 0L)
            res["min.p.value"] <- res["max.p.value"]
        } else if (p.value < res["min.p.value"]){
          res["min.p.value"] <- p.value
        }
        if (p.value > alpha || (alpha == 1 && p.value >= alpha)){
          if (debug >= 3)
            cat(sprintf("   > node %s is independent from %s given %s (p-value: %g)\n", 
                        x, y, paste(sx[S], collapse = " "), as.numeric(p.value)))
          break
        }
        if (debug >= 3)
          cat(sprintf("   > node %s is dependent on %s given %s (p-value: %g)\n", 
                      x, y, paste(sx[S], collapse = " "), as.numeric(p.value)))
        
        if (dsep.size == 0) break
      }
      
      ## get next conditioning set
      nextS <- pcalg::getNextSet(n = length(sx), k = dsep.size, set = S)
      if (nextS$wasLast)
        break
      S <- nextS$nextS
    }
    if (p.value > alpha || (alpha == 1 && p.value >= alpha)) break
  }
  res["p.value"] <- res["max.p.value"]
  return(res)
}



# Compute distances from statistics (mutual information, correlation).

get_distances <- function(nodes, nnodes, key, skeleton, node.pairs, 
                          x, test, B, alpha, complete){
  
  ## TODO: test and support more tests from bnlearn:::test.labels
  
  if (test %in% bnlearn:::available.discrete.tests){
    
    ## based on normalized mutual information for discrete data
    
    distances <- Matrix::sparseMatrix(
      i = key[c(skeleton[,2], nodes)], 
      j = key[c(skeleton[,1], nodes)], 
      x = c(sapply(node.pairs, `[[`, "statistic"), rep(1, nnodes)), 
      dims = rep(nnodes, 2), dimnames = list(nodes, nodes)
    )
    Matrix::diag(distances) <- sapply(nodes, function(node){
      bnlearn:::indep.test(x = node, y = node, sx = character(0), 
                           data = x, test = "mi", 
                           B = B, alpha = alpha, learning = FALSE, 
                           complete = complete)$statistic
    })
    diag <- Matrix::diag(distances)
    distances <- normalize_muti(muti = distances, entropy = diag)
    attr(distances, "Dimnames") <- list(nodes, nodes)
  } 
  else if (test %in% bnlearn:::available.continuous.tests){
    
    ## based on correlation for continuous data
    
    distances <- Matrix::sparseMatrix(
      i = key[c(skeleton[,2], nodes)], 
      j = key[c(skeleton[,1], nodes)], 
      x = c(1 - abs(sapply(node.pairs, `[[`, "statistic")), rep(0, nnodes)), 
      dims = rep(nnodes, 2), dimnames = list(nodes, nodes)
    )
  }
  return(distances)
}



# Algorithm 1 in Gu, J., & Zhou, Q. (2020). Learning Big Gaussian Bayesian Networks: Partition, Estimation and Fusion. J. Mach. Learn. Res., 21, 158-1.
# Get group given distances using modified hierarchical clustering.

get_group <- function(distances, nnodes, max_groups = 20){
  
  clust <- stats::hclust(d = stats::as.dist(distances),  # hierarchical clustering based on distance
                         method = "average")
  
  large <- max(3, ceiling(nnodes * min(0.05, fraction <- 1 / max_groups)))
  k <- cut <- 0
  
  for (i in seq_len(nnodes-1)){
    
    group <- stats::cutree(clust, k = i)
    n_large <- sum(table(group) >= large)
    if (n_large > max_groups ||  # if exceed maximum number of groups
        n_large == 0)      # or if no groups are big enough
      break
    if (n_large > k){
      k <- n_large
      cut <- i
    }
  }
  group <- stats::cutree(clust, cut)
  
  sgroup <- table(group)  # size of groups
  ugroup <- as.integer(names(sgroup))  # unique groups
  temp_distances <- as.matrix(distances + Matrix::t(distances))
  group <- assign_small(temp_distances,
                        group, ugroup, sgroup, 
                        fraction = min(0.05, fraction),
                        linkage = 0)  # 0 for average linkage, 1 for single linkage
  group <- c(group)
  names(group) <- colnames(distances)
  
  return(group)
}



# Convert node.pairs to a parent edge list.

node.pairs2parentList <- function(node.pairs, nodes, alpha){
  
  ## helper helper function for converting between node.pairs and parentList
  ## TODO: optimize 
  
  arcs.still.present = lapply(node.pairs, function(x) {
    if (x$p.value < alpha) 
      return(x$arc)
    else return(NULL)
  })
  arcs = do.call(rbind, arcs.still.present)
  
  amat <- Matrix::Matrix(bnlearn:::arcs2amat(arcs, nodes), sparse = TRUE)
  amat <- amat + Matrix::t(amat)
  parentList <- sparsebnUtils::edgeList(x = amat)
  parentList <- lapply(parentList, function(x) nodes[x])
  
  return(parentList)
}



# Get logical vector indicating which arcs in node.pairs connect
# nodes in different groups according to the vector group.

get_bool_btwn <- function(group){
  
  ## based on the group, creates a vector of length p * (1 - p) / 2 
  ## that indicates whether or not the corresponding nodes in 
  ## node.pairs is between clusters (i.e., whether arc[1] and arc[2]
  ## in node.pairs[[i]]$arc are in different clusters)
  
  nnodes <- length(group)
  nodes <- names(group)
  
  btwn <- do.call(base::c, lapply(seq(nnodes), function(x){
    after <- seq(from = x+1, to = nnodes)
    after <- after[group[after] != group[x]]
    (x > 1) * sum((nnodes-1):(nnodes-x+1)) + (after-x)
  }))
  
  bool_btwn <- rep(FALSE, nnodes*(nnodes-1)/2)
  bool_btwn[btwn] <- TRUE
  rm(btwn)
  
  return(bool_btwn)
}



# Get logical vector indicating which arcs in node.pairs have not
# been found to be conditionally independent.

get_bool_pairs <- function(min_size, node.pairs, alpha, parentList){
  
  sapply(node.pairs, function(x){
    x$p.value <= alpha &&  # not separated
      any(sapply(parentList[x$arc], 
                 length) >= min_size)  # nonzero neighborhoods
  })
}



# Indicate which node pairs are affected in that their neighbor
# set is updated by adding the edges between clusters.

get_affected <- function(affected_indices, nnodes){
  
  affected <- rep(FALSE, nnodes*(nnodes-1)/2)
  if (length(affected_indices) == 0)
    return(affected)
  
  ## from affected indices to all other nodes
  from <- do.call(base::c, lapply(affected_indices, function(x){
    if (x < nnodes) 
      (x > 1) * sum((nnodes-1):(nnodes-x+1)) + seq_len(nnodes-x) 
    else 
      NULL
  }))
  
  ## from all other nodes to affected indices
  to <- do.call(base::c, lapply(affected_indices, function(y){
    if (y > 1) sapply(seq_len(y-1), function(x) (x > 1) * sum((nnodes-1):(nnodes-x+1)) + (y-x)) else NULL
  }))
  
  affected[c(from, to)] <- TRUE
  return(affected)
}



# Find column index corresponding to an arc in a matrix generated by 
# combn(nnodes, 2). Drastically improves scalability as compared to 
# which() for looking up arcs in node.pairs.

find_arc_index <- function(arc_ind, nnodes){
  
  x <- min(arc_ind)
  y <- max(arc_ind)
  
  return((x > 1) * sum((nnodes-1):(nnodes-x+1)) + (y-x))
}



# Print out node.pairs, for debugging

print_node.pairs <- function(node.pairs, alpha){
  
  arcs.still.present = lapply(node.pairs, function(x) {
    if (x$p.value < alpha) 
      return(x$arc)
    else return(NULL)
  })
  skeleton = do.call(rbind, arcs.still.present)
  
  cat("----------------------------------------------------------------\n")
  cat("* remaining arcs:\n")
  print(arcs.rbind(skeleton, skeleton, reverse2 = TRUE))
}