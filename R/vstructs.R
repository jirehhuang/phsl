# Faster version of bnlearn:::vstruct.detect() that uses Rcpp function
# vstruct_centered_on for detecting v-structures using a scalable
# look-up solution for finding the index corresponding to an arc in
# node.pairs.

vstruct.detect_ <- function (nodes, arcs, mb, data, alpha, B = NULL, test, blacklist, 
                             max.sx = ncol(data), complete, true_bn = NULL, debug = FALSE) 
{
  key_cpp <- build_key(nodes = nodes, cpp = TRUE)
  vstruct.list <- NULL
  ## TODO: check why previously required !is.null(blacklist)
  if (!is.null(attr(mb, "dsep.set"))){
    vstruct.list <- tryCatch({
      vstruct.list <- lapply(key_cpp[nodes], vstruct_centered_on, 
                             arcs = apply(arcs, 2, function(x) key_cpp[x]), 
                             dsep_set = attr(mb, "dsep.set"), alpha = alpha, 
                             nodes = nodes, debug = debug)
      lapply(vstruct.list, function(x){
        if (nrow(x)){
          return(data.frame(max_a = x[,1],
                            y = nodes[x[,2]+1],
                            x = nodes[x[,3]+1],
                            z = nodes[x[,4]+1],
                            stringsAsFactors = FALSE))
        } else return(NULL)
      })
    }, error = function(x){
      cat("Error running vstruct_centered_on. Reattempting with vstruct.centered.on. \n")
      return(NULL)
    })
  }
  ## if cpp not successful, reattempt with vstruct.centered.on(),
  ## which is in general much slower and less scalable
  if (is.null(vstruct.list)){
    ## TODO: rewrite to use key
    vstruct.centered.on = function(x, mb, data, dsep.set) {
      if (debug) {
        cat("----------------------------------------------------------------\n")
        cat("* v-structures centered on", x, ".\n")
      }
      tos = arcs[(arcs[, "to"] == x), "from"]
      if (length(tos) < 2) 
        return(NULL)
      tos.combs = bnlearn:::subsets(tos, 2)
      vs = NULL
      for (j in 1:nrow(tos.combs)) {
        y = tos.combs[j, 1]
        z = tos.combs[j, 2]
        if (debug) 
          cat("  * checking", y, "->", x, "<-", z, "\n")
        if (bnlearn:::is.listed(arcs, c(y, z), either = TRUE)) 
          next
        if (is.null(blacklist)) 
          blacklisted.arc = FALSE
        else blacklisted.arc = bnlearn:::is.listed(blacklist, c(y, 
                                                                z), both = TRUE)
        if (!is.null(dsep.set) && !blacklisted.arc) {
          el = dsep.set[[which(sapply(dsep.set, function(x) setequal(x$arc, 
                                                                     c(y, z))))]]
          if (!x %in% el$dsep.set) {
            if (debug) 
              cat("    @ detected v-structure", y, "->", 
                  x, "<-", z, "from d-separating set.\n")
            vs = rbind(vs, data.frame(max_a = el$p.value, 
                                      y, x, z, stringsAsFactors = FALSE))
          }
        }
        else {
          sx = bnlearn:::smaller(setdiff(mb[[y]][["mb"]], c(x, z)),
                                 setdiff(mb[[z]][["mb"]], c(x, y)))
          if (debug)
            cat("    > chosen d-separating set: '", sx,
                "'\n")
          ## allow for true_bn
          a <- allsubs.test_(x = y, y = z, fixed = x, sx = sx,
                             data = data, test = test, B = B, alpha = alpha,
                             max = min(max.sx, length(sx)), complete = complete, 
                             true_bn = true_bn, debug = debug)
          ## TODO: don't remember what this was for
          # if ((a["p.value"] < alpha || (alpha < 1 && a["p.value"] <= alpha)) && 
          #     (is.null(true_vs) || !is.na(prodlim::row.match(data.frame(matrix(c(y, x, z), nrow=1)), 
          #                                                    data.frame(true_vs))))
          # ){
          ## TODO: don't remember why didn't just keep this original condition; dsep?
          # if (a["p.value"] <= alpha){
          if (a["p.value"] < alpha || (alpha < 1 && a["p.value"] <= alpha)){
            if (debug)
              cat("    @ detected v-structure", y, "->",
                  x, "<-", z, "\n")
            vs = rbind(vs, data.frame(max_a = a["max.p.value"],
                                      y, x, z, stringsAsFactors = FALSE))
          }
        }
      }
      return(vs)
    }
    vstruct.list <- sapply(nodes, vstruct.centered.on, mb = mb, data = data, 
                           dsep.set = attr(mb, "dsep.set"), simplify = FALSE)
  }
  ## no v-structures
  if (is.null(vstruct.list))
    vstruct.list <- list()
  return(vstruct.list)
}



# Extended version of bnlearn:::vstruct.apply() with the following
# modifications.
# 1) Sorts v-structures in decreasing order of p-values.
# 2) Sorts v-structures according to how many are in agreement.
# 3) Uses vstruct_apply_cpp() for default applications.
# 4) Uses greedy vstruct_apply_hgi() for HGI.

vstruct.apply_ <- function(arcs, vs, nodes, data = NULL, 
                           rs = NULL, lambda = 0.5, maxp = 8, 
                           debug = FALSE){
  
  start_time <- Sys.time()
  
  ## reorder consideration priority according to
  ## 1) number of agreeing v-structures; 2) p-value
  vs <- vs[order(vs[, "max_a"], decreasing = TRUE), ]
  key_cpp <- build_key(nodes = nodes, cpp = TRUE)
  vs_cpp <- vs2vs_cpp(vs = vs, 
                      nodes = nodes)
  ## for each v-structure, count the total number
  ## of v-structures that share a directed edge
  count <- integer(nrow(vs_cpp))
  vstruct_count_cpp(p = length(nodes), count = count, 
                    vs = vs_cpp, debug = debug)
  vs <- vs[order(count, decreasing = TRUE), , drop = FALSE]
  vs_cpp <- vs2vs_cpp(vs = vs, 
                      nodes = nodes)
  if (is.null(data)){
    
    ## no hgi
    undirMat <- bnlearn:::arcs2amat(arcs = arcs, nodes = nodes)
    start_time_apply <- Sys.time()
    dirMat <- vstruct_apply_cpp(undirMat = undirMat,
                                vs = vs_cpp,
                                nodes = nodes,
                                maxp = maxp,
                                debug = debug)
    dirMat <- apply(dirMat, 2, as.integer)
    rownames(dirMat) <- colnames(dirMat) <- colnames(undirMat)
    arcs <- bnlearn:::amat2arcs(dirMat, nodes)
    attr(arcs, "nscores") <- 0
    attr(arcs, "nvs") <- sum(vs_cpp[,1] < 0)
    vs$max_a <- vs_cpp[,1]
    
  } else{
    
    ## hgi
    undirMat <- bnlearn:::arcs2amat(arcs = arcs, nodes = nodes)
    is_discrete <- class(data[[1]]) %in% c("factor", "integer")
    if (is.null(rs))
      rs <- sapply(nodes, R_loglik_dnode, parents = character(0), 
                   data = data, k = lambda * log(nrow(data)), debug = debug)
    delta <- vs_cpp[,1]
    
    ## greedily apply v-structures
    out <- vstruct_apply_hgi(undirMat = undirMat,
                             reference = rs,
                             nodes = nodes,
                             vs = vs_cpp,
                             data = data,
                             delta = delta,
                             reverse = 0L,
                             is_discrete = is_discrete,
                             k = lambda * log(nrow(data)),
                             maxp = maxp,
                             debug = debug)
    ## check
    # bn <- bnlearn::empty.graph(nodes)
    # bnlearn::amat(bn) <- apply(out$graph, 2, as.integer)
    # bnlearn::score(bn, data)
    # sum(rs)
    ## prepare PDAG
    dirMat <- apply(out$graph, 2, as.integer)
    undirMat <- undirMat - (dirMat | t(dirMat))
    dirMat <- dirMat + undirMat
    rownames(dirMat) <- colnames(dirMat) <- colnames(data)
    arcs <- bnlearn:::amat2arcs(dirMat, nodes)
    attr(arcs, "rs") <- rs
    attr(arcs, "nscores") <- out$nscores
    attr(arcs, "nvs") <- sum(out$ignore == 2)
  }
  end_time <- Sys.time()
  attr(arcs, "time") <- as.numeric(end_time - start_time, unit = 'secs')
  debug_sprintf(debug, "Applied %s v-structures in %s seconds with %s scores",
                attr(arcs, "nvs"), attr(arcs, "time"), attr(arcs, "nscores"))
  return(arcs)
}



# Rcpp implementation of bnlearn::vstructs() due to memory allocation
# issues with bnlearn::vstructs() in 4.7. Example errors include:
#
# double free or  corruption (!prev)
# realloc(): invalid next size
# Aborted (core dumped)
#
# Consider the following example (warning: may crash your R):
# 
# true_bn <- bnrepository("hepar2")
# bnlearn::vstructs(true_bn)

vstructs_ <- function(x, arcs = FALSE, debug = FALSE){
  
  bnlearn:::check.bn.or.fit(x)
  bnlearn:::check.logical(debug)
  
  amat <- bnlearn::amat(x)
  vsmat <- amat * 0
  vs <- vstructs_cpp(amat, vsmat, debug)
  vs <- apply(vs, 2, function(y) bnlearn::nodes(x)[y+1])
  colnames(vs) <- c("X", "Z", "Y")
  
  if (arcs){
    mode(vsmat) <- "integer"
    return(bnlearn:::amat2arcs(a = vsmat, nodes = bnlearn::nodes(x)))
  } else{
    return(vs)
  }
}



# Convert v-structures to Rcpp indices.

vs2vs_cpp <- function(vs, nodes){
  
  key_cpp <- build_key(nodes = nodes, cpp = TRUE)
  vs[,2:4] <- apply(vs[,2:4, drop=FALSE], 2, function(x) key_cpp[x])
  
  return(as.matrix(vs))
}