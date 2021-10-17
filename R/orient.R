# Orient edges
# 
# Extended version of bnlearn:::learn.arc.directions() that includes
# p-value adjacency thresholding (PATH) and hybrid greedy initialization
# (HGI). 

orient <- function(x, cluster = NULL, local.structure, whitelist, blacklist, 
                   test, score, alpha, complete, B = NULL, max.sx = ncol(x)-2, 
                   maxp = ncol(x)-1, path = 1, min_alpha = 1e-5, hgi = FALSE, 
                   true_bn = NULL, debug = FALSE)
{
  times <- tests <- 
    c(detect = bnlearn::test.counter(), apply = 0, rules = 0, select = 0)
  
  test <- bnlearn:::check.test(test, x)
  score <- bnlearn:::check.score(score, x)
  alpha <- bnlearn:::check.alpha(alpha)
  B <- bnlearn:::check.B(B, test)
  max.sx <- bnlearn:::check.largest.sx.set(max.sx, 
                                           x)
  maxp <- bnlearn:::check.maxp(maxp,
                               x)
  nodes = names(x)
  arcs = bnlearn:::mb2arcs(local.structure, nodes)
  to.drop = !apply(arcs, 1, function(x) {
    bnlearn:::is.blacklisted(blacklist, x)
  })
  arcs = arcs[to.drop, , 
              drop = FALSE]
  
  is_discrete <- class(x[[1]]) %in% c("factor", "integer")
  if (path <= 1){
    path <- 1
    min_alpha <- alpha
    alpha_path <- alpha
  } else{
    stopifnot(min_alpha < alpha)
    min_alpha <- bnlearn:::check.alpha(min_alpha)
  }
  best_alpha <- alpha
  
  
  ## detect v-structures
  start_time_d <- Sys.time()
  if (!is.null(true_bn)){  
    ## if true_bn is supplied, local.structure is the true skeleton,
    ## so vstruct.detect_ will recover all and only true v-structures
    ## thus, for speed purposes, simply extract v-structures from true_bn
    vs <- cbind(data.frame(max_a = 0), 
                as.data.frame(vstructs_(x = true_bn), 
                              stringsAsFactors = FALSE))
    names(vs) <- c("max_a", "y", "x", "z")
  } else{
    vs <- do.call("rbind", vstruct.detect_(nodes = nodes, arcs = arcs, mb = local.structure, 
                                           data = x, alpha = alpha, B = B, test = test, 
                                           blacklist = blacklist, max.sx = max.sx, 
                                           complete = complete, true_bn = true_bn, 
                                           debug = debug > 1))
    if (is.null(vs))
      vs <- matrix(nrow = 0, ncol = 4)
  }
  rownames(vs) = NULL
  end_time_d <- Sys.time()
  times[1] <- ifelse(is.null(true_bn), 
                     as.numeric(end_time_d - start_time_d, unit = 'secs'), 0)
  tests[1] <- bnlearn::test.counter() - tests[1]
  debug_sprintf(debug, "Detected %g v-structures in %g seconds with %g tests", 
                nrow(vs), times[1], tests[1])
  
  ## initialize
  lambda <- 0.5  # TODO: lambda parameter
  attr(arcs, "alpha") <- alpha
  
  pdag <- list(learning = list(whitelist = whitelist, blacklist = blacklist, test = test, 
                               args = list(alpha = alpha, min_alpha = alpha, best_alpha = alpha),
                               ntests = 0, max.sx = max.sx, 
                               group = attr(local.structure, "learning")$group,
                               path = 1, hgi = hgi, valid = TRUE),
               nodes = NULL, arcs = arcs, 
               arcs_path = list(arcs))
  if (!is.null(B)) 
    pdag$learning$args$B = B
  names(pdag$arcs_path) <- alpha
  
  ## if there are any estimated arcs
  if (nrow(arcs)){
    
    key <- build_key(nodes = nodes)
    
    if (path > 1){
      
      ## build alpha_path
      arc_indices <- apply(arcs, 1, function(arc) find_arc_index(arc_ind = key[arc], 
                                                                 nnodes = length(nodes)))
      p_vec <- unname(sapply(arc_indices, function(x) attr(local.structure, 
                                                           "dsep.set")[[x]]$p.value))
      path <- min(path, sum((unique_p_vec <- unique(p_vec)) <= alpha) + 1)
      alpha_path <- build_alpha_path(p_vec = p_vec, alpha = alpha, 
                                     min_alpha = min_alpha, n_alpha = path)
      min_alpha <- max(min_alpha, min(alpha_path))
      path <- length(alpha_path)
      
      ## list of v-structure indicators for each alpha
      bool_vs_path <- lapply(alpha_path, function(a){
        bool_arcs <- (p_vec <= a)
        if (!is.null(vs) && nrow(vs) && any(bool_arcs)){
          bool_vs <- apply(vs, 1, function(v){
            i <- find_arc_index(arc_ind = sort(key[v[c("x", "y")]]), 
                                nnodes = length(nodes))
            j <- find_arc_index(arc_ind = sort(key[v[c("x", "z")]]), 
                                nnodes = length(nodes))
            if (!setequal(v[c("x", "y")],
                          attr(local.structure, "dsep.set")[[i]]$arc) ||
                !setequal(v[c("x", "z")],
                          attr(local.structure, "dsep.set")[[j]]$arc)){
              cat("Problem\n")
            }
            max(attr(local.structure, "dsep.set")[[i]]$p.value,
                attr(local.structure, "dsep.set")[[j]]$p.value) <= a
          })
          bool_vs <- vs[, "max_a"] > a & bool_vs
        } else{
          bool_vs <- logical(0)
        }
        attr(bool_vs, "alpha") <- a
        return(bool_vs)
      })
    } else {
      p_vec <- numeric(nrow(arcs))
      bool_vs_path <- list(sapply(seq_len(nrow(vs)), function(x) TRUE))
      attr(bool_vs_path[[1]], "alpha") <- alpha
    }
    ## cache initial reference scores and score deltas for hgs
    if (hgi){
      rs_delta <- get_rs_delta(arcs = arcs, vs = vs, nodes = nodes, data = x, debug = FALSE)
    } else{
      rs_delta <- NULL
    }
    ## determine arcs_path
    pdag$arcs_path <- lapply(seq_len(length(bool_vs_path)), function(i){
      
      bool_vs <- bool_vs_path[[i]]
      
      ## TODO: check if all necessary
      nundir <- 0
      arcs <- arcs[p_vec <= attr(bool_vs, "alpha"), , drop=FALSE]
      vs_temp <- vs[as.logical(bool_vs), , drop=FALSE]
      rs <- sapply(nodes, function(x) 0)
      
      ## apply v-structures, if any
      if (nrow(vs_temp)){
        
        ## TODO: check if still having memory reference issues
        if (hgi){
          deep_copy_NumericVector(original = rs_delta$rs, copy = rs)
          deep_copy_NumericVector(original = rs_delta$delta[bool_vs], 
                                  copy = vs_temp$max_a)
        }
        arcs <- vstruct.apply_(arcs = arcs, vs = vs_temp, nodes = nodes, 
                               data = if (hgi) x else NULL, rs = rs,
                               lambda = lambda, maxp = maxp, debug = debug > 1)
        times[2] <- attr(arcs, "time")
        tests[2] <- attr(arcs, "nscores")
        
        ## TODO: check if still having memory reference issues
        if (!is.null(attr(arcs, "rs")))
          deep_copy_NumericVector(original = attr(arcs, "rs"), copy = rs)
      }
      ## apply orientation rules (R1-4 or extend)
      if (nrow(arcs)){
        
        start_time <- Sys.time()
        if (!hgi){
          nscores <- 0
          pdag = list(learning = list(), nodes = structure(rep(0, length(nodes)), 
                                                           names = nodes), arcs = arcs)
          # arcs <- bnlearn:::cpdag.backend(pdag, moral = TRUE, 
          #                                 fix = TRUE, debug = debug > 1)$arcs
          amat <- bnlearn:::arcs2amat(arcs = pdag$arcs, nodes = nodes)
          amat <- apply_cpdag_rules(pdag = amat, nodes = nodes, 
                                    remove_invalid = TRUE, debug = debug)
          arcs <- bnlearn:::amat2arcs(a = amat, nodes = nodes)
          
        } else{
          
          ## TODO: check why would this be necessary
          # if (sum(rs) == 0){
          #   if (sum(rs_delta$rs) != 0){
          #     rs <- sapply(nodes, function(x) 0)
          #     deep_copy_NumericVector(original = rs_delta$rs, copy = rs)
          #   } else{
          #     rs <- sapply(nodes, R_loglik_dnode, parents = character(0),
          #                  data = data, k = lambda * log(nrow(x)), debug = debug)
          #   }
          # }
          
          ## initialize
          pdag <- bnlearn:::arcs2amat(arcs, nodes)
          dag <- pdag - pdag * t(pdag)
          mode(dag) <- "integer"
          nscores <- numeric(1)
          
          ## greedy extension
          a <- pdag2dag_greedy(a = pdag,
                               d = dag,
                               reference = rs,
                               nodes = nodes,
                               data = x,
                               nscores = nscores,
                               is_discrete = is_discrete,
                               k = lambda * log(nrow(x)), 
                               maxp = maxp,
                               verbose = FALSE)
          
          ## transfer to arcs, including deleted arcs as undirected edges
          ## TODO: check if undirecting cycles necessary
          # dag <- undirect_cycles(dag)
          dag <- dag - dag * t(dag)
          pdag <- bnlearn:::arcs2amat(arcs, nodes)
          pdag[dag | t(dag)] <- 0
          mode(pdag) <- "integer"
          arcs <- bnlearn:::amat2arcs(pdag + dag, nodes)
          attr(arcs, "rs") <- rs
          # attr(arcs, "nundir") <- (sum(a * t(a)) / 2)
        }
        end_time <- Sys.time()
        times[3] <- as.numeric(end_time - start_time, unit = 'secs')
        tests[3] <- nscores
        debug_sprintf(debug, "Applied rules for estimate %g in %g seconds with %g scores", 
                      i, times[3], tests[3])
      }
      attr(arcs, "alpha") <- attr(bool_vs, "alpha")
      attr(arcs, "times") <- times
      attr(arcs, "tests") <- tests
      
      return(arcs)
    })
    names(pdag$arcs_path) <- alpha_path
    for (i in pdag$arcs_path){
      times[2:3] <- times[2:3] + attr(i, "times")[2:3]
      tests[2:3] <- tests[2:3] + attr(i, "tests")[2:3]
    }
    if (length(pdag$arcs_path) > 1){
      
      start_time_s <- Sys.time()
      
      if (!hgi){
        
        pdag$dag_path <- score_arcs_path(arcs_path = pdag$arcs_path, nodes = nodes, x = x, 
                                         score = score, lambda = lambda, debug = debug)
        tests[4] <- attr(pdag$dag_path, "nscores")
        which_best <- attr(pdag$dag_path, "which_best")
        pdag$learning$valid <- attr(pdag$dag_path, "success")
        pdag$arcs <- pdag$arcs_path[[which_best]]
        
      } else{
        
        which_best <- which.max(sapply(pdag$arcs_path, 
                                       function(arcs) sum(attr(arcs, "rs"))))
        ## keep all undirected edges from the original alpha
        ## but keep only directed edges in the best hgi estimate
        a0 <- bnlearn:::arcs2amat(pdag$arcs_path[[1]], nodes)
        a1 <- bnlearn:::arcs2amat(pdag$arcs_path[[which_best]], nodes)
        a1 <- (a0 | t(a0)) - (a1 | t(a1)) + a1
        pdag$arcs <- bnlearn:::amat2arcs(a1,
                                         nodes)
        attr(pdag$arcs, "alpha") <- attr(pdag$arcs_path[[which_best]], "alpha")
      }
      end_time_s <- Sys.time()
      times[4] <- as.numeric(end_time_s - start_time_s, unit = "secs")
    } else{
      
      pdag$arcs <- pdag$arcs_path[[1]]
      pdag$learning$valid <- attr(arcs2dag(arcs = pdag$arcs, 
                                           nodes = nodes), "success")
    }
    pdag$learning[c("tests", "path")] <- list(sum(tests),
                                              length(alpha_path))
    pdag$learning$arcs1 <- pdag$arcs_path[[1]]
    pdag$learning$args[c("min_alpha", "best_alpha")] <- list(min(alpha_path),
                                                             attr(pdag$arcs, "alpha"))
    pdag$nodes <- bnlearn:::cache.structure(nodes, 
                                            arcs = pdag$arcs)
  }
  debug_sprintf(debug, "Oriented edges for %g estimate(s) in %g seconds with %g calls", 
                length(pdag$arcs_path), sum(times), sum(tests))
  return(pdag)
}



# Score solution path generated by PATH and selects highest-scoring estimate.

score_arcs_path <- function(arcs_path, nodes, x, score, 
                            lambda = 0.5, debug = FALSE){
  
  start_time <- Sys.time()
  nscores <- bnlearn::test.counter()
  ns0 <- sapply(nodes, function(x) 0)
  
  if (length(arcs_path) > 1){
    
    bn_i <- bn_im1 <- bnlearn::empty.graph(nodes = nodes)
    # extra.args <- list(k = lambda * log(nrow(x)))
    dag_path <- list()
    
    for (i in seq_len(length(arcs_path))){
      
      arcs <- arcs_path[[i]]
      if (is.null(attr(arcs, "dag"))) 
        attr(arcs, "dag") <- arcs2dag(arcs = arcs, nodes = nodes, 
                                      direct_all = TRUE)
      dag_path[[i]] <- list(dag = attr(arcs, "dag"),
                            delta = 0, ns = ns0,
                            alpha = attr(arcs, 'alpha'))
      dag_path[[i]]$eL <- sparsebnUtils::edgeList(dag_path[[i]]$dag)
      
      if (i > 1){
        differences <- Matrix::colSums(!(dag_path[[i-1]]$dag == 
                                           dag_path[[i]]$dag))
        different <- names(differences[differences > 0])
        dag_path[[i]]$ns <- dag_path[[i-1]]$ns
        if (length(different)){
          
          bnlearn::amat(bn_i) <- as.matrix(dag_path[[i]]$dag)
          bnlearn::amat(bn_im1) <- as.matrix(dag_path[[i-1]]$dag)
          
          for (node in different){
            if (dag_path[[i-1]]$ns[node] >= 0){
              ## score previous node if not yet scored
              # dag_path[[i-1]]$ns[node] <- R_loglik_dnode(node, nodes[dag_path[[i-1]]$eL[[node]]],
              #                                            x, extra.args$k, debug > 1)
              extra.args <- bnlearn:::check.score.args(score = score, network = bn_im1, data = x, 
                                                       extra.args = list(), learning = TRUE)
              dag_path[[i-1]]$ns[node] <- bnlearn:::per.node.score(bn_im1, x, score, node,
                                                                   extra.args, debug = debug > 1)
            }
            ## score current node
            # dag_path[[i]]$ns[node] <- R_loglik_dnode(node, nodes[dag_path[[i]]$eL[[node]]],
            #                                                   x, extra.args$k, debug > 1)
            increment.test.counter_(1)  # one per score difference
            extra.args <- bnlearn:::check.score.args(score = score, network = bn_i, data = x, 
                                                     extra.args = list(), learning = TRUE)
            dag_path[[i]]$ns[node] <- bnlearn:::per.node.score(bn_i, x, score, node,
                                                               extra.args, debug = debug > 1)
          }
        }
        ## compute score difference between estimates
        dag_path[[i]]$delta <- sum(dag_path[[i]]$ns[different] - 
                                     dag_path[[i-1]]$ns[different])
        debug_sprintf(debug, "Score delta of %s between estimates %s and %s", 
                      dag_path[[i]]$delta, i-1, i)
      }
    }
  } else{
    
    dag_path <- list(list(dag = arcs2dag(arcs = arcs_path[[1]], nodes = nodes),
                          delta = 0, ns = ns0,
                          alpha = attr(arcs_path[[1]], 'alpha')))
  }
  successes <- sapply(dag_path, function(x) attr(x$dag, "success"))
  if (success <- any(sapply(successes, function(x) is.logical(x) && x))){
    which_best <- which(successes)[which.max(cumsum(sapply(dag_path, 
                                                           `[[`, 'delta'))[successes])]
    debug_sprintf(debug, "%g admissible PDAG(s) found out of %g", 
                  sum(successes), length(successes))
  } else{
    which_best <- which.max(cumsum(sapply(dag_path, `[[`, 'delta')))
    debug_sprintf(debug, "No admissible PDAG found out of %g", 
                  length(successes))
  }
  names(dag_path) <- sapply(arcs_path, function(x) attr(x, 'alpha'))
  attr(dag_path, 'which_best') <- which_best
  attr(dag_path, 'success') <- success
  attr(dag_path, 'nscores') <- bnlearn::test.counter() - nscores
  end_time <- Sys.time()
  debug_sprintf(debug, 
                "Chose alpha = %g out of %g different significance levels with %g scores in %g seconds", 
                dag_path[[which_best]]$alpha, length(dag_path), attr(dag_path, 'nscores'), 
                as.numeric(end_time - start_time, units = "secs"))
  return(dag_path)
}



# Generate evenly spaced out sequence of alpha values for PATH.

build_alpha_path <- function(p_vec, alpha, 
                             min_alpha = 1e-5, n_alpha = 10){
  
  if (min_alpha >= alpha || 
      abs(alpha - min_alpha) < alpha^2){
    if (n_alpha == 1)
      return(alpha)
    min_alpha <- alpha * 1e-1
  }
  
  # n_alpha <- n_alpha - (n_alpha > 1)
  alpha_path <- unique(sort(p_vec[p_vec >= min_alpha & p_vec <= 1], 
                            decreasing = TRUE))[-1]  # first element corresponds to alpha
  alpha_path <- union(alpha, alpha_path)
  
  if (length(alpha_path) <= 1)
    return(alpha)
  
  ind <- unique(floor(seq(1, length(alpha_path), length = n_alpha)))
  return(alpha_path[ind])
}



# Compute reference scores (rs) and score differences (delta). 

get_rs_delta <- function(arcs, vs, nodes, data, lambda = 0.5, debug = FALSE){
  
  ## currently, only discrete implementation of hgi
  is_discrete <- TRUE
  rs <- sapply(nodes, R_loglik_dnode, parents = character(0),
               data = data, k = lambda * log(nrow(data)), debug = debug, USE.NAMES = TRUE)
  
  ## TODO: eventually generalizable
  # network <- bnlearn::empty.graph(nodes = nodes)
  # score <- bnlearn:::check.score(NULL, data = data)
  # extra.args <- bnlearn:::check.score.args(score = score, network = network, 
  #                                          data = x, extra.args = list(), learning = TRUE)
  # rs <- bnlearn:::per.node.score(network = network, data = data, score = score, 
  #                                targets = nodes, extra.args = extra.args, debug = debug)
  
  if (is.null(vs) || nrow(vs) == 0)
    return(list(rs = rs, delta = numeric(0), 
                nscores = integer(1)))
  
  vs_cpp <- vs2vs_cpp(vs, nodes)
  undirMat <- bnlearn:::arcs2amat(arcs = arcs, 
                                  nodes = nodes)
  out <- vstruct_apply_hgi(undirMat = undirMat,
                           reference = rs,
                           nodes = nodes,
                           vs = vs_cpp,
                           data = data,
                           delta = numeric(0),
                           just_delta = TRUE,
                           is_discrete = is_discrete,
                           k = lambda * log(nrow(data)),
                           debug = debug)
  return(list(rs = rs, 
              delta = c(out$delta), nscores = out$nscores))
}



# Extend a PDAG in the form of arcs to a DAG in the form
# of an adjacency matrix.

arcs2dag <- function(arcs, nodes, direct_all = FALSE){
  
  ## TODO: add maxp argument
  
  a <- bnlearn:::arcs2amat(arcs, nodes)
  ## TODO: check if necessary
  # if (undirect)
  #   a <- undirect_cycles(a)
  ## if any undirected edges
  if (any(a == 1 & a == t(a))){
    extend <- pdag2dag_cpp(g = a, nodes = nodes,
                           direct_all = direct_all)
    dag <- Matrix::Matrix(extend$graph, sparse = TRUE)
    attr(dag, "success") <- extend$success
    attr(dag, "undirected") <- extend$undirected
    attr(dag, "nundir") <- (sum(a | t(a)) - sum(extend$graph | t(extend$graph))) / 2
  } else{
    dag <- Matrix::Matrix(a, sparse = TRUE)
    attr(dag, "success") <- TRUE
    attr(dag, "undirected") <- 0
    attr(dag, "nundir") <- 0
  }
  rownames(dag) <- colnames(dag) <- nodes
  return(dag)
}