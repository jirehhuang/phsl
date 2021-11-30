# Learn skeleton of a Bayesian network structure
# 
# Returns \code{local.structure} from bnlearn:::bnlearn().

skeleton <- function(x, cluster = NULL, whitelist = NULL, blacklist = NULL, test = NULL, 
                     alpha = NULL, B = NULL, method = "gs", max.sx = ncol(x)-2,
                     max_wthn_sx = max.sx, max_btwn_sx = max.sx, max_btwn_nbr = ncol(x)-2, 
                     sort_pval = TRUE, max_groups = 20, true_bn = NULL, debug = FALSE){
  start_time <- Sys.time()
  bnlearn::reset.test.counter()
  res = NULL
  parallel = FALSE
  data.info = bnlearn:::check.data(x, allow.missing = TRUE)
  if (!method %in% c("ppc", "true", "cig"))
    bnlearn:::check.learning.algorithm(method, class = "constraint")
  test = bnlearn:::check.test(test, x)
  debug_cli(!(is.numeric(debug) || is.logical(debug)), cli::cli_abort, 
            "debug must be logical or numeric")  # numeric for ppc_skeleton()
  alpha = bnlearn:::check.alpha(alpha)
  B = bnlearn:::check.B(B, test)
  max.sx = bnlearn:::check.largest.sx.set(max.sx, x)
  if (!is.null(cluster)) {
    bnlearn:::check.cluster(cluster)
    parallel = TRUE
    bnlearn:::slaves.setup(cluster)
    if (debug) {
      debug_cli(TRUE, cli::cli_alert_warning,
                "disabled debugging output for parallel computing")
      debug = FALSE
    }
  }
  whitelist = bnlearn:::build.whitelist(whitelist, nodes = names(x), 
                                        data = x, algo = method, criterion = test)
  blacklist = bnlearn:::build.blacklist(blacklist, whitelist, names(x), 
                                        algo = method)
  full.blacklist = bnlearn:::arcs.rbind(blacklist, bnlearn:::check.arcs.against.assumptions(NULL, 
                                                                                            x, test))
  full.blacklist = bnlearn:::unique.arcs(full.blacklist, 
                                         names(x)) 
  if (method == "true"){
    bnlearn:::check.bn(true_bn)
    trueSKEL <- bnlearn::amat(true_bn)
    trueSKEL <- trueSKEL + t(trueSKEL)
    skeleton <- bnlearn:::amat2arcs(trueSKEL, names(x))
    local.structure <- 
      bnlearn:::cache.structure(names(x), bnlearn:::arcs.rbind(skeleton, skeleton, 
                                                               reverse2 = TRUE))
  } 
  else if (method == "cig"){
    bnlearn:::check.bn(true_bn)
    nodes <- bnlearn::nodes(true_bn)
    vs <- vstructs_(true_bn)
    trueCIG <- bnlearn::amat(true_bn)
    trueCIG[vs[, c(1, 3), drop=FALSE]] <- 1L
    trueCIG <- 1L * (trueCIG | t(trueCIG))
    skeleton <- bnlearn:::amat2arcs(trueCIG, nodes)
    local.structure <- 
      bnlearn:::cache.structure(nodes, bnlearn:::arcs.rbind(skeleton, skeleton,
                                                            reverse2 = TRUE))
  } 
  else if (method == "ppc"){
    max_wthn_sx <- bnlearn:::check.largest.sx.set(max_wthn_sx, x)
    max_btwn_sx <- bnlearn:::check.largest.sx.set(max_btwn_sx, x)
    max_btwn_nbr <- bnlearn:::check.largest.sx.set(max_btwn_nbr, x)
    local.structure <- 
      ppc_skeleton(x = x, cluster = cluster, whitelist = whitelist, 
                   blacklist = full.blacklist, test = test, alpha = alpha, 
                   B = B, complete = data.info$complete.nodes,
                   max.sx = max.sx, max_wthn_sx = max_wthn_sx,
                   max_btwn_sx = max_btwn_sx, max_btwn_nbr = max_btwn_nbr,
                   sort_pval = sort_pval, max_groups = max_groups, 
                   true_bn = true_bn, debug = debug)
  }
  else if (method == "pc.stable") {
    local.structure = 
      bnlearn:::pc.stable.backend(x = x, whitelist = whitelist, 
                                  blacklist = full.blacklist, test = test, alpha = alpha, 
                                  B = B, max.sx = max.sx, debug = debug >= 3, cluster = cluster, 
                                  complete = data.info$complete.nodes)
  }
  else if (method == "gs") {
    local.structure = 
      bnlearn:::grow.shrink(x = x, whitelist = whitelist, 
                            blacklist = full.blacklist, test = test, alpha = alpha, 
                            B = B, max.sx = max.sx, debug = debug >= 3, cluster = cluster, 
                            complete = data.info$complete.nodes)
  }
  else if (method == "iamb") {
    local.structure = 
      bnlearn:::incremental.association(x = x, whitelist = whitelist, 
                                        blacklist = full.blacklist, test = test, alpha = alpha, 
                                        B = B, max.sx = max.sx, debug = debug >= 3, cluster = cluster, 
                                        complete = data.info$complete.nodes)
  }
  else if (method == "fast.iamb") {
    local.structure = 
      bnlearn:::fast.incremental.association(x = x, 
                                             whitelist = whitelist, blacklist = full.blacklist, 
                                             test = test, alpha = alpha, B = B, max.sx = max.sx, 
                                             debug = debug >= 3, cluster = cluster, complete = data.info$complete.nodes)
  }
  else if (method == "inter.iamb") {
    local.structure = 
      bnlearn:::inter.incremental.association(x = x, 
                                              whitelist = whitelist, blacklist = full.blacklist, 
                                              test = test, alpha = alpha, B = B, max.sx = max.sx, 
                                              debug = debug >= 3, cluster = cluster, complete = data.info$complete.nodes)
  }
  else if (method == "iamb.fdr") {
    local.structure = 
      bnlearn:::incremental.association.fdr(x = x, 
                                            whitelist = whitelist, blacklist = full.blacklist, 
                                            test = test, alpha = alpha, B = B, max.sx = max.sx, 
                                            debug = debug >= 3, cluster = cluster, complete = data.info$complete.nodes)
  }
  else if (method == "mmpc") {
    local.structure = 
      bnlearn:::maxmin.pc(x = x, whitelist = whitelist, 
                          blacklist = full.blacklist, test = test, alpha = alpha, 
                          B = B, debug = debug >= 3, max.sx = max.sx, cluster = cluster, 
                          complete = data.info$complete.nodes)
  }
  else if (method == "si.hiton.pc") {
    local.structure = 
      bnlearn:::si.hiton.pc.backend(x = x, whitelist = whitelist, 
                                    blacklist = full.blacklist, test = test, alpha = alpha, 
                                    B = B, max.sx = max.sx, debug = debug >= 3, cluster = cluster, 
                                    complete = data.info$complete.nodes)
  }
  else if (method == "hpc") {
    local.structure = 
      bnlearn:::hybrid.pc.backend(x = x, whitelist = whitelist, 
                                  blacklist = full.blacklist, test = test, alpha = alpha, 
                                  B = B, max.sx = max.sx, debug = debug >= 3, cluster = cluster, 
                                  complete = data.info$complete.nodes)
  }
  attr(local.structure, 
       "learning") <- list(whitelist = whitelist, blacklist = blacklist, test = test, 
                           args = list(alpha = alpha, max.sx = max.sx, max_wthn_sx = max_wthn_sx,
                                       max_btwn_sx = max_btwn_sx, max_btwn_nbr = max_btwn_nbr,
                                       max_groups = max_groups), 
                           group = attr(local.structure, "learning")$group,
                           ntests = bnlearn::test.counter(), 
                           complete = data.info$complete.nodes)
  if (!is.null(B)) 
    attr(local.structure, "learning")$args$B = B
  
  end_time <- Sys.time()
  debug_cli(debug, cli::cli_alert_success,
            c("completed skeleton estimation with {method} ",
              "in {prettyunits::pretty_sec(as.numeric(end_time - start_time, unit = 'secs'))} ",
              "with {bnlearn::test.counter()} calls"))
  
  return(invisible(local.structure))
}