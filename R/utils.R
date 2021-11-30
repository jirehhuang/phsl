#' Download from Bayesian network repository
#'
#' Read bn.fit objects from the \code{\link[bnlearn]{bnlearn}} Bayesian network repository directly from \href{www.bnlearn.com}{www.bnlearn.com} using the RDS links.
#'
#' @param x character value indicating desired Bayesian network.
#' @return A Bayesian network as an object of class \code{bn.fit}.
#' @author Jireh Huang (\email{jirehhuang@@ucla.edu})
#' @examples
#' ## Read Bayesian network object 
#' true_bn <- bnrepository("child")
#' @export

bnrepository <- function(x){
  
  avail_bn <- c(
    # Discrete Bayesian networks
    "asia", "cancer", "earthquake", "sachs", "survey",
    "alarm", "barley", "child", "insurance", "mildew", "water",
    "hailfinder", "hepar2", "win95pts",
    "munin", "munin1", "munin2", "munin3", "munin4",
    "andes", "diabetes", "link", "pathfinder", "pigs",
    # Gaussian Bayesian networks
    "ecoli70", "magic-niab", "magic-irri", "arth150",
    # Conditional linear Gaussian Bayesian networks
    "healthcare", "sangiovese", "mehra"
  )
  if (!x %in% avail_bn){
    
    stop(sprintf("x must be one of %s",
                 paste(avail_bn, collapse = ", ")))
  }
  
  x1 <- ifelse(x %in% sprintf("munin%s", seq_len(4)), "munin4", x)
  x2 <- ifelse(x == "mehra", "mehra-complete", x)
  
  bn.fit <- readRDS(
    file = url(sprintf("https://www.bnlearn.com/bnrepository/%s/%s.rds",
                       x1, x2))
  )
  return(bn.fit)
}



# Width of function portion of debugging output

DEBUG_WIDTH <- 18



# Determine minimum debug width

debug_width <- function(name = "phsl", add = 2){
  
  ns <- ls(getNamespace(name = name))
  ns <- ns[!grepl(sprintf("_%s", name), ns)]
  max(nchar(ns)) + add
}



# Print debugging output with cli

debug_cli <- function(debug,
                      fun = cli::cli_alert,
                      text = "",
                      ...){
  
  if (debug){
    
    ## identify calling function in namespace
    ns <- ls(getNamespace(name = "phsl"))
    which <- -1
    repeat{
      fn <- sys.call(which = which)[1]
      fn <- gsub("\\(.*", "", as.character(fn))
      fn <- gsub(".*:", "", fn)
      if (length(fn) == 0 || fn %in% ns) break
      which <- which - 1
    }
    if (length(fn) == 0)
      fn <- "[UNKNOWN]"
    
    fn <- sprintf("{.field {.strong %s}}:", fn)
    fn <- stringr::str_pad(fn, width = max(DEBUG_WIDTH + 10 + 9,
                                           nchar(fn) + 2), side = "right")
    
    ## text message
    text <- c(fn, text)  # glue
    
    ## prepare and execute function
    if (!is.function(fun)){
      
      ## TODO: replace with cli::cli_text
      fun <- cli::cli_alert
    }
    if (identical(fun, cli::cli_progress_bar)){
      
      text <- c(cli::symbol$arrow_right, " ", text)
    }
    
    args <- list(...)
    if (is.null(args[[".envir"]]))
      args[[".envir"]] <- sys.frame(which = which)
    
    ## add text
    formals_nms <- names(formals(fun))
    if ("text" %in% formals_nms){
      
      args$text <- text
      
    } else if ("msg" %in% formals_nms){
      
      args$msg <- text
      
    } else if ("message" %in% formals_nms){
      
      ## TODO: check glue behavior of cli::cli_abort()
      args$message <- text
      
    } else if ("format" %in% formals_nms){
      
      args$format <- text
      
    }
    ## modify other arguments
    if ("format_done" %in% names(args)){
      
      args$format_done <- c(green_tick, " ", fn, args$format_done)
    }
    if ("format_failed" %in% names(args)){
      
      args$format_failed <- c(red_cross, " ", fn, args$format_failed)
    }
    do.call(what = fun, args = args)
  }
}



# Print debugging output
# 
# Convenience function for printing debugging output.
# 
# @param debug logical value that activates printing debugging output.
# @param fmt character value input to \code{\link[base]{sprintf}}.
# @param ... additional arguments passed into \code{\link[base]{sprintf}}.
# @return None.
# @author Jireh Huang (\email{jirehhuang@@ucla.edu})
# @examples
# fn <- function(debug = FALSE){
# 
#   set.seed(1)
#   number <- rnorm(1)
#   string = "error"
# 
#   debug_sprintf(debug, "number = %g, string = %s",
#                 number, string)
# }
# fn(debug = TRUE)

debug_sprintf <- function(debug, fmt, ...){
  
  # deprecated; replaced by debug_cli()
  
  if (debug){
    
    ## version 1: can be slow
    # fn <- gsub("\\(.*",
    #            "", as.character(sys.calls()))
    # fn <- sprintf("%s:", fn[length(fn)-1])
    ## version 2: still slow
    # fn <- sys.calls()
    # fn <- fn[length(fn)-1]
    # fn <- gsub("\\(.*", "", as.character(fn))
    # fn <- sprintf("%s:", fn)
    ## version 3: so far so good
    fn <- sys.call(-1)[1]
    fn <- gsub("\\(.*", "", as.character(fn))
    fn <- sprintf("%s:", fn)
    ## version n: gave up
    # fn <- character(0)
    if (length(fn) == 0)
      fn <- "[UNKNOWN]:"
    
    cat(sprintf("%s%s\n",
                stringr::str_pad(fn, width = DEBUG_WIDTH, side = "right"),
                sprintf(fmt, ...)))
  }
}



# Increment test counter
# 
# Wrapper for \code{\link[bnlearn]{increment.test.counter}}. Sometimes,
# \code{\link[bnlearn]{increment.test.counter}} doesn't register, hence
# the hack with \code{\link[bnlearn]{test.counter}} to interact with it.

increment.test.counter_ <- function(i = 1){
  bnlearn::increment.test.counter(1)
  interact <- bnlearn::test.counter()
}



# Build key to avoid using match().

build_key <- function(nodes, cpp = FALSE){
  key <- seq(length(nodes))
  names(key) <- nodes
  if (cpp)
    key <- key - 1L
  return(key)
}



# Check if there is a directed path 
# from i to j in amat

has_path <- function(i, j, amat, nodes){
  
  if ((bool_i <- is.character(i)) ||
      (bool_j <- is.character(j))){
    
    key <- build_key(nodes = nodes, 
                     cpp = FALSE)
    if (bool_i)
      i <- key[i]
    if (bool_j)
      j <- key[j]
  }
  
  cpp_has_path(i = i-1L, j = j-1L, 
               amat = amat, nodes = nodes)
}