debug <- TRUE

true_bn <- bnlearn::empty.graph(nodes = sprintf("X%s", seq_len(15)))

amat <- bnlearn::amat(true_bn)
amat[1, 2] <-
  amat[3, 2] <-
  amat[2, 4] <-
  amat[2, 5] <-
  amat[5, 6] <-
  amat[4, 8] <-
  amat[5, 9] <-
  amat[6, 12] <-
  amat[7, 11] <-
  amat[8, 11] <-
  amat[9, 12] <-
  amat[11, 10] <-
  amat[11, 13] <-
  amat[12, 14] <-
  amat[12, 15] <- 
  amat[14, 15] <- 1L
bnlearn::amat(true_bn) <- amat

set.seed(1)
x <- do.call(cbind, 
  lapply(list(1:5, c(7:8, 10:11, 13), c(6, 9, 12, 14:15)), function(x){
    
    z <- rbinom(n = 100, size = 1, prob = 0.5)
    data <- as.data.frame(sapply(x, function(y){
      
      z[seq_len(10)] <- sample(z[seq_len(10)])
      
      return(as.factor(z))
    }), stringsAsFactors = TRUE)
    names(data) <- sprintf("X%s", x)
    return(data)
  })
)
x <- x[,bnlearn::nodes(true_bn)]

result <- phsl::ppc(x = x, max.sx = 3, true_bn = true_bn, debug = debug)