% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{bnrepository}
\alias{bnrepository}
\title{Download from Bayesian network repository}
\usage{
bnrepository(x)
}
\arguments{
\item{x}{character value indicating desired Bayesian network.}
}
\value{
A Bayesian network as an object of class \code{bn.fit}.
}
\description{
Read bn.fit objects from the \code{\link[bnlearn]{bnlearn}} Bayesian network repository directly from \href{www.bnlearn.com}{www.bnlearn.com} using the RDS links.
}
\examples{
## Read Bayesian network object 
true_bn <- bnrepository("child")
}
\author{
Jireh Huang (\email{jirehhuang@ucla.edu})
}
