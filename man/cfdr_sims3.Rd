% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cfdr_sims3.R
\name{cfdr_sims3}
\alias{cfdr_sims3}
\title{Run Simulations}
\usage{
cfdr_sims3(x, pk.ht.funcs, type.sequence, seed = NULL, n.perms = 500,
  s0 = c(0, 0, 0), n.seg = c(), auto.min.length = c(), level = c(0.02,
  0.05, 0.1, 0.2), huber.maxit = 50, save.data = FALSE, file.name = NULL,
  data.only = FALSE, random.peak.loc = FALSE, min.peak.sep = 2)
}
\arguments{
\item{x}{Trait}

\item{pk.ht.funcs}{List of peak height functions. Each should take in a value of x
and return a a list of "height"=peak height (or vector of peak heights) and "assoc" as a
vector indicating if each peak is associated with x or not. "assoc" should be fixed for
all values of x.}

\item{type.sequence}{Vector of integers stating which peak height function to use.}

\item{seed}{Optional - set a seed.}

\item{n.perms}{Number of permutations to use to estimate lambda.}

\item{s0}{Vector of length 3 giving additional variance to add to Poisson, Huber and t-test statistics.}
}
\description{
Run simulations
}

