% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CDFtestingSuite.R
\name{diffat25}
\alias{diffat25}
\title{Determine the distance between CDFs at the .25 quantile.}
\usage{
diffat25(Y, est, ...)
}
\arguments{
\item{Y}{The vector output of a non-differentially private CDF 
computation (cumulative count bins)}

\item{est}{The vector output of a differentially private CDF
computation (cumulative count bins)}

\item{...}{Optionally add additional parameters. This is primarily used to allow automated
execution of varied diagnostic functions.}
}
\value{
The error at the .25 quantile
}
\description{
Find the error (between 0 and 1) introduced by DP-Noise at 
the .25 quantile.
}
\examples{
diffat25(c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1), c(.1,.2,.3,.3,.3,.3,.3,.3,.4,1))
}

