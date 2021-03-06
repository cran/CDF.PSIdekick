% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CDFtestingSuite.R
\name{getMean}
\alias{getMean}
\title{Calculate the private mean from the DP-CDF}
\usage{
getMean(est, range, gran, ...)
}
\arguments{
\item{est}{The vector output of a differentially private CDF 
computation (cumulative count bins)}

\item{range}{A vector length 2 containing user-specified min and max 
Note that the gran and range must be the same as used to make the DP-CDF!}

\item{gran}{The smallest unit of measurement in the data (one [year]
for a list of ages)}

\item{...}{Optionally add additional parameters. This is primarily used to allow automated
execution of varied diagnostic functions.}
}
\description{
Calculates the mean value from a CDF plot.
}
\examples{
getMean(c(.1,.2,.3,.3,.3,.3,.3,.3,.4,1), c(1,10), 1)
}

