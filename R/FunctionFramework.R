#'@importFrom Rcpp evalCpp
#'@useDynLib CDF.PSIdekick
#install.packages("Rcpp")
#library("Rcpp")
# This script is made to translate foreign (cpp) DP-CDF algorithm implentations into 
#      an R format compatible with the CDF Visualization and Evaluation Suite.
#      Authored by Dan Muise and George Kellaris, Harvard CRCS
#      dmuise@stanford.edu

##########################################################################


### here's the function platform with all of the inputs I normally use
### (works with cdf.cpp by Georgios Kellaris, available with this R script):

#' @title Create a DP-CDF by creating a K-degree noisy tree
#' @description This function creates a storage tree of degree K using gran and range,
#'  adds independent noise to each node proportional to epsilon, and then searches
#'  the tree to create a DP-CDF.
#'
#' @param eps Epsilon value for Differential privacy control
#' @param cdfstep The step sized used in outputting the approximate CDF; the values output are [min, min + cdfstep], [min, min + 2 * cdfstep], etc. 
#'   Setting cdfstep equal to 0 (default) will set cdfstep = granularity
#' @param data A vector of the data (single variable to compute CDFs from) 
#' @param range A vector length 2 containing user-specified min and max to truncate the universe to
#' @param gran The smallest unit of measurement in the data (one [year] for a list of ages)
#' @param K This sets the degree of the underlying  tree.
#' @param ... Optionally add additional parameters. 
#' @export
#' @return A list with 2 vectors: one is the y coordinates of the DP-CDF, the other is the
#'    abs values of the anlytically expected bounds on it at 95 percent probability.
#'
#' @examples 
#' functionH(eps = .01, cdfstep = .1, data = rexp(10000,.4), range= c(1,10), gran = .1, K= 2)

functionH <- function(eps, cdfstep, data, range, gran,K=2,...) {

  #sourceCpp("cdf.cpp")
  # here the function will redefine the range
          min <- range[1]
          max <- range[2]
    # here the function will redefine the bins and gran from sizing into division
          bin    <- floor((max - min)/gran +1)

          
             
            theCDFandCI            <- TreeCDF(eps,bin,K,"H",min,max,gran,data)
           GoodOutput           <- list()
  GoodOutput$theCDF      <-  theCDFandCI[1:(length(theCDFandCI)/2)]
  GoodOutput$bound       <-  theCDFandCI[ ((length(theCDFandCI)/2)+1):(length(theCDFandCI))] 

          return(GoodOutput) 
        }

####################################################

#' @title Create a monotonically increasing DP-CDF by creating a K-degree noisy tree
#' @description This function creates a storage tree of degree K using gran and range,
#'  adds independent noise to each node proportional to epsilon, and then searches
#'  the tree to create a DP-CDF. It then enforces monotonicity on the resuling 
#'  dpCDF.
#'
#' @param eps Epsilon value for Differential privacy control
#' @param cdfstep The step sized used in outputting the approximate CDF; the values output are [min, min + cdfstep], [min, min + 2 * cdfstep], etc. 
#' @param data A vector of the data (single variable to compute CDFs from) 
#' @param range A vector length 2 containing user-specified min and max to truncate the universe to
#' @param gran The smallest unit of measurement in the data (one [year] for a list of ages)
#' @param K This sets the degree of the underlying tree.
#' @param ... Optionally add additional parameters. 
#' @export
#' @return A list with 2 vectors: one is the y coordinates of the DP-CDF, the other is the
#'    abs values of the anlytically expected bounds for a similarly-constructed
#'    non-monotonized DP-CDF, at 95 percent probability.
#'
#' @examples 
#' functionHmono(eps = .01, cdfstep = .1, data = rexp(10000,.4), range= c(1,10), gran = .1, K= 2)
functionHmono <- function(eps, cdfstep, data, range, gran,K=2,...) {

  #sourceCpp("cdf.cpp")
  # here the function will redefine the range
          min <- range[1]
          max <- range[2]
    # here the function will redefine the bins and gran from sizing into division
          bin    <- floor((max - min)/gran +1)

  
            theCDFandCI            <- TreeCDF(eps,bin,K,"H",min,max,gran,data)
                GoodOutput           <- list()
  GoodOutput$theCDF      <- theCDFandCI[1:(length(theCDFandCI)/2)]
  GoodOutput$theCDF      <- smoothVector2(GoodOutput$theCDF)
  GoodOutput$bound       <- theCDFandCI[((length(theCDFandCI)/2)+1):(length(theCDFandCI))]  

          return(GoodOutput) 
        }

####################################################

#' @title Build dpCDFs through Histogram smoothing and minimized expected L2 per bin
#' @description The function seperates the epsilon value in two.
#' The first epsilon component is used to privately discover
#'  the best way to merge contiguous histogram bins in order to reduce the L2 error 
#'  due to the noise addition. It then applies the discovered bin merging to the original histogram,
#'  and outputs it by utilizing epsilon2. 
#'  Finally, it utilizes this output to compute and release the private CDF.
#'
#' @param eps Epsilon value for Differential privacy control
#' @param cdfstep The step sized used in outputting the approximate CDF; 
#' the values output are [min, min + cdfstep], [min, min + 2 * cdfstep], etc. 
#' @param data A vector of the data (single variable to compute CDFs from) 
#' @param range A vector length 2 containing user-specified min and max to truncate the universe to
#' @param gran The smallest unit of measurement in the data (one [year] for a list of ages)
#' @param K This sets the degree of the underlying tree
#' @param ... Optionally add additional parameters 
#' @export
#' @return A list with 2 vectors: one is the y coordinates of the DP-CDF, the other is the
#'    abs values of the anlytically expected bounds for a similarly-constructed
#'    non-monotonized DP-CDF made without merging of bins, at 95 percent probability.
#'
#' @examples 
#' functionS2(eps = .01, cdfstep = .1, data = rexp(10000,.4), range= c(1,10), gran = .1, K= 2)
functionS2 <- function(eps, cdfstep, data, range, gran,K=16,...) {

  #sourceCpp("cdf.cpp")
  # here the function will redefine the range
  min <- range[1]
  max <- range[2]
  # here the function will redefine the bins and gran from sizing into division
  bin    <- floor((max - min)/gran +1)
 K<-bin

  theCDFandCI            <- TreeCDF(eps,bin,K,"S2",min,max,gran,data)
  GoodOutput           <- list()
  GoodOutput$theCDF      <- theCDFandCI[1:(length(theCDFandCI)/2)]
  GoodOutput$bound       <- theCDFandCI[ ((length(theCDFandCI)/2)+1):(length(theCDFandCI))] 

  return(GoodOutput) 
}
#' @title Build dpCDFs through use of a noisy tree with bin merging.
#' @description The function first creates a k-ary aggregate tree on the histogram bins.
#'  It then utilizes epsilon1 in order to privately discover the best way to
#'  prune sub-trees in order to reduce the L2 error due to the noise addition.
#'  It then prunes the sub-trees of the original tree, and outputs it by utilizing epsilon2.
#'  Finally, it utilizes this output to compute and release the private CDF.
#'
#' @param eps Epsilon value for Differential privacy control
#' @param cdfstep The step sized used in outputting the approximate CDF; 
#'  the values output are [min, min + cdfstep], [min, min + 2 * cdfstep], etc. 
#' @param data A vector of the data (single variable to compute CDFs from) 
#' @param range A vector length 2 containing user-specified min and max to truncate the universe to
#' @param gran The smallest unit of measurement in the data (one [year] for a list of ages)
#' @param K This sets the degree of the underlying tree.
#' @param ... Optionally add additional parameters. 
#' @export
#' @return A list with 2 vectors: one is the y coordinates of the DP-CDF, the other is the
#'    abs values of the anlytically expected bounds for a similarly-constructed
#'    DP-CDF, at 95 percent probability made without merging.
#'
#' @examples 
#' functionSUB(eps = .01, cdfstep = .1, data = rexp(10000,.4), range= c(1,10), gran = .1, K= 2)
functionSUB <- function(eps, cdfstep, data, range, gran,K=2,...) {

  #sourceCpp("cdf.cpp")
  # here the function will redefine the range
  min <- range[1]
  max <- range[2]
  # here the function will redefine the bins and gran from sizing into division
  bin    <- floor((max - min)/gran +1)
  
  theCDFandCI            <- TreeCDF(eps,bin,K,"SUB",min,max,gran,data)
  GoodOutput           <- list()
  GoodOutput$theCDF      <- theCDFandCI[1:(length(theCDFandCI)/2)]
  GoodOutput$bound       <- theCDFandCI[ ((length(theCDFandCI)/2)+1):(length(theCDFandCI))]  

  return(GoodOutput) 
}