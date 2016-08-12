
#'@importFrom Rcpp evalCpp
#'@useDynLib CDF.PSIdekick

#' @import grDevices
#' @import graphics
#' @import stats
#' @import utils

##########################################################
# DP_CDF Testing Suite: system for determining the utility
# and privacy of DP-CDF creation methods. 
#
# Scribe: Daniel Muise (2015)
# Harvard University

# Unless otherwise specified, the following processes were authored by
# Daniel Muise, Kobbi Nissim, Georgios Kellaris, Mark Bun, Victor Balcer  
##########################################################

## The UI for this testing suite is pasted and commented out at the end

# NULLing to rectify global variable binding
Abounds <- functionname <- AnalyticBounds <- smoothvector <- smoothVector <- NULL

### Functions Available
##########################################################
   # Diagnostic and minor functions defined
   #   MovetoRange        (val, range)
   #   getMaxError        (Y, est)
   #      MaxError_CDF    (Y, est)  
   #      MaxError_PDF    (Y, est) 
   #   findMaxError       (Y, est, databins)
   #       MaxErrorAt_CDF (Y, est, databins)
   #       MaxErrorAt_PDF (Y, est, databins)
   #   diffatQuantile     (Y, est, quantile = 0.50)
   #       diffat25       (Y, est)
   #       diffatMed      (Y, est)
   #       diffat75       (Y, est)
   #   horzdiffatQuantile (Y, est, range, gran, quantile)
   #       horzdiffat25   (Y, est, range, gran)
   #       horzdiffatMed  (Y, est, range, gran)
   #       horzdiffat75   (Y, est, range, gran)
   #   QuantileFromCDF    (est, range, gran, quantile)
   #       Medians        (est, range, gran)
   #   smoothVector       (vector)
   #   MSEanalytic         (eps, range, gran, data)
   #       nodes          (height, k, l) (helper to MSEanalytic)
   #   MAE                (Y, est)
   #   MSE                (Y, est)
   #   SDempiric          (Y, est)
   #   getMean            (est, range, gran)
   #   DerivDiff          (Y, est)
   #   smoothVector       (vector) 
   #      partialMeans    (vector, i, j) (helper to smoothVector)
   #   smoothVector2      (cdf)
   #   Abbrev             (value)
   ############################################


##########################################################
# The following segment defines internal helper functions
#########################################################################
##################################################################

#' @title Clamp a value to a specified range.
#'
#' @description  Returns a vector of elements clamped to the specified minimum
#'    and maximum
#' @param val A value to clamp.
#' @param range A vector of length 2 in the form \code{c(min, max)}
#' @return A single value that is either unchanged or clamped upward to minimum 
#'    or clamped downward to the maximum
#' @export
#' @examples 
#' MovetoRange(11, c(1,10))

MovetoRange <- function(val, range) {
 
  if(val < range[1]) {
    return(range[1])
  }
  else if(val > range[2]) {
    return(range[2])
  }
  else
    return(val)
  }
#################################################
#' @title Determine an approximate CDF's maximum error.
#'
#' @description Find the maximum direct error between a non-private CDF and a 
#'    DP approximation of that CDF.
#'
#' @param Y The vector output of a non-differentially private CDF 
#'    computation (cumulative count bins)
#' @param est The vector output of a differentially private CDF 
#'    computation (cumulative count bins) 
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#' @return A single value, the largest absolute vertical difference between 
#'    parallel observations in the private- and true-CDF vectors.
#' @export
#' @examples 
#' getMaxError(c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1), c(.1,.2,.3,.3,.3,.3,.3,.3,.4,1))
getMaxError <- function(Y, est,...) {
  diff <- (Y-est) #store the point-by-point deviations
  return <-(max(abs(diff))) #output the max
}

#################################################
#' @title Determine an approximate CDF's maximum error.
#'
#' @description Find the maximum direct error between a non-private CDF and a
#'    DP approximation of that CDF.
#' @param Y The vector output of a non-differentially private CDF 
#'    computation (cumulative count bins)
#' @param est The vector output of a differentially private CDF
#'    computation (cumulative count bins) 
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#' @return A single value, the largest absolute vertical difference between 
#'    parallel observations in the private- and true-CDF vectors.
#' @export
#' @examples 
#' MaxError_CDF(c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1), c(.1,.2,.3,.3,.3,.3,.3,.3,.4,1))
MaxError_CDF <- function(Y,est,...){
  return(getMaxError(Y, est,...))
}

#################################################
#' @title Determine an approximate PDF's maximum error.
#'
#' @description Find the maximum direct error between a non-private PDF and a 
#'    DP approximation of that PDF.
#'
#' @param Y The vector output of a non-differentially private PDF 
#'    computation (heights of bins)
#' @param est The vector output of a differentially private PDF
#'    computation (heights of bins) 
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#' @return A single value, the largest absolute vertical difference between 
#'    parallel observations in the private- and true-PDF vectors.
#' @export
#' @examples 
#' MaxError_PDF(c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1), c(.1,.2,.3,.3,.3,.3,.3,.3,.4,1))
MaxError_PDF <- function(Y,est,...){
  return(getMaxError(Y, est,...))
}

#################################################
#' @title Locate where the maximum error occurs between two CDFs
#'
#' @description Find the location of the maximum direct error between a 
#'     non-private CDF and a DP approximation of that CDF.
#'
#' @param Y The vector output of a non-differentially private CDF 
#'     computation (cumulative count bins)
#' @param est The vector output of a differentially private CDF
#'     computation (cumulative count bins) 
#' @param range A vector length 2 containing user-specified min and max to
#'     truncate the universe to
#' @param gran The smallest unit of measurement in the data (one [year] for 
#'     a list of ages)
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#' @return A single value, the value at which the  largest absolute vertical
#'     difference between
#'   parallel observations in the private- and true-CDF vectors occurs.
#' @export
#' @examples 
#' findMaxError(c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1), c(.1,.2,.3,.3,.3,.3,.3,.3,.4,1), c(1,10),1)
findMaxError <- function(Y, est, range, gran,...) {
  databins <- seq(from=range[1], to=range[2], by=gran)
  diff <- (Y-est) #store the point-by-point deviations
  maxx <- (max(abs(diff))) #output the max
  for (i in 1: length(diff)){
    if (abs(diff[i]) == maxx){
      MaxErrorAt <- databins[i]
    }
  }
  return <-(MaxErrorAt) #output the max error location
}

#################################################
# define specific versions to help with syntax in CDFtestTrack

#' @title Locate where the maximum error occurs between two CDFs
#' @description Find the location of the maximum direct error between a 
#'    non-private CDF and a DP approximation of that CDF.
#'
#' @param Y The vector output of a non-differentially private CDF 
#'    computation (cumulative count bins)
#' @param est The vector output of a differentially private CDF
#'    computation (cumulative count bins) 
#' @param range A vector length 2 containing user-specified min and max to 
#'    truncate the universe to
#' @param gran The smallest unit of measurement in the data (one [year] 
#'    for a list of ages)
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#' @return A single value, the value at which the  largest absolute vertical
#'    difference between
#'   parallel observations in the private- and true-CDF vectors occurs.
#' @export
#' @examples 
#' MaxErrorAt_CDF(c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1), c(.1,.2,.3,.3,.3,.3,.3,.3,.4,1),
#'     range= c(1,10), gran =1)
MaxErrorAt_CDF <- function(Y, est,range,gran,...){
  return(findMaxError(Y, est, range,gran,...))
}
#################################################
#' @title Locate where the maximum error occurs between two PDFs
#' @description Find the location of the maximum direct error between a
#'    non-private PDF and a DP approximation of that PDF.
#'
#' @param Y The vector output of a non-differentially private PDF
#'    computation (values within bins)
#' @param est The vector output of a differentially private PDF 
#'    computation (values within bins) 
#' @param range A vector length 2 containing user-specified min and max to
#'    truncate the universe to
#' @param gran The smallest unit of measurement in the data (one [year]
#'    for a list of ages)
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#' @return A single value, the value at which the  largest absolute vertical
#'    difference between
#'   parallel observations in the private- and true-PDF vectors occurs.
#' @export
#' @examples 
#' MaxErrorAt_PDF(c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1), c(.1,.2,.3,.3,.3,.3,.3,.3,.4,1),
#'    range= c(1,10), gran =1)
MaxErrorAt_PDF <- function(Y, est, range, gran,...){
  return(findMaxError(Y, est, range, gran,...))
}

#################################################
#' @title Determine the distance between CDFs at key quantiles.
#' @description Find the error (between 0 and 1) introduced by DP-Noise at a
#'    given quantile in the CDF
#'
#' @param Y The vector output of a non-differentially private CDF
#'    computation (cumulative count bins).
#' @param est The vector output of a differentially private CDF 
#'    computation (cumulative count bins).
#' @param quantile A quantile value between 0 and 1, defaults to 0.5
#'    for the median.
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#' @return The error at the quantile specified by \code{quantile}
#' @export
#' @examples 
#' diffatQuantile(c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1),
#'     c(.1,.2,.3,.3,.3,.3,.3,.3,.4,1), .05)

#################################################

diffatQuantile <- function(Y, est, quantile=0.5,...) {
xY <- c()
xest <- c()
#for storing the real quantile value and the corresponding estimate value
tempvector <- c(0,0)  
#locate the quantile in the CDF vector using the index, and store that value
# the median of quantile values must be found in the case of non-monotonic CDFs
if (min(Y-quantile)<0){
  for (i in 2:(length(Y)-1)) {  
    if ((Y[i]-quantile) >=0)  {
      if ((Y[i-1] -quantile) <0 ){
      xY[length(xY)+1] <- Y[i] 
      tempvector[1] <- median(xY)
      }
    }
  }
}
else{
  tempvector[1]<-Y[1]
}
  # do the same for the private-CDF vector
  if (min(est-quantile)<0){
  for (i in 2:(length(est)-1)) {  
    if ((est[i]-quantile) >=0)  {
      if ((est[i-1] -quantile) <0 ){
     xest[length(xest)+1] <- est[i] 
     tempvector[2] <- median(xest) 
      }
    }
  }
}
else{
  tempvector[2] <-est[1]
}
#find the absolute difference
diff <- abs(tempvector[2] - tempvector[1])
return(diff)
}
#################################################
####################
#' @title Determine the distance between CDFs at the .25 quantile.
#'
#' @description Find the error (between 0 and 1) introduced by DP-Noise at 
#' the .25 quantile.
#'
#' @param Y The vector output of a non-differentially private CDF 
#'     computation (cumulative count bins)
#' @param est The vector output of a differentially private CDF
#'     computation (cumulative count bins) 
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#' @return The error at the .25 quantile
#'@export
#' @examples 
#' diffat25(c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1), c(.1,.2,.3,.3,.3,.3,.3,.3,.4,1))
# define specific versions to help with syntax in CDFtestTrack
diffat25 <- function(Y, est,...){
return(diffatQuantile(Y, est, .25,...)) 
}

#################################################
#' @title Determine the distance between CDFs at the median.
#'
#' @description Find the error (between 0 and 1) introduced by DP-Noise at 
#'    the median
#'
#' @param Y The vector output of a non-differentially private CDF 
#'    computation (cumulative count bins)
#' @param est The vector output of a differentially private CDF 
#'    computation (cumulative count bins) 
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#' @return The error at the .5 quantile
#'@export
#' @examples 
#' diffatMedian(c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1), c(.1,.2,.3,.3,.3,.3,.3,.3,.4,1))
diffatMedian <- function(Y, est,...){
return(diffatQuantile(Y, est, .50,...)) 
}

#################################################
#' @title Determine the distance between CDFs at the .75 quantile.
#'
#' @description Find the error (between 0 and 1) introduced by DP-Noise 
#'   at the .75 quantile.
#'
#' @param Y The vector output of a non-differentially private 
#'     CDF computation (cumulative count bins)
#' @param est The vector output of a differentially private
#'     CDF computation (cumulative count bins) 
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#' @return The error at the .75 quantile
#'@export
#' @examples 
#' diffat75(c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1), c(.1,.2,.3,.3,.3,.3,.3,.3,.4,1))
diffat75 <- function(Y, est,...){
return(diffatQuantile(Y, est, .75,...)) 
}

###########################################
#' @title Determine the distance between the quantile values 
#'   returned by two CDFs.
#' @description Find the distance between the quantile value and that returned
#'   by the dpCDF at a given quantile.
#'
#' @param Y The vector output of a non-differentially private CDF
#'   computation (cumulative count bins)
#' @param est The vector output of a differentially private CDF 
#'   computation (cumulative count bins) 
#' @param range A vector length 2 containing user-specified min and max 
#'   to truncate the universe to
#' @param gran The smallest unit of measurement in the data (one [year]
#'   for a list of ages)
#' @param quantile A quantile value between 0 and 1, defaults to
#'   0.5 for the median
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#' @return The horizontal error at the quantile specified by \code{quantile}
#'@export
#' @examples 
#' diffatQuantile(c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1), 
#'    c(.1,.2,.3,.3,.3,.3,.3,.3,.4,1),c(1,10), 1, .05)
horzdiffatQuantile <- function(Y, est, range, gran, quantile,...) {
databins <- seq(from=range[1], to=range[2], by=gran)

xY <- c()
xest <- c()
tempvector <- c(0,0) #for storing the real and estimated quantile values 
#locate the quantile in the true-CDF vector using the index, then store 
#  the corresponding databins obeservation.
# the median of quantile values must be found in the case of non-monotonic CDFs.
if (min(Y-quantile)<0){
for (i in 2:(length(Y)-1)) {  
  if ((Y[i]-quantile) >=0)  {
    if ((Y[i-1] -quantile) <0 ){
     xY[length(xY)+1] <- databins[i] 
     tempvector[1] <- median(xY)  
   }
  }
 }
}
else{
  tempvector[1]<-databins[1]
}
# do the same for the private-CDF
if (min(est-quantile)<0){
for (i in 2:(length(est)-1)) {  
  if ((est[i]-quantile) >=0)  {
    if ((est[i-1] -quantile) <0 ){
     xest[length(xest)+1] <- databins[i] 
     tempvector[2] <- median(xest)  
   }
  } 
 }
}
else{
  tempvector[2]<-databins[1]
}
# store the absolute difference
diff <- abs(tempvector[2] - tempvector[1])
return(diff)
}
#################################################
########################
# define specific versions to help with syntax in CDFtestTrack
#' @title Determine the distance between the .25 quantile values returned 
#'   by two CDFs.
#' @description Find the distance between the .25 quantile value and that
#'   returned by the dpCDF.
#'
#' @param Y The vector output of a non-differentially private CDF
#'   computation (cumulative count bins)
#' @param est The vector output of a differentially private CDF
#'   computation (cumulative count bins) 
#' @param range A vector length 2 containing user-specified min and max to 
#'   truncate the universe to
#' @param gran The smallest unit of measurement in the data (one [year]
#'   for a list of ages)
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#' @return The horizontal error at the .25 quantile
#' @export
#' @examples 
#' horzdiffat25(c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1), 
#'  c(.1,.2,.3,.3,.3,.3,.3,.3,.4,1),c(1,10), 1)
#################################################
horzdiffat25 <- function(Y, est, range, gran,...){
return(horzdiffatQuantile(Y, est, range, gran, .25,...)) 
}
#' @title Determine the distance between the median values returned by two CDFs.
#' @description Find the distance between the median value and that returned
#'    by the DP CDF.
#'
#' @param Y The vector output of a non-differentially private CDF
#'    computation (cumulative count bins)
#' @param est The vector output of a differentially private CDF
#'    computation (cumulative count bins) 
#' @param range A vector length 2 containing user-specified min and max to
#'    truncate the universe to
#' @param gran The smallest unit of measurement in the data (one [year] 
#'    for a list of ages)
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#' @return The horizontal error at the median
#'@export
#' @examples 
#' horzdiffatMed(c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1), 
#'    c(.1,.2,.3,.3,.3,.3,.3,.3,.4,1),c(1,10), 1)
#################################################
horzdiffatMed <- function(Y, est, range, gran,...){
return(horzdiffatQuantile(Y, est, range, gran, .50,...)) 
}
#' @title Determine the distance between the .75 quantile values
#'    returned by two CDFs.
#' @description Find the distance between the .75 quantile value
#'    and that returned by the DP CDF.
#'
#' @param Y The vector output of a non-differentially private CDF
#'   computation (cumulative count bins)
#' @param est The vector output of a differentially private CDF 
#'   computation (cumulative count bins) 
#' @param range A vector length 2 containing user-specified min
#'   and max to truncate the universe to
#' @param gran The smallest unit of measurement in the data (one [year]
#'   for a list of ages)
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#' @return The horizontal error at the .75 quantile
#'@export
#' @examples 
#' horzdiffat75(c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1),
#'             c(.1,.2,.3,.3,.3,.3,.3,.3,.4,1),c(1,10), 1)
horzdiffat75 <- function(Y, est, range, gran,...){
return(horzdiffatQuantile(Y, est, range, gran, .75,...)) 
}
###########################################

#' @title Retrieve a private quantile estimate from the dpCDF
#' @description Determines a quantile value from a CDF vector.
#'
#' @param est The vector output of a differentially private CDF 
#'    computation (cumulative count bins) 
#' @param range A vector length 2 containing user-specified min and max to 
#'    truncate the universe to
#' @param gran The smallest unit of measurement in the data (one [year] 
#'    for a list of ages),
#'    the Domain (ie gran and range) should be identical to those used to 
#'    create the CDF!
#' @param quantile the quantile score in question (for testing the median,
#'    use quantile = 0.5) 
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#' @return A quantile value obtained from a (differentially private) CDF vector,
#'    not using any extra privacy budget
#'@export
#' @examples 
#' QuantileFromCDF(c(.1,.2,.3,.3,.3,.3,.3,.3,.4,1),c(1,10), 1, .05)
QuantileFromCDF <- function(est, range, gran, quantile,...){
  Domain <- seq(from=range[1], to=range[2], by=gran)
  outputvalue <- 0
  DifferenceVector <- est-quantile

  for (z in 1: length(DifferenceVector)){
    if (DifferenceVector[z] >0) {
      DifferenceVector[z] <- -10
    }
  }
  maxDiff <- max(DifferenceVector)

  for (x in (length(DifferenceVector)):1){
    if (DifferenceVector[x] == maxDiff){
      outputvalue <- Domain[x]
    }
  }
  return(outputvalue)
}
#################################################
###
#' @title Retrieve a median estimate from the dpCDF
#' @description Determines a median value from a CDF vector.
#'
#' @param est The vector output of a differentially private CDF 
#'    computation (cumulative count bins) 
#' @param range A vector length 2 containing user-specified min and max to 
#'    truncate the universe to
#' @param gran The smallest unit of measurement in the data (one [year]
#'    for a list of ages),
#'    the Domain (ie gran and range) should be identical to those used to 
#'    create the CDF!
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#' @return A vector of medians obtained from a (differentially private) 
#'    CDF vector, not using any extra privacy budget,
#'    there may be more than one due to random noise causing the DPCDF doubling
#'    back over the .5 probablity latitude
#'@export
#' @examples 
#' Medians(c(.1,.2,.3,.3,.3,.3,.3,.3,.4,1),c(1,10), 1)
Medians <- function(est, range, gran,...){
  return(QuantileFromCDF(est, range, gran, 0.5,...))
}

################################################
#'@title Node parser.
#'@description Runs through tree nodes (assists MSE analytic)
#'
#' @param height The height of the tree
#' @param k The tree degree
#' @param l The leaf length
#'
#' @return A nodesum containing information for MSEanalytic
#'@export
#' @examples 
#' nodes(10,4,2)
nodes <- function (height, k, l){
  #(FOR MSEanalytic)
  v<-l
  j <- k^height
  sum <- 0
  while (j >= k){
    sum <- sum + v %/% j
    v <- v %% j
    height  <- height - 1
    j <- k^height
  }
  sum <- sum + v %% k
  
  # var[i]=2.12*var[i]/(double)n;

  return(sum)
}
##################
#' @title Determine the expected MSE of a simple DPCDF from its parameters.
#' @description Generates the analytically expected 
#'    Mean Squared Error of a dpCDF.
#'    introduced by random noise, SUPPOSING that the DP-CDF is through the use
#'    of a noisy binary tree.
#'
#' @param eps Epsilon value for differential privacy control
#' @param range A vector length 2 containing user-specified min and max to
#'    truncate the universe to
#' @param gran The smallest unit of measurement in the data (one [year]
#'    for a list of ages)
#' @param data The vector of data from which the DP CDF was/is computed
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#' @return The MSE guaranteed by the given parameter combination assuming 
#'    it's built from 
#'    the min and max inward from a DP-Histogram, with 95% probability
#'@export
#' @examples 
#' MSEanalytic(.01, c(1,10),1, rexp(10000,.4))
MSEanalytic <- function(eps, range, gran, data,...) {

 k <- 2
  universe_size <- floor((range[2] - range[1]) / gran) + 1
  height <- ceiling(log(universe_size)/log(k))

  avg <- 0
  for (i in 1:(universe_size-1) ){
    num1 <- nodes(height, k, i)
    num2 <- nodes(height, k, universe_size-i)
    avg <- avg + (num1*(num2/(num1+num2))*(num2/(num1+num2))*8.0/(eps^2)*(height^2)+num2*(num1/(num1+num2))*(num1/(num1+num2))*8.0/(eps^2)*(height^2))/(universe_size);
  }

  return ((avg)/length(data))
}

############################################################
#' @title Calculate the empirical L2norm between two CDFs.
#' @description Calculates the L2 (squared error) area between the 
#'    non-private CDF and the dpCDF
#' @param Y The vector output of a non-differentially private CDF 
#'    computation (cumulative count bins)
#' @param est The vector output of a differentially private CDF
#'    computation (cumulative count bins) 
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#' @return The empirical L2 norm
#'@export
#' @examples 
#' L2empiric(c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1), c(.1,.2,.3,.3,.3,.3,.3,.3,.4,1))
L2empiric <- function(Y, est, ...){
   sqvect <- c()
   sqvect <- (Y-est)^2
  
L2 <- sqrt(sum(sqvect))
return(L2)
}

################################
#' @title Calculate the area between two CDFs.
#' @description Calculates the L1 (distance error) area between the non-private 
#'    CDF and the dpCDF
#' @param Y The vector output of a non-differentially private CDF
#'    computation (cumulative count bins)
#' @param est The vector output of a differentially private CDF
#'    computation (cumulative count bins) 
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#' @return The empirical L1 norm
#'@export
#' @examples 
#' L1empiric(c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1), c(.1,.2,.3,.3,.3,.3,.3,.3,.4,1))
L1empiric <- function(Y,est,...){

  diffs <- c()    
  diffs <- abs(Y-est)
  L1 <- sum(diffs)
return(L1)
}

################################
#' @title Calculate the MAE of a dpCDF relative to that of the non-private CDF.
#' @description Calculates the Mean Absolute Error area between the
#'    non-private CDF and the dpCDF
#' @param Y The vector output of a non-differentially private 
#'    CDF computation (cumulative count bins)
#' @param est The vector output of a differentially private
#'    CDF computation (cumulative count bins) 
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#' @return The MAE
#'@export
#' @examples 
#' MAE(c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1), c(.1,.2,.3,.3,.3,.3,.3,.3,.4,1))
MAE <- function(Y,est,...){

  diffs <- c()
  diffs <- abs(Y-est)
  L1 <- sum(diffs)
  MAE <- L1/length(Y)

return(MAE)
}
################################
#' @title Calculate the MSE of a DP-CDF relative to the non-private CDF.
#' @description Calculates the Mean Squared Error area between the non-private 
#'     CDF and the DP-CDF
#' @param Y The vector output of a non-differentially private CDF
#'     computation (cumulative count bins)
#' @param est The vector output of a differentially private CDF
#'     computation (cumulative count bins)
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#' @return The MSE
#' 
#' @export
#' @examples 
#' MSE(c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1), c(.1,.2,.3,.3,.3,.3,.3,.3,.4,1))
MSE <- function(Y,est,...){
  sqvect <- c()    
  sqvect <- (Y-est)^2 
  L2 <- sqrt(sum(sqvect))
  MSE <- (L2^2)/length(Y)
return(MSE)
}

#######################################
#' @title Calculate the std. dev. on a DPCDF.
#' @description Calculates the standard deviation across bins between the 
#'    non-private CDF and the DP-CDF
#' @param Y The vector output of a non-differentially private CDF
#'    computation (cumulative count bins)
#' @param est The vector output of a differentially private CDF
#'    computation (cumulative count bins) 
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#' @return The standard deviation
#'@export
#' @examples 
#' SDempiric(c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1), c(.1,.2,.3,.3,.3,.3,.3,.3,.4,1))

SDempiric <- function(Y, est,...){
  sqvect <- c()
  sqvect <- (Y - est)^2
  L2 <- sqrt(sum(sqvect))
  SD <- L2/(sqrt(length(Y)))
return(SD)
}

################################
#' @title Calculate the private mean from the DP-CDF
#' @description Calculates the mean value from a CDF plot.
#'
#' @param est The vector output of a differentially private CDF 
#'    computation (cumulative count bins) 
#' @param gran The smallest unit of measurement in the data (one [year]
#'    for a list of ages)
#' @param range A vector length 2 containing user-specified min and max 
#'    Note that the gran and range must be the same as used to make the DP-CDF!
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#'@export
#' @examples 
#' getMean(c(.1,.2,.3,.3,.3,.3,.3,.3,.4,1), c(1,10), 1)
  getMean <- function(est, range, gran,...){
  Domain <- seq(from=range[1], to=range[2], by=gran)
  PDFvector <- c(est[1])
    for (x in (length(est)):2){
    PDFvector[x] <- est[x] - est[x-1]
    }
  weights <-c()
    for (x in 1:length(Domain)){
      weights[x] <- ( PDFvector[x] * Domain[x] )
    }

  mean <- sum( weights )
return(mean)
}

#####################################################

############################################################
# use it on PDF output to save time

#' @title Error in mean from CDF
#' @description Calculate difference between the private mean and the original mean (from CDFs)
#'
#' @param Y The vector output of a non-differentially private CDF
#'    computation (cumulative count bins)
#' @param est The vector output of a differentially private CDF 
#'    computation (cumulative count bins) 
#' @param gran The smallest unit of measurement in the data (one [year]
#'    for a list of ages)
#' @param range A vector length 2 containing user-specified min and max 
#'    Note that the gran and range must be the same as used to make the DP-CDF!
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#' @return A single difference value
#'
MeanDiffpdf <- function(Y, est, range, gran){

Domain <- seq(from=range[1], to=range[2], by=gran)

trueweights <-c()
dpweights   <-c()

  trueweights<- ( Y * Domain )
    dpweights <- ( est * Domain )

truemean <- sum( trueweights )
dpmean <-sum(dpweights)
return(truemean-dpmean)
}
#############################################
#' @title Error in Skewness from CDF (under development)
#' @description Calculate difference between the private Skewness and the original Skewness (from CDFs)
#'
#' @param Y The vector output of a non-differentially private CDF
#'    computation (cumulative count bins)
#' @param est The vector output of a differentially private CDF 
#'    computation (cumulative count bins) 
#' @param gran The smallest unit of measurement in the data (one [year]
#'    for a list of ages)
#' @param range A vector length 2 containing user-specified min and max 
#'    Note that the gran and range must be the same as used to make the DP-CDF!
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#' @return A single difference value
#'
SkewDiffpdf <- function(Y, est, range, gran){

  #determines the skewness of a univariate distribution taken from a DP-CDF vector

Domain <- seq(from=range[1], to=range[2], by=gran)

trueweights <-c()
dpweights   <-c()

  trueweights <- ( Y * Domain)
    dpweights <- ( est * Domain )

truemean <- sum( trueweights )
dpmean <-sum(dpweights)


truedistances <- c()
  dpdistance  <- c()  

  truedistances <- ((truemean - Domain)^2) * (Y)
    dpdistances <- ((  dpmean - Domain)^2) * (est)


truestddev <- sqrt(sum(truedistances))
dpstddev <- sqrt(sum(dpdistances))

trueskewness <- mean(((trueweights-truemean)/ truestddev)^3)
  dpskewness <- mean(((  dpweights-  dpmean)/   dpstddev)^3)
return(trueskewness-dpskewness)
}

#############################################

#' @title Error in Variance from CDF
#' @description Calculate difference between the private Variance and the original Variance (from CDFs)
#'
#' @param Y The vector output of a non-differentially private CDF
#'    computation (cumulative count bins)
#' @param est The vector output of a differentially private CDF 
#'    computation (cumulative count bins) 
#' @param gran The smallest unit of measurement in the data (one [year]
#'    for a list of ages)
#' @param range A vector length 2 containing user-specified min and max 
#'    Note that the gran and range must be the same as used to make the DP-CDF!
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#' @return A single difference value
#'
VarDiffpdf <- function(Y, est, range, gran){

Domain <- seq(from=range[1], to=range[2], by=gran)

trueweights <-c()
dpweights   <-c()

  trueweights <- ( Y * Domain )
    dpweights <- ( est * Domain )


truemean <- sum( trueweights )
dpmean <-sum(dpweights)


truedistances <- c()
  dpdistance  <- c()  

  truedistances <- ((truemean - Domain)^2) * (Y)
    dpdistances <- ((  dpmean - Domain)^2) * (est)



truevariance <- sum(truedistances)
dpvariance <- sum(dpdistances)
return(truevariance-dpvariance)
}

################################################
#' @title Error in Mode from CDF
#' @description Calculate difference between the private Mode and the original Mode (from CDFs)
#'
#' @param Y The vector output of a non-differentially private CDF
#'    computation (cumulative count bins)
#' @param est The vector output of a differentially private CDF 
#'    computation (cumulative count bins) 
#' @param gran The smallest unit of measurement in the data (one [year]
#'    for a list of ages)
#' @param range A vector length 2 containing user-specified min and max 
#'    Note that the gran and range must be the same as used to make the DP-CDF!
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#' @return A single difference value
#'

ModeDiffpdf <- function(Y,est, range, gran,...){
  Domain <- seq(from=range[1], to=range[2], by=gran)
  
  PDF1k <- Y*10000000000000000
  PDFblank <- PDF1k - (max(PDF1k)-1)
  PDFblank <- (PDFblank + abs(PDFblank)) / 2
  valuevect<- PDFblank * Domain
  
  truemode <- (sum(valuevect)/ sum(PDFblank))
  
  PDF1k <- est*1000000000000000
  PDFblank <- PDF1k - (max(PDF1k)-1)
  
  PDFblank <- (PDFblank + abs(PDFblank)) / 2
  valuevect<- PDFblank * Domain
  
  dpmode <- (sum(valuevect)/ sum(PDFblank))
  
  
  return(truemode-dpmode)
}


#############################################
#' @title Error in Standard Deviation from CDF 
#' @description Calculate difference between the private Standard Deviation and the original Standard Deviation (from CDFs)
#'
#' @param Y The vector output of a non-differentially private CDF
#'    computation (cumulative count bins)
#' @param est The vector output of a differentially private CDF 
#'    computation (cumulative count bins) 
#' @param gran The smallest unit of measurement in the data (one [year]
#'    for a list of ages)
#' @param range A vector length 2 containing user-specified min and max 
#'    Note that the gran and range must be the same as used to make the DP-CDF!
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#' @return A single difference value
#'
StdDiffpdf <- function(Y, est, range, gran){

Domain <- seq(from=range[1], to=range[2], by=gran)

trueweights <-c()
dpweights   <-c()

  trueweights <- ( Y * Domain )
    dpweights <- ( est * Domain )


truemean <- sum( trueweights )
dpmean <-sum(dpweights)


truedistances <- c()
  dpdistance  <- c()  

  truedistances <- ((truemean - Domain)^2) * (Y)
    dpdistances <- ((  dpmean - Domain)^2) * (est)



trueSTD <- sqrt(sum(truedistances))
dpSTD <- sqrt(sum(dpdistances))
return(trueSTD-dpSTD)
}

################################################
#' @title Error in Kurtosis from CDF (under development)
#' @description Calculate difference between the private Kurtosis and the original Kurtosis (from CDFs)
#'
#' @param Y The vector output of a non-differentially private CDF
#'    computation (cumulative count bins)
#' @param est The vector output of a differentially private CDF 
#'    computation (cumulative count bins) 
#' @param gran The smallest unit of measurement in the data (one [year]
#'    for a list of ages)
#' @param range A vector length 2 containing user-specified min and max
#'    Note that the gran and range must be the same as used to make the DP-CDF!
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#' @return A single difference value
#'
KurtDiffpdf <- function(Y, est, gran, range){

Domain <- seq(from=range[1], to=range[2], by=gran)

trueweights <-c()
dpweights   <-c()

  trueweights <- ( Y * Domain )
    dpweights <- ( est * Domain )

truemean <- sum( trueweights )
dpmean <-sum(dpweights)

truedistances <- c()
  dpdistance  <- c()  

  truedistances <- ((truemean - Domain)^2) * (Y)
    dpdistances <- ((  dpmean - Domain)^2) * (est)

truevariance <- sum(truedistances)
truemoment4 <- (mean(truedistances))
truesigma4  <-  truevariance^2
#normalize the 4th moment
truekurtosis <- truemoment4/truesigma4

dpvariance <- sum(dpdistances)
dpmoment4 <- (mean(dpdistances))
dpsigma4  <-  dpvariance^2
#normalize the 4th moment
dpkurtosis <- dpmoment4/dpsigma4


return(truekurtosis-dpkurtosis)

}

###############################################
#' @title Determine how well a single DPCDF matches the shape of its data.
#' @description Calculates a score for how much the DP-CDF's slope varies from
#'    the true CDF's slope at various resolutions.
#'
#' @param Y The vector output of a non-differentially private CDF
#'    computation (cumulative count bins)
#' @param est The vector output of a differentially private CDF
#'    computation (cumulative count bins) 
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#' @return A single so-called derivative score; lower scores suggest
#'    better performance 
#'@export
#' @examples 
#' DerivDiff(c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1),c(.1,.2,.3,.3,.3,.3,.3,.3,.4,1))
DerivDiff <- function(Y, est,...) {

workingY <- Y
workingE <- est

crunch <- c()
Ypdf <- c()
Epdf <- c()
Ypdf[1]<-workingY[1]
Epdf[1]<-workingE[1]
  for(v in length(workingY):2){
                Ypdf[v] <- workingY[v]-  workingY[v-1]
                Epdf[v] <- workingE[v]-  workingE[v-1]
  }
  diff      <- abs(Ypdf - Epdf)
  L1area    <- sum(diff)
  MAE <- L1area/length(Y)

  crunch[length(crunch)+1] <-MAE


for (k in 1:(floor(log2(length(Y))))){
Yx <- c()
Ex <- c()

  for (i in 1:floor(length(workingY)/2)){
    Yx[i] <- (workingY[2*i]+workingY[(2*i)-1])/2
  }

  if ((floor(length(Yx))*2) != length(workingY) ){
    remainder <- c(workingY[length(workingY)-(2*length(Yx))]:workingY[length(workingY)])
    Yx[length(Yx)+1] <- mean(remainder)
  }
####
  for (i in 1:floor(length(workingE)/2)){
    Ex[i] <- (workingE[2*i]+workingE[(2*i)-1])/2
  }
  if ((floor(length(Ex))*2) != length(workingE) ){
    remainder <- c(workingE[length(workingE)-(2*length(Ex))]:workingE[length(workingE)])
    Ex[length(Ex)+1] <- mean(remainder)
  }
    workingE <- Ex
    workingY <- Yx
Ypdf <-c()
Epdf <-c()
Ypdf[1]<-workingY[1]
Epdf[1]<-workingE[1]
if(length(workingY) >1){
  for(v in length(workingY):2){
                Ypdf[v] <- workingY[v]-  workingY[v-1]
                Epdf[v] <- workingE[v]-  workingE[v-1]
  }
  diff      <- abs(Ypdf - Epdf)
  L1area    <- sum(diff)
  MAE <- L1area/length(workingY)
   crunch[length(crunch)+1] <-MAE
}
}
output <- mean(crunch)

return(output)
}
########################################################
# Smooth <- function() {
#   .Call('dpCDFtester_Smooth', PACKAGE = 'dpCDFtester')
# }
############################################################
#' @title Enforce monotnocity on a vector.

#' @description Forces DP-CDFs into the nearest monotonic vector
#'    (by euclidean distance minimization).
#'
#' @param cdf The vector output of a differentially private CDF 
#'    computation (cumulative count bins) 
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#' @return A single monotonically increasing vector which is the post-processed
#'    DP-CDF's Y coordinates
#'@export
#' @examples 
#' smoothVector2(c(.1,.2,.3,.2,.3,.3,.3,.3,1))
smoothVector2 <- function(cdf){
  #Rcpp::sourceCpp("smooth.cpp")
  returnVector  <- Smooth(cdf)
return(returnVector)
}

##########################################

#' @title Tranforms long numbers into short strings.
#' @description Abbreviates long numeric values into visually shorter strings
#'
#' @param value A single numeric value
#' @export
#' @return A string value such as 1k for 1000
#' @export
#' @examples 
#' Abbrev(1700000)
Abbrev <- function(value){
     returnValue <- value

        if(value >= 1000000000000){
            returnValue <- paste(as.numeric(value)/1000000000000,"t", sep="")
        }
        else if(value >= 1000000000){
            returnValue <- paste(as.numeric(value)/1000000000,"b", sep="")
        }
         else if(value >= 1000000){
            returnValue <- paste(as.numeric(value)/1000000,"m", sep="")
        }
         else if (value >= 1000){
            returnValue <- paste(as.numeric(value)/1000,"k", sep="")
              }
     return(returnValue)
} 

############################################################################################################################
#' @title Make a straight-line faux CDF.
#' @description Creates a placeholder CDF (a uniform straight line) 
#'    for demonatration.
#'
#' @param range A vector length 2 containing user-specified min and max to
#'    truncate the universe to
#' @param gran The smallest unit of measurement in the data (one [year]
#'    for a list of ages)
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#' @return A fake CDF for demonstration only.
#' @export
#' @examples 
#' badCDF(c(1,50), 1)
badCDF <- function(range, gran,...){
  Domain <- seq(range[1], range[2], gran)
  theCDF <- Domain/range[2]
return(theCDF)
}

############################################################################################################################
############################################################################################################################
############################################################################################################################
#' @title Test a single CDF implementation with one set of parameters.
#' @description Generates mean/median empirical error measurements, complete
#'    results, single iterations of DP CDFs at each combination of parameters,
#'    and diagnostic functions used.
#'
#' @param funct The differentially-private CDF-generating function to be tested
#' @param eps Epsilon value for Differential privacy control
#' @param cdfstep The step sized used in outputting the approximate CDF; the 
#'   values output are [min, min + cdfstep], [min, min + 2 * cdfstep], etc. 
#'   Setting cdfstep equal to 0 (default) will set cdfstep = granularity
#' @param data A vector of the data (single variable to compute CDFs from) 
#' @param range A vector length 2 containing user-specified min and max to
#'   truncate the universe to.
#' @param gran The smallest unit of measurement in the data (one [year] for a
#'   list of ages). The Domain (ie gran and range) should be identical to those 
#'   used to create the CDF!
#' @param reps The number of times the combination of CDFfunction, dataset, and 
#'   epsilon will be tested
#' @param SmoothAll Applies L2 monotnocity post-processing to every DP-CDF
#' @param ABounds This is a flag and should be set to "true" if the functions 
#'   being tested are expected to output 
#'   analytical variance bounds. The proper output form is 
#'   output = list(DPCDFvector, LowerBoundVector, UpperBoundVector)
#' @param EmpiricBounds When TRUE, outputted graphs
#'   depict the minimum and maximum values taken by each bin across reps
#' @param ExtraTests_CDF If a user wishes to add extra diagnostics, the proper 
#'   syntax would be:
#'   ExtraTests_CDF = list( functionName1 = function1, functionName2 = function2) 
#' @param  ExtraTests_PDF See above
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#' @return A list in the form of:
#' 
#' ...\code{$meanscores} Contains mean diagnostic
#'  results for each diagnostic across reps iterations;
#'  
#' ...\code{$medianscores} Contains median diagnostic 
#'  results for each diagnostic across reps iterations;
#'  
#'...\code{$yourCDFoutput} Containing a single dpCDF iteration;
#'
#'...\code{$yourPDFoutput}  Containing a single dpPDF iteration;
#'
#' ...\code{$realCDFoutput} Containing the real (non-DP) CDF output;
#' 
#' ...\code{$realPDFoutput} Containing the real (non-DP) PDF output;
#' 
#' ...\code{$databins} Containing the domain used to construct the CDFs;
#' 
#' ...\code{$TestPack_CDF} Containing the definitions of diagnostic functions used on dpCDFs;
#' 
#' ...\code{$TestPack_PDF} Containing the definitions of diagnostic functions used on dpPDFs;
#' 
#' ...\code{$allscores} Containing all raw diagnostic output.
#' @export
#' @examples 
#' CDFtestTrack(badCDF, eps = .01, cdfstep = 1, data = rexp(10000,.4),
#' range= c(1,10), gran = .1, reps = 20)
CDFtestTrack <- function(funct, eps, cdfstep = 1, data, range, gran, reps, SmoothAll = FALSE, ABounds = FALSE, EmpiricBounds = FALSE, ExtraTests_CDF = list(), ExtraTests_PDF = list(), ...) {


# determine range, universe_size, and create vectors to hold outputs of
# the non-private and private-function (for reps-many iterations)

range         <- as.vector(range)
universe_size <- floor(((range[2] - range[1]) / gran) +1)
output_size   <- floor(((range[2] - range[1]) / cdfstep) +1)

yourCDFoutput <- matrix(nrow=output_size,ncol= reps)
yourPDFoutput <- matrix(nrow=output_size,ncol= reps)

# set up the vector for the x axis
databins      <- seq(from=(range[1]), to=(range[2]), by=cdfstep) 

# this will all output scores
allscores     <- list()

#this is how we'll check later if bounds should be used
LBoundAnalytic <- 100

CDF           <- ecdf(data)     # use the ecdf function to define my CDF function, for an ordindary CDF
realCDFoutput <- CDF (databins) # fill the realCDF vector
realPDFoutput <-c()
realPDFoutput[1]<-realCDFoutput[1]
  for(v in length(realCDFoutput):2){
                realPDFoutput[v] <- realCDFoutput[v]-  realCDFoutput[v-1]
  }


for(i in 1:reps) {

              #run the function as many times as specified by "reps"
              results<- funct(eps=eps, data=data, range=range,gran= gran, 
                        varbounds=Abounds)
              # if AnalyticBounds is true, the inputted function is expected to output 3 vectors:
              # the DP-CDF, a lower bound and an upper bound

if (SmoothAll == FALSE){

if (cdfstep==gran){
                        if (ABounds == TRUE){ 
                             
                                    yourCDFoutput[,i]  <-results[[1]]                                                      
                                if (i ==1){

                                    LBoundAnalytic<- results[[1]] - results[[2]]
                                    UBoundAnalytic<- results[[1]] + results[[2]]
                                }
                        }
                             else{
                                    yourCDFoutput[,i]  <-results[[1]]                                 
                             }
             }else{
                    output_seq <- (seq(range[1],range[2], by=cdfstep)/gran)-(range[1]/gran)+range[1]
                       if (ABounds == TRUE){                               
                                    WorkingVector <-results[[1]]
                              if (i ==1){
                                    WorkingL<- results[[1]] - results[[2]]
                                    WorkingU<- results[[1]] + results[[2]]

                                   for(c in 1:length(output_seq)){
                                     yourCDFoutput[c,i] <- WorkingVector[output_seq[c]]
                                     LBoundAnalytic<- WorkingL[output_seq[c]]
                                     UBoundAnalytic<- WorkingU[output_seq[c]]
                                   }     
                              }
                       }
                       else if (length(results[[1]]) > 1){
                               WorkingVector <-results[[1]]                     
                               for(c in 1:length(output_seq)){                
                               yourCDFoutput[c,i] <- WorkingVector[output_seq[c]]
                               }
                       }
                       else{
                               for(c in 1:length(output_seq)){
                                 yourCDFoutput[,i] <-results[output_seq[c]]
                               }

                       }
             }

                     yourPDFoutput[1,i] <-yourCDFoutput[1,i]
                   for(k in length(yourCDFoutput[,1]):2){
                       yourPDFoutput[k,i] <- yourCDFoutput[k,i] -  yourCDFoutput[k-1,i]
                   }
                    yourPDFoutput[,i] <- round(yourPDFoutput[,i],3)  
                   for (j in 1:length(yourPDFoutput[,i])){
                    yourPDFoutput[j,i] <- MovetoRange(yourPDFoutput[j,i], c(0,1))
                   }
           
}else{
nonMonoCDFs <- matrix(nrow=output_size,ncol= reps)

 if (cdfstep==gran){
                        if (ABounds == TRUE){
                            nonMonoCDFs[,i]   <- results[[1]]   
                            yourCDFoutput[,i] <- smoothVector2(results[[1]])                                                            
                                if (i ==1){
                                    LBoundAnalytic<- results[[1]] - results[[2]]
                                    UBoundAnalytic<- results[[1]] + results[[2]]
                                }
                        }
                        else if (length(results[[1]]) > 1){
                             nonMonoCDFs[,i]   <- results[[1]]                                  
                             yourCDFoutput[,i] <- smoothVector2(results[[1]])                                 
                             }
                        else{                                 
                            yourCDFoutput[,i] <- smoothVector2(results[[1]])                                 
                            }
             }else{
                    output_seq <- (seq(range[1],range[2], by=cdfstep)/gran)              
                       if (ABounds == TRUE){                               
                           WorkingVector <- smoothVector2(results[[1]])
                           NonMonoWV   <- results[[1]]                               
                              if (i ==1){
                                    WorkingL<- WorkingVector - results[[2]]
                                    WorkingU<- WorkingVector + results[[2]]
                                    LBoundAnalytic <-c()
                                    UBoundAnalytic <-c()
                                   for(c in 1:length(output_seq)){
                                     yourCDFoutput[c,i] <- WorkingVector[output_seq[c]]
                                     LBoundAnalytic[c]<- WorkingL[output_seq[c]]
                                     UBoundAnalytic[c]<- WorkingU[output_seq[c]]
                                     nonMonoCDFs[c,i] <- NonMonoWV[output_seq[c]]
                                   }     
                              }
                       }
                       else if (length(results[[1]]) > 1){                               
                                    WorkingVector <- smoothVector2(results[[1]])   
                                    NonMonoWV     <- results[[1]]                                         
                               for(c in 1:length(output_seq)){                
                               yourCDFoutput[c,i] <- WorkingVector[output_seq[c]]
                               nonMonoCDFs[c,i] <- NonMonoWV[output_seq[c]]
                               }
                       }
                       else{
                               for(c in 1:length(output_seq)){
                               yourCDFoutput[,i] <-results[output_seq[c]]
                               nonMonoCDFs[c,i]  <- NonMonoWV[output_seq[c]]
                               }
                       }
             }

                     yourPDFoutput[1,i] <-nonMonoCDFs[1,i]
                   for(k in length(nonMonoCDFs[,1]):2){
                       yourPDFoutput[k,i] <- nonMonoCDFs[k,i] -  nonMonoCDFs[k-1,i]
                   }
                    yourPDFoutput[,i] <- round(yourPDFoutput[,i],3)     
                   for (j in 1:length(yourPDFoutput[,i])){
                    yourPDFoutput[j,i] <- MovetoRange(yourPDFoutput[j,i], c(0,1))
                   }
      }
}         

################################################
TestPack_CDF <-list(MaxError_CDF   = MaxError_CDF,
                    MaxErrorAt_CDF = MaxErrorAt_CDF,
                    diffat25       = diffat25,
                    diffatMedian   = diffatMedian,
                    diffat75       = diffat75,
                    horzdiffat25   = horzdiffat25,
                    horzdiffatMed  = horzdiffatMed,
                    horzdiffat75   = horzdiffat75,
                    Medians        = Medians,
                    MAE_CDF        = MAE,
                    MSE            = MSE,
                  #  L1emp          = L1empiric,
                  #  L2emp          = L2empiric,
                    DerivScore     = DerivDiff             
              )

TestPack_PDF <-list(MaxError_PDF   = MaxError_PDF,
                    MaxErrorAt_PDF = MaxErrorAt_PDF,
                    MAE_PDF        = MAE,
                    MSE_PDF        = MSE,
                    MeanDiff = MeanDiffpdf,
                     ModeDiff = ModeDiffpdf,
                      StdDiff = StdDiffpdf,
                       VarDiff = VarDiffpdf,
                        SkewDiff = SkewDiffpdf,
                         KurtDiff= KurtDiffpdf)


## if the user has any custom built DP-CDF or DP-PDF empirical diagnostic functions, they're introduced here.

if (length(ExtraTests_CDF) >0){
assist <- names(ExtraTests_CDF)
  
    for(q in 1:length(ExtraTests_CDF)){
      TestPack_CDF[length(TestPack_CDF)+q] <- ExtraTests_CDF[q]
       names(TestPack_CDF)[length(TestPack_CDF)] <- assist[q]
    }
  }

if (length(ExtraTests_PDF) >0){
assist <- names(ExtraTests_PDF)
    for(q in 1:length(ExtraTests_PDF)){
      TestPack_PDF[length(TestPack_PDF)+q] <- ExtraTests_PDF[q]
       names(TestPack_PDF)[length(TestPack_PDF)] <- assist[q]
    }
  }

############################################################
# initialize dataframes to hold the single output of each diagnostic.
TheNames <- names(TestPack_CDF)
      MeanScores   <- data.frame( TheNames = 0)
      MedianScores <- data.frame( TheNames = 0)

for (u in 1:length(TestPack_CDF)){
      TempVector <- rep(0, times = reps)

      #and fill those frames
      for (i in 1:reps){
           TempVector[i] <- TestPack_CDF[[u]](Y = realCDFoutput, est = yourCDFoutput[,i], range = range, gran = gran, eps = eps, data= data)
      }

      MeanScores  [[TheNames[u]]] <- mean  (TempVector)
      MedianScores[[TheNames[u]]] <- median(TempVector)
      allscores   [[TheNames[u]]] <-        TempVector
      }


#####
TheNames <- names(TestPack_PDF)
for (u in 1:length(TestPack_PDF)){
      TempVector <- rep(0, times = reps)
      #and fill those frames

         for (i in 1:reps){

         TempVector[i] <- TestPack_PDF[[u]](Y = realPDFoutput, est = yourPDFoutput[,i], range = range, gran = gran)
         }

      MeanScores  [[TheNames[u]]] <- mean  (TempVector)
      MedianScores[[TheNames[u]]] <- median(TempVector)
      allscores   [[TheNames[u]]] <-        TempVector
}


  #make empirical bounds for the CDF
  if (EmpiricBounds){

        differences <- matrix(0,nrow = length(databins), ncol = reps)

        for (i in 1: reps){
             differences[,i]<- yourCDFoutput[,i] - realCDFoutput
        }
        UBoundEmpiric <- rep(0,length(databins)) 
        LBoundEmpiric <- rep(0,length(databins)) 

        for (b in 1:length(databins)){
             UBoundEmpiric[b] <- max(differences[b,]) + realCDFoutput[b]
             LBoundEmpiric[b] <- min(differences[b,]) + realCDFoutput[b]
        } 
}
############################################################
# store everything! Including two arbitrary iterations of the dp-CDF function,
# used for plotting in the the CDFtest graphical ouput
scores <- list()
scores$meanscores      <-  MeanScores
scores$medianscores    <-  MedianScores
scores$yourCDFoutput   <-  yourCDFoutput[,1]
scores$yourPDFoutput   <-  yourPDFoutput[,1]
scores$realCDFoutput   <-  realCDFoutput
scores$realPDFoutput   <-  realPDFoutput
scores$databins        <-  databins
scores$TestPack_CDF    <-  TestPack_CDF
scores$TestPack_PDF    <-  TestPack_PDF
scores$allscores       <-  allscores

if(EmpiricBounds == TRUE){
    scores$yourCDF_LboundE <-  LBoundEmpiric
    scores$yourCDF_UboundE <-  UBoundEmpiric
}
if(ABounds == TRUE){
    scores$yourCDF_LboundA <-  LBoundAnalytic
    scores$yourCDF_UboundA <-  UBoundAnalytic
}


return(scores)
}

############################################################################
#This is a streamlined version of \code{CDFtestTrack} (see above), used in Data Collection mode. 
#' @title Test a single CDF implementation with one set of parameters.
#' @description Applies diagnostic functions to a single dpCDF, and only 
#'   releases a complete set of diagnostic results (called within\code{CDFtest}
#'   in Data Collection mode --- e.g., when \code{Visualization = FALSE})
#'
#' @param funct The differentially-private CDF-generating function to be tested
#' @param eps Epsilon value for Differential privacy control
#' @param data A vector of the data (single variable to compute CDFs from) 
#' @param range A vector length 2 containing user-specified min and max to
#'   truncate the universe to
#' @param gran The smallest unit of measurement in the data (one [year] for a 
#'   list of ages)
#' @param reps The number of times the combination of CDFfunction, dataset,
#'   and epsilon will be tested
#' @param samplesize The specified sample size is randomly selected from each dataset
#'   without replacement. 
#' @param SmoothAll Applies L2 monotonicity post-processing to every DP-CDF
#' @param ExtraTests_CDF If a user wishes to add extra diagnostics, the proper 
#'   syntax would be:   
#'   \code{ExtraTests_CDF = list( functionName1 = function1, functionName2 = function2)}
#' @param  ExtraTests_PDF See above
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.

#' @return A complete set of diagnostic results in the form of
#'       \code{...$allscores}, which holds out a row of output for each of \code{reps} results.
#' @export
#' @examples 
#' CDFtestTrackx(badCDF, eps = .01, cdfstep = 0, data = rexp(10000,.4),
#'   range= c(1,10), gran = .1, reps = 20, samplesize = 10000)
CDFtestTrackx <- function(funct, eps, data, range=range, gran, reps,samplesize, SmoothAll = FALSE,ExtraTests_CDF = list(), ExtraTests_PDF = list(),...) {
# source relevant scripts and libraries
# determine range, universe_size, and create vectors to hold outputs of the non-private and private-function (for reps-many iterations)
  smoothVector <- NULL
range         <- as.vector(range)
universe_size <- floor(((range[2] - range[1]) / gran) +1)
output_size   <- universe_size

data <- sample(data, samplesize, replace = TRUE)

yourCDFoutput <- matrix(nrow=output_size,ncol= reps)
yourPDFoutput <- matrix(nrow=output_size,ncol= reps)

# set up the vector for the x axis
databins      <- seq(from=(range[1]), to=(range[2]), by=gran) 

# this will all output scores
allscores     <- list()

#this is how we'll check later if bounds should be used
LBoundAnalytic <- 100

CDF           <- ecdf(data)     # use the ecdf function to define my CDF function, for an ordindary CDF
realCDFoutput <- CDF (databins) # fill the realCDF vector
realPDFoutput <-c()
realPDFoutput[1]<-realCDFoutput[1]
  for(v in length(realCDFoutput):2){
                realPDFoutput[v] <- realCDFoutput[v]-  realCDFoutput[v-1]
  }
MSEanalytic <- c()
for(i in 1:reps) {
              #run the function as many times as specified by "reps"
              results<- funct(eps=eps, data=data, range=range,gran= gran, 
                        Analyticbounds=AnalyticBounds)

              # if AnalyticBounds is true, the inputted function is expected to output 2 vectors:
              # the DP-CDF, and MSEanalytic abs diffs (analytic bounds)

                     yourCDFoutput[,i]  <-results[[1]]
                    MSEanalytic[i]            <-mean(results[[2]])                                  
           
                   CDFshift<- c(0)
                   CDFshift[2:length(yourCDFoutput[,i])]<-yourCDFoutput[(1:length(yourCDFoutput[,i])-1),i]
                   yourPDFoutput[,i] <-yourCDFoutput[,i] - CDFshift
                   yourPDFoutput[,i] <- round(yourPDFoutput[,i],3)  


if (SmoothAll){                               
                     yourCDFoutput[,i] <- smoothVector2(results[[1]])      
      }
}         

################################################
TestPack_CDF <-list(MaxError_CDF   = MaxError_CDF,
                    MaxErrorAt_CDF = MaxErrorAt_CDF,
                    diffat25       = diffat25,
                    diffatMedian   = diffatMedian,
                    diffat75       = diffat75,
                    horzdiffat25   = horzdiffat25,
                    horzdiffatMed  = horzdiffatMed,
                    horzdiffat75   = horzdiffat75,
                    Medians        = Medians,
                    MAE_CDF        = MAE,
                    MSE            = MSE
                  #  L1emp          = L1empiric,
                  #  L2emp          = L2empiric,
                  #  DerivScore     = DerivDiff,
                  # MSEanalytic     = MSEanalytic                
              )

TestPack_PDF <-list(MaxError_PDF   = MaxError_PDF,
                    MaxErrorAt_PDF = MaxErrorAt_PDF,
                    MAE_PDF        = MAE,
                    MSE_PDF        = MSE,
                    MeanDiff = MeanDiffpdf,
                     ModeDiff = ModeDiffpdf,
                      StdDiff = StdDiffpdf,
                       VarDiff = VarDiffpdf,
                        SkewDiff = SkewDiffpdf,
                         KurtDiff= KurtDiffpdf)

## if the user has any custom built DP-CDF or DP-PDF empirical diagnostic functions, they're introduced here.

if (length(ExtraTests_CDF) >0){

assist <- names(ExtraTests_CDF)
  
    for(q in 1:length(ExtraTests_CDF)){
      TestPack_CDF[length(TestPack_CDF)+1] <- ExtraTests_CDF[q]
       names(TestPack_CDF)[length(TestPack_CDF)] <- assist[q]
    }
  }

if (length(ExtraTests_PDF) >0){

assist <- names(ExtraTests_PDF)

    for(q in 1:length(ExtraTests_PDF)){
      TestPack_PDF[length(TestPack_PDF)+1] <- ExtraTests_PDF[q]
       names(TestPack_PDF)[length(TestPack_PDF)] <- assist[q]

    }
  }

############################################################
# initialize dataframes to hold the single output of each diagnostic.

TheNames <- names(TestPack_CDF)
      MeanScores   <- data.frame( TheNames = 0)
      MedianScores <- data.frame( TheNames = 0)
for (u in 1:length(TestPack_CDF)){
      TempVector <- rep(0, times = reps)

      #and fill those frames
      for (i in 1:reps){
           TempVector[i] <- TestPack_CDF[[u]](Y = realCDFoutput, est = yourCDFoutput[,i], range = range, gran = gran, eps = eps, data= data)
      }

      MeanScores  [[TheNames[u]]] <- mean  (TempVector)
      MedianScores[[TheNames[u]]] <- median(TempVector)
      allscores   [[TheNames[u]]] <-        TempVector
      }

#####
TheNames <- names(TestPack_PDF)
for (u in 1:length(TestPack_PDF)){
      TempVector <- rep(0, times = reps)
      #and fill those frames
         for (i in 1:reps){

         TempVector[i] <- TestPack_PDF[[u]](Y = realPDFoutput, est = yourPDFoutput[,i], range = range, gran = gran)

         }

      MeanScores  [[TheNames[u]]] <- mean  (TempVector)
      MedianScores[[TheNames[u]]] <- median(TempVector)
      allscores   [[TheNames[u]]] <-        TempVector
}

      allscores$MSEanalytic        <- MSEanalytic

############################################################
# store everything! Including two arbitrary iterations of the dp-CDF function, used for plotting in the the CDFtest graphical ouput
scores <- list()
scores$allscores       <-  allscores


return(scores)
}

###############################################################################
###############################################################################
###############################################################################
#' @title Comprehensively evaluate and visualize the utility of 
#'     CDF-generating implementations.
#' @description The suite is a system for determining the utility
#'     of differentially private cumulative distribution function (DP-CDF)
#'     algorithm implementations. The system can empirically evaluate and 
#'     provide visualizations for several DP-CDF algorithms simultaneously,
#'     under various parameters. It can also be set to focus strictly on data
#'     collection, rather than spending time on visualization. 
#'     
#'     It comes with several pre-loaded adjustable synthetic datasets, and can
#'     also analyze functions on user-defined datasets.
#'     
#'     dpCDF implementations to test must take the following as arguments:
#'     \code{data, epsilon, granularity, range}, and any number of other inputs.
#'     Use "?functionH" for an example of an implementation drawing on C++ files
#'      through Rcpp.
#'     
#'       USERS SHOULD NOTE: the
#'     following included diagnostic functions are under development:
#'     \code{SkewDiffpdf,KurtDiffpdf, StdDiffpdf}, corresponding to error measurements of
#'     skewness, kurtoses, and standard deviations generated from dpCDFs.
#'     This is evident through the occasional result of \code{NA}.
#' @param Visualization Sets the testing suite into Visualization mode (default, 
#'        \code{Visualization = TRUE}) 
#'     or Data Collection mode \code{(Visualization = FALSE)}
  #'  In Visualization mode (default):
  #'    A \code{.csv} file conatining the mean and median results (across \code{reps}
  #'        iterations) of diagnostic functions on DP-CDF algorithms per each 
  #'        combination of data, function, and epsilon. 
  #'    A \code{.pdf} file containing one graphical example DP CDF for each combination
  #'        of dataset, function, and epsilon, as well as a set of boxplots 
  #'        showing the distribution of all diagnostic results for all 
  #'        combinations of parameters.
  #'  In Data Collection mode (set \code{Visualization = FALSE}):
  #'    A \code{.csv} file containing the entire (raw) results (across \code{reps}
  #'    iterations) of diagnostic functions on DP-CDF algorithms
  #'    per each combination of dataset, and function, seperately looped over
  #'    all epsilons, then all granularities, and all samplesizes. 
#' @param OutputDirectory The location of the folder which will hold the
#'    output (\code{.csv} and \code{.pdf} files). This defaults to the
#'    \code{tempdir()} directory.
#' @param functlist A list of CDF-computing functions to be tested on the
#'    \code{CDFtestTrack} (if \code{visualization = TRUE}) or \code{CDFtestTrackx}
#'    (if \code{Visualization =FALSE})) 
#' @param Fnameslist A vector of function names corresponding to the functions
#' @param epslist  A vector of epsilon values for differential privacy
#' @param datalist A list containing vectors of data, each to be used in a test
#' @param Dnameslist A list of dataset names corresponding to the data/variables 
#'   being tested; used for labelling the output
#' @param SmoothAll Applies L2 monotnocity post-processing to every DP-CDF
#' @param synthsets This script generates pre-defined synthetic datasets upon 
#'  request, and fully incorporates them into testing. To call them, users
#'  should input a string vector containing the names of the sets they desire.
#'    For example, \code{synthsets = list(list(type,size,shape),list(type,size,shape))}.
#'    There are no limits on the amounts of datasets included.
#'    Sets available include:
#'       type:  \code{"age"} (which ranges from about 0 to 100, \code{gran =1})
#'       and \code{"wage"} which ranges from 0 to 500k);
#'       size:  Any positive integer. Type in exact numerical representation
#'         (eg, for ten thousand use 10000 not 10k and not 10^4);
#'       shape: gaussian, sparse, uniform, bimodal;
#'       It is assumed that the data input is rounded to the granularity
#' @param range The range of the domain as a vector \code{c(min, max)}. Defined based 
#'   on user intuition. to preserve differential privacy, the domain is constructed
#'   using this range. Setting the min too high will bias output upward. 
#'   Same in reverse for a low max. However, setting min too low and max 
#'   too high could reveal the true limits of your data,
#'   compromising some privacy.
#' @param gran FOR Visualization MODE ONLY. refer to \code{granlist} for setting 
#'   granularities (thus domain sizes) in Data Collection mode.
#'   This command is irrelevant in Data Collection mode.          
#'   The granularity of the domain between the min and max. ie, if age is
#'   measureds per 1 year of age, \code{gran =1}.
#'   The same granularity is applied to all datasets, so using comparable 
#'   (or scaled) data is necessary.
#' @param granlist FOR Data Collection MODE ONLY. refer to \code{gran} for selecting
#'   samplesizes in Data Collection mode. This command is irrelevant in
#'   Visualization mode. A list of granularities of the domain between the
#'   min and max. ie, if age is measure per 1 year of age, \code{gran =1}.
#' @param samplesize FOR Visualization MODE ONLY. refer to \code{nlist} for selecting
#'   samplesizes in Data Collection mode. This command is irrelevant in 
#'   Data Collection mode. when set to zero, the entire dataset is used.
#'   Otherwise, the specified sample size is randomly selected from each dataset
#'   without replacement. 
#' @param nlist FOR Data Collection MODE ONLY. refer to \code{samplesize} for 
#'   selecting samplesizes in visualization mode. This command is irrelevant in
#'   Visualization mode. Sets the absolute sample sizes to draw from each
#'   dataset, with replacement. Any vector of integer values is appropriate.
#' @param cdfstep The step size used in outputting the approximate CDF;
#' @param reps  The number of times to repeat each diagnostic. higher \code{reps} lends
#'   greater accuracy, but comsumes time and power. Author recommends \code{reps = 10}
#'   for quick examples and \code{reps = 100} for more robust examinations.
#' @param ExtraTests_CDF If a user wishes to add extra diagnostics, the proper
#    syntax would be:
#'     \code{ExtraTests_CDF = list(functionName1=function1, functionName2=function2)}.  
#'    Diagnostic Functions should have inputs such as \code{Y} for a public CDF,
#'     \code{est} for a DP-representation of that CDF,
#'     \code{range} and \code{gran}, and the output should be just one value.
#' @param ExtraTests_PDF See above
#' @param setseed In the function, each combination of data, epsilon, and
#'   function is executed with a separate seed, which by default is randomly
#'   generated and reported. Users interested in replicating specific results
#'   can locate the reported seed and parameter combination to replicate tests.
#' @param comments "Comments written here print to a log in excel" 
#' @param EmpiricBounds FOR Visualization MODE ONLY. When TRUE, outputted graphs
#'   depict the minimum and maximum values taken by each bin across reps
#' @param AnalyticBounds FOR Visualization MODE ONLY. This is a flag and should
#'   be set to \code{TRUE} if the functions being tested are expected to output 
#'  analytical variance bounds. The proper output form for such a function is 
#'   \code{output = list(DPCDFvector, LowerBoundVector, UpperBoundVector)}.
#' @param AnalyticProbSleeve FOR Visualization MODE ONLY. When \code{TRUE}, outputted
#'   DP-CDFs will have a 'fuzzy' analytic sleeve around them, approximating
#'   probabalitity density for each point given by DP. This also requires the
#'   function format specified above in the description for \code{AnalyticBounds}.
#' @param SuppressRealCDF FOR Visualization MODE ONLY. When \code{TRUE}, outputted
#'   graphs will not include real (non-private) CDFs.
#' @param SuppressDPCDF   FOR Visualization MODE ONLY. When \code{TRUE}, outputted
#'   graphs will not include DP-CDFs (but if \code{SmoothAll = TRUE}, monotonized
#'   DP CDFs still appear).
#' @param SuppressLegends FOR Visualization MODE ONLY. When \code{TRUE}, outputted
#'   graphs will not include legends
#' @param ... Optionally add additional parameters. This is primarily used to allow automated
#'   execution of varied diagnostic functions.
#' @export
#' @return If \code{Visualization = TRUE}, a list containing:                        
#'
#'  ...\code{$means} Contains mean diagnostic
#'  results for each diagnostic across reps iterations for each parameter combination;
#'  
#' ..\code{$medians} Contains median diagnostic 
#'  results for each diagnostic across reps iterations for each parameter combination;\\
#'
#'...\code{$yourCDFoutput} Containing a single dpCDF iteration for each parameter combination;\\
#'
#'...\code{$yourPDFoutput}  Containing a single dpPDF iteration for each parameter combination;\\
#'
#'...\code{$realCDFoutput} Containing the real (non-DP) CDF output for each relevant parameter combination;
#'
#'...\code{$realPDFoutput} Containing the real (non-DP) PDF output for each relevant parameter combination;
#'
#'...\code{$databins} Containing the domain used to construct the CDFs;
#'
#'...\code{$TestPack_CDF} Containing the definitions of diagnostic functions used on dpCDFs;
#'
#'...\code{$TestPack_PDF} Containing the definitions of diagnostic functions used on dpPDFs;
#'
#'...\code{$allscores} Containing all raw diagnostic output.
#'
#'...\code{$seed} Containing the list of seeds used in the test
#'
#'...\code{$permetric} holding a rearranged dataframe (ordered by parameter
#'                   combinations) useful for plotting.
#'  
#'  A \code{.pdf} file:
#'      with boxplots showing the distributions of diagnostic outputs, 
#'      and categorized plots of dp-CDF function output. Each such graph with
#'      show one arbitrary CDF iterations and empirical boundaries.
#'      the empirical boundaries are the max and min values reached by that
#'      function (and parameters) during the test.
#'              
#'   A \code{.csv} file:
#'      containing the mean and median scores of each diagnostic on each 
#'      combination of data, eps, function, and the seedlist for reproduction.
#'
#'  Notes on Visualization mode: Both the \code{.pdf} and \code{.csv} components are named with a time stamp index,
#'  in the form of \code{YearMonthDayHourMinuteSecond}. To locate particular tests,
#'  look at the \code{CDFtestindexchart.csv}, which automatically records the
#'  parameters and index of each test. These can be found in the file specified by 
#'  \code{OutputDirectory}, which defaults to the R temp files \code{tempdir()}.
#'  
#'  Alternatively in Data Collection mode (\code{Visualization = FALSE}), a list containing:
#'
#'      ...\code{$allscores} holding the output of each combination of parameters,
#'      which is that each \code{eps} in epslist is varied across the first value
#'      specified in \code{granlist} and \code{nlist}. The same is true for varying
#'      granularity and sample size. In that way, only one variable is varied at a time
#'      while the other two are held fixed. All such combinations of parameters are
#'      executed on all combinations of data and function (specified within       
#'      ...\code{datalist} and \code{functlist});
#'      
#'      ...\code{$seed} holding the list of seeds used in the test.
#'      
#'    A \code{.csv} file conatining the entire (raw) results (across reps iterations) 
#'    of diagnostic functions on DP-CDF algorithms per each combination of 
#'    dataset, and function, looped over epsilon, granularity, and sample size values
#'    as described directly above.\\
#'    This mode was designed for collecting metric data for subsequent supervised 
#'    learning modelling. 
#'
#' @examples 
#' CDFtest( Visualization = TRUE,OutputDirectory = 0, functlist = c(functionH),
#' Fnameslist = c("H"), epslist  = c(.1, .01), datalist = list(),
#' Dnameslist = c(), synthsets= list(list("wage", 100000, "uniform"), 
#'  list("wage",100000,"sparse"), list("wage",100000,"bimodal")),
#'  range    = c(1,500000),gran =1000,granlist =c(2500,1250,1000,500), 
#'  samplesize = 0,nlist = c(100,1000,10000,100000,1000000),
#'  cdfstep  =0, reps = 5,  ExtraTests_CDF = list(),ExtraTests_PDF = list(),
#'  setseed = c(-100),
#'  comments = "x",SmoothAll = FALSE,EmpiricBounds = FALSE,
#'  AnalyticBounds = FALSE,AnalyticProbSleeve = FALSE,
#'  SuppressRealCDF = FALSE,SuppressDPCDF = FALSE,SuppressLegends = FALSE)

CDFtest <- function(
 Visualization = TRUE,  OutputDirectory= 0, 
 functlist, Fnameslist, epslist = c(.05,.1,1),
 datalist,  Dnameslist, synthsets =NULL, range, gran = 1, granlist = c(1), samplesize =0, nlist = (10000), cdfstep =1, 
 reps = 5,  ExtraTests_CDF = list(), ExtraTests_PDF = list(), setseed = -100,  comments="none",
 SmoothAll = FALSE,       EmpiricBounds = FALSE, AnalyticBounds  = FALSE,  AnalyticProbSleeve = FALSE, 
 SuppressRealCDF = FALSE, SuppressDPCDF = FALSE, SuppressLegends = FALSE,... ){


if(OutputDirectory == 0){
  print("No output directory specified; plotting to R temp directory:")
  OutputDirectory <- tempdir()
}else{
  print("Plotting to the following directory:")
}
print(OutputDirectory)

if(Visualization){
print("Using Visualization Mode")
}else{
  print("Using Data Collection Mode")
}

# stops
if (length(epslist)<1){
  print("STOP: At least one epsilon value required")
  return("function stopped")
}
for(w in 1:length(epslist)){
  if (epslist[w] <= 0| is.numeric(epslist[w]) == FALSE){
    print("STOP: Epsilon values must be positive.")
    return("function stopped")
  } 
}

if(range[1] > range[2]){
    print("STOP: range[1] value must be the minimum observation value; range[2] is the maximum.")
    return("function stopped")
}


if (cdfstep == 0){
print("cdfstep set to 0, using maximum resolution (cdfstep = gran)")
cdfstep <- gran
}
else if (Visualization){
 if(cdfstep < gran){
    print("STOP: cdfstep must be greater than gran, or set to 0 to equal gran.")
    return("function stopped")
}}
else{
  for(z in 1:length(granlist)){
    if (cdfstep < granlist[z]){
    print("STOP: cdfstep must always be greater than gran, or set to 0 to equal gran.")
    return("function stopped")
    }
  }
}

if (reps<1 | reps%%1 != 0){
    print("STOP: reps must be a positive integer")
    return("function stopped")
}

if(Visualization){
if (samplesize < 0){
    print("STOP: samplesize must be a positive integer, or 0 to use the full data")
    return("function stopped")
}
else if (samplesize == 0){
  print("samplesize set to 0, using entire dataset")
}
}else{

  for (z in 1:length(nlist)){
    if (nlist[z] <1){
          print("STOP: samplesize must be a positive integer, or 0 to use the full data")
    return("function stopped")
   }
  }

}
################## SAMPLE THE DATASETS #######################

if(length(datalist)!=0){

  if(samplesize != 0){
      for(d in 1:length(datalist)){
        datalist[[d]]<- sample(datalist[[d]],samplesize)
      }
    }
 
 dataSizeNaming<- list()
for(d in 1:length(datalist)){
  dataSizeNaming[d] <- Abbrev(length(datalist[[d]]))
         Dnameslist[d]     <- paste(Dnameslist[d],dataSizeNaming[d], sep = ".")
  }
 }
######################### CREATE SYNTHETIC DATASETS ###############################

if(length(synthsets)!=0){

 for (s in 1:length(synthsets)){

   #########GAUSSIAN-SHAPED##########

   if(!is.null(synthsets[[s]]) && synthsets[[s]][3] == "gaussian"){
      if (synthsets[[s]][1] == "age"){
        set.seed=100
        NewData <- round(rnorm(as.numeric(synthsets[[s]][2]),38,25),0)
        }

      if (synthsets[[s]][1] == "wage"){
        set.seed=100
        NewData <- round(rnorm(as.numeric(synthsets[[s]][2]),220000,100000),-3)
        }
      datalist[[length(datalist)+1]]       <-  NewData

      synthsets[[s]][2] <-Abbrev(synthsets[[s]][2])

      Dnameslist[length(Dnameslist)+1]     <- paste(synthsets[[s]][1],synthsets[[s]][2],synthsets[[s]][3], sep = ".")
   }

   #########UNIFORMLY-SHAPED##########

   if(!is.null(synthsets[[s]]) && synthsets[[s]][3] == "uniform"){
      if (synthsets[[s]][1] == "age"){  
           databins       <- seq(from=(range[1]), to=(range[2]), by=gran)
           tool           <- floor(as.numeric(synthsets[[s]][2])/(range[2]-range[1]+1))
           NewData <- rep(databins, times=tool)
        }
      if (synthsets[[s]][1] == "wage"){
           databins       <- seq(from=(range[1]), to=(range[2]), by=gran)
           tool           <- floor(as.numeric(synthsets[[s]][2])/((range[2]/gran)-((range[1]/gran)+1)))
           NewData <- rep(databins, times=tool)
        }
      datalist[[length(datalist)+1]]       <- NewData

               synthsets[[s]][2] <-Abbrev(synthsets[[s]][2])
                
      Dnameslist[length(Dnameslist)+1]     <- paste(synthsets[[s]][1],synthsets[[s]][2],synthsets[[s]][3], sep=".")
   }

    #########BIMODALLY-SHAPED##########

   if(!is.null(synthsets[[s]]) && synthsets[[s]][3] == "bimodal"){
           if (synthsets[[s]][1] == "age"){
           set.seed=100
           A <- round(rnorm(as.numeric(synthsets[[s]][2])/2,20,10),0)
           set.seed=100
           B <- round(rnorm(as.numeric(synthsets[[s]][2])/2,60,10),0)
           NewData <- c(A,B)
        }
      if (synthsets[[s]][1] == "wage"){

           set.seed=100
           A <- round(rnorm(as.numeric(synthsets[[s]][2]),100000,35000),-3)
           set.seed=100
           B <- round(rnorm(as.numeric(synthsets[[s]][2]),300000,35000),-3)
           NewData <- c(A,B)
        }
      datalist [[length(datalist)  +1]]    <- NewData

             synthsets[[s]][2] <-Abbrev(synthsets[[s]][2])
      Dnameslist[length(Dnameslist)+1]     <- paste(synthsets[[s]][1],synthsets[[s]][2],synthsets[[s]][3], sep=".")
   }
   #########SPARSE##########

   if(!is.null(synthsets[[s]]) && synthsets[[s]][3] == "sparse"){
      if (synthsets[[s]][1] == "age"){

        agerep   <- c(10,11,13,20,21,35,60,61,65)
        timesrep <- round(c(2500,2000,5000,7000,3500,3000,5000,10000,8000)*(35/46),0)
        rand     <- c(15000,30)
        set.seed=100
        NewData <- rpois(round( rand[1]*(as.numeric(synthsets[[s]][2])/50000) ,0),rand[2])
           for(w in 1:9){
           NewData <- c(NewData, rep(agerep[w], times= round((timesrep[w] * (as.numeric(synthsets[[s]][2])/50000)),0)))
          }
         }

      if (synthsets[[s]][1] == "wage"){
       
        wagerep <- c(0,6000,10000,12000,14000,19000,20000,23000,25000,26000,30000,32000,35000,40000,42000,46000,47000,49000,
             50000,53000,60000,65000, 67000, 70000, 75000,80000,86000,87000,92000,95000, 97000,100000,105000, 
             110000,120000,125000,130000,135000,140000,145000,150000,160000, 170000,172000,175000, 180000,190000,200000,211000,220000,225000,230000,240000,250000,260000, 265000, 275000,300000,310000,
             330000,350000,360000,375000, 3850000, 400000)
        timesrep<- c(1000,200,50,2600,300,100,50,50,300,150,150,25,1400,3000,1050,50,100,1750,5600,50,50,1500,900,100,900,250,
              1250,50,50,100,2000,500,25,50,3500,300,200,25,50,1050,50,4000,100,50,150,125,100,1250,2100, 25, 1000, 50,
              1350, 2925, 100, 25, 850, 2050, 100, 75, 25, 25, 2600, 25, 25)
        rand    <- c(5000,200000)
        set.seed=100
       NewData <- rpois(round( rand[1]*(as.numeric(synthsets[[s]][2])/50000) ,0),rand[2])
           for(w in 1:length(wagerep)){
           NewData <- c(NewData, rep(wagerep[w], times= round((timesrep[w] * (as.numeric(synthsets[[s]][2])/50000)),0)))
          }
         }
      datalist[[length(datalist)+1]]       <- NewData

              synthsets[[s]][2] <-Abbrev(synthsets[[s]][2])
      Dnameslist[length(Dnameslist)+1]     <- paste(synthsets[[s]][1],synthsets[[s]][2],synthsets[[s]][3])
   }
}
}

##################################################################################################################################
if (Visualization == FALSE){

  # this takes the given seedlist
seedlist <- setseed  

# if you're looking for a random seed, this will overwrite the 
suppressWarnings

  if (setseed[1] == -100){
          seedlist <- round(runif((length(epslist) * length(functlist)* length(datalist)) + 
                                  (length(granlist)* length(functlist)* length(datalist)) +
                                  (length(nlist)   * length(functlist)* length(datalist)), min=0, max=1000),0)
            }
set.seed(seedlist)

MainOutput <- matrix(nrow = (((length(datalist))*length(granlist)*length(functlist)) +
                             ((length(datalist))*length(epslist)* length(functlist)) + 
                             ((length(datalist))*length(nlist)*   length(functlist)))  *reps,
                     ncol = 0)
MainOutput <- data.frame(MainOutput)

  index <- 1-reps
  # run the function for every CDF function
  for (f in 1:length(functlist)){
  print(paste("__this is functlist",f))

    # run the function for every dataset
    for(d in 1:length(datalist)){
    datalist[[d]] <- as.numeric(lapply(datalist[[d]], function(x) MovetoRange(x, range)))
    print(paste("____this is datalist",d))

      # run the function for every epsilon, hold dSize and n fixed
      for(i in 1:length(epslist)){
      print(paste("_______this is epslist",i))  
   #   index <-  ((f-1)*length(epslist)*length(datalist) + (d-1)*length(epslist)+i + counter)
index <- index+(1*reps)

     dSize <- ceiling(((range[2]-range[1])/1)+1)
     MainOutput$imp  [index:(index+reps-1)] <- Fnameslist[f] 
     MainOutput$eps  [index:(index+reps-1)] <- epslist[i]
     MainOutput$n    [index:(index+reps-1)] <- nlist[1]
     MainOutput$dSize[index:(index+reps-1)] <- dSize
     MainOutput$data [index:(index+reps-1)] <- Dnameslist[d]
     MainOutput$gran [index:(index+reps-1)] <- granlist[1]
     MainOutput$min  [index:(index+reps-1)] <- range[1]
     MainOutput$max  [index:(index+reps-1)] <- range[2]
     MainOutput$seed [index:(index+reps-1)] <- seedlist[index]

          scores1<- CDFtestTrackx(funct = functlist[[f]],
                                eps= epslist[i],
                                data=datalist[[d]],
                                range,
                                gran = granlist[1],
                                samplesize = nlist[1],
                                reps,
                                SmoothAll,
                                ExtraTests_CDF = ExtraTests_CDF,
                                ExtraTests_PDF = ExtraTests_PDF
                                )
       #   MainOutput$Mean          [index:(index+reps-1)] <-scores1$allscores$Mean
          MainOutput$MSEanalytic    [index:(index+reps-1)] <-scores1$allscores$MSEanalytic
          MainOutput$SDempiric     [index:(index+reps-1)] <-scores1$allscores$SDempiric
          MainOutput$MaxError_CDF  [index:(index+reps-1)] <-scores1$allscores$MaxError_CDF 
          MainOutput$MaxErrorAt_CDF[index:(index+reps-1)] <-scores1$allscores$MaxErrorAt_CDF
          MainOutput$diffat25      [index:(index+reps-1)] <-scores1$allscores$diffat25
          MainOutput$diffatMedian  [index:(index+reps-1)] <-scores1$allscores$diffatMedian 
          MainOutput$diffat75      [index:(index+reps-1)] <-scores1$allscores$diffat75
          MainOutput$horzdiffat25  [index:(index+reps-1)] <-scores1$allscores$horzdiffat25
          MainOutput$horzdiffatMed [index:(index+reps-1)] <-scores1$allscores$horzdiffatMed 
          MainOutput$horzdiffat75  [index:(index+reps-1)] <-scores1$allscores$horzdiffat75
          MainOutput$Medians       [index:(index+reps-1)] <-scores1$allscores$Medians
          MainOutput$MAE_CDF       [index:(index+reps-1)] <-scores1$allscores$MAE_CDF
          MainOutput$MSE_CDF       [index:(index+reps-1)] <-scores1$allscores$MSE
       #   MainOutput$DerivScore    [index:(index+reps-1)] <-scores1$allscores$DerivScore
          MainOutput$MaxError_PDF  [index:(index+reps-1)] <-scores1$allscores$MaxError_PDF 
          MainOutput$MaxErrorAt_PDF[index:(index+reps-1)] <-scores1$allscores$MaxErrorAt_PDF
          MainOutput$MAE_PDF       [index:(index+reps-1)] <-scores1$allscores$L1_area_PDF
          MainOutput$MSEanalytic   [index:(index+reps-1)] <-scores1$allscores$MSEanalytic
          MainOutput$MeanDiff      [index:(index+reps-1)] <-scores1$allscores$MeanDiff
          MainOutput$ModeDiff      [index:(index+reps-1)] <-scores1$allscores$ModeDiff                     
          MainOutput$StdDiff       [index:(index+reps-1)] <-scores1$allscores$StdDiff                     
          MainOutput$VarDiff       [index:(index+reps-1)] <-scores1$allscores$VarDiff                     
          MainOutput$SkewDiff      [index:(index+reps-1)] <-scores1$allscores$SkewDiff                     
          MainOutput$KurtDiff      [index:(index+reps-1)] <-scores1$allscores$KurtDiff                                    
}

#     # run the function for every n, hold eps and dSize fixed
    for(i in 1:length(nlist)){
    print(paste("_______this is nlist",i))

index <- index+(1*reps)
     dSize <- ceiling(((range[2]-range[1])/1)+1)
     MainOutput$imp  [index:(index+reps-1)] <- Fnameslist[f] 
     MainOutput$eps  [index:(index+reps-1)] <- epslist[1]  
     MainOutput$n    [index:(index+reps-1)] <- nlist[i]
     MainOutput$dSize[index:(index+reps-1)] <- dSize
     MainOutput$data [index:(index+reps-1)] <- Dnameslist[d] 
     MainOutput$gran [index:(index+reps-1)] <- granlist[1] 
     MainOutput$min  [index:(index+reps-1)] <- range[1]
     MainOutput$max  [index:(index+reps-1)] <- range[2]
     MainOutput$seed [index:(index+reps-1)] <- seedlist[index]

          scores2<- CDFtestTrackx(funct = functlist[[f]],
                                eps= epslist[1],
                                data=datalist[[d]], 
                                range, 
                                gran <- granlist[1],
                                samplesize = nlist[i],
                                reps,
                                SmoothAll,
                                ExtraTests_CDF = ExtraTests_CDF,
                                ExtraTests_PDF = ExtraTests_PDF
                                )
       #   MainOutput$Mean          [index:(index+reps-1)] <-scores2$allscores$Mean
          MainOutput$MSEanalytic    [index:(index+reps-1)] <-scores2$allscores$MSEanalytic
          MainOutput$SDempiric     [index:(index+reps-1)] <-scores2$allscores$SDempiric
          MainOutput$MaxError_CDF  [index:(index+reps-1)] <-scores2$allscores$MaxError_CDF 
          MainOutput$MaxErrorAt_CDF[index:(index+reps-1)] <-scores2$allscores$MaxErrorAt_CDF
          MainOutput$diffat25      [index:(index+reps-1)] <-scores2$allscores$diffat25
          MainOutput$diffatMedian  [index:(index+reps-1)] <-scores2$allscores$diffatMedian 
          MainOutput$diffat75      [index:(index+reps-1)] <-scores2$allscores$diffat75
          MainOutput$horzdiffat25  [index:(index+reps-1)] <-scores2$allscores$horzdiffat25
          MainOutput$horzdiffatMed [index:(index+reps-1)] <-scores2$allscores$horzdiffatMed 
          MainOutput$horzdiffat75  [index:(index+reps-1)] <-scores2$allscores$horzdiffat75
          MainOutput$Medians       [index:(index+reps-1)] <-scores2$allscores$Medians
          MainOutput$MAE_CDF       [index:(index+reps-1)] <-scores2$allscores$MAE_CDF
          MainOutput$MSE_CDF       [index:(index+reps-1)] <-scores2$allscores$MSE
       #   MainOutput$DerivScore    [index:(index+reps-1)] <-scores2$allscores$DerivScore
          MainOutput$MaxError_PDF  [index:(index+reps-1)] <-scores2$allscores$MaxError_PDF 
          MainOutput$MaxErrorAt_PDF[index:(index+reps-1)] <-scores2$allscores$MaxErrorAt_PDF
          MainOutput$MAE_PDF       [index:(index+reps-1)] <-scores2$allscores$L1_area_PDF
          MainOutput$MSEanalytic   [index:(index+reps-1)] <-scores2$allscores$MSEanalytic
          MainOutput$MeanDiff      [index:(index+reps-1)] <-scores2$allscores$MeanDiff
          MainOutput$ModeDiff      [index:(index+reps-1)] <-scores2$allscores$ModeDiff                     
          MainOutput$StdDiff       [index:(index+reps-1)] <-scores2$allscores$StdDiff                     
          MainOutput$VarDiff       [index:(index+reps-1)] <-scores2$allscores$VarDiff                     
          MainOutput$SkewDiff      [index:(index+reps-1)] <-scores2$allscores$SkewDiff                     
          MainOutput$KurtDiff      [index:(index+reps-1)] <-scores2$allscores$KurtDiff      
           
}

    # run the function for every dSize, hold eps and n fixed
    for(i in 1:length(granlist)){
    print(paste("_______this is granlist",i))


index <- index+(1*reps)
     dSize <- ceiling(((range[2]-range[1])/granlist[i])+1)
     MainOutput$imp  [index:(index+reps-1)] <- Fnameslist[f] 
     MainOutput$eps  [index:(index+reps-1)] <- epslist[1]  
     MainOutput$n    [index:(index+reps-1)] <- nlist[1]
     MainOutput$dSize[index:(index+reps-1)] <- dSize
     MainOutput$data [index:(index+reps-1)] <- Dnameslist[d] 
     MainOutput$gran [index:(index+reps-1)] <- granlist[i]
     MainOutput$min  [index:(index+reps-1)] <- range[1]
     MainOutput$max  [index:(index+reps-1)] <- range[2]
     MainOutput$seed [index:(index+reps-1)] <- seedlist[index]

          scores3<- CDFtestTrackx(funct = functlist[[f]],
                                eps= epslist[1],
                                data=datalist[[d]], 
                                range, 
                                gran= granlist[i], 
                                samplesize = nlist[1],
                                reps,
                                SmoothAll,
                                ExtraTests_CDF = ExtraTests_CDF,
                                ExtraTests_PDF = ExtraTests_PDF
                                )
       #   MainOutput$Mean          [index:(index+reps-1)] <-scores3$allscores$Mean
          MainOutput$MSEanalytic    [index:(index+reps-1)] <-scores3$allscores$MSEanalytic
          MainOutput$SDempiric     [index:(index+reps-1)] <-scores3$allscores$SDempiric
          MainOutput$MaxError_CDF  [index:(index+reps-1)] <-scores3$allscores$MaxError_CDF 
          MainOutput$MaxErrorAt_CDF[index:(index+reps-1)] <-scores3$allscores$MaxErrorAt_CDF
          MainOutput$diffat25      [index:(index+reps-1)] <-scores3$allscores$diffat25
          MainOutput$diffatMedian  [index:(index+reps-1)] <-scores3$allscores$diffatMedian 
          MainOutput$diffat75      [index:(index+reps-1)] <-scores3$allscores$diffat75
          MainOutput$horzdiffat25  [index:(index+reps-1)] <-scores3$allscores$horzdiffat25
          MainOutput$horzdiffatMed [index:(index+reps-1)] <-scores3$allscores$horzdiffatMed 
          MainOutput$horzdiffat75  [index:(index+reps-1)] <-scores3$allscores$horzdiffat75
          MainOutput$Medians       [index:(index+reps-1)] <-scores3$allscores$Medians
          MainOutput$MAE_CDF       [index:(index+reps-1)] <-scores3$allscores$MAE_CDF
          MainOutput$MSE_CDF       [index:(index+reps-1)] <-scores3$allscores$MSE
        #  MainOutput$L1empiric     [index:(index+reps-1)] <-scores3$allscores$L1emp
        #  MainOutput$L2empiric     [index:(index+reps-1)] <-scores3$allscores$L2emp
        #  MainOutput$DerivScore    [index:(index+reps-1)] <-scores3$allscores$DerivScore
          MainOutput$MaxError_PDF  [index:(index+reps-1)] <-scores3$allscores$MaxError_PDF 
          MainOutput$MaxErrorAt_PDF[index:(index+reps-1)] <-scores3$allscores$MaxErrorAt_PDF
          MainOutput$MAE_PDF       [index:(index+reps-1)] <-scores3$allscores$L1_area_PDF
          MainOutput$MSEanalytic   [index:(index+reps-1)] <-scores3$allscores$MSEanalytic
          MainOutput$MeanDiff      [index:(index+reps-1)] <-scores3$allscores$MeanDiff
          MainOutput$ModeDiff      [index:(index+reps-1)] <-scores3$allscores$ModeDiff                     
          MainOutput$StdDiff       [index:(index+reps-1)] <-scores3$allscores$StdDiff                     
          MainOutput$VarDiff       [index:(index+reps-1)] <-scores3$allscores$VarDiff                     
          MainOutput$SkewDiff      [index:(index+reps-1)] <-scores3$allscores$SkewDiff                     
          MainOutput$KurtDiff      [index:(index+reps-1)] <-scores3$allscores$KurtDiff      
         
}
}}
MainOutput$comments<- comments

########################## BEGIN PRINTING ##############################
 
  ##### Excel Files w/ means and medians of diagnostic results #####

            t <- Sys.time()
            timestamp <-as.POSIXlt(t)
            timeindex <-strftime(t, "%Y%m%d%H%M%S")

            write.table(MainOutput,
                  paste( OutputDirectory,"/DataCollector",timeindex,".csv", sep =""), row.names=FALSE, col.names=TRUE,  sep=",")

            maxlength <- max(length(datalist), length(epslist), length(functlist), length(range), length(granlist), length(nlist), length(setseed))

            forchart= data.frame( "index"         = c(timeindex, rep(" ",times= maxlength-1)),
                                  "functions"     = c(Fnameslist,rep(" ",times= maxlength-(length(Fnameslist)))),
                                  "datasets"      = c(Dnameslist,rep(" ",times= maxlength-(length(Dnameslist)))),
                                  "epsilons"      = c(epslist,   rep(" ",times= maxlength-(length(epslist)))),
                                  "range"         = c(range,     rep(" ",times= maxlength-2)),
                                  "granularities" = c(granlist,  rep(" ",times= maxlength-length(granlist))),
                                  "nlist"         = c(nlist,     rep(" ",times= maxlength-length(nlist))),
                                  "cdfstep"       = c(cdfstep,   rep(" ",times= maxlength-1)),               
                                  "repetitions"   = c(reps,      rep(" ",times= maxlength-1)),
                                  "SmoothAll"     = c(SmoothAll, rep(" ",times= maxlength-1)),
                                  "setseed"       = c(setseed,   rep(" ",times= maxlength-length(setseed))), 
                                  "comments"      = c(comments,  rep(" ",times= maxlength-1))
                                 )          
            suppressWarnings(
            write.table(forchart,   
                 paste(OutputDirectory,"/CDFtestindexchart.csv", sep =""),    row.names=FALSE, col.names=TRUE,  sep=",", append=TRUE)
            )


return(MainOutput)
}else{

##################################################################################################################################

# this takes the given seedlist
seedlist <- setseed  

# if you're looking for a random seed, this will overwrite the 
suppressWarnings

  if (setseed[1] == -100){
          seedlist <- round(runif((length(epslist)*length(functlist)* length(datalist)), min=0, max=1000),0)
            }


set.seed(seedlist)



universe         <- seq(from=(range[1]), to=(range[2]), by=gran)
databins         <- seq(from=(range[1]), to=(range[2]), by=cdfstep)
yourCDF          <- matrix(nrow=length(databins), ncol= (length(epslist)*length(functlist)* length(datalist)))
realCDF          <- matrix(nrow=length(databins), ncol=length(datalist))

 if (EmpiricBounds){
        yourCDF_LboundE  <- matrix(nrow=length(databins), ncol= (length(epslist)*length(functlist)* length(datalist)))
        yourCDF_UboundE  <- matrix(nrow=length(databins), ncol= (length(epslist)*length(functlist)* length(datalist)))
    }
 if (AnalyticBounds|AnalyticProbSleeve){
        yourCDF_LboundA  <- matrix(nrow=length(databins), ncol= (length(epslist)*length(functlist)* length(datalist)))
        yourCDF_UboundA  <- matrix(nrow=length(databins), ncol= (length(epslist)*length(functlist)* length(datalist)))
    }

yourPDF        <- matrix(nrow=length(databins), ncol= (length(epslist)*length(functlist)* length(datalist)))
realPDF        <- matrix(nrow=length(databins), ncol=length(datalist))

# run the function for every dataset
for(d in 1:length(datalist)){


datalist[[d]] <- as.numeric(lapply(datalist[[d]], function(x) MovetoRange(x, range)))
print(paste("this is datalist",d))

    # run the function for every epsilon
    for(i in 1:length(epslist)){
    print(paste("this is epslist",i))

        # run the function for every CDF function
        for (f in 1:length(functlist)){
        print(paste("this is functlist",f))

            Z <- FALSE
               if(AnalyticBounds|AnalyticProbSleeve){
            Z <- TRUE
   } 
          scores<- CDFtestTrack(funct = functlist[[f]],
                                eps= epslist[i],
                                cdfstep, 
                                data=datalist[[d]], 
                                range, 
                                gran, 
                                samplesize = samplesize,
                                reps,
                                SmoothAll,
                                ABounds = Z,
                                EmpiricBounds   = EmpiricBounds,
                                ExtraTests_CDF, 
                                ExtraTests_PDF
                                )
            # extract everything from the output list. 
                meansBase <- data.frame (  blank    = 0, 
                                           blank    = 0, 
                                           blank    = 0, 
                                           blank    = 0,
                                           blank    = 0,
                                           blank    = 0,
                                           blank    = 0,
                                           blank    = 0
                                     )
                mediansBase <- data.frame (blank    = 0, 
                                           blank    = 0, 
                                           blank    = 0, 
                                           blank    = 0,
                                           blank    = 0,
                                           blank    = 0,
                                           blank    = 0,
                                           blank    = 0
                                     ) 

                scores$meanscores   <- c(meansBase,   scores$meanscores  )
                scores$medianscores <- c(mediansBase, scores$medianscores)
          if ( (f + d + i) == 3){   

                means <- data.frame (  stat          = 0, 
                                       seed          = 0, 
                                       dataset       = 0, 
                                       epsilon       = 0, 
                                       functionName  = 0,
                                       n             = 0,
                                       domainSize    = 0,
                                       minValue      = 0,
                                       maxValue      = 0
                                     )
                medians <- data.frame (stat          = 0, 
                                       seed          = 0, 
                                       dataset       = 0, 
                                       epsilon       = 0, 
                                       functionName  = 0,
                                       n             = 0,
                                       domainSize    = 0,
                                       minValue      = 0,
                                       maxValue      = 0
                                     )

                    # create objects to store outputs
                    for (p in 1:length(scores$TestPack_CDF)){
                    means[p+9] <- 0
                    colnames(means)[p+9]   <- names(scores$TestPack_CDF[p])
                    medians[p+9] <- 0
                    colnames(medians)[p+9] <- names(scores$TestPack_CDF[p])
                    }

                    # create objects to store outputs 
                    for (p in 1:length(scores$TestPack_PDF)){
                    means[p+9+length(scores$TestPack_CDF)]   <- 0
                    colnames(means)[p+9+length(scores$TestPack_CDF)]   <- names(scores$TestPack_PDF[ p ])
                    medians[p+9+length(scores$TestPack_CDF)] <- 0
                    colnames(medians)[p+9+length(scores$TestPack_CDF)]   <- names(scores$TestPack_PDF[ p ])
                    }

          }

      index <-  ((d-1)*length(epslist)*length(functlist) + (i-1)*length(functlist)+f)

            means  [index, ]     <- scores$meanscores
            means  [index,9]     <- range[2]
            means  [index,8]     <- range[1]
            means  [index,7]     <- length(seq(range[1], range[2],gran))
            means  [index,6]     <- length(datalist[[d]])
            means  [index,5]     <- Fnameslist[f]
            means  [index,4]     <- epslist   [i]
            means  [index,3]     <- Dnameslist[d]
            means  [index,2]     <- seedlist[index]
            means  [index,1]     <- "means"

            medians  [index, ]  <- scores$medianscores
            medians  [index,9]     <- range[2]
            medians  [index,8]     <- range[1]
            medians  [index,7]     <- length(seq(range[1], range[2],gran))
            medians  [index,6]     <- length(datalist[[d]])
            medians  [index,5]     <- Fnameslist[f]
            medians  [index,4]     <- epslist   [i]
            medians  [index,3]     <- Dnameslist[d]
            medians  [index,2]     <- seedlist[index]
            medians  [index,1]     <- "medians"
      
           # round each numeric value in the vector (for asthetics, mostly)
             for (k in 10:(length(means[i,]))){             
            means[,k] <- round(means[,k],3)
            medians[,k] <- round(medians[,k],3)
            }

            yourCDF        [,index] <- scores$yourCDFoutput   # this one will be used for plotting
            if (EmpiricBounds){ 
                yourCDF_LboundE[,index] <- scores$yourCDF_LboundE # this one will be used for plotting (empirical lower bound)
                yourCDF_UboundE[,index] <- scores$yourCDF_UboundE # this one will be used for plotting
            }
            if (AnalyticBounds|AnalyticProbSleeve){    
                yourCDF_LboundA[,index] <- scores$yourCDF_LboundA # this one will be used for plotting (analytic Lower bound)
                yourCDF_UboundA[,index] <- scores$yourCDF_UboundA # this one will be used for plotting
            }

            realCDF [,d]           <- scores$realCDFoutput
            yourPDF [,index]       <- scores$yourPDFoutput   # this one will be used for plotting
            realPDF [,d]           <- scores$realPDFoutput
            databins               <- scores$databins

            # store for out. Add the seed as well, for reproducibility.       
            masterlist <- list()
            masterlist$allscores <-scores$allscores
            masterlist$seed      <-seedlist
            masterlist$means     <- means
            masterlist$medians   <- medians

            foroutput          <- scores$allscores
            dataset            <- rep(Dnameslist[d], reps)
            epsilon            <- rep(epslist   [i], reps)
            functName          <- rep(Fnameslist[f], reps)
            foroutput$dataset      <- dataset
            foroutput$epsilon      <- epsilon
            foroutput$functionname <- functName
            foroutput <- data.frame(foroutput)
            if(d+i+f ==3){
              container  <-foroutput
              }

            else{                                 
            container <- rbind(container, foroutput, make.row.names = FALSE)
            }
         } # ends the function forloop
        } # ends the epsilon forloop
       } # ends the dataset forloop

########################## BEGIN PRINTING ##############################
 
  ##### Excel Files w/ means and medians of diagnostic results #####

            t <- Sys.time()
            timestamp <-as.POSIXlt(t)
            timeindex <-strftime(t, "%Y%m%d%H%M%S")
            maxlength <- max(length(datalist), length(epslist), length(functlist), length(range), length(setseed))

             forchart= data.frame("index"          = c(timeindex, rep(" ",times= maxlength-1)),
                                  "functions"      = c(Fnameslist,rep(" ",times= maxlength-(length(Fnameslist)))),
                                  "datasets"       = c(Dnameslist,rep(" ",times= maxlength-(length(Dnameslist)))),
                                  "epsilons"       = c(epslist,   rep(" ",times= maxlength-(length(epslist)))),
                                  "range"          = c(range,     rep(" ",times= maxlength-2)),
                                  "granularity"    = c(gran,      rep(" ",times= maxlength-1)),
                                  "cdfstep"        = c(cdfstep,   rep(" ",times= maxlength-1)),
                                  "samplesize"     = c(samplesize,rep(" ",times= maxlength-1)),
                                  "repetitions"    = c(reps,      rep(" ",times= maxlength-1)),
                                  "SmoothAll"         = c(SmoothAll,            rep(" ",times= maxlength-1)),
                                  "EmpiricBounds"     = c(EmpiricBounds,        rep(" ",times= maxlength-1)),
                                  "AnalyticBounds"    = c(AnalyticBounds,       rep(" ",times= maxlength-1)),
                                  "AnalyticProbSleeve"= c(AnalyticProbSleeve,   rep(" ",times= maxlength-1)), 
                                  "SuppressRealCDF"   = c(SuppressRealCDF,      rep(" ",times= maxlength-1)), 
                                  "ExtraTests_CDF"    = c(names(ExtraTests_CDF),rep(" ",times= maxlength-length(names(ExtraTests_CDF)))), 
                                  "ExtraTests_PDF"    = c(names(ExtraTests_PDF),rep(" ",times= maxlength-length(names(ExtraTests_PDF)))), 
                                  "setseed"           = c(setseed,              rep(" ",times= maxlength-length(setseed))), 
                                  "comments"          = c(comments,             rep(" ",times= maxlength-1))
                                 )          

            write.table(masterlist[[3]],
                  paste( OutputDirectory,"/CDFtesting",timeindex,".csv", sep=""), row.names=FALSE, col.names=TRUE,  sep=",")

            suppressWarnings(
            write.table(masterlist[[4]],
                  paste( OutputDirectory,"/CDFtesting",timeindex,".csv", sep=""), row.names=FALSE, col.names=TRUE,  sep=",", append=TRUE)
            )

            suppressWarnings(
            write.table(forchart,   
                 paste(OutputDirectory,"/CDFtestindexchart.csv", sep=""),    row.names=FALSE, col.names=TRUE,  sep=",", append=TRUE)
            )
########################################################################################
##############################################################
##### PDF file with graphs each showing:
######## REAL CDF,
######## 1 DP CDF,
######## Optional Analytical Bounds and Empirical Bounds ####

## make the PDF file
g <- paste(OutputDirectory,"/CDFtesting",timeindex,".pdf", sep="")
pdf(file =g, onefile = TRUE, paper = "a4r", width = 0, height = 0)   

for(d in 1:length(datalist)){
 par(mfrow =c(length(functlist),length(epslist))) 

  if (length(functlist) == 1 ){
   par(mar= c(2,2,4,1))
  }else{
     par(mar= c(0,2,2,1))
  }

  for (f in 1:length(functlist)){ 

    for (i in 1:length(epslist)){
        
         index <-  ((d-1)*length(epslist)*length(functlist) + (i-1)*length(functlist)+f)
            
            #plot the original CDF for each dataset
          if(SuppressRealCDF){
            fakeline <- rep(1.1, times = (length(databins)))
            if (f==1 ) {
                plot(databins,
                fakeline,
                col  = "white",
                lwd  = 6, 
                type = "l", 
                ylim = c(0, 1.05),
                ylab = "cumulative proportion", 
                xlab = " "
                )  
              title(main = " ")  
                }
            else{
              if (f== length(functlist)){
                plot(databins,
                fakeline,
                col  = "white",
                lwd  = 6, 
                type = "l", 
                ylim = c(0, 1.05),
                ylab = "cumulative proportion",
                xlab = " " 
                )
                 title(main =".",
                sub  = "Functions are shown by line thickness, and epsilon is shown by color.",
                )  
              }
              else{
               plot(databins,
                fakeline,
                col  = "white",
                lwd  = 6, 
                type = "l", 
                ylim = c(0, 1.05),
                ylab = "cumulative proportion",
                xlab= " "
                )
              }
              }   
          }else {

            if (f==1 ) {
                plot(databins,
                realCDF[,d],
                col  = "black",
                lwd  = 6, 
                type = "l", 
                ylim = c(0, 1.05),
                ylab = "cumulative proportion", 
                xlab = " "
                )  
              title(main = " ")  
                }
            else{
              if (f== length(functlist)){
                plot(databins,
                realCDF[,d],
                col  = "black",
                lwd  = 6, 
                type = "l", 
                ylim = c(0, 1.05),
                ylab = "cumulative proportion",
                xlab = "binned data" 
                )
                 title(main =".",
                sub  = "_",
                )  
              }
              else{
               plot(databins,
                realCDF[,d],
                col  = "black",
                lwd  = 6, 
                type = "l", 
                ylim = c(0, 1.05),
                ylab = "cumulative proportion",
                xlab= " "
                )
              }
              }   
            }
# lines(databins, yourCDF        [,index],  col=colors()[26 + 10*((i-1)%%4)],  type="l", lwd=2)

      lines(databins, rep(.25,length(databins)), col="gray", lwd =(1), type="l")
      lines(databins, rep(.50,length(databins)), col="gray", lwd =(1), type="l")
      lines(databins, rep(.75,length(databins)), col="gray", lwd =(1), type="l")


 if (EmpiricBounds){
      lines(databins, yourCDF_LboundE[,index],  col=colors()[26 + 10*((i-1)%%4)],  type="l", lwd=1)
      lines(databins, yourCDF_UboundE[,index],  col=colors()[26 + 10*((i-1)%%4)],  type="l", lwd=1)
 }

## plot the curves
 if (AnalyticBounds){
      lines(databins, yourCDF_LboundA[,index],  col="black", type="l", lwd=1)
      lines(databins, yourCDF_UboundA[,index],  col="black", type="l", lwd=1)
 }

if (AnalyticProbSleeve){
  if(Z){
    SleeveVector <- seq(1,100, 1)
    for(v in 1:100){
       probsleeveLine <- c()

          for(w in 1: length(yourCDF[,index])){
            probsleeveLine[w] <- abs((yourCDF[,index][w])-(yourCDF_LboundA[,index][w]))  /max(SleeveVector) * SleeveVector[v]
          }
      lines(databins, yourCDF        [,index]+ probsleeveLine, col= adjustcolor(colors()[26 + 10*((i-1)%%4)], alpha.f=0),  type="n", lwd=0)
      lines(databins, yourCDF        [,index]- probsleeveLine, col= adjustcolor(colors()[26 + 10*((i-1)%%4)], alpha.f=0),  type="n", lwd=0)
      
      polygon(c(databins, rev(databins)), c((yourCDF [,index]+ probsleeveLine), rev(yourCDF [,index]- probsleeveLine)),
      col = adjustcolor(colors()[26 + 10*((i-1)%%4)], alpha.f=0.01), border = NA)
      # polygon(c(databins, rev(databins)), c((yourCDF [,index]+ probsleeveLine), rev(yourCDF [,index]- probsleeveLine)),
      # col = adjustcolor(colors()[555], alpha=0.01), border = NA)

    }
  
  }else{
  SleeveVector <- rexp(100,40)
  sort(SleeveVector)
  for(v in 1:40){
      lines(databins, yourCDF        [,index]+SleeveVector[v], col= adjustcolor(colors()[26 + 10*((i-1)%%4)], alpha.f=0.017),  type="n", lwd=15)
      lines(databins, yourCDF        [,index]-SleeveVector[v], col= adjustcolor(colors()[26 + 10*((i-1)%%4)], alpha.f=0.017),  type="n", lwd=15)
        }
      }      
}

if(SuppressDPCDF == FALSE){
  lines(databins, yourCDF        [,index],  col=colors()[26 + 10*((i-1)%%4)],  type="l", lwd=2)
#lines(databins, yourCDF        [,index],  col=colors()[555],  type="l", lwd=2)
#lines(databins, yourCDF        [,index],  col="black",  type="l", lwd=6)
}

if(SuppressLegends == FALSE){
     # add the legend

          #shrink ranges if necessary
          printRange <- c()
          printRange[1] <- Abbrev(range[1])
          printRange[2] <- Abbrev(range[2])
          printgran     <- Abbrev(gran)
          printcdfstep  <- Abbrev(cdfstep)

         # Determine plot boundaries, in units of the data
              xmin <- par("usr")[1]
              xmax <- par("usr")[2]
              ymin <- par("usr")[3]
              ymax <- par("usr")[4]
               
         # Now determine the size of the legend you would like to plot.  
              lgd <- legend(x = mean(c(xmin,xmax)), y =  mean(c(ymin,ymax)),
              c(paste("d:",Dnameslist[d]), 
                             paste("f: ",Fnameslist[f]),
                             paste("e:",epslist[i]),
                             paste("r: ",range[1],",",range[2]),
                             paste("gran: ", printgran),
                             paste("binSize: ",printcdfstep)),  
              pch = rep(1, times=6),
              col = rep(colors()[26 + 10*((i-1)%%4)], times=6),
              bg ="white",
              plot = F)

              # Add legend in the lower right corner:
              legend(x = xmax - lgd$rect$w, y =  ymin + lgd$rect$h,
              c(paste("d:",Dnameslist[d]), 
                paste("f: ",Fnameslist[f]),
                paste("e:",epslist[i]),
                paste("r: ",printRange[1],",",printRange[2]),
                paste("gran: ", printgran),
                paste("binSize: ",printcdfstep)), 
              pch = rep(1, times=6),
              col = rep(colors()[26 + 10*((i-1)%%4)], times=6),
              bg ="white",
              plot = T)

       text(c(xmin*1.1,ymax*.95), "CDF", font=4)

         }
        }
       }
      }

#####################################################################
#####  NOW DO THE SAME FOR PROB. DENSE. FUNCTIONS ###
for(d in 1:length(datalist)){
 par(mfrow =c(length(functlist),length(epslist))) 
 par(mar= c(0,2,2,2))

  for (f in 1:length(functlist)){ 

    for (i in 1:length(epslist)){
        
       index <-  ((d-1)*length(epslist)*length(functlist) + (i-1)*length(functlist)+f)
           mx <- max(realPDF)*2


            #plot the original curve for each dataset
            if (f==1 ) {
              
                plot(databins,
                realPDF[,d],
                col  = 1,
                lwd  = 6, 
                type = "l", 
                ylim = c(0,.35),
                ylab = "proportion", 
                xlab = " "
                )  
              title(main = " ")  
                }
            else{
              if (f== length(functlist)){
                plot(databins,
                realPDF[,d],
                col  = 1,
                lwd  = 6, 
                type = "l", 
                ylim = c(0,.35),
                ylab = "proportion",
                xlab = "binned data" 
                )
                 title(main =" ",
                sub  = "Functions are shown by line thickness, and epsilon is shown by color.",
                )  
              }
              else{
               plot(databins,
                realPDF[,d],
                col  = 1,
                lwd  = 6, 
                type = "l", 
                ylim = c(0,.35),
                ylab = "proportion",
                xlab= " "
                )
              }
              }   
if (SuppressDPCDF == FALSE){

tool <- ceiling(length(databins)/50)
LowResDPPDF <- c()
z <-1

  for (w in 1:50){
    if ((z-1+tool)<length(databins)){
  LowResDPPDF[z:(z-1+tool)] <- mean(yourPDF[w:(w+1),index])
  z<- z+tool
  }
}
if(length(LowResDPPDF)<length(databins)){
  LowResDPPDF[(length(LowResDPPDF)+1):length(databins)] <- LowResDPPDF[length(LowResDPPDF)]
}

             lines(databins, LowResDPPDF,  col=colors()[26 + 10*((i-1)%%4)],  type="l", lwd=2)
     }
      # add the legend


if (SuppressLegends==FALSE){
              # Determine plot boundaries, in units of the data
              xmin <- par("usr")[1]
              xmax <- par("usr")[2]
              ymin <- par("usr")[3]
              ymax <- par("usr")[4]
               
              #Now determine the size of the legend you would like to plot.  Right now the exact
              #location is not important, we just want to know the dimension!  Note that we are
              # treating the lengend as a variable and we are NOT plotting the legend on the figure!
              lgd <- legend(x = mean(c(xmin,xmax)), y =  mean(c(ymin,ymax)),
              c(paste("d:",Dnameslist[d]), 
                             paste("f: ",Fnameslist[f]),
                             paste("e:",epslist[i]),
                             paste("r: ",printRange[1],",",printRange[2]),
                             paste("gran: ",printgran),
                             paste("binSize: ",printcdfstep)), 
              pch = rep(1, times=6),
              col = rep(colors()[26 + 10*((i-1)%%4)], times=6),
              bg ="white",
              plot = F)

              # Add legend in the lower right corner:
              legend(x = xmax - lgd$rect$w, y =  ymax,
              c(paste("d:",Dnameslist[d]), 
                paste("f: ",Fnameslist[f]),
                paste("e:",epslist[i]),
                paste("r: ",printRange[1],",",printRange[2]),
                             paste("gran: ",printgran),
                             paste("binSize: ",printcdfstep)),  
              pch = rep(1, times=6),
              col = rep(colors()[26 + 10*((i-1)%%4)], times=6),
              bg ="white",
              plot = T)

       text(c(xmin*1.1,ymax*.95), "PDF", font=4)

            }
           }
          }
         }

################################################
# lastly, create some boxplots representing errors
################################################
par(mfrow= c(1,1))
permetric <- list()
for (p in 1:(length(container)-3)){
    forplotting <- data.frame("1"= reps, stringsAsFactors=FALSE)

    nextstep <- data.frame("epsilon" = container$epsilon, "dataset" = container$dataset, "functionName" = container$functionname, "metric"=container[[p]])
    nextstep$dataset      <- factor(nextstep$dataset)
    nextstep$functionname <- factor(nextstep$functionName)
     
    forplots<- list()
    for(i in 1:length(epslist)){
        forplotting <- data.frame("1"= c(seq(from=1, to=reps,by=1)), stringsAsFactors=FALSE)
        for (d in 1:length(Dnameslist)){
          for (f in 1:length(Fnameslist)){

                A   <- subset(nextstep, epsilon == epslist[i])
                AA  <- subset(A,        dataset == Dnameslist[d])
                AAA <- subset(AA,  functionname == Fnameslist[f])
                forplotting <- cbind(forplotting,AAA$metric)
        }}
        forplotting[,1] <-NULL
        forplots[[i]]<-forplotting
    }

    permetric[[p]] <- forplots
}
masterlist$permetric <- permetric


#######################################################33

for (p in 1:length(permetric)){

    forplots <- permetric[[p]]
    for (i in 1:length(epslist)){
        forplotting <- forplots[[i]]

    #set up the colors of the boxplot columns
        cols <- numeric(0)
        for (f in 1:length(functlist)){
             cols[length(cols)+1] <- f+1
          }
        cols <- rep(cols, times=length(datalist))

    #set up the locations of the boxplot columns (function sequence grouped by dataset)
    # if 3 functions and 3 datasets, creates |F1_D1| |F2_D1| |F3_D1| _____ |F1_D2| |F2_D2| |F3_D2|
        ats <- numeric(0)
        for (d in 1:length(datalist)){
             for (f in 1:length(functlist)){
                  ats[length(ats)+1] <- (d-1)*(length(functlist)+2)+f+1
          }}

        namer <- rep("", times= length(datalist)*length(functlist))
        mx <- 1.25*max(forplotting)

################################################
        labellist<-list()
        labellist[1]  <- paste("Maximum Vertical Error ('Linf norm')\n epsilon=",         epslist[i],",",reps,"repetitions")
        labellist[2]  <- paste("Location of Maximum Vertical Error\n epsilon=",           epslist[i],",",reps,"repetitions")
        labellist[3]  <- paste("Vertical Error at the True 25th Percentile\n epsilon=",   epslist[i],",",reps,"repetitions")
        labellist[4]  <- paste("Vertical Error at the True Median\n epsilon=",            epslist[i],",",reps,"repetitions")
        labellist[5]  <- paste("Vertical Error at the True 75th Percentile\n epsilon=",   epslist[i],",",reps,"repetitions")
        labellist[6]  <- paste("Horizontal Error at the True 25th Percentile\n epsilon=", epslist[i],",",reps,"repetitions")
        labellist[7]  <- paste("Horizontal Error at the True Median\n epsilon=",          epslist[i],",",reps,"repetitions")
        labellist[8]  <- paste("Horizontal Error at the True 75th Percentile\n epsilon=", epslist[i],",",reps,"repetitions")
        labellist[9]  <- paste("Medians Returned\n epsilon=",                             epslist[i],",",reps,"repetitions")
        labellist[10] <- paste("Mean Absolute Error (emp, L1/DomainSize)\n epsilon=",     epslist[i],",",reps,"repetitions")
        labellist[11] <- paste("Mean Squared Error (emp (L2^2)/DomainSize)\n epsilon=",   epslist[i],",",reps,"repetitions")
        labellist[12] <- paste("Derivative Score (higher values are 'better')\n epsilon=", epslist[i],",",reps,"repetitions")


        if (length(ExtraTests_CDF) >0){
        assist <- names(ExtraTests_CDF)              
             TestNames <- names(scores$TestPack_CDF)
             for (L in (length(labellist)+1):(length(labellist)+length(assist))){
                  labellist[L] <- paste("Custom:",assist[L-length(labellist)], "\n epsilon=", epslist[i],",",reps,"repetitions" )
             }
        }
TestPack_PDF <-list(MaxError_PDF   = MaxError_PDF,
                    MaxErrorAt_PDF = MaxErrorAt_PDF,
                    MAE_PDF        = MAE,
                    MSE_PDF        = MSE,
                    MeanDiff = MeanDiffpdf,
                     ModeDiff = ModeDiffpdf,
                      StdDiff = StdDiffpdf,
                       VarDiff = VarDiffpdf,
                        SkewDiff = SkewDiffpdf,
                         KurtDiff= KurtDiffpdf)


        labellist[length(labellist)+1] <- paste("Maximum Vertical Error (PDF)\n epsilon=",            epslist[i],",",reps,"repetitions")
        labellist[length(labellist)+1] <- paste("Location of Maximum Vertical Error (PDF)\n epsilon=",epslist[i],",",reps,"repetitions")
        labellist[length(labellist)+1] <- paste("Mean Absolute Error (PDF)\n epsilon=",               epslist[i],",",reps,"repetitions")
        labellist[length(labellist)+1] <- paste("Mean Squared Error (PDF)\n epsilon=",                epslist[i],",",reps,"repetitions")
        labellist[length(labellist)+1] <- paste("Difference between Means given by CDFs\n epsilon=",                epslist[i],",",reps,"repetitions")
        labellist[length(labellist)+1] <- paste("Difference between Modes given by CDFs\n epsilon=",                epslist[i],",",reps,"repetitions")
        labellist[length(labellist)+1] <- paste("Difference between Standard Deviation given by CDFs\n (under development, see ?CDFtest) epsilon=",                epslist[i],",",reps,"repetitions")
        labellist[length(labellist)+1] <- paste("Difference between Variances given by CDFs\n epsilon=",                epslist[i],",",reps,"repetitions")
        labellist[length(labellist)+1] <- paste("Difference between Skewnesses given by CDFs \n (under development, see ?CDFtest) epsilon=",                epslist[i],",",reps,"repetitions")
        labellist[length(labellist)+1] <- paste("Difference between Kurtoses given by CDFs \n (under development, see ?CDFtest) epsilon=",                epslist[i],",",reps,"repetitions")
 
        if (length(ExtraTests_CDF) >0){
        assist <- names(ExtraTests_PDF)
             TestNames <- names(scores$TestPack_PDF)
             for (L in (length(labellist)+1):(length(labellist)+length(assist))){
                   labellist[L] <- paste(assist[L-length(labellist)], "\n epsilon=", epslist[i],",",reps,"repetitions" )
             }
        }

        ylablist <-list()
        ylablist[1]   <- paste("Vertical Distance between CDFs at maximum distance point")
        ylablist[2]   <- paste("Databin (X-axis) Values")
        ylablist[3]   <- paste("Distance (Error) Values")
        ylablist[4]   <- paste("Distance (Error) Values")
        ylablist[5]   <- paste("Distance (Error) Values")
        ylablist[6]   <- paste("Distance between 25th percentile values")
        ylablist[7]   <- paste("Distance between median values")
        ylablist[8]   <- paste("Distance between 75th percentile values")
        ylablist[9]   <- paste("Median values") 
        ylablist[10]  <- paste("Average distance between real and estimated CDF")
        ylablist[11]  <- paste("SqRt of average squared distance between real and estimated CDF")
        ylablist[12]  <- paste("Normalized difference in derivatives across gran resolutions")


       if (length(ExtraTests_CDF) >0){
        assist <- names(ExtraTests_CDF)
   
             TestNames <- names(scores$TestPack_CDF)
             for (L in (length(ylablist)+1):(length(ylablist)+length(assist))){
                  ylablist[L] <- paste("Custom: user-defined units (CDF)" )
             }
        }




        ylablist[length(ylablist)+1]  <- paste("Vertical Distance between PDFs at maximum distance point")
        ylablist[length(ylablist)+1]  <- paste("Databin (X-axis) Values")
        ylablist[length(ylablist)+1]  <- paste("Average distance between real and estimated CDF")
        ylablist[length(ylablist)+1]  <- paste("Average squared distance between real and estimated CDF")
        ylablist[length(ylablist)+1]  <- paste("Error Introduced by DP noise")
        ylablist[length(ylablist)+1]  <- paste("Error Introduced by DP noise")
        ylablist[length(ylablist)+1]  <- paste("Error Introduced by DP noise")
        ylablist[length(ylablist)+1]  <- paste("Error Introduced by DP noise")
        ylablist[length(ylablist)+1]  <- paste("Error Introduced by DP noise")
        ylablist[length(ylablist)+1]  <- paste("Error Introduced by DP noise")

        if (length(labellist) > length(ylablist)){
            for (Y in length(ylablist):length(labellist)){
              ylablist[Y] <- paste("Custom: user-defined units (PDF)")
            }
        }
        boxplot(forplotting, 
                las = 2,
                col = cols, 
                at = ats,
                names = namer,
                par(mar = c(12, 5, 4, 2) + 0.1),
                ylim = c(0,mx),        
                main = labellist[p],
                ylab = " "
                )    
        mtext(ylablist[p],side=2, line =4)

        for (d in 1:length(datalist)){
        for (f in 1:length(functlist)){
          labelindex <- (d-1)*(length(functlist)+2)+f+1
        mtext(Fnameslist[f], side=1, line=1, at=labelindex, las=2, font=2, col=(f+1))
        }}

        datatitlesAt <- numeric(0)
        for (d in 1:length(datalist)) {
          if(length(functlist)>2){
        datatitlesAt[d] <- (d-1)*(length(functlist)+2)+ mean(1,length(functlist)) +2.5
        }
          if(length(functlist)<=2){
        datatitlesAt[d] <- (d-1)*(length(functlist)+2)+ mean(1,length(functlist)) +1
        }
        text(datatitlesAt[d],.9*mx, Dnameslist[d], font=2)
         if(p==10){
            text(datatitlesAt[d],median(datalist[[d]]),"----------", font=2, col=1)
                  }
                  }                   
    }
}

dev.off()
return(masterlist)
}}


##begin UI:
## To use, highlight from here downwards and press ctrl+/
## some arbitrary example values have been input.
## The 3 functions here are available from the Harvard DP group in privateZelig on Github
#   results <- CDFtest(
# 
#   Visualization = TRUE,
# #           #Sets the testing suite into Visualization mode (default) or Data Collection mode (visualization = FALSE)
# #           #   In Visualization mode (default):
# #           #       An excel sheet conatining the mean and median results (across reps iterations) of diagnostic functions on DP-CDF algorithms
# #           #         per each combination of dataset, function, and epsilon.
# #           #       A .pdf file containing one graphical example DP CDF for each combination of dataset, function, and epsilon,
# #           #        as well as a set of boxplots showing the distribution of all diagnostic results for all combinations of parameters.
# #           #   In Data Collection mode (set Visualization = FALSE):
# #           #       An excel sheet conatining the entire (raw) results (across reps iterations) of diagnostic functions on DP-CDF algorithms
# #           #        per each combination of dataset, and function, sepearately looped over all epsilons, then all granularities, and all samplesizes
# 
#            OutputDirectory = "C:\\Users\\Dan\\Dropbox\\Harv\\CDFtesting\\",
# #           # the location of the folder which will hold the output (csv and pdf files)
# 
#            functlist = c(functionH),
# #           # A list of CDF-computing functions to be tested on the CDFtestTrack
# 
#            Fnameslist = c("H"),
# #           # A list of function names corresponding to the functions being tested; used for labelling the output.
# epslist = c(.2),
# #           epslist   = c(0.01,1,0.1,0.001,0.0001),
# #           # A vector of epsilon values for differential privacy
# 
#            datalist  = list(),
# #           # A list containing vectors of data, each to be used in a test
# 
#            Dnameslist     = c(),
# #           # A list of dataset  names corresponding to the functions being tested; used for labelling the output.
# synthsets = list(list("wage", 10000, "uniform")),
# #           synthsets      = list(list("wage", 100000, "uniform"), list("wage",100000,"sparse"), list("wage",100000,"bimodal"), list("wage", 100000, "gaussian")),
# #           # This script generates pre-defined synthetic datasets upon request, and fulling incorporates them into testing.
# #           # To call them, users should input a character vector containing the names of the sets they desire.
# #           # For example, synthsets = list( list(type,size,shape) , list(type,size,shape) ). There are no limits on the amounts of datasets included.
# #           # sets available include:
# #               # type: ("age" which ranges from about 0 to 100, gran =1) ("wage" which ranges from 0 to 450k, gran = 1000)
# #               # size: Any positive integer. Type in exact numerical representation (eg, for ten thousand use 10000 not 10k and not 10^4)
# #               # shape: gaussian, sparse, uniform, bimodal
# 
#            range     = c(1,500000),
# #           # The range of the universe as a vector (min, max). Defined based on user intuition.
# #           # Setting min too high will bias output upward. Same in reverse for a low max.
# #           # However, setting min too low and max too high could reveal the true limits of your data, compromising some privacy.
# 
#            gran    =1000,
# #           # FOR Visualization MODE ONLY. refer to "granlist" for setting granularities (thus domain sizes) in Data Collection mode.
# #           #     This command is irrelevant in Data Collection mode.
# #           # The granularity of the domain between the min and max. ie, if age is measure per 1 year of age, gran =1
# #           # The same granularity is applied to all datasets, so using comparable (or scaled) data is necessary
#     #  granlist = c(1000,100000,50000,25000,12500,10000,5000,2500,1250,500,250,125,100,50,25,10),  
#           granlist =c(5000),
#         #    granlist =c(100000,50000,25000,12500,10000,5000,2500,1250,1000,500,250,125,100,50,25),
# #           # FOR Data Collection MODE ONLY. refer to "gran" for selecting granularities in Data Collection mode.
# #           #     This command is irrelevant in Visualization mode.
# #           # A list of granularities of the universe between the min and max. ie, if age is measure per 1 year of age, gran =1
# 
#            samplesize = 0,
# #           # FOR Visualization MODE ONLY. refer to "nlist" for selecting samplesizes in Data Collection mode.
# #           #     This command is irrelevant in Data Collection mode.
# #           #     when set to zero, the entire dataset is used.
# #           # Otherwise, the specified sample size is randomly selected from each dataset without replacement.
# #           # the samples are used.
# nlist = c(10000),
# #           nlist = c(100000,100,1000,10000,1000000),
# #           # FOR Data Collection MODE ONLY. refer to "samplesize" for selecting samplesizes in visualization mode.
# #           #     This command is irrelevant in Visualization mode.
# #           #sets the absolute sample sizes to draw from each dataset, with replacement. Any vector of integer values is appropriate.
# 
#            cdfstep   =0,
# #           # The step size used in outputting the approximate CDF.
# #           # The values output are [min, min + cdfstep], [min, min + 2 * cdfstep], etc.
# #           # For best results (maximum CDF resolution), set cdfstep = 0, which gives cdfstep = gran
# 
#            reps    = 50,
# #           # The number of times to repeat each diagnostic. higher reps lends greater accuracy,
# #           # but comsumes time and power. Author recommends reps = 10 for quick examples and 300 for robust examinations
# 
#            ExtraTests_CDF = list(),
# #           # If a user wishes to add extra diagnostics, the proper syntax would be:
# #           # ExtraTests_CDF = list( functionName1 = function1, functionName2 = function2)
# #           # Diagnostic Functions should have similar inputs as diagnostics defined above, and output just one value
# 
#            ExtraTests_PDF = list(),
# #           # see above
# 
#            setseed    = c(-100),
# #           # In the function, each combination of data, epsilon, and function is executed with a separate seed, which
# #           # by default is randomly generated and reported. Users interested in replicating specific results
# #           # can locate the reported seed and parameter combination to replicate tests.
# 
#            comments       = ("x"),
# #           # In the function, each combination of data, epsilon, and function is executed with a separate seed, which
# #           # by default is randomly generated and reported. Users interested in replicating specific results
# #           # can locate the reported seed and parameter combination to replicate tests.
# 
#            SmoothAll      = FALSE,
# #           # Applies a monotonizing post-processing function to all dpCDFs
# #           # Monotonization does not affect prob dense functions
# 
#            EmpiricBounds  = FALSE,
# #           # When TRUE, outputted graphs depict the minimum and maximum values taken by each bin across reps
# 
#            AnalyticBounds = FALSE,
# #           # AnalyticBounds: This is a flag and should be set to "true" if the functions being tested are expected to output
# #           # analytical variance bounds. The proper output form is output = list(DPCDFvector, LowerBoundVector, UpperBoundVector)
# #           # This will not switch on an Analytical bounds flag within the inputted function if the function's flag is set/defaulted to "FALSE".
# #           # However, if this flag here is set FALSE, it will suppress bounds from the inputted function.
# #           # Best practice is to have an analytical bounds flag defaulted to "TRUE" within the inputted function, such that it can be turned on and off
# #           # with this flag.
# 
#            AnalyticProbSleeve = TRUE,
# #           # When TRUE, outputted DP-CDFs will have a fuzzy analytic sleeve aroudn them, approximating probabalitity density for each point.
# 
# #           SuppressRealCDF = FALSE,
# #           # When TRUE, outputted graphs will not include Real (non-private) CDFs
# 
#            SuppressDPCDF = FALSE,
# #           # When TRUE, outputted graphs will not include DP-CDFs (but if SmoothAll = TRUE, monotonized DP CDFs still appear)
# 
#            SuppressLegends = FALSE,
# #           # When TRUE, outputted graphs will not include legends
#           )