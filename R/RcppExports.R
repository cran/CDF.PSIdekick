# This file was generated by Rcpp::compileAttributes
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @title Creates a Tree then a CDF
#' @description This thing sure does make a fine CDF
#' @param eps An epsilon value for Differential Privacy
#' @param ds The data or something
#' @param Ks the degree of the tree
#' @param methods Either H or S2 or SUB
#' @param mins the minimum of the domain's range
#' @param maxs the maximum of the domain's range
#' @param grans The granularity
#' @param datas The data to be CDFd
#' @return A dpCDF
#' @export
TreeCDF <- function(eps, ds, Ks, methods, mins, maxs, grans, datas) {
    .Call('CDF_PSIdekick_TreeCDF', PACKAGE = 'CDF.PSIdekick', eps, ds, Ks, methods, mins, maxs, grans, datas)
}

#' @title Monotonicity enforcement
#' @description When CDFs get out of line, we call the enforcer 
#' @param x A numeric vector to be enforced 
#' @export 
#' @return A monotonized vector 
Smooth <- function(x) {
    .Call('CDF_PSIdekick_Smooth', PACKAGE = 'CDF.PSIdekick', x)
}

