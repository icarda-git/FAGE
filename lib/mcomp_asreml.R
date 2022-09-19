#' Calculates multiple comparisons for a given factor using ASReml-R output
#'
#' Generates a table with all possible pairwise comparisons between levels
#' of a given factor. The test corresponds to a simple t-test without any
#' control for multiple testing, hence, it is equivalent to a LSD test. The
#' input corresponds to elements obtained from the predict.asreml() function
#' from ASReml-R. It also produces, if requested, and overal F-test of
#' the null hypothesis of all means equal.
#'
#' @param mu A vector with the predicted means of the factor to be evaluated.
#' @param V The variance-covariance matrix of the predicted means as part of
#' a predict statement. It needs to have the same number of rows and columns
#' as in mu.
#' @param nedf Denominator degrees of freedom, as provided by the LM object
#' or from the ANOVA table (default = 9999).
#' @param test The type of test to be requested, where \code{'overall'} 
#' corresponds to the F-test comparing all menas, and \code{'single'} is all
#' the pairwise comparisons (default = \code{'overall'}).
#' @param labels A vector of character names corresponding to each of the 
#' levels of the factor of interest (default = NULL).
#'
#' @return A list with ONLY one of these two elements: \cr
#' \itemize{
#' \item \code{p.value.F}: the p-value associated with the overall F-test. \cr
#' \item \code{mcomp}: a data frame with the result of all pairwise comparisons. \cr
#' }
#' 
#' @author 
#' Salvador A. Gezan
#' VSN International
#' 
#' @example
#' 
#' 
#' @export

mcomp.asreml <- function(mu=NULL, V=NULL, nedf=9999, test='overall', labels=NULL){
  
  V <- as.matrix(V)
  # Check if the class of V is matrix
  if (length(V)==0 || is.null(V)) {
    stop("Matix V is null or of dimension zero.")
  }
  if (length(mu)==0 || is.null(mu)) {
    stop("Vector mu is null or of dimension zero.")
  }
  if (length(labels)==0 || is.null(mu)) {
    message("Input of labels is null or of dimension zero. They will be made numbers 1:t.")
    labels <- as.character(1:length(mu))
  }

  if (sum(is.na(mu))>0) {
    message('Some missing values in mu, these will be droped form mu and V.')
    V <- V[!is.na(mu),!is.na(mu)]
    mu <- mu[!is.na(mu)]
    labels <-labels[!is.na(mu)]
  }
  pvalF <- NULL
  pvalT <- NULL
  t <- length(mu)
  if (t==1){
    stop('Only a single term in mu provided, no need for testing.')
  }
  
  # Internal function to generate L matrix
  LmatrixF <- function(ntrt=NULL) {
    L <- matrix(0, nrow=(ntrt-1), ncol=ntrt)
    for (i in 1:(ntrt-1)){
      L[i,i] <-1
      L[i,i+1] <- -1
    }
    return(L)
  }
  
  # Internal function to generate L matrix
  LmatrixT <- function(ntrt=NULL, labels=NULL) {
    ntest <- ntrt*(ntrt-1)/2
    L <- matrix(0, nrow=ntest, ncol=ntrt)
    df.comp <- data.frame(matrix(0, nrow=ntest, ncol=2))
    m <- 1
    for (i in 1:(ntrt-1)){
      for (j in (i+1):ntrt) {
        L[m,i] <- 1
        L[m,j] <- -1
        df.comp[m,1] <- labels[i]
        df.comp[m,2] <- labels[j]
        m <- m + 1
      }
    }
    return(list(L=L, df.comp=df.comp))
  }
    
  # Doing overall F test
  if (test=='overall') {
    L <- LmatrixF(ntrt=t)
    Z <- L %*% mu
    V.Z <- L %*% V %*% t(L)
    numdf <- rankMatrix(L)[1]
    Fobs <- as.numeric( t(Z) %*% solve(V.Z) %*% Z / numdf )
    pvalF <- 1-pf(abs(Fobs), numdf, nedf, log=FALSE)
  }
  
  # Doing single t-tests
  if (test=='single') {
    Lobj <- LmatrixT(ntrt=t, labels=labels)
    L <- Lobj$L
    m <- nrow(L)
    pvalT <- matrix(0, ncol=4, nrow=m)
    for (k in 1:m) {
      Z <- L[k,] %*% mu
      V.Z <- L[k,] %*% V %*% (L[k,])
      Tobs <- as.numeric( Z/sqrt(V.Z) )
      pvalT[k,1] <- Z
      pvalT[k,2] <- sqrt(V.Z)
      pvalT[k,3] <- Tobs
      pvalT[k,4] <- 1-pt(abs(Tobs), nedf, log=FALSE)
    }
    # Building data.frame
    pvalT <- data.frame(Label1=Lobj$df.comp[,1], Label2=Lobj$df.comp[,2], 
                        diff=pvalT[,1], std.diff=pvalT[,2], 
                        t.value=pvalT[,3], p.value=pvalT[,4])
  }
    
  return(list(p.value.F=pvalF, mcomp=pvalT))
  
}
