#' Verifies data for MET analyses
#'
#' \code{auditMET} Evaluates and verifies the data originating from several trials with the 
#' aim of determining levels of: connectivity, variability, etc., which are reported as 
#' statistics. Only non-NA observations are considered in the reports. 
#
#' @param data dataframe with all relevant columns for MET analyses.
#' @param gen factor name for genotypes (or treatments)
#' @param trial factor name for trial (or environment) 
#' @param resp column name for the response variable to evaluate
#' @param single if TRUE it is assume there is a single observation per genotype-trial (default = TRUE)
#'
#' @return Several objects with reports of the MET data. 
#' incidence:   Matrix of incidence with counts (non-NA) of genotypes by trial (on diagonal), and common genotypes
#'               between pairs of trials (off-diagonal).
#' gt:          Statistics for the MET data for the response variable of interest, including columns with: 
#'               names of trials, number of non-NA genotypes (n), mean, minimum, maximum, standard deviation
#'               and coefficient of variation (CVp, %).
#'               
#' @author 
#' Salvador A. Gezan. VSN International
#' 
#' @examples
#' # Example 1: 
#' 

auditMET <- function(data=NULL, gen=NULL, trial=NULL, resp=NULL, single=TRUE) {
  
  n<-nrow(data)
  if (n==0) { stop('No information in data.frame provided.')}
  # Defining factors
  df <- data.frame(IDSORT=c(1:n)) 
  if (is.null(gen)) { 
    stop('No genotype column provided.')
  } else {
    df$gen <- as.factor(data[,gen])
  }
  if (is.null(trial)) { 
    stop('No trial column provided.')
  } else {
    df$trial <- as.factor(data[,trial])
  }
  if (is.null(resp)) { 
    stop('No response column provided.')
  } else {
    df$resp <- data[,resp]
  }
  
  # df: IDSORT, gen, trial, resp
  # Eliminating NA
  df <- df[!is.na(df$resp),]
  # Checking Connectivity INC.matrix
  s <- length(levels(df$trial))
  inc <- as.matrix(table(df$gen,df$trial))
  for (i in 1:s) {
    inc[inc[,i]>0,i] <- 1
  }
  INC.matrix <- t(inc) %*% inc
  if (min(INC.matrix) < 5) {
    message('You have on some pairs of trials less than 5 genotypes on common')
    message('hence, your analysis might be unstable.')
  }
  
  # Obtaining - summary statitcs by trial: mean, min, max, SD, CV%   
  agg.mean <- aggregate(df$resp,by=list(df$trial),FUN=mean,na.rm=TRUE)
  agg.min  <- aggregate(df$resp,by=list(df$trial),FUN=min,na.rm=TRUE)
  agg.max  <- aggregate(df$resp,by=list(df$trial),FUN=max,na.rm=TRUE)
  agg.n    <- table(df$trial)
  agg.sd   <- aggregate(df$resp,by=list(df$trial),FUN=sd,na.rm=TRUE)
  
  gf <- data.frame(agg.n, agg.mean[,2], agg.min[,2], agg.max[,2], agg.sd[,2])
  gf$CVp <- 100*gf[,6]/gf[,3]
  colnames(gf)<-c(trial,'n','mean','min','max','sd','CVp')

  ratio <- max(gf$sd,na.rm=TRUE)^2/min(gf$sd,na.rm=TRUE)^2
  if (ratio > 5) {
    message('The ratio of max to min variance is greater than 5. Consider scaling your data.')
  }
  
  return(list(incidence=INC.matrix,stats=gf))
  
}
