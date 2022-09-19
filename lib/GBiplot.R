#' Generates a Biplot from multivariate data or MET data
#'
#' \code{GBiplot} Generates a Biplot for multivariate data with vectors and points based on 
#' Principal Component Analyses (PCA). It reads data with column for vectors (e.g. traits or 
#' environments) and column for units (e.g. genotypes) and one or more column for response 
#' variables (e.g. predicted means from an MET analysis). If the matrix of variance-covariance
#' (or correlations) between vectors is provided, this is used to define principal components.
#' Missing values are allowed and imputed using the fucntion imputePCA from package missMDA. 
#'
#' @param data dataframe with column for vectors and units together with response variables
#' @param vcovM matrix with variance-covariance estimates of the vectors (traits or environments) 
#' @param corrM matrix with correlation estimates of the vectors (traits or environments) 
#' @param vector factor name for vectors (traits or environments)
#' @param unit factor name for units (often genotypes)
#' @param resp column name for the response variable to evaluate
#' @param scale logical, if TRUE PCA analysis will scale the response variable (and correlations 
#'              are used), otherwise original responses are used (and variance-covaraince matrix)
#'              (default = TRUE)
#' @param groups column number that identifies the factor (or variate) to separate units according 
#'               to groups. Each group will be plotted with a different color.
#' @param vector.label logical, if TRUE vector labels will be included in the plot (default = TRUE). 
#' @param unit.label logical, if TRUE unit labels will be included in the plot (default = FALSE). 
#'
#' @return A graphical output with a biplot
#'
#' @author 
#' Salvador A. Gezan. VSN International
#' Amanda Avelar de Oliveira. VSN International
#' 
#' @examples
#' # Example 1: Default Biplot
#' head(GxE_Preds)
#' GBiplot(X=GxE_Preds, vector=2, unit=3, response=4, scale=TRUE)
#' GBiplot(X=GxE_Preds, vector=2, unit=3, response=4, legend=TRUE)
#' 

GBiplot <- function(data=NULL, vcovM=NULL, corrM=NULL, vector=NULL, unit=NULL, resp=NULL, 
                    groups=NULL, scale=TRUE, vector.label=TRUE, unit.label=FALSE) {

  #library(ggfortify)
  #library(missMDA)
  n <- nrow(data)
  if (n==0) { stop('No information in data.frame provided.') }
  if (is.null(vector)) { stop('No vector column provided.') } 
  if (is.null(unit))   { stop('No unit column provided.') }
  if (is.null(resp))   { stop('No response column provided.') }

  # Reformating data Matrix 
  n <- length(unique(data[,vector]))
  m <- length(unique(data[,unit]))
  Xmatrix <- data.frame(matrix(data[,resp], nrow=m, byrow=TRUE))
  colnames(Xmatrix) <- unique(data[,vector])
  rownames(Xmatrix) <- unique(data[,unit])
  #print(head(Xmatrix))
  # Some additional checks
  if (nrow(Xmatrix) < 3) { stop("Matrix needs at least 3 units.") }
  if (ncol(Xmatrix) < 2) { stop("Matrix needs at least 2 vectors (trait or environment).") }
  
  # If the dataset contains missing data impute it before generating Biplot
  pctMiss <- sum(is.na(Xmatrix))/(nrow(Xmatrix)*ncol(Xmatrix))
  if (pctMiss > 0) {
    cat("Missing values detected (", round(100*pctMiss,0), "%)\n", sep = "")
    Data_imputed <- missMDA::imputePCA(Xmatrix, ncp = 2)
    Data_imputed <- Data_imputed$completeObs
    # Check if the missing were imputed
    pctMiss <- sum(is.na(Data_imputed))/(nrow(Data_imputed)*ncol(Data_imputed))
    if (pctMiss == 0){ message("Missing values were imputed.") }
    Xmatrix <- as.matrix(Data_imputed)
  }

  # Older code in case we need it
  #PCA1 <- prcomp(Xmatrix, scale.=TRUE)
  #GGEBiplot <- autoplot(PCA1, data=NULL, colour=NULL, loadings=TRUE,  
  #                      loadings.colour='blue', loadings.label=TRUE, 
  #                      loadings.label.size=3, loadings.label.colour='blue', 
  #                      frame=FALSE) #, scale=1)
  #GGEBiplot   # This is a bit different in terms of scale of the axis, a bit different same points
  
  # If we want scale - default
  if (scale) { 
    Xscale <- scale(Xmatrix, center=TRUE, scale=TRUE)

    if (!is.null(corrM)) {
      message('Provided correlation matrix, corrM, used.')
    }     
    if (is.null(corrM) & !is.null(vcovM)) {
      message('Provided var-covariance matrix, vcovM, used.')
      corrM <- cov2cor(vcovM)
    }
    # No corrM or vcovM provided
    if (is.null(vcovM) & is.null(corrM)) {
       corrM <- t(Xscale) %*% Xscale / (nrow(Xscale)-1)
    }

    # PCA-eigenvalues
    Peig <- eigen(corrM)
    #sqrt(Peig$values)  # same as s$sdev
    prop.eig <- Peig$values/sum(Peig$values) 
    #print(prop.eig)
    loadings <- Peig$vectors  # Loadings
    PCAcoord <- as.matrix(Xscale) %*% as.matrix(loadings)  # New PCA coordinates
    colnames(PCAcoord) <- paste('PC',c(1:ncol(Xscale)),sep='')
    colnames(loadings) <- paste('PC',c(1:ncol(Xscale)),sep='')
    rownames(loadings) <- paste('PC',c(1:ncol(Xscale)),sep='')
    xlab <- paste('PC1 (', round(100*prop.eig[1],2), '%)', sep='')
    ylab <- paste('PC2 (', round(100*prop.eig[2],2), '%)', sep='')
    Biplot <- ggfortify::ggbiplot(plot.data=data.frame(PCAcoord), loadings.data=data.frame(loadings),
                                  loadings=TRUE, loadings.label=vector.label, loadings.label.colour='blue', 
                                  xlab=xlab, ylab=ylab, asp=1,
                                  loadings.label.label=colnames(Xmatrix),
                                  label=unit.label, label.label=rownames(Xmatrix))
  }
  # No scaling
  if (!scale) { 
    
    Xcent <- scale(Xmatrix, center=TRUE, scale=FALSE)
    if (!is.null(vcovM)) {
      message('Provided var-covariance matrix, vcovM, used.')
    }
    # No corrM or vcovM provided
    if (is.null(vcovM)) {
      vcovM <- t(Xcent) %*% Xcent / (nrow(Xcent)-1)
      #print(vcovM)
    }
    # PCA-eigenvalues
    Peig <- eigen(vcovM)
    #sqrt(Peig$values)  # same as s$sdev
    prop.eig <- Peig$values/sum(Peig$values) 
    #print(prop.eig)
    loadings <- Peig$vectors  # Loadings
    PCAcoord <- as.matrix(Xcent) %*% as.matrix(loadings)  # New PCA coordinates
    colnames(PCAcoord) <- paste('PC',c(1:ncol(Xcent)),sep='')
    colnames(loadings) <- paste('PC',c(1:ncol(Xcent)),sep='')
    rownames(loadings) <- paste('PC',c(1:ncol(Xcent)),sep='')
    xlab <- paste('PC1 (', round(100*prop.eig[1],2), '%)', sep='')
    ylab <- paste('PC2 (', round(100*prop.eig[2],2), '%)', sep='')
    Biplot <- ggfortify::ggbiplot(plot.data=data.frame(PCAcoord), loadings.data=data.frame(loadings),
                                  loadings=TRUE, loadings.label=vector.label, loadings.label.colour='blue', 
                                  xlab=xlab, ylab=ylab, asp=1,
                                  loadings.label.label=colnames(Xmatrix),
                                  label=unit.label, label.label=rownames(Xmatrix))
  }
  
  return(Biplot)
  
}
