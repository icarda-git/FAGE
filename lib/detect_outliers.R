#' Title line
#' 
#' Description paragraph
#' 
#' @param df              description 
#' @param mod.ref         description
#' @param threshold.start description
#' @param threshold.stop  description
#' @param threshold.step  description
#' 
#' @return description
#' 
#' @author 
#' Khaled Al-Shamaa. ICARDA
#' 
#' @examples 
#' # R code
detect_outliers <- function(df=NULL, mod.ref=NULL, threshold.start=4, threshold.stop=3, threshold.step=0.25){
  
  model.formula <- mod.ref$call
  
  # create an empty vectors to save the outliers entry id and associated residual value
  outlier.id <- c();   outlier.scres <- c()
  
  for(outlier.threshold in seq(threshold.start, threshold.stop, -1*threshold.step)){
    modr <- paste('trial.diag <-', c(model.formula), sep = ' ')
    eval(parse(text=modr))
    
    # update and re-fit an asreml model to return standardized conditional residuals 
    trial.diag     <- update(trial.diag, aom = TRUE)
    trial.diag.sum <- summary(trial.diag, nice = TRUE)$nice[[1]]
    
    # identify outliers
    df$scres <- trial.diag$aom$R[, 2]
    trial.outliers <- df[!is.na(df$scres) & abs(df$scres) > outlier.threshold, ]
    
    # list records to check by breeders
    check.records <- as.vector(unlist(apply(trial.outliers[, c("trial", "gen")], 1, 
                                            function(x) which(df$trial == x[1] & df$gen == x[2]))))
    
    # report outliers records associated with response values for the same genotypes 
    # in other replications for ease of compare 
    if(!empty(trial.outliers)) {
      outlier.id <- c(outlier.id, trial.outliers$id)
      outlier.scres <- c(outlier.scres, trial.outliers$scres)
      
      if(exists("check.data")){
        check.data <- rbind(check.data, df[check.records, ])
      }else{
        check.data <- df[check.records, ]
      }
    }
    
    # make outliers as NA
    df$resp[!is.na(df$scres) & abs(df$scres) > outlier.threshold] <- NA
  }
  
  if(!exists("check.data")) {
    check.data <- NULL
  }else{
    check.data$outlier <- as.numeric(abs(check.data$scres) > threshold.stop)
  }
  
  return(list(data=df, verify=check.data, id=outlier.id, scres=outlier.scres))
}