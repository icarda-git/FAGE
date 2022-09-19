#' Routine to obtain several goodness-of-fit statistics for spatial and non-spatial models 
#'
#' \code{gral.stats} Internal routine that calculates goodness-of-fit statisitcs for a
#' specific fitted model provided as input. Statistics obtained are log-likelihood, AIC, 
#' BIC, heritability based on predictor error variance (heritPEV, only for genotypes 'random'),
#' A-optimality value, and logarithm of D-optimality value,
#' 
#' @param object model fit output for object class asreml. Genotype factor is called 'gen'
#' @param solution data frame with the random solutions from asreml (BLUPs)
#' @param preds data frame with the predictions for genotypes (required to calculate Aopt, Dopt)
#' 
#' @return A table with goodness-of-fit statistics for models evaluated. This includes columns: 
#' number of variance components in the model (n.VC), log-likelihood (logL), Akaike information 
#' criteria (AIC), Bayesian information criteria (BIC), heritability based on variance components,
#' heritability based on predictor error variance (heritPEV, only for genotypes 'random'), 
#' A-optimality value, and logarithm of D-optimality value. In addition a heritability based on
#' predictor error variance (herit.PEV, only for genotypes random), and heritability calculated 
#' based on variance componens (herit.VC = var.gen/var.total), the genetic variance estimate (v.gen),
#' the residual variance estimate (v.res), the coefficient of variation (%) based on the mean of 
#' predictions and the residual variance, and the 5% lsd based on the average of the standard error
#' of the difference and infinity degrees of freedom (i.e., t0 = 1.96)
#'
#' @author
#' Salvador A. Gezan. VSN International
#' 
#' @examples
#' # Example: Pending
#' 


gral.stats <- function(object=NULL, solution=NULL, preds=NULL, ctable=NULL){
  
  if (is.null(object)) { stop('No fitted model provided.')}
  ss <- summary(object)
  vc <- summary(object)$varcomp
  
  # Extracting the solutions and calculating H2PEV
  H2PEV <- NA
  H2VC <- NA
  var.res <- NA
  var.gen <- vc[grep('gen',rownames(vc)),1]

  if (length(var.gen) > 0){
    if (!is.null(solution)) {
      #sol.gen <- as.data.frame(solution[grep('gen',rownames(solution)),])
      H2PEV <- 1-mean(solution$std.error^2)/(2*var.gen)  
    }
    # Calculating H2VC
    vctemp <- vc[vc$bound!='F' & vc$bound!='U',1]
    var.res <- vctemp[length(vctemp)]
    H2VC <- var.gen/sum(vctemp)
  }

  if (length(var.gen) == 0){ var.gen = NA }
    
  # Dealing with the predictions: Aopt / Dopt
  # Note that these can be restricted to only comparisons against controls
  if (!is.null(preds)) {
      pvals <- preds$pvals
      vcov <- as.matrix(preds$vcov)
      sel <- matrix(1, ncol=1, nrow=length(pvals$predicted.value))
      sel[is.na(pvals$predicted.value),] <- 0
      vcov <- vcov[sel==1,sel==1]
      Aopt <- mean(diag(vcov))
      Dopt <- log(det(vcov))
      # Calculating cv
      cv.VC <- 100*sqrt(var.res)/mean(pvals$predicted.value[sel==1])
      # Calculating LSD 5% 
      mcomp <- mcomp.asreml(mu=pvals$predicted.value[sel==1], V=vcov, test='single', labels=1:nrow(vcov))
      mean.sed <- sqrt( mean(mcomp$mcomp$std.diff^2)  )
      lsd <- mean.sed*1.96
  }
  if (is.null(preds)) {
    Aopt <- NA; Dopt <- NA
    cv.VC <- NA; lsd <- NA
  }
  
  return(list(n.vc=nrow(vc), logL=ss$loglik, aic=ss$aic, bic=ss$bic, 
              Aopt=Aopt, Dopt=Dopt, herit.PEV=H2PEV, herit.VC=H2VC,
              v.gen=var.gen, v.res=var.res, cv.VC=cv.VC, lsd=lsd))
}