#' Routine to fit several spatial models and suggest the best model
#'
#' \code{spatial.selector} Generates a table with the goodness-of-fit statitics to select the
#' best spatial model for an analysis of different analyses (unreplicated and replicated). 
#' Models that can not fit are reported as NA. 
#' The models require some input conditions from user: 1) Genotypes: fixed/random; 2)	Blocks: 
#' TRUE/FALSE (if TRUE then fixed/random); 3) Covariates (always fixed, max 2); and 4) 
#' Incomplete Blocks: TRUE/FALSE (always random).
#' The model selector will run over a series of models with the conditions: 1) Rows within 
#' replicate: TRUE/FALSE (always random); 2) Columns within replicate: TRUE/FALSE (always
#' random); 3) Spline on rows: TRUE/FALSE (fixed linear covariate + random spline); 4) Spline
#' on columns: TRUE/FALSE (fixed linear covariate + random spline); and 5) Residual structure:
#' 'indep' or 'ar1'. In addition the best spatial model (with ar1) and selected according to 
#' A-optimality is fitted with an additional nugget effect. List of 32 models is found on data
#' MODLIST.Rda.
#' 
#' @param data dataframe with all relevant columns for spatial model and response variables.
#' @param gen factor name for genotypes (or treatments)
#' @param check column name for variable identifying checks (categories 0:TEST, 1:CHECK) (optional)
#' @param block factor name for full block (or replicates) (optional) 
#' @param ibk factor name for incomplete block (optional) (optional)
#' @param row column name for row coordinates of each experimental unit
#' @param col column name for column coordinates of each experimental unit
#' @param cov1 column name with additional covariate 1 (optional)
#' @param cov2 column name with additional covariate 2 (optional)
#' @param resp column name for the response variable to analyze
#' @param type.gen model assumption for genotypes: 'random' or 'fixed' (default: 'random')
#' @param type.block model assumption for full blocks: 'random' or 'fixed' (default: 'fixed')
#' @param nugget logical to add nugget to any spatial model fitted (default: FALSE)
#' @param model model number to be fitted (optional)
#' @param strict number for level of selection for models to fit (0:none, 1:moderate, 2:strict) (default = 1)
#' 
#' @return A table with goodness-of-fit statistics for models evaluated. This includes columns: 
#'    number of variance components in the model (n.VC), log-likelihood (logL), Akaike information 
#'    criteria (AIC), Bayesian information criteria (BIC), and heritability based on predictor 
#'    error variance (heritPEV, only for genotypes 'random').
#' 
#'    If parameter 'model' is specified the object with fitted model number is provided, with 
#'    several objects with details of the fitted model. 
#' 
#'    aov:         Wald-test (mixed model ANOVA-like table)
#'    call:        String with the ASReml-R call used to fit the requested model
#'    gt:          Goodness-of-fit statistics for the model evaluated are reported as summary. This includes columns: 
#'                  number of variance components in the model (n.VC), log-likelihood (logL), Akaike information 
#'                  criteria (AIC), Bayesian information criteria (BIC), A-optimality value, logarithm of the 
#'                  D-optimality value, and heritability based on predictor error variance (heritPEV, only for 
#'                  genotypes 'random')
#'    mod:         ASReml-R object with all information from the fitted model
#'    predictions: Predictions for all genotypes, with an additional column of 'weight' is to be used on a 
#'                  on a second-stage analysis. These weights are the diagonal of the inverse of the variance-
#'                  covariance matrix of predictions
#'
#' @author 
#' Salvador A. Gezan. VSN International
#' 
#' @examples
#' # Example 1: Selecting best model from Replicated Trial
#' library(agridat)
#' testREP <- durban.rowcol
#' testREP$bed <- as.factor(testREP$bed) 
#' testREP$row <- as.factor(testREP$row) 
#' testREP$gen <- as.factor(testREP$gen) 
#' head(testREP)
#' test.sel <- spatial.selector(data=testREP, gen='gen', row='row', col='bed', 
#'                              resp='yield', type.gen='random', strict=2)
#' ls(test.sel)
#' View(test.sel$parms)
#' test.sel$mod
#' 
#' # Example 2: Fitting only selected model 17
#' mod.sel <- spatial.selector(data=testREP, gen='gen', row='row', col='bed', 
#'                              resp='yield', type.gen='random', model=17)
#' ls(mod.sel)
#' mod.sel$call                  # Call for model fitted
#' plot(mod.sel$mod)             # Residual plots
#' plot(varioGram(mod.sel$mod))  # Variogram fitted model
#' summary(mod.sel$mod)$varcomp  # Variance Components
#' mod.sel$aov                   # Wald-test Table
#' mod.sel$gt$herit.PEV          # Heritability (other statistics on gt)
#' head(mod.sel$predictions)     # Predictions for fitted model with additional weights

spatial.selector <- function(data=NULL, gen=NULL, check=NULL, block=NULL, ibk=NULL, row=NULL, col=NULL, 
                        cov1=NULL, cov2=NULL, resp=NULL, 
                        type.gen='random', type.block='fixed', nugget=FALSE, model=NULL, strict=1){

  load(file='./lib/MODLIST.Rda')
  #asreml.options(maxiter=20, extra=3)
  # Generic parameters for all models
  if (!is.null(block)) { add.block <- TRUE; type.block <- type.block }
  if (is.null(block))  { add.block <- FALSE; type.block <- NULL } ### Fixed
  if (!is.null(ibk)) { add.ibk <- TRUE }
  if (is.null(ibk)) { add.ibk <- FALSE }
  if (!is.null(cov1)) { add.cov1 <- TRUE }
  if (is.null(cov1)) { add.cov1 <- FALSE }
  if (!is.null(cov2)) { add.cov2 <- TRUE }
  if (is.null(cov2)) { add.cov2 <- FALSE }

  if (is.null(model)) { 
  
    # Selecting models
    set <- mod.list$Model[mod.list$List >= strict]
    gfit <- matrix(NA, ncol=9, nrow=length(set))

    for (i in 1:length(set)) {
      add.row <- mod.list$RepRow[set[i]]
      add.col <- mod.list$RepCol[set[i]]
      add.spl.row <- mod.list$splineRow[set[i]]
      add.spl.col <- mod.list$splineCol[set[i]]
      type.residual <- mod.list$Resid[set[i]]
      if (type.residual=='ar1') {
        if (nugget) { add.nugget <- TRUE }
        if (!nugget) { add.nugget <- mod.list$Nugget[set[i]] }
      }
      mod.ref <-  gral.single(data=data, gen=gen, check=check, block=block, ibk=ibk, row=row, 
                             col=col,  cov1=cov1, cov2=cov2, resp=resp,
                             add.block=add.block,     add.ibk=add.ibk,  
                             add.row=add.row,         add.col=add.col,
                             add.spl.row=add.spl.row, add.spl.col=add.spl.col,
                             add.cov1=add.cov1,       add.cov2=add.cov2,
                             add.nugget=add.nugget,   type.gen=type.gen,
                             type.block=type.block,   type.residual=type.residual)
      # Storing the statistics
      gfit[i,1] <- mod.list$Model[set[i]]
      if(!is.na(mod.ref$call)) {
        gt <- mod.ref$gt
        gfit[i,2] <- gt$n.vc; gfit[i,3] <- gt$logL; gfit[i,4] <- gt$aic; gfit[i,5] <- gt$bic; 
        gfit[i,6] <- gt$herit.PEV; gfit[i,7] <- gt$herit.VC; gfit[i,8] <- gt$Aopt; gfit[i,9] <- gt$Dopt
      }
      # create progress bar
      #pb <- txtProgressBar(min=1, max=length(set), style=3)
      #setTxtProgressBar(pb, i)
    }
    #close(pb)
    
    colnames(gfit) <- c('MODEL','n.VC','logL','AIC','BIC','heritPEV','heritVC','A-opt','D-opt')
    gfit <- data.frame(gfit)
    gfit <- data.frame(gfit,mod.list[set,-1])
    
    return(list(mod=NULL, parms=gfit))
  }
  
  if (!is.null(model)) { 
    if (!(model %in% mod.list$Model)) { 
      stop('Model selected is not on the list of pre-defined models.') 
    } 
    gfit <- matrix(NA, ncol=8, nrow=1)
    add.row <- mod.list$RepRow[model]
    add.col <- mod.list$RepCol[model]
    add.spl.row <- mod.list$splineRow[model]
    add.spl.col <- mod.list$splineCol[model]
    add.nugget <- mod.list$Nugget[model]
    type.residual <- mod.list$Resid[model]
    mod <-  gral.single(data=data, gen=gen, block=block, ibk=ibk, row=row, 
                               col=col,  cov1=cov1, cov2=cov2, resp=resp,
                               add.block=add.block,     add.ibk=add.ibk,  
                               add.row=add.row,         add.col=add.col,
                               add.spl.row=add.spl.row, add.spl.col=add.spl.col,
                               add.cov1=add.cov1,       add.cov2=add.cov2,
                               add.nugget=add.nugget,   type.gen=type.gen,
                               type.block=type.block,   type.residual=type.residual)

    #print('Details Model Fitted')
    #print(mod.list[model,])
    #return(list(mod=mod.ref, parms=gfit))
    return(mod)
    
  }
   
}