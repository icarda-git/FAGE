#' Routine to fit a given spatial or non-spatial models for field trials
#'
#' \code{gral.single} Fit the spatial or non-spatial model according to the specification of conditions
#' given a provided dataset. It also calculates goodness-of-fit statisitcs to later compare
#' with other spatial model. This routine works best for replicated genotypes.
#' Some checks on data are incorporated, but data for spatial it is assumed to come in a full grid.
#' Accept a variable to identify checks (optional) that are fitted as fixed effects.
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
#' @param add.block logical to add to model block effects (default = FALSE)
#' @param add.ibk logical to fit incomplete block random effects (default = FALSE)
#' @param add.row logical to fit row within block random effects (default = FALSE)
#' @param add.col logical to fit column within block random effects (default = FALSE)
#' @param add.spl.row logical to fit splines across rows (default = FALSE)
#' @param add.spl.col logical to fit splines across columns (default = FALSE)
#' @param add.cov1 logical to fit additional covariate 1 (default = FALSE)
#' @param add.cov2 logical to fit additional covariate 2 (default = FALSE)
#' @param add.nugget logical to fit nugget random effects (only for autoregressive errors) (default = FALSE)
#' @param type.gen model assumption for genotypes: 'random' or 'fixed' (default = 'random')
#' @param type.block model assumption for full blocks: 'random' or 'fixed' (default = 'fixed')
#' @param type.residual model assumption for residual terms: 'indep' or 'ar1' (default = 'ar1')
#' @param fix.vc logical to fix variance components based on model with genotype random (default = FALSE)
#' @param vc.table dataframe with table of fixed variance components (as with start.values = TRUE)
#'  
#' @return Several objects with details of the fitted model. 
#' aov:         Wald-test (mixed model ANOVA-like table)
#' call:        String with the ASReml-R call used to fit the requested model
#' gt:          Goodness-of-fit statistics for the model evaluated are reported as summary. This includes columns: 
#'               number of variance components in the model (n.VC), log-likelihood (logL), Akaike information 
#'               criteria (AIC), Bayesian information criteria (BIC), A-optimality value, logarithm of the 
#'               D-optimality value, a heritability based on predictor error variance (herit.PEV, only for
#'               genotypes random), and heritability calculated based on variance componens (herit.VC = 
#'               var.gen/var.total), the genetic variance estimate (v.gen), the residual variance estimate 
#'               (v.res), the coefficient of variation (%) based on the mean of predictions and the residual
#'               variance, and the 5% lsd based on the average of the standard error of the difference and 
#'               infinity degrees of freedom (i.e., t0 = 1.96)
#' mod:         ASReml-R object with all information from the fitted model
#' predictions: Predictions for all genotypes, with an additional column of 'weight' is to be used on a 
#'               on a second-stage analysis. These weights are the diagonal of the inverse of the variance-
#'               covariance matrix of predictions
#' 
#' @author 
#' Salvador A. Gezan. VSN International
#' 
#' @examples
#' # Example 1: Replicated Trial - Genotype random + spatial terms 
#' library(agridat)
#' testREP <- durban.rowcol
#' testREP$bed <- as.factor(testREP$bed) 
#' testREP$row <- as.factor(testREP$row) 
#' testREP$gen <- as.factor(testREP$gen) 
#' head(testREP)
#' output.REP <- gral.single(data=testREP, gen='gen', row='row', col='bed', resp='yield',
#'                             add.block=TRUE,  add.row=TRUE,  add.col=TRUE, 
#'                             type.gen='random', type.block='fixed', type.residual='ar1')
#' output.REP$gt$herit.PEV           # Heritability-PEV
#' output.REP$aov                    # Wald Test
#' summary(output.REP$mod)$varcomp   # Variance Components
#' output.REP$gt
#' head(output.REP$predictions)
#' 
#' # Example 2: Same data but with Checks as fixed effects
#' head(testREP)
#' testREP$Ck <- 0
#' testREP$Ck[testREP$gen=='G001'] <- 1
#' testREP$Ck[testREP$gen=='G002'] <- 1
#' testREP$Ck[testREP$gen=='G003'] <- 1
#' output.Ck <- gral.single(data=testREP, gen='gen', check='Ck', row='row', col='bed', resp='yield',
#'                             add.block=TRUE,  add.row=TRUE,  add.col=TRUE, 
#'                             type.gen='random', type.block='fixed', type.residual='ar1')
#' output.Ck$gt$herit.PEV           # Heritability-PEV
#' output.Ck$aov                    # Wald Test
#' summary(output.Ck$mod)$varcomp   # Variance Components
#' output.Ck$gt
#' head(output.Ck$predictions)
#' 
#' # Example 3: Case with Fixed VC for generating Weights
#' output.REPr <- gral.single(data=testREP, gen='gen', row='row', col='bed', resp='yield',
#'                             add.block=TRUE,  add.row=TRUE,  add.col=TRUE, 
#'                             type.gen='random', type.block='fixed', type.residual='ar1',
#'                             fix.vc=TRUE)
#' summary(output.REPr$mod)$varcomp
#' head(output.REPr$predictions)

gral.single <- function(data=NULL, gen=NULL, check=NULL, block=NULL, ibk=NULL, row=NULL, 
                           col=NULL,  cov1=NULL, cov2=NULL, resp=NULL,
                           add.block=FALSE,   add.ibk=FALSE,     add.row=FALSE,  add.col=FALSE,
                           add.spl.row=FALSE, add.spl.col=FALSE, add.cov1=FALSE, add.cov2=FALSE,
                           add.nugget=FALSE,  type.gen='random', type.block='fixed', 
                           type.residual='indep', fix.vc=FALSE, vc.table=NULL){

  asreml.options(trace=FALSE, extra=3)
  n<-nrow(data)
  if (n==0) { stop('No information in data.frame provided.')}
  # Defining factors
  df <- data.frame(IDSORT=c(1:n)) 
  if (is.null(gen)) { 
     stop('No genotype column provided.')
  } else {
     df$gen <- as.factor(data[,gen])
  }
  if (!is.null(check)) { 
    df$check <- as.factor(data[,check])
    l <- length(levels(df$check))
    if (l!=2) { stop('Number of levels for check column is not equal to 2.')}
    ck.lev<-levels(df$check)
    if (ck.lev[1]!='0' & ck.lev[1]!='1') { stop('Your levels for check are not 0 or 1.') }
    if (ck.lev[2]!='0' & ck.lev[2]!='1') { stop('Your levels for check are not 0 or 1.') }
  }
  if (is.null(row) || is.null(col)) { 
    if (type.residual == 'ar1') { stop('No spatial coordinates provided.') }
  } else {
       df$row <- as.factor(data[,row])
       df$col <- as.factor(data[,col])
       nr <- length(unique(df$row))
       nc <- length(unique(df$col))
       if (n != nr*nc) { warning('Number of observations does not conform with #rows x #columns.') }
  }
  if (!is.null(block)) { 
    df$block <- as.factor(data[,block]) 
    if (type.block!='fixed' & type.block!='random') { 
      stop('Specification of blocks (random/fixed) not indicated or incorrect.') 
    }
  }
  if (!is.null(ibk)) { df$ibk <- as.factor(data[,ibk]) }
  if (!is.null(cov1)) { df$cov1 <- as.numeric(data[,cov1]) }
  if (!is.null(cov2)) { df$cov2 <- as.numeric(data[,cov2]) }
  if (is.null(resp)) { 
    stop('No response variable indicated.')
  } else {
    df$resp <- data[,resp]
  }
  if (type.gen!='fixed' & type.gen!='random') { 
    stop('Specification of genotypes (random/fixed) not indicated or incorrect.') 
  } 
  if (type.residual!='indep' & type.residual!='ar1') { 
    stop('Specification of type of residuals (indep/ar1) not indicated or incorrect.') 
  } 
  
  if (!is.null(row) && !is.null(col)) {
    df <- df[order(df$row, df$col),] # order data by row then column
  }
  
  # Code Strings for ASReml-R
  code.asr <- as.character()
  code.asr[1] <- 'asreml(fixed=resp~1'
  code.asr[2] <- 'random=~'
  code.asr[3] <- 'residual=~' 
  code.asr[4] <- 'na.action=list(x=\"include\",y=\"include\"),data=df)'
  nrand <- 0  # Number of random terms
  #mod.sp<-paste(code.asr[1],code.asr[2],code.asr[3],code.asr[4],sep=',')
  #mod.sp 
  
  # Adding cov1 (fixed)
  if (!is.null(cov1) & add.cov1) {
    code.asr[1] <- paste(code.asr[1], 'cov1', sep='+')
  }
  # Adding cov2 (fixed)
  if (!is.null(cov2) & add.cov2) {
    code.asr[1] <- paste(code.asr[1], 'cov2', sep='+')
  }
  # Adding blocks (random or fixed)
  if (!is.null(block) & add.block) {  # 
     if (type.block=='random') {
       code.asr[2] <- paste(code.asr[2], 'block', sep='+')
       nrand <- nrand + 1
     }
    if (type.block=='fixed') {
      code.asr[1] <- paste(code.asr[1], 'block', sep='+')
    }
  }
  # Adding block:ibk (random)
  if (!is.null(ibk) & add.ibk) {
    if (!add.block) { # No blocks
      code.asr[2] <- paste(code.asr[2], 'ibk', sep='+')
      nrand <- nrand + 1  
    }
    if (!is.null(block) & add.block) { # with blocks
      code.asr[2] <- paste(code.asr[2], 'block:ibk', sep='+')
      nrand <- nrand + 1  
    }
  }
  # Adding block:row (random)
  if (!is.null(row) & add.row) {
    if (!add.block) { # No blocks
      code.asr[2] <- paste(code.asr[2], 'row', sep='+')
      nrand <- nrand + 1  
    }
    if (!is.null(block) & add.block) { # with blocks
      code.asr[2] <- paste(code.asr[2], 'block:row', sep='+')
      nrand <- nrand + 1  
    }
  }
  # Adding block:col (random)
  if (!is.null(col) & add.col) {
    if (!add.block) { # No blocks
      code.asr[2] <- paste(code.asr[2], 'col', sep='+')
      nrand <- nrand + 1  
    }
    if (!is.null(block) & add.block) { # with blocks
      code.asr[2] <- paste(code.asr[2], 'block:col', sep='+')
      nrand <- nrand + 1  
    }
  } 
  # Adding spl.row (random and fixed lin)
  if (add.spl.row) {
    code.asr[1] <- paste(code.asr[1], 'lin(row)', sep='+')
    code.asr[2] <- paste(code.asr[2], 'spl(row)', sep='+')
    nrand <- nrand + 1  
  }
  # Adding spl.col (random and fixed lin)
  if (add.spl.col) {
    code.asr[1] <- paste(code.asr[1], 'lin(col)', sep='+')
    code.asr[2] <- paste(code.asr[2], 'spl(col)', sep='+')
    nrand <- nrand + 1  
  }
  # Adding gen (fixed or random)
  if (type.gen=='fixed') {
    code.asr[1] <- paste(code.asr[1], 'gen', sep='+')
  }
  if (type.gen=='random') {
    if (is.null(check)) {
      code.asr[2] <- paste(code.asr[2], 'gen', sep='+')
      nrand <- nrand + 1
    }
    if (!is.null(check)) {
      code.asr[1] <- paste(code.asr[1], "at(check,'1'):gen", sep='+')
      code.asr[2] <- paste(code.asr[2], "at(check,'0'):gen", sep='+')
      nrand <- nrand + 1
    }
  }
  # Specifying residual type and adding nugget (only ar1)
  if (type.residual=='indep') {
    code.asr[3] <- 'residual=~idv(units)'
    #code.asr[3] <- 'residual=~id(row):idv(col)'
  }
  if (type.residual=='ar1') {
    if (!is.null(row) & !is.null(col)) {
    #if (add.row & add.col) {
       if (!add.nugget) {
         code.asr[3] <- 'residual=~ar1(row):ar1v(col)'
       }
       if (add.nugget) {
         code.asr[3] <- 'residual=~ar1(row):ar1v(col)'
         code.asr[2] <- paste(code.asr[2], 'units', sep='+')
         nrand <- nrand + 1
       }
    } else {
      stop('Structure ar1 can not be fitted as add.row or add.col are FALSE.')
    }
  }
  code.asr[1] <- paste('mod.ref<-', code.asr[1], sep='') 
  if (!is.null(vc.table)) {
     code.asr[4] <- paste('G.param=vc.table,R.param=vc.table', code.asr[4], sep=',')
  }

  # Running final spatial models in ASReml-R
  code.asr[1] <- paste('mod.ref<-', code.asr[1], sep='') 
  if (nrand==0) { str.mod <- paste(code.asr[1],code.asr[3],code.asr[4],sep=',') }
  if (nrand!=0) { str.mod <- paste(code.asr[1],code.asr[2],code.asr[3],code.asr[4],sep=',') }
  
  #print(str.mod)
  mod.ref <- tryCatch({ eval(parse(text=str.mod)) }, error=function(msg){return(NA)})

  if(class(mod.ref)=='asreml') {
    if (!mod.ref$converge) { eval(parse(text='mod.ref<-update.asreml(mod.ref)')) }
    if (fix.vc) {
       vc.old <- summary(mod.ref)$varcomp
       vc.table <- data.frame(Component=rownames(vc.old),
                              Value=vc.old$component,Constraint=as.character(vc.old$bound))
       vc.table$Constraint <- rep('F',length(vc.table$Constraint))
       vc.table <- vc.table[vc.table$Component!='gen',]

       # Running same model but gen fixed 
       mod.ref <- gral.single(data=data, gen=gen, check=check, block=block, ibk=ibk, row=row, 
                           col=col, cov1=cov1, cov2=cov2, resp=resp,
                           add.block=add.block, add.ibk=add.ibk, add.row=add.row, add.col=add.col,
                           add.spl.row=add.spl.row, add.spl.col=add.spl.col, add.cov1=add.cov1, add.cov2=add.cov2,
                           add.nugget=add.nugget, type.gen='fixed', type.block=type.block, 
                           type.residual=type.residual, fix.vc=FALSE, vc.table=vc.table)$mod
    }
    
    # Obtaining predictions for models (and ANOVA)
    aov <- wald.asreml(mod.ref, denDF='algebraic', ssType='incremental')$Wald
    
    # Obtaining solutions (only random)
    if (type.gen=='fixed') { 
      sols <- NULL
      preds <- predict(mod.ref,classify='gen',vcov=TRUE)
    }
    if (type.gen=='random') {
      sols <- as.data.frame(summary(mod.ref,coef=TRUE)$coef.random) 
      if (is.null(check)) {
        sols <- sols[grep("gen",rownames(sols)),]
        preds <- predict(mod.ref,classify='gen',vcov=TRUE)
      }
      if (!is.null(check)) {
        sols <- sols[grep("gen",rownames(sols)),]
        preds <- predict(mod.ref,classify='check:gen',levels=list(check='0'),vcov=TRUE)
      }
    }
    
    # Calling the statistics
    gt <- gral.stats(object=mod.ref, solution=sols, preds=preds)

    # Adding Weights for 2-Stage
    pvals <- preds$pvals
    vcov <- as.matrix(preds$vcov)
    sel <- matrix(1, ncol=1, nrow=length(pvals$predicted.value))
    sel[is.na(pvals$predicted.value),] <- 0
    vcov <- vcov[sel==1,sel==1]
    pvals$weight[sel==1] <- diag(solve(vcov))  
    
    return(list(call=str.mod, mod=mod.ref, gt=gt, predictions=pvals, aov=aov))
  }
  if(class(mod.ref)!='asreml') {
    return(list(call=NA, mod=NA, gt=NA, predictions=NA, aov=NA))
  }
      
}