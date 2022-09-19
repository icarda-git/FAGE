#' Routine to perform final step of a two-stage MET analyses
#'
#' \code{stageMET} Performs the genetic analysis of an MET dataset corresponding
#' to the final step of a two-stage analysis, where input corresponds to meand (or 
#' adjusted values) of Evaluates and verifies the data originating from several trials with the 
#' aim of determining levels of: connectivity, variability, etc., which are reported as 
#' statistics. Only non-NA observations are considered in the reports. Note that trial is 
#' always considered a fixed effect.
#
#' @param data dataframe with all relevant columns for MET analyses.
#' @param gen factor name for genotypes (or treatments)
#' @param trial factor name for trial (or environment) 
#' @param resp column name for the response variable to evaluate
#' @param weight column name for the weight of response (default = 1)
#' @param type.gen model assumption for genotype effects: 'random' or 'fixed' (default = 'random')
#' @param type.trial model assumption for trial effects: 'random' or 'fixed' (default = 'fixed')
#' @param vc.model variance-covariance model to fit: 'diag', 'corv', 'corh', 'fa1', 'fa2', 
#'                 'fa3', 'fa4', 'fa5', 'fa6', 'fa7', 'fa8', 'fa9', 'fa10', 'corgh' 
#'                 (default = 'corh') (only for type.gen =' random')
#' @param rotation type of rotation implemented: 'svd' or 'varimax' (default = 'svd')
#' @param vcovp logical if TRUE the variance-covariance matrix of predictions is reported
#'              (default = FALSE)   
#'
#' @return Several objects with reports of the MET analysis. 
#' call:        String with the ASReml-R call used to fit the requested model
#' mod:         ASReml-R object with all information from the fitted model
#' predictions: Predictions for all genotypes across all sites, with their standard
#'                error and reliability
#' vcov.M:      GxE variance-covariance matrix 
#' corr.M:      GxE correlation matrix
#' vcov.P:      Variance-covariance matrix of predictions
#'               
#' @author 
#' Salvador A. Gezan. VSN International
#' 
#' @examples
#' # Example 1: 
#' 

stageMET <- function(data=NULL, gen=NULL, trial=NULL, resp=NULL, weight=NULL,
                     type.gen='random', type.trial='fixed', vc.model='corh',
                     rotation='svd', vcovp='FALSE'){

  asreml.options(trace=FALSE)
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
    s <- length(levels(df$trial))
    if (s <= 1) { stop('Only 1 trial on the data, this is not an MET.')}
  }
  if (is.null(resp)) { 
    stop('No response column provided.')
  } else {
    df$resp <- data[,resp]
  }
  if (is.null(weight)) { 
    message('No weight column provided, all weights = 1.')
    df$weight <- NA
    df$weight[!is.na(df$resp)] <- 1
    if (sum(is.na(df$weight))>0) { 
      message('Weight column has some missing values, respective records were eliminated.')
      # Eliminating NA on weights
      df <- df[!is.na(df$weight),]
    }
  } else {
    df$weight <- data[,weight]
    if (sum(is.na(df$weight))>0) { 
      message('Weight column has some missing values, respective records were eliminated.')
      # Eliminating NA on weights
      df <- df[!is.na(df$weight),]
    }
  }
  
  # df: IDSORT, gen, trial, resp, weight
  # Code Strings for ASReml-R
  code.asr <- as.character()
  code.asr[1] <- 'asreml::asreml(fixed=resp~1'
  code.asr[2] <- 'random=~'
  code.asr[3] <- 'weights=weight'
  code.asr[4] <- 'family=asr_gaussian(dispersion=1)'
  code.asr[5] <- 'na.action=list(x="include",y="include"),workspace=128e06,data=df)'
  nrand <- 0  # Number of random terms
  
  # Adding gen (fixed or random)
  if (type.gen=='fixed') {
    if (type.trial=='fixed') { 
      code.asr[1] <- paste(code.asr[1], 'gen+trial+trial:gen', sep='+') 
    }
    if (type.trial=='random') { 
      code.asr[1] <- paste(code.asr[1], 'gen', sep='+')
      code.asr[2] <- paste(code.asr[2], 'trial+trial:gen', sep='+')
      nrand <- nrand + 1
    }
  }
  # gen random
  if (type.gen=='random') {
    nrand <- nrand + 1
    if (type.trial=='fixed') { 
      code.asr[1] <- paste(code.asr[1], 'trial', sep='+') 
    }
    if (type.trial=='random') { 
      code.asr[2] <- paste(code.asr[2], 'trial', sep='+')
    }
    if (vc.model=='diag') {
      code.asr[2] <- paste(code.asr[2], 'diag(trial):id(gen)', sep='+')
    }
    if (vc.model=='corv') {
      code.asr[2] <- paste(code.asr[2], 'corv(trial):id(gen)', sep='+')
    }
    if (vc.model=='corh') {
      code.asr[2] <- paste(code.asr[2], 'corh(trial):id(gen)', sep='+')
    }
    if (vc.model=='fa1') {
      if (s<=2) { message('This Factor Analytic 1 analysis is over-parametrized.')}
      code.asr[2] <- paste(code.asr[2], 'fa(trial,1):id(gen)', sep='+')
    }
    if (vc.model=='fa2') {
      if (s<=4) { message('This Factor Analytic 2 analysis is over-parametrized.')}
      code.asr[2] <- paste(code.asr[2], 'fa(trial,2):id(gen)', sep='+')
    }
    if (vc.model=='fa3') {
      if (s<=6) { message('This Factor Analytic 3 analysis is over-parametrized.')}
      code.asr[2] <- paste(code.asr[2], 'fa(trial,3):id(gen)', sep='+')
    }
    if (vc.model=='fa4') {
      if (s<=8) { message('This Factor Analytic 4 analysis is over-parametrized.')}
      code.asr[2] <- paste(code.asr[2], 'fa(trial,4):id(gen)', sep='+')
    }
    if (vc.model=='fa5') {
      if (s<=10) { message('This Factor Analytic 5 analysis is over-parametrized.')}
      code.asr[2] <- paste(code.asr[2], 'fa(trial,5):id(gen)', sep='+')
    }
    if (vc.model=='fa6') {
      if (s<=12) { message('This Factor Analytic 6 analysis is over-parametrized.')}
      code.asr[2] <- paste(code.asr[2], 'fa(trial,6):id(gen)', sep='+')
    }
    if (vc.model=='fa7') {
      if (s<=14) { message('This Factor Analytic 7 analysis is over-parametrized.')}
      code.asr[2] <- paste(code.asr[2], 'fa(trial,7):id(gen)', sep='+')
    }
    if (vc.model=='fa8') {
      if (s<=16) { message('This Factor Analytic 8 analysis is over-parametrized.')}
      code.asr[2] <- paste(code.asr[2], 'fa(trial,8):id(gen)', sep='+')
    }    
    if (vc.model=='fa9') {
      if (s<=18) { message('This Factor Analytic 9 analysis is over-parametrized.')}
      code.asr[2] <- paste(code.asr[2], 'fa(trial,9):id(gen)', sep='+')
    }
    if (vc.model=='fa10') {
      if (s<=20) { message('This Factor Analytic 10 analysis is over-parametrized.')}
      code.asr[2] <- paste(code.asr[2], 'fa(trial,10):id(gen)', sep='+')
    }
    if (vc.model=='corgh') {
      code.asr[2] <- paste(code.asr[2], 'corgh(trial):id(gen)', sep='+')
    }
  }
  
  # Running final MET models in ASReml-R
  code.asr[1] <- paste('mod.ref<-', code.asr[1], sep='') 
  if (nrand==0) { str.mod <- paste(code.asr[1],code.asr[3],code.asr[4],code.asr[5],sep=',') }
  if (nrand!=0) { str.mod <- paste(code.asr[1],code.asr[2],code.asr[3],code.asr[4],code.asr[5],sep=',') } 
  #print(str.mod)
  eval(parse(text=str.mod) )
  if (!mod.ref$converge) { eval(parse(text='mod.ref<-asreml::update.asreml(mod.ref)')) }
  
  # Obtaining predictions (and ebvs) for models
  if (vcovp) {
    opred <- predict(mod.ref,classify='trial:gen', pworkspace=1e08, sed=FALSE, vcov=TRUE)
    pvals <- opred$pvals
    vcovp <- opred$vcov
    #opred <- predict(mod.ref,classify='trial', pworkspace=1e08, sed=FALSE, vcov=FALSE)
    #pvals.trial <- as.data.frame(opred$pvals)
  } else {
    pvals <- predict(mod.ref,classify='trial:gen', pworkspace=1e08, sed=FALSE, vcov=FALSE)$pvals
    vcovp <- NULL
    #opred <- predict(mod.ref,classify='trial', pworkspace=1e08, sed=FALSE, vcov=FALSE)
    #pvals.trial <-  as.data.frame(opred$pvals)
  }
  # Obtaining Extrapolation TRUE/FALSE
  pvals <- as.data.frame(pvals)
  pvals$extrap <- as.logical(1-as.vector(table(df$gen, df$trial)))
  
  # Obtaining some statistics
  gfit <- matrix(NA, ncol=4, nrow=1)
  gfit[1] <- nrow(summary(mod.ref)$varcomp)
  gfit[2] <- summary(mod.ref)$loglik; gfit[3] <- summary(mod.ref)$aic; gfit[4] <- summary(mod.ref)$bic
  colnames(gfit) <-  c('n.VC','logL','AIC','BIC')

  # Extracting Variance-Covariance\Correlation Matrix
  if (type.gen=='fixed') { VCOV <- NA; CORR <- NA }
  if (type.gen=='random') {
    vc <- summary(mod.ref)$varcomp
    VCOV <- matrix(0, ncol=s, nrow=s)
    CORR <- matrix(0, ncol=s, nrow=s)
    diag(CORR) <- rep(1,s)
    if (vc.model=='diag') {
      vc <-  vc[grep('trial:gen',rownames(vc)),]
      diag(VCOV) <- vc[,1]
    }
    if (vc.model=='corv') {
      vc <-  vc[grep('trial:gen',rownames(vc)),]
      CORR <- matrix(1, ncol=s, nrow=s)
      CORR <- vc[1,1]*CORR
      diag(CORR) <- rep(1,s)
      D <- rep(vc[2,1],s)
      VCOV <- diag(sqrt(D)) %*% CORR %*% diag(sqrt(D))
    }
    if (vc.model=='corh') {
      vc <-  vc[grep('trial:gen',rownames(vc)),]
      CORR <- matrix(1, ncol=s, nrow=s)
      CORR <- vc[1,1]*CORR
      diag(CORR) <- rep(1,s)
      D <- vc[2:(s+1),1]
      VCOV <- diag(sqrt(D)) %*% CORR %*% diag(sqrt(D))
    }
    if (vc.model=='fa1') {
      vc.var <-  vc[grep('!var',rownames(vc)),]
      vc.fa1 <-  vc[grep('!fa1',rownames(vc)),]
      R <- as.matrix(vc.var[,1]); L <- as.matrix(vc.fa1[,1])
      VCOV <- L %*% t(L) + diag(R) # L is a row-vector
      CORR <- cov2cor(VCOV)
    }
    if (vc.model=='fa2') {
      vc.var <-  vc[grep('!var',rownames(vc)),]
      vc.fa1 <-  vc[grep('!fa1',rownames(vc)),]
      vc.fa2 <-  vc[grep('!fa2',rownames(vc)),]
      R <- vc.var[,1]; L1 <- vc.fa1[,1]; L2 <- vc.fa2[,1]
      L <- cbind(L1, L2)
      VCOV <- L %*% t(L) + diag(R)
      CORR <- cov2cor(VCOV)
    }
    if (vc.model=='fa3') {
      vc.var <-  vc[grep('!var',rownames(vc)),]
      vc.fa1 <-  vc[grep('!fa1',rownames(vc)),]
      vc.fa2 <-  vc[grep('!fa2',rownames(vc)),]
      vc.fa3 <-  vc[grep('!fa3',rownames(vc)),]
      R <- vc.var[,1]
      L1 <- vc.fa1[,1]; L2 <- vc.fa2[,1]; L3 <- vc.fa3[,1]
      L <- cbind(L1, L2, L3)
      VCOV <- L %*% t(L) + diag(R)
      CORR <- cov2cor(VCOV)
    }
    if (vc.model=='fa4') {
      vc.var <-  vc[grep('!var',rownames(vc)),]
      vc.fa1 <-  vc[grep('!fa1',rownames(vc)),]
      vc.fa2 <-  vc[grep('!fa2',rownames(vc)),]
      vc.fa3 <-  vc[grep('!fa3',rownames(vc)),]
      vc.fa4 <-  vc[grep('!fa4',rownames(vc)),]
      R <- vc.var[,1]
      L1 <- vc.fa1[,1]; L2 <- vc.fa2[,1]; L3 <- vc.fa3[,1]; L4 <- vc.fa4[,1]
      L <- cbind(L1, L2, L3, L4)
      VCOV <- L %*% t(L) + diag(R)
      CORR <- cov2cor(VCOV)
    }
    if (vc.model=='fa5') {
      vc.var <-  vc[grep('!var',rownames(vc)),]
      vc.fa1 <-  vc[grep('!fa1',rownames(vc)),]
      vc.fa2 <-  vc[grep('!fa2',rownames(vc)),]
      vc.fa3 <-  vc[grep('!fa3',rownames(vc)),]
      vc.fa4 <-  vc[grep('!fa4',rownames(vc)),]
      vc.fa5 <-  vc[grep('!fa5',rownames(vc)),]
      R <- vc.var[,1]
      L1 <- vc.fa1[,1]; L2 <- vc.fa2[,1]; L3 <- vc.fa3[,1]; L4 <- vc.fa4[,1]
      L5 <- vc.fa5[,1]
      L <- cbind(L1, L2, L3, L4, L5)
      VCOV <- L %*% t(L) + diag(R)
      CORR <- cov2cor(VCOV)
    }
    if (vc.model=='fa6') {
      vc.var <-  vc[grep('!var',rownames(vc)),]
      vc.fa1 <-  vc[grep('!fa1',rownames(vc)),]
      vc.fa2 <-  vc[grep('!fa2',rownames(vc)),]
      vc.fa3 <-  vc[grep('!fa3',rownames(vc)),]
      vc.fa4 <-  vc[grep('!fa4',rownames(vc)),]
      vc.fa5 <-  vc[grep('!fa5',rownames(vc)),]
      vc.fa6 <-  vc[grep('!fa6',rownames(vc)),]
      R <- vc.var[,1]
      L1 <- vc.fa1[,1]; L2 <- vc.fa2[,1]; L3 <- vc.fa3[,1]; L4 <- vc.fa4[,1]
      L5 <- vc.fa5[,1]; L6 <- vc.fa6[,1]
      L <- cbind(L1, L2, L3, L4, L5, L6)
      VCOV <- L %*% t(L) + diag(R)
      CORR <- cov2cor(VCOV)
    }
    if (vc.model=='fa7') {
      vc.var <-  vc[grep('!var',rownames(vc)),]
      vc.fa1 <-  vc[grep('!fa1',rownames(vc)),]
      vc.fa2 <-  vc[grep('!fa2',rownames(vc)),]
      vc.fa3 <-  vc[grep('!fa3',rownames(vc)),]
      vc.fa4 <-  vc[grep('!fa4',rownames(vc)),]
      vc.fa5 <-  vc[grep('!fa5',rownames(vc)),]
      vc.fa6 <-  vc[grep('!fa6',rownames(vc)),]
      vc.fa7 <-  vc[grep('!fa7',rownames(vc)),]
      R <- vc.var[,1]
      L1 <- vc.fa1[,1]; L2 <- vc.fa2[,1]; L3 <- vc.fa3[,1]; L4 <- vc.fa4[,1]
      L5 <- vc.fa5[,1]; L6 <- vc.fa6[,1]; L7 <- vc.fa7[,1]
      L <- cbind(L1, L2, L3, L4, L5, L6, L7)
      VCOV <- L %*% t(L) + diag(R)
      CORR <- cov2cor(VCOV)
    }
    if (vc.model=='fa8') {
      vc.var <-  vc[grep('!var',rownames(vc)),]
      vc.fa1 <-  vc[grep('!fa1',rownames(vc)),]
      vc.fa2 <-  vc[grep('!fa2',rownames(vc)),]
      vc.fa3 <-  vc[grep('!fa3',rownames(vc)),]
      vc.fa4 <-  vc[grep('!fa4',rownames(vc)),]
      vc.fa5 <-  vc[grep('!fa5',rownames(vc)),]
      vc.fa6 <-  vc[grep('!fa6',rownames(vc)),]
      vc.fa7 <-  vc[grep('!fa7',rownames(vc)),]
      vc.fa8 <-  vc[grep('!fa8',rownames(vc)),]
      R <- vc.var[,1]
      L1 <- vc.fa1[,1]; L2 <- vc.fa2[,1]; L3 <- vc.fa3[,1]; L4 <- vc.fa4[,1]
      L5 <- vc.fa5[,1]; L6 <- vc.fa6[,1]; L7 <- vc.fa7[,1]; L8 <- vc.fa8[,1]
      L <- cbind(L1, L2, L3, L4, L5, L6, L7, L8)
      VCOV <- L %*% t(L) + diag(R)
      CORR <- cov2cor(VCOV)
    }
    if (vc.model=='fa9') {
      vc.var <-  vc[grep('!var',rownames(vc)),]
      vc.fa1 <-  vc[grep('!fa1',rownames(vc)),]
      vc.fa2 <-  vc[grep('!fa2',rownames(vc)),]
      vc.fa3 <-  vc[grep('!fa3',rownames(vc)),]
      vc.fa4 <-  vc[grep('!fa4',rownames(vc)),]
      vc.fa5 <-  vc[grep('!fa5',rownames(vc)),]
      vc.fa6 <-  vc[grep('!fa6',rownames(vc)),]
      vc.fa7 <-  vc[grep('!fa7',rownames(vc)),]
      vc.fa8 <-  vc[grep('!fa8',rownames(vc)),]
      vc.fa9 <-  vc[grep('!fa9',rownames(vc)),]
      R <- vc.var[,1]
      L1 <- vc.fa1[,1]; L2 <- vc.fa2[,1]; L3 <- vc.fa3[,1]; L4 <- vc.fa4[,1]
      L5 <- vc.fa5[,1]; L6 <- vc.fa6[,1]; L7 <- vc.fa7[,1]; L8 <- vc.fa8[,1]
      L9 <- vc.fa9[,1]
      L <- cbind(L1, L2, L3, L4, L5, L6, L7, L8, L9)
      VCOV <- L %*% t(L) + diag(R)
      CORR <- cov2cor(VCOV)
    }
    if (vc.model=='fa10') {
      vc.var <-  vc[grep('!var',rownames(vc)),]
      vc.fa1 <-  vc[grep('!fa1 ',rownames(vc)),]
      vc.fa2 <-  vc[grep('!fa2',rownames(vc)),]
      vc.fa3 <-  vc[grep('!fa3',rownames(vc)),]
      vc.fa4 <-  vc[grep('!fa4',rownames(vc)),]
      vc.fa5 <-  vc[grep('!fa5',rownames(vc)),]
      vc.fa6 <-  vc[grep('!fa6',rownames(vc)),]
      vc.fa7 <-  vc[grep('!fa7',rownames(vc)),]
      vc.fa8 <-  vc[grep('!fa8',rownames(vc)),]
      vc.fa9 <-  vc[grep('!fa9',rownames(vc)),]
      vc.fa10 <- vc[grep('!fa10',rownames(vc)),]
      R <- vc.var[,1]
      L1 <- vc.fa1[,1]; L2 <- vc.fa2[,1]; L3 <- vc.fa3[,1]; L4 <- vc.fa4[,1]
      L5 <- vc.fa5[,1]; L6 <- vc.fa6[,1]; L7 <- vc.fa7[,1]; L8 <- vc.fa8[,1]
      L9 <- vc.fa9[,1]; L10 <- vc.fa10[,1]
      L <- cbind(L1, L2, L3, L4, L5, L6, L7, L8, L9, L10)
      VCOV <- L %*% t(L) + diag(R)
      CORR <- cov2cor(VCOV)
    }
    if (vc.model=='corgh') {
      vc.corr <-  vc[grep('.cor',rownames(vc)),]
      vc.var <-  vc[-grep('.cor',rownames(vc)),]
      k <- 1
      for (i in 1:s) {
        for (j in 1:i) {
          if (i != j) { 
            CORR[i,j] <- vc.corr[k,1]; CORR[j,i] <- vc.corr[k,1]
            k <- k+1
          }
        }
      }
      D <- vc.var[1:s,1]
      VCOV <- diag(sqrt(D)) %*% CORR %*% diag(sqrt(D))
    }
    # Calculating Variance Explained for fa1 to fa10 & report loadings
    fa.loadings <- NA
    if (vc.model=='fa1' | vc.model=='fa2' | vc.model=='fa3' | vc.model=='fa4' | vc.model=='fa5' |
        vc.model=='fa6' | vc.model=='fa7' | vc.model=='fa8' | vc.model=='fa9' | vc.model=='fa10'){
      DEN <- sum(diag(VCOV))
      k <- ncol(L)
      var.exp <- matrix(NA, ncol=2, nrow=k)
      for (j in 1:k) {
        var.exp[j,1] <- 100*sum(diag(L[,j] %*% t(L[,j])))/DEN
        var.exp[j,2] <- 100*sum(diag(L[,1:j] %*% t(L[,1:j])))/DEN
      }
      fa.loadings <- L
      colnames(fa.loadings) <- paste0('FA',c(1:k))
      rownames(var.exp) <- colnames(fa.loadings)
      colnames(var.exp) <- c('var.exp%', 'cum.var.exp%')
      if (vc.model=='fa1') {
        fa.rot.loadings <- fa.loadings  # no need for FA1
      } else {
        if (rotation=='varimax') {
          fa.rot.loadings <- varimax(fa.loadings)$loadings[]
        }
        if (rotation=='svd') {
          ss <- svd(fa.loadings)
          fa.rot.loadings <- -fa.loadings %*% ss$v
        }
      }
      rownames(fa.loadings) <- levels(df$trial)
      rownames(fa.rot.loadings) <- levels(df$trial)
      colnames(fa.rot.loadings) <- colnames(fa.loadings)
      var.exp.site <- matrix(0, nrow=nrow(L), ncol=k+1)
      for (i in 1:k) {
        var.exp.site[,i] <- 100*diag(fa.rot.loadings[,i] %*% t(fa.rot.loadings[,i]))/diag(VCOV)
      }
      if (vc.model=='fa1') {
        var.exp.site[,(k+1)] <- var.exp.site[,k]
      } else {
        var.exp.site[,(k+1)] <- rowSums(var.exp.site[,1:k])
      }
      
      colnames(var.exp.site) <- c(colnames(fa.rot.loadings),'cum.var.exp%')
      rownames(var.exp.site) <- levels(df$trial)
    } else {
      var.exp <- NA
      var.exp.site <- NA
      fa.loadings <- NA
      fa.rot.loadings <- NA
    }
    
    #if (vc.model=='fa1' | vc.model=='fa2' | vc.model=='fa3' | vc.model=='fa4' | vc.model=='fa5' |
    #    vc.model=='fa6' | vc.model=='fa7' | vc.model=='fa8' | vc.model=='fa9' | vc.model=='fa10'){
    #    gfit[5] <- 100*sum(diag(L %*% t(L)))/sum(diag(VCOV))  # var.exp
    #    fa.loadings <- L
    #}
    colnames(VCOV) <- levels(df$trial)
    colnames(CORR) <- levels(df$trial)
    rownames(VCOV) <- levels(df$trial)
    rownames(CORR) <- levels(df$trial)
  }
  
  return(list(call=str.mod, mod=mod.ref, predictions=pvals, # breeding.values=ebvs, 
              gfit=gfit, vcov.M=VCOV, corr.M=CORR, vcov.P=vcovp, 
              fa.loadings=fa.loadings, fa.rot.loadings=fa.rot.loadings, 
              cum.var=round(var.exp,3), cum.var.site=round(var.exp.site,3)))
  
}