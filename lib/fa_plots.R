#' Routine to obtain plots from a factor analytic (FA) analysis
#'
#' Takes output from an multi-environment trial (MET) analyses fitted with an FA
#' of any order k in ASReml-R and generates a single plot associated with its output.
#' These are the multiple regressions implicit in the FA model, and the two
#' types of plots generated are: a) regression against given factor loadings
#' or b) added variable plot against given factor loadings
#' There is an option to plot prediction or BLUP value. 
#' Note that a single plot is generated, but these can be later blunded.
#' 
#' @param object output object from fitting an MET analysis for a FA model with function stageMAT()
#' @param gen levels of factor genotype to be plotted
#' @param n.fa dimension or factors to plot (maximum is k) (default = 1) 
#' @param type.resp type of response to consider for plots: 'prediction' or 'blup' (default = 'prediction') 
#' @param type.plot type of plot to generate: 'regression' or 'added.variable' (default = 'regression')
#' @param fa.rotated logical if rotated fa loadings are used (default = FALSE)
#' @param flag.extrap identifies extrapolated genotypes from MET analysis (default = FALSE)
#' @param save.plot logical to save plots to .tiff file (default = FALSE)
#'
#' @return The following objects are returned. 
#' plot.fa:  Plot of loading against response requested
#' df.gen:   Data frame with relevant information from corresponding genotype   
#' comp.gen: Data frame with the component (slope) associated with the genotyte
#' sort.fa:
#'  
#' @author 
#' Salvador A. Gezan. VSN International
#' 
#' @examples
#' # Example 1: 


fa.plots <- function(object=NULL, gen=NULL, n.fa=1,
                     type.resp='prediction', type.plot='regression', fa.rotated=FALSE,
                     flag.extrap=FALSE, save.plot=FALSE){

  require(ggplot2)
  if (is.null(object)) { 
    stop('No MET object provided.')
  } 
  if (is.null(gen)) { 
    warning('No genotype levels indicated, the first level from predictions will be plotted.')
    gen <- levels(object$predictions$gen)[1] 
  } 
  
  # Getting fa loadings
  if (fa.rotated) {
    fa.loadings <- object$fa.rot.loadings  # rotated
  } else {
    fa.loadings <- object$fa.loadings  # original
  }

  # Calculating the Genetic Values
  k <- ncol(fa.loadings)  # number of factors (k)
  s <- nrow(fa.loadings)  # number of sites/env. (k)
  
  # Obtaining predictions and BLUPs
  pvals <- as.data.frame(object$predictions)
  BLUP <- as.data.frame(summary(object$mod,coef=TRUE)$coef.random)
  
  # Obtaining overall y axis limits to unified generated plots for ease of comparible
  pred.uni.limits <- c(min(pvals$predicted.value), max(pvals$predicted.value))
  BLUP.uni.limits <- c(min(BLUP[,1]), max(BLUP[,1]))
  
  # Searching over each requested genotype
  # pred.gen <- pvals[grep(gen, pvals$gen), ] 
  # BLUP.gen <- BLUP[grep(gen, rownames(BLUP)), ] 
  # sel.comp <- BLUP.gen[grep('Comp', rownames(BLUP.gen)), ] # slopes
  pred.gen <- pvals[pvals$gen == gen, ] 
  BLUP.gen <- BLUP[grep(paste0(':gen_', gen, '$'), rownames(BLUP)), ] 
  sel.comp <- BLUP.gen[grep('_Comp\\d+:gen_', rownames(BLUP.gen)), ] # slopes
  
  sel.gen.blup <- BLUP.gen[1:(nrow(fa.loadings)),] # genetic values
  # Creating data matrix of requested genotype
  dat_gen <- data.frame(pred.gen[,1:3], BLUP.value=sel.gen.blup[,1])
  
  # Obtaining the delta (contribution) by each k
  for (i in 1:k) {
    dat_gen[, 4+i] <- fa.loadings[,i] * sel.comp$solution[i]
    colnames(dat_gen)[4+i] <- paste0('DFA', i, sep='')
  }
  
  # Generating the requested plot
  xlab <- paste0('Factor ', n.fa, ' (', round(object$cum.var[n.fa,1],2), '%)' )
  dat_gen$xresp <- fa.loadings[,n.fa] 
  if (type.resp == 'prediction') {
    ylab <- 'Predicted Values'
    uni.limits <- pred.uni.limits
    dat_gen$yresp <- dat_gen$predicted.value
  } 
  if (type.resp == 'blup') {
    ylab <- 'BLUP Values'
    uni.limits <- BLUP.uni.limits
    dat_gen$yresp <- dat_gen$BLUP.value
  } 
  
  dat_gen$extrap <- pred.gen$extrap

  if (type.plot == 'regression') {
    plot.fa <- ggplot(dat_gen) +
       theme_bw() + 
       geom_smooth(aes(y=yresp, x=xresp), formula='y ~ x', 
                   col="red", method="lm", size=0.5, se=FALSE) +
       #geom_point(aes(y=yresp, x=xresp, colour=factor(trial)), cex=2) +
       geom_text(aes(y=yresp, x=xresp, label=as.numeric(factor(trial)), colour=extrap)) +
       scale_colour_manual(values=c("black", "red")) + 
       theme(legend.position="none")+
       #guides(colour=guide_legend(nrow=4), shape="none")+
       #labs(colour='Trials')+
       labs(y=ylab, x=xlab) +
       ylim(uni.limits)
  }
  if (type.plot == 'added.variable') {
    if (n.fa == 1) {
      message('Added variable plot of dimension 1 is equivalent to regression plot.')
      plot.fa <- ggplot(dat_gen) +
         theme_bw() +
         geom_smooth(aes(y=yresp, x=xresp), formula='y ~ x',
                     col="red", method="lm", size=0.5, se=FALSE) +
         #geom_point(aes(y=yresp, x=xresp, colour=factor(trial)), cex=2) +
         geom_text(aes(y=yresp, x=xresp, label=as.numeric(factor(trial)), colour=extrap)) +
        scale_colour_manual(values=c("black", "red")) + 
        theme(legend.position="none")+
         #guides(colour=guide_legend(nrow=4), shape="none")+
         #labs(colour='Trials')+
         labs(y=ylab, x=xlab) +
         ylim(uni.limits)
    } else {
      ylab <- paste0(ylab, ' - Added Variable', sep='')
      dat_gen$yresp <- dat_gen$yresp - dat_gen[,3+n.fa]
      plot.fa <- ggplot(dat_gen) +
         theme_bw() + 
         geom_smooth(aes(y=yresp, x=xresp), formula='y ~ x', 
                     col="red", method="lm", size=0.5, se=FALSE) +
         #geom_point(aes(y=yresp, x=xresp, colour=factor(trial)), cex=2) +
         geom_text(aes(y=yresp, x=xresp, label=as.numeric(factor(trial)), colour=extrap)) +
        scale_colour_manual(values=c("black", "red")) + 
        theme(legend.position="none")+
         #guides(colour=guide_legend(nrow=4), shape="none")+
         #labs(colour='Trials')+
         labs(y=ylab, x=xlab) +
         ylim(uni.limits)
    }
    
    dat_gen$grp <- factor(ifelse(dat_gen$yresp > 0, "blue", "red"))
    dat_gen$trial2 <- paste0(levels(factor(dat_gen$trial)), " (", as.numeric(factor(dat_gen$trial)), ")")
    sort.fa <- ggbarplot(dat_gen, x = "trial2", y = "yresp",
                         fill = "grp",               # change fill color by grp
                         color = "white",            # Set bar border colors to white
                         palette = "jco",            # jco journal color palett. see ?ggpar
                         sort.val = "desc",          # Sort the value in descending order
                         sort.by.groups = FALSE,     # Don't sort inside each group
                         x.text.angle = 90,          # Rotate vertically x axis texts
                         ylab = "BLUE Values - Added Variable",
                         xlab = "",
                         legend = "none",
                         rotate = TRUE,
                         ggtheme = theme_minimal())
  }

  # Save Plot    
  if(save.plot == TRUE){
    fname <- paste0('Plot_fa', n.fa, '.tiff', sep='')
    ggplot2::ggsave(filename=fname, plot=plot.fa, width=8, height=6, units="in", dpi=100)
  }
  
  return(list(plot.fa=plot.fa, df.gen=dat_gen, comp.gen=sel.comp, sort.fa=sort.fa))
  
}