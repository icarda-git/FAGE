



fieldmap <- function (x, column, colorscale = NULL, textColumn = NULL, borderColumns = c("EXPT", 
                                                                             "BLOCK"), asFactor = FALSE, colors = NULL, rows = "ROW", 
          ranges = "RANGE", ...) 
{
  if (length(borderColumns) > 3) 
    stop("The 'borderColumn' argument can only take up to three factors. \n", 
         "Please check the argument.")
  if (!all(borderColumns %in% colnames(x))) {
    missingColumns <- setdiff(borderColumns, colnames(x))
    warning("The following borders column(s) is/are not found in 'x': ", 
            paste(missingColumns, collapse = ", "))
    borderColumns <- borderColumns[-which(borderColumns %in% 
                                            missingColumns)]
  }
  if (!is.null(colorscale) && (length(colorscale) != 2 | !is.numeric(colorscale))) {
    warning("The 'colorscale' is provided in the wrong format, it won't be used.")
    colorscale <- NULL
  }
  if (!rows %in% colnames(x)) 
    stop("The row-indicator 'rows' is not found as column of 'x'.")
  if (!ranges %in% colnames(x)) 
    stop("The range-indicator 'ranges' is not found as column of 'x'.")
  x$ROW <- x[, rows]
  if (!is.numeric(x$ROW)) 
    x$ROW <- as.numeric(as.character(x$ROW))
  x$RANGE <- x[, ranges]
 
  
  ##comment no ---need to develop the fillDesign function
   #if (!is.numeric(x$RANGE)) 
    #x$RANGE <- as.numeric(as.character(x$RANGE))
  #x$tmpColumn <- "FIELD1"
 # x <- fillDesign(x, fields = "tmpColumn")
 ###comment end
  
  
  
   data <- x
  data$Y <- data[, column]
  if (asFactor) 
    data$Y <- factor(data$Y)
  if (is.numeric(data$Y)) {
    if (is.null(colors)) 
      colors <- rev(c("#A50026", "#D73027", "#F46D43", 
                      "#FDAE61", "#FEE090", "#FFFFBF", "#E0F3F8", "#ABD9E9", 
                      "#74ADD1", "#4575B4", "#313695"))
    if (!is.null(colorscale)) {
      breakpoints <- seq(colorscale[1], colorscale[2], 
                         length = 11)
    }
    else {
      if (all(c(-1, 1) %in% sign(data$Y))) {
        M <- max(abs(data$Y), na.rm = TRUE)
        breakpoints <- seq(-M, M, length = 11)
      }
      else {
        breakpoints <- seq(min(data$Y, na.rm = TRUE), 
                           max(data$Y, na.rm = TRUE), length = 11)
      }
    }
    colorkey <- TRUE
    border <- "white"
  }
  else {
    levels <- levels(factor(data$Y))
    if (is.null(colors)) 
      colors <- colorRampPalette(c("cyan", "green", "blue", 
                                   "red", "yellow"))
    breakpoints <- 0:length(levels)
    colorkey <- list(labels = list(labels = levels, at = 0:length(levels) + 
                                     0.5))
    border <- "white"
  }
  if (!is.null(textColumn)) {
    data$textColumn <- data[, textColumn]
    data$textColumn <- ifelse(is.na(data$textColumn), "", 
                              as.character(data$textColumn))
  }
  nBorderColumns <- length(borderColumns)
  if (nBorderColumns >= 1) {
    data$border1 <- factor(data[, borderColumns[1]])
  }
  if (nBorderColumns >= 2) {
    data$border2 <- factor(data[, borderColumns[2]])
  }
  if (nBorderColumns == 3) {
    data$border3 <- factor(data[, borderColumns[3]])
  }
  require(lattice)
  levelplot(Y ~ ROW * RANGE, data = data, contour = FALSE, 
            colorkey = colorkey, at = breakpoints, col.regions = colors, 
            border = border, panel = function(x, y, z, ...) {
              panel.levelplot(x, y, z, ...)
              if (!is.null(textColumn)) {
                panel.text(x = x, y = y, data$textColumn)
              }
              if ("border1" %in% colnames(data)) {
                for (iBorder1 in levels(data$border1)) {
                  sub1 <- subset(data, border1 == iBorder1)
                  panel.rect(min(sub1$ROW - 0.5), min(sub1$RANGE - 
                                                        0.5), max(sub1$ROW + 0.5), max(sub1$RANGE + 
                                                                                         0.5), lwd = 2, lty = 1, border = "black")
                }
              }
              if (all(c("border1", "border2") %in% colnames(data))) {
                for (iBorder1 in levels(data$border1)) {
                  for (jBorder2 in levels(data$border2)) {
                    sub2 <- subset(data, border1 == iBorder1 & 
                                     border2 == jBorder2)
                    if (nrow(sub2) > 0) {
                      panel.rect(min(sub2$ROW - 0.5), min(sub2$RANGE - 
                                                            0.5), max(sub2$ROW + 0.5), max(sub2$RANGE + 
                                                                                             0.5), lwd = 1, lty = 1, border = "black")
                    }
                  }
                }
              }
              if (all(c("border1", "border2", "border3") %in% colnames(data))) {
                for (iBorder1 in levels(data$border1)) {
                  for (jBorder2 in levels(data$border2)) {
                    for (kBorder3 in levels(data$border3)) {
                      sub3 <- subset(data, border1 == iBorder1 & 
                                       border2 == jBorder2 & border3 == kBorder3)
                      if (nrow(sub3) > 0) {
                        panel.rect(min(sub3$ROW - 0.5), min(sub3$RANGE - 
                                                              0.5), max(sub3$ROW + 0.5), max(sub3$RANGE + 
                                                                                               0.5), lwd = 1, lty = 2, border = "black")
                      }
                    }
                  }
                }
              }
            }, ...)
}