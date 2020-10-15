# Function: Embedding Heatmap, 2D histogram of expression

embeddingHeatmap <- function(SeuratObject, GOI, z.range = NULL, use.internal.col = TRUE, col=list(low='lightgrey', mid='blue', high='green'), n = 20, x.axis.symrange = 40, y.axis.symrange = 40, ...){
  
  # SeuratObject: an object created using the Seurat package.
  # GOI: gene of interest (string)
  # use.internal.col: escape-option to plot full z-range using internally specified color scheme
  # col: list containing low, mid and high color values to be passed on to marray::maPalette()
  # n: number of shades in color palette (approximate value, maPalette does some trimming here and then, according to internal aesthetics)
  # speclim: specify limits of x- and y-axis labels
  # ...: additional options passed on to image()
  
  # specify color palette
  require(marray)
  palette = maPalette(low = col$low, mid = col$mid, high = col$high, k = n)
  
  # extract data from Seurat object
  x <- SeuratObject@dr$tsne@cell.embeddings[,1]
  y <- SeuratObject@dr$tsne@cell.embeddings[,2]
  if (is.character(GOI)){
   z <- SeuratObject@data[SeuratObject@data@Dimnames[[1]] %in% GOI, ]
  } else { # GOI can also be a vector of expression values, in this case, replace z with GOI
    z = GOI
    GOI = ''
  }
  
  # determine x-y range
  x.range = 1.04*range(x)
  y.range = 1.04*range(y)
  
  # define intervals
  x.breaks <- seq.int(x.range[1], x.range[2], length.out = round(diff(x.range))/3)
  y.breaks <- seq.int(y.range[1], y.range[2], length.out = round(diff(y.range))/3)
  
  # determine x-y-bin for each cell
  x.int <- cut(x = x, breaks = x.breaks, ordered_result = T)
  y.int <- cut(x = y, breaks = y.breaks, ordered_result = T)
  
  # average across area
  emb.hm <- tapply(z, list(x.int, y.int), mean, na.rm=T)
  
  # plot
  if (use.internal.col) { image(emb.hm, axes=F, col=palette, main=GOI, useRaster=F, ...); box();
  } else {
    # shrink color-range to 1st-99th percentile
    if(is.null(z.range)) z.range <- quantile(z, c(.01, .99))
    # define intervals
    z.breaks <- seq(z.range[1], z.range[2], length.out = n+1)
    # determine category for each cell
    g <- matrix(findInterval(x = emb.hm, vec = z.breaks, all.inside = TRUE), nrow = nrow(emb.hm))
    image(g, axes=F, col=palette, main=GOI, useRaster=F, ...); box();
    rm(z.breaks)
  }
  
    # add x-axis 
    
    # find symmetric limits rounded to the next decade of the original range and transform to 0-1 range 
    if (is.null(x.axis.symrange)) x.axis.symrange <- min(abs(signif(x.range*0.96, -1)))
    # position of zero if original range is not symmetric
    x.axis.zero <- abs(x.range[1]) / diff(x.range)
    # position of the lower label
    x.axis.low <- x.axis.zero * (1 - x.axis.symrange / abs(x.range[1]))
    # position of the upper axis label
    x.axis.high <- x.axis.zero + (1 - x.axis.zero) * x.axis.symrange / x.range[2]
    # plot axis
    axis(side = 1, las=1, line = 0, at = c(x.axis.low,  x.axis.zero, x.axis.high), labels = c(-x.axis.symrange, 0, x.axis.symrange));
    mtext('Dimension 1', 1, line = 2.5, cex = 0.75) # x
    
    rm(x.axis.symrange, x.axis.zero, x.axis.low, x.axis.high)
    
    # add y-axis 
    
    # find symmetric limits rounded to the next decade of the original range and transform to 0-1 range 
    if (is.null(y.axis.symrange)) y.axis.symrange <- min(abs(signif(y.range*0.96, -1)))
    # position of zero if original range is not symmetric
    y.axis.zero <- abs(y.range[1]) / diff(y.range)
    # position of the lower label
    y.axis.low <- y.axis.zero * (1 - y.axis.symrange / abs(y.range[1]))
    # position of the upper axis label
    y.axis.high <- y.axis.zero + (1 - y.axis.zero) * y.axis.symrange / y.range[2]
    # plot axis
    axis(side = 2, las=1, line = 0, at = c(y.axis.low,  y.axis.zero, y.axis.high), labels = c(-y.axis.symrange, 0, y.axis.symrange));
    mtext('Dimension 2', 2, line = 2.5, cex = 0.75) # y
    
    rm(y.axis.symrange, y.axis.zero, y.axis.low, y.axis.high)
    
  rm(x, y, z, x.range, y.range, z.range, x.breaks, y.breaks, x.int, y.int, palette, emb.hm)
}