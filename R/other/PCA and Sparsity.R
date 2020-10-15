# Viz PCA 

col <- sapply(levels(mIdent),
              function(x) { as.vector(mAnnot.res1.1$DPcolor[mAnnot.res1.1$CellType == x][1])})

col.h <- sapply(levels(hIdent),
                function(x) { as.vector(mAnnot.res1.1$DPcolor[mAnnot.res1.1$CellType == x][1])})

mIdent.col <- as.vector(factor(mIdent, labels=col))
hIdent.col <- as.vector(factor(hIdent, labels=col.h))

png(filename = "PCA_mouse.png", width=1024, height=1024, bg = "transparent")
mBMMNC@dr$pca@cell.embeddings[,1:11] %>% pairs(col=mIdent.col, lower.panel=NULL)
dev.off()

png(filename = "PCA_human.png", width=1024, height=1024, bg = "transparent")
hBMMNC@dr$pca@cell.embeddings[,1:11] %>% pairs(col=hIdent.col, upper.panel=NULL)
dev.off()





# Sparsity

sparsity <- list( mouse = mBMMNC@data %>%
                    Matrix::t() %>%
                    by(INDICES = mBMMNC@meta.data$Specimen,
                       FUN= function(x) { 1 - (Matrix::nnzero(x) / prod(dim(x)))}) %>% unlist,
                  human = hBMMNC@data %>% Matrix::t() %>%
                    by(INDICES = hBMMNC@meta.data$Specimen,
                       FUN= function(x) { 1 - (Matrix::nnzero(x) / prod(dim(x)))}) %>% unlist)

at <- sparsity %>% lapply(FUN=mean) %>% unlist %>%  barplot(beside=T, las=1, ylim=c(0,1), ylab='Sparsity')
sparsity %>% stripchart(add=T, at=at, method='jitter', vertical = T, pch=1)

# mean of x mean of y 
# 0.9325577 0.9594367 
wilcox.test(x=as.vector(sparsity[[1]]), y=as.vector(sparsity[[2]]))