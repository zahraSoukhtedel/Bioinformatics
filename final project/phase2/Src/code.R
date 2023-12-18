
setwd("~/Desktop/bio_proj/")
library(GEOquery)
library(limma)
library(umap)
library(pheatmap)
library(gplots)
library(ggplot2)
library(reshape2)
library(plyr)
library(Rtsne)

res <- read.delim("Results/results.txt")
head(res)

# load series and platform data from GEO
series <- "GSE48558"
gset <- getGEO(series, GSEMatrix =TRUE, AnnotGPL=TRUE, destdir = "Data/")
if (length(gset) > 1) idx <- grep("GPL6244", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

gr <- c(rep("AML",13), rep("X",27), "healthy", rep("X",3), "healthy", rep("X",23),  "healthy", "X","healthy",rep("X",3),"healthy","X", rep("healthy",4),"X","healthy",rep("X",2), rep("healthy",2), rep("X",2),rep("healthy",2), "X", "healthy","X", "healthy", "X", "healthy","X", "healthy","X", "healthy", rep("X",3), "healthy", rep("X",3), "healthy",rep("X",29),rep("healthy",7), rep("AML",2), "healthy", rep("AML",3),rep("healthy",20))

ex <- exprs(gset)
# ex <- log2(ex+1) if the values were large it means it needs log
# exprs(gset) <- ex

#log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

```
### Quality Control
In each step we might have some errors, hence we should always check the quality of data before proceeding to next steps.

#### Normalization
Raw gene expression data from the high-throughput technologies must be normalized to remove technical variation so that meaningful biological comparisons can be made. Data normalization is a crucial step in the gene expression analysis as it ensures the validity of its downstream analyses.
\\
First we check if data is normalized or not

```{r, message=FALSE}
#  the first step is to obtain a box plot to check if values are normalized or not

pdf("Results/boxplot.pdf", width = 64)
boxplot(ex)
dev.off()

# if it needs normalization, we can use the code below:
# ex <-normalizeQuantiles(ex) 
# exprs(gset) <- ex


ex <-normalizeQuantiles(ex) 
exprs(gset) <- ex

pdf("Results/CorHeatmap.pdf", width=15, height=15)
pheatmap(cor(ex), color=bluered(256), labels_col = gr, labels_row=gr, border_color = NA) 
dev.off()

dim(ex)
pc <- prcomp(ex)
pdf("Results/PC.pdf")
plot(pc)
plot(pc$x[,1:2]) # this plots the first two principal components
dev.off()

names(pc)
dim(pc$x)
# colnames(pc$x)


ex.scale <- t(scale(t(ex), scale = FALSE)) #we don't want to eliminate variance
# mean(ex.scale[,1])
pc <- prcomp(ex.scale)
pdf("Results/PC_scaled.pdf")
plot(pc)
plot(pc$x[,1:2])
dev.off()

pcr <- data.frame(pc$r[,1:2], Group=gr)
head(pcr)
pdf("Results/PCA_samples.pdf")
ggplot(pcr, aes(PC1, PC2, color=Group)) + geom_point(size=3) + geom_point(size=3) + theme_bw()
dev.off()

tsne <- Rtsne(t(ex), perplexity=7, check_duplicates = FALSE)
pdf("Results/t-SNE.pdf")
plot(tsne$Y)
dev.off()

tsne_df <- data.frame(tsne$Y, Group=gr)
pdf("Results/t-SNE_samples.pdf")
ggplot(tsne_df, aes(X1, X2, color=Group)) + geom_point(size=3)+ geom_point(size=3)  + theme_bw()
dev.off()


ex_new <- cmdscale(dist(t(ex)))
pdf("Results/MDS.pdf")
plot(ex_new[,1], ex_new[,2])
dev.off()

pdf("Results/MDS_samples.pdf")
ggplot(data.frame(ex_new[,1:2], Group=gr), aes(X1, X2, color=Group)) + geom_point(size=3) + geom_point(size=3) +  theme_bw()
dev.off()

gr <- c(rep("AML",13), rep("X",27), "Granulocytes", rep("X",3), "Granulocytes", rep("X",23),  "Bcell", "X","Tcell",rep("X",3),"Granulocytes","X", "Granulocytes","Monocytes","Monocytes","Bcell","X","Tcell",rep("X",2), rep("Tcell",2), rep("X",2),rep("Tcell",2), "X", "Bcell","X", "Tcell", "X", "Bcell","X", "Tcell","X", "CD34", rep("X",3), "CD34", rep("X",3), "CD34",rep("X",29),rep("Granulocytes",7), rep("AML",2), "Tcell", rep("AML",3),rep("Bcell",7),"Tcell",rep("Monocytes",4),"Granulocytes",rep("Tcell",7))

pdf("Results/CorHeatmap_q4.pdf", width=15, height=15)
pheatmap(cor(ex), color=bluered(256), labels_col = gr, labels_row=gr, border_color = NA) 
dev.off()