setwd("~/Desktop/bio_proj/")
library(GEOquery)
library(limma)
library(umap)
library(pheatmap)
library(gplots)
library(ggplot2)
library(reshape2)
library(plyr)

series <- "GSE48558" 
gset <- getGEO(series, GSEMatrix =TRUE, getGPL=FALSE, destdir = "Data/")

if (length(gset) > 1) idx <- grep("GPL6244", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

#gr <- c("CD34", rep("BM",10), rep("CD34",7), rep("AML",26), rep("PB",10), rep("CD34",10))

ex <- exprs(gset)
# ex <- log2(ex+1) if the values were large it means it needs log
# exprs(gset) <- ex

# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
ex <- log2(ex) }

### quality control
# the first step is to obtain a box plot to check if values are normalized or not

pdf("Results/boxplot.pdf", width = 64)
boxplot(ex)
dev.off()

# ex <-normalizeQuantiles(ex) 
# exprs(gset) <- ex
#### Correlation heatmap
pdf("Results/CorHeatmap.pdf", width=15, height=15)
# pheatmap(cor(ex), labels_row = gr, labels_col = gr) 
pheatmap(cor(ex), color=bluered(256), border_color = NA) 
dev.off()

# excor <- cor(ex)
# excor[1:5,1:5]

#PCA
pc <- prcomp(ex)
pdf("Results/PC.pdf")
plot(pc)
plot(pc$x[,1:2])
dev.off()
# names(pc)
# dim(pc$x)
# colnames(pc$x)


ex.scale <- t(scale(t(ex), scale = FALSE)) #we don't want to eliminate variance
# mean(ex.scale[,1])
pc <- prcomp(ex.scale)
pdf("Results/PC_scaled.pdf")
plot(pc)
plot(pc$x[,1:2])
dev.off()

pcr <- data.frame(pc$r[,1:3])
# pcr <- data.frame(pc$r[,1:3], Group=gr)
# head(pcr)
pdf("Results/PCA_samples.pdf")
ggplot(pcr, aes(PC1, PC2)) + geom_point(size=3) + theme_bw()
# ggplot(pcr, aes(PC1, PC2, color=Group)) + geom_point(size=3)
dev.off()


#### Differential expression analysis
gr <- factor(gr) 
gset$description <- gr
design <- model.matrix(~ description + 0 , gset)
colnames(design) <- levels(gr)
fit <- lmFit(gset, design)
cont.matrix  <- makeContrasts(AML-CD34, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by = "B", number = Inf)

tT <- subset(tT, select = c("Gene.symbol",  "Gene.ID", "adj.P.Val", "logFC"))
write.table(tT, "Results/AML_CD34.txt", row.names = F, sep = "\t", quote=Fq)


fit <- lmFit(gset, design)
cont_matrix <- makeContrasts(aml - healthy, levels=design)
fit2 <- contrasts.fit(fit, cont_matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
tT <- subset(tT, select=c("Gene.symbol", "Gene.ID", "adj.P.Val", "logFC"))
write.table(tT, "Results/AML_vs_Healthy.txt", row.names=FALSE, sep='\t', quote=FALSE)


# for phase 2
#aml.up <- subset(tT, logFC > 1 & adj.P.Val < 0.05)
#aml.up.genes <- unique(as.character(strsplit2(aml.up$Gene.symbol, "///")))
#write.table(aml.up.genes, file="Results/AML_vs_Healthy_UP.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
#aml.down <- subset(tT, logFC < -1 & adj.P.Val < 0.05)
#aml.down.genes <- unique(as.character(strsplit2(aml.down$Gene.symbol, "///")))
#write.table(aml.down.genes, file="Results/AML_vs_Health_Down.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)


