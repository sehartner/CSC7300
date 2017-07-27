# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Tue Jul 25 18:25:08 EDT 2017

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO
basedir <- "/home/rstudio/Dropbox/CSC7300/Project"

setwd(basedir)
resfolder <- "GEO2R.GSE1751.results2"

# if (! file.exists(resfolder) ) {
#   dir.create(resfolder, showWarnings = FALSE, recursive = FALSE, mode = "0777")
#   Sys.chmod(resfolder, mode = "0777", use_umask = TRUE)
# }
# 
# if (! exists("gset") ){
#   gset <- getGEO("GSE1751", GSEMatrix =TRUE, AnnotGPL=TRUE)
# }

# if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
# gset <- gset[[idx]]
# 
# save(gset, file = "gset.RData")
load(file = "gset.RData")

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- "000000000000XXXXX11111111111111"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
tTall <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(fit2))

filename <- "GSE1751.data.tsv"
outfile <- paste(basedir, resfolder, filename, sep="/")
write.table(tT, file = outfile, quote = FALSE, dec = ",", sep = "\t", col.names = NA, row.names = T)

filenameAll <- "GSE1751.dataAll.tsv"
outfileAll <- paste(basedir, resfolder, filenameAll, sep="/")
write.table(tTall, file = outfileAll, quote = FALSE, dec = ",", sep = "\t", col.names = NA, row.names = T)


tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

tTall <- subset(tTall, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
#write.table(tTall, file=stdout(), row.names=F, sep="\t")


# ################################################################
# #   Boxplot for selected GEO samples
# library(Biobase)
# library(GEOquery)
# 
# # load series and platform data from GEO
# 
# gset <- getGEO("GSE1751", GSEMatrix =TRUE, getGPL=FALSE)
# if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
# gset <- gset[[idx]]
# 
# # group names for all samples in a series
# gsms <- "000000000000XXXXX11111111111111"
# sml <- c()
# for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
# sml <- paste("G", sml, sep="")  #set group names
# 
# # eliminate samples marked as "X"
# sel <- which(sml != "X")
# sml <- sml[sel]
# gset <- gset[ ,sel]
# 
# # order samples by group
# ex <- exprs(gset)[ , order(sml)]
# sml <- sml[order(sml)]
# fl <- as.factor(sml)
# labels <- c("HD","Control")
# 
# # set parameters and draw the plot
# palette(c("#dfeaf4","#f2cb98", "#AABBCC"))
# dev.new(width=4+dim(gset)[[2]]/5, height=6)
# par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
# title <- paste ("GSE1751", '/', annotation(gset), " selected samples", sep ='')
# boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
# legend("topleft", labels, fill=palette(), bty="n")

############################################
# Simple Volcano Plot

plot(tTall$logFC, 1-tTall$adj.P.Val, xlim = c(-6,6),
     main="Effect of Huntington's Disease on Transcriptional Regulation",
     xlab = "log2Ratio", ylab = "1-adj.P.Val")
abline(h=0.95, col = "red")

############################################
# Principal Component Analysis

pca.data <- ex
grps <- c(rep("HD",12), rep("Control",14)) # rename groups for clarity
grpcol <- c(rep("red",12), rep("blue", 14)) # vector of colors for groups
colnames(pca.data) <- paste(grps, colnames(pca.data), sep="-")

pca.data <- na.omit(as.matrix(pca.data)) # remove NAs
pca.data <- t(pca.data) # transpose data to get rows = samples and columns = probes

pca.data[,1:5] # preview first 5 probe sets

pca <- prcomp(pca.data, scale = TRUE) # compute principal component analysis
summary(pca)

# components 1 and 2
plot(pca$x[,1], pca$x[,2], xlab="PCA1", ylab="PCA2", main="PCA for components 1&2", type="p", pch=10, col=grpcol)
#text(pca$x[,1], pca$x[,2], rownames(pca.data), cex=0.75)
legend(75, 90, legend=c("HD", "Control"),
       col=c("red", "blue"), pch = c(10, 10), bty="y")

# components #1 and #3
plot(pca$x[,1], pca$x[,3], xlab="PCA1", ylab="PCA3", main="PCA for components 1&3", type="p", pch=10, col=grpcol)
#text(pca$x[,1], pca$x[,3], rownames(pca.data), cex=0.75)
legend(75, 125, legend=c("HD", "Control"),
       col=c("red", "blue"), pch = c(10, 10), bty="y")

# components #2 and #3
plot(pca$x[,2], pca$x[,3], xlab="PCA2", ylab="PCA3", main="PCA for components 2&3", type="p", pch=10, col=grpcol)
#text(pca$x[,2], pca$x[,3], rownames(pca.data), cex=0.75)
legend(75, 125, legend=c("HD", "Control"),
       col=c("red", "blue"), pch = c(10, 10), bty="y")

barplot(summary(pca)$importance["Proportion of Variance", ], ylim = c(0,1))

############################################
# Simple Heatmap

smpls <- as.character(pData(gset)$title)
colnames(gset) = pData(gset)$title
cells <- factor(sub("2_(\\w+)_[ab]", "\\1", smpls))
replicate <- factor(sub("2_\\w+_([ab])", "\\1", smpls))
pData(gset) <- cbind(pData(gset), cells = cells, replicate = replicate)
rownames(gset) = make.names(rownames(gset), unique = TRUE)
heatmap(cor(exprs(gset)), col = rev(brewer.pal(11, "RdBu")), labRow = smpls, labCol = smpls, scale = "none",
        margins = c(8,8), ColSideColors = c("darkgreen", "orange", "darkviolet")[cells])

############################################
# Clustering/Basic Machine Learning

e <- ex
d <- dist( t(e) )

#install.packages("rafalib")
library(rafalib)
mypar()
hc <- hclust(d)
samples <- colnames(gset)
plot(hc, labels = sml, cex = 0.5)

myplclust(hc, labels = sml, lab.col = as.fumeric(sml), cex = 0.5)
hclusters <-cutree(hc, h=80000)
table(true = sml, cluster = hclusters)

hclusters <- cutree(hc, k=3)
table(true = sml, cluster = hclusters)

set.seed(1)
km <- kmeans(t(e[1:2,]), centers = 3)
names(km)
mypar(1,2)
plot(e[1,], e[2,], col = as.fumeric(sml), pch = 16)
plot(e[1,], e[2,], col = km$cluster, pch = 16)

km <- kmeans(t(e), centers = 3)
mds <- cmdscale(d)
mypar(1,2)
plot(mds[,1], mds[,2])
plot(mds[,1], mds[,2], col = km$cluster, pch = 16)


###########################################
# Heatmap

library(RColorBrewer)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
library(genefilter)
rv <- rowVars(e)
idx <- order(-rv[1:40])
library(gplots)
cols <- palette(brewer.pal(8, "Dark2"))[as.fumeric(sml)]
head(cbind(colnames(e),cols))
heatmap.2(e[idx,], labCol = sml, trace = "none", ColSideColors = cols, col = hmcol)

#############################################
# Complex volcano plot

library(ggplot2)

# Make a basic volcano plot
with(tTall, plot(logFC, -log10(P.Value), pch=20, main="Volcano plot", xlim=c(-4,4), ylim=c(-1,15), cex = 0.73))

# Add colored points: red if padj<0.05, orange of log2FC>2, green if both)
with(subset(tTall, P.Value<.05 ), points(logFC, -log10(P.Value), pch=20, col="gray", cex = 0.73))
# with(subset(tTall, abs(logFC)>2), points(logFC, -log10(P.Value), pch=20, col="orange", cex = 0.75))
with(subset(tTall, P.Value<.05 & (logFC)>2), points(logFC, -log10(P.Value), pch=20, col="red", cex = 0.75))
with(subset(tTall, P.Value<.05 & (logFC)<(-2)), points(logFC, -log10(P.Value), pch=20, col="green", cex = 0.75))
abline(v=2, col="black")
abline(v=-2, col="black")
abline(h=-log10(.05), col="black")


gray = length(which(tTall$P.Value<.05 & tTall$logFC<=2 & tTall$logFC>=(-2)))
red = length(which(tTall$P.Value<.05 & tTall$logFC > 2))
green = length(which(tTall$P.Value<.05 & tTall$logFC <(-2)))
black = dim(tTall)[1] - (gray+red+green)

# legend(1, 95, legend=c("Line 1", "Line 2"), col=c("red", "blue"), lty=1:2, cex=0.8)
legend("topleft", legend="N = 171", cex = 0.8)
legend("top", legend="N = 2384", cex = 0.8)
legend("topright", legend="N = 353", cex = 0.8)
legend("bottom", legend="N = 42193", cex = 0.8)

legend(-1,13, legend=c("Sig. & Negative FC", "Significant", "Sig. & Postive FC", "Not Significant"),
       col=c("green", "gray", "red", "black"), pch = c(20, 20, 20, 20), bty="n", cex=0.6)


