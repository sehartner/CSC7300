           # Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Thu Jul 27 22:15:57 EDT 2017

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

basedir <- "/home/rstudio/Dropbox/CSC7300/Project/Final"

setwd(basedir)
resfolder <- "GEO2R.GSE1751.results2"

 if (! file.exists(resfolder) ) {
   dir.create(resfolder, showWarnings = FALSE, recursive = FALSE, mode = "0777")
   Sys.chmod(resfolder, mode = "0777", use_umask = TRUE)
 }

 if (! exists("gset") ){
   gset <- getGEO("GSE1751", GSEMatrix =TRUE, AnnotGPL=TRUE)
 }

 if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
 gset <- gset[[idx]]

save(gset, file = "gset.RData")

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
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=322)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

# Save that table
filename <- "GSE1751.data.tsv"
outfile <- paste(basedir, resfolder, filename, sep="/")
write.table(tT, file = outfile, quote = FALSE, dec = ",", sep = "\t", col.names = NA, row.names = T)

geneID = tT$Gene.symbol
write.csv(geneID, file = "geneID.csv", row.names = FALSE, quote = FALSE)

################################################################
# Volcano Plot

tTall <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(fit2))
tTall <- subset(tTall, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
#write.table(tTall, file=stdout(), row.names=F, sep="\t")

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


################################################################
# Heat Map

library(ComplexHeatmap)
dataHeatmap = exprs(gset)  #22.283 probesets versus HD and N samples
mostDEgenes <- tT$ID # use 322 most differentially expressed genes
dataHeatmap = dataHeatmap[mostDEgenes, ] #filter most differentially expressed genes
colnames(dataHeatmap) = pData(gset)$title
Heatmap(dataHeatmap, show_row_names = F)

################################################################
# Principal Component Analysis

pca.data <- ex
pca.data <- pca.data[mostDEgenes, ] #filter most DE genes
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
text(pca$x[,1], pca$x[,2], rownames(pca.data), cex=0.5)
#legend(10, -5, legend=c("HD", "Control"),
       #col=c("red", "blue"), pch = c(10, 10), bty="y")

# components #1 and #3
plot(pca$x[,1], pca$x[,3], xlab="PCA1", ylab="PCA3", main="PCA for components 1&3", type="p", pch=10, col=grpcol)
text(pca$x[,1], pca$x[,3], rownames(pca.data), cex=0.5)
#legend(10, -2, legend=c("HD", "Control"),
       #col=c("red", "blue"), pch = c(10, 10), bty="y")

# components #2 and #3
plot(pca$x[,2], pca$x[,3], xlab="PCA2", ylab="PCA3", main="PCA for components 2&3", type="p", pch=10, col=grpcol)
text(pca$x[,2], pca$x[,3], rownames(pca.data), cex=0.5)
#legend(10, -2, legend=c("HD", "Control"),
       #col=c("red", "blue"), pch = c(10, 10), bty="y")

barplot(summary(pca)$importance["Proportion of Variance", ], ylim = c(0,1))


################################################################
# Kmeans Clustering

d <- dist( t(dataHeatmap) )

#install.packages("rafalib")
library(rafalib)
mypar()
hc <- hclust(d)
samples <- colnames(gset)
plot(hc, labels = sml, cex = 0.5)

myplclust(hc, labels = sml, lab.col = as.fumeric(sml), cex = 0.5)
hclusters <-cutree(hc, h=25)
table(true = sml, cluster = hclusters)

hclusters <- cutree(hc, k=2)
table(true = sml, cluster = hclusters)

set.seed(1)
km <- kmeans(t(dataHeatmap[1:2,]), centers = 2)
names(km)
mypar(1,2)
plot(dataHeatmap[1,], dataHeatmap[2,], col = as.fumeric(sml), pch = 16)
plot(dataHeatmap[1,], dataHeatmap[2,], col = km$cluster, pch = 16)

km <- kmeans(t(dataHeatmap), centers = 2)
mds <- cmdscale(d)
mypar(1,2)
plot(mds[,1], mds[,2])
plot(mds[,1], mds[,2], col = km$cluster, pch = 16)

################################################################
# Top Genes

topGenes = subset(tT, (adj.P.Val<.05) & (abs(logFC)>2))

