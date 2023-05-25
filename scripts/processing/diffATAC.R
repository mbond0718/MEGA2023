## DIFFERENTIAL ATAC-SEQ ANALYIS 
## Identify number of differential ATAC peaks across timecourse and cluster temporally 
## This script (1) performs differential analysis and (2) generates heatmaps for normalized counts and lineplots showing clusters 

## Load libraries -----------------------------------------
library(DESeq2)
library(apeglm)
library(RColorBrewer)
library(ggplot2)
library(Sushi)
library(Sushi2)

## Define functions  --------------------------------------
# Set up DESeq dds
setupDeSeq <- function(matrixPath){
  
  peakMatrix <- read.csv(matrixPath, header = T, sep = "\t")
  rownames(peakMatrix) <- paste0("peak", 1:nrow(peakMatrix))
  assign("peakMatrix", peakMatrix, envir = globalenv())
  
  countMatrix <- peakMatrix[,4:ncol(peakMatrix)]
  assign("countMatrix", countMatrix, envir = globalenv())
  
  sampleNames <- colnames(countMatrix)
  
  colData <- data.frame(time, biorep)
  rownames(colData) <- sampleNames
  assign("colData", colData, envir = globalenv())
  
  dds <- DESeqDataSetFromMatrix(countMatrix, colData, ~biorep + time)
  dds <- DESeq(dds, test = "LRT", full = ~biorep + time, reduced = ~biorep)
  assign("dds", dds, envir = globalenv())
  
  samples <- unique(time)[2:length(unique(time))]
  assign("samples", samples, envir = globalenv())
}

# Convert samples to Coef objects for use with calling differential peaks
coefConvert <- function(samples){
  value <- paste0("time_", samples, "_vs_0")
  assign("coefs", value, envir = globalenv())
}

# Convert samples to Results objects for use with calling differential peaks
resultConvert <- function(samples){
  value <- paste0("atacRes", samples)
  assign("results", value, envir = globalenv())
}

# Convert samples to Sig objects for use with calling differential peaks 
sigConvert <- function(samples){
  value <- paste0("atacRes", samples, "sig")
  assign("sig", value, envir = globalenv())
}

# LFC Shrink
simpleShrink <- function(coefs){
  lfcShrink(dds, coef = coefs, type = "apeglm")
}

# Summarize findings of differential peak calls
sigSummary <- function(samples){
  inputName <- sigConvert(samples)
  input <- get(inputName)
  paste0("There are ", nrow(input), " differential peaks after ", samples, " min.")
}

# Generate and plot PCA
makePCA <- function(dds){
  vsd <- vst(dds, blind = FALSE)
  pca <- plotPCA(vsd, intgroup = c("time", "biorep"), returnData = TRUE)
  percentVar <- round(100 * attr(pca, "percentVar"))
  ggplot(pca, aes(PC1, PC2, color = time, shape=biorep)) + 
    geom_point(size = 3) + 
    xlab(paste0("PC1: ", percentVar[1], "%variance")) + 
    ylab(paste0("PC2: ", percentVar[2], "%variance")) + 
    labs(title = "Peak Calls") + 
    coord_fixed()
}


#Normalize Matrix function that combines replicates and normalizes the matrix 
normalizeMatrix <- function(countMatrix){
  #Subset for only the loops that are in allsig
  subMatrix <- countMatrix[rownames(countMatrix) %in% rownames(atacAllsig), ]
  #Get the treatment conditions to determine how to merge samples 
  treatment <- levels(samples)
  x <- ncol(subMatrix)/length(treatment)
  #Merge by the number determined above 
  mergedMatrix <- c()
  for (i in seq(1,ncol(subMatrix), by = x)){
    tempMatrix <- matrix(apply((subMatrix[,i:(i+1)]), 1, mean))
    mergedMatrix <- cbind(mergedMatrix,tempMatrix)
  }
  #Set correct row and col names
  colnames(mergedMatrix) <- treatment
  rownames(mergedMatrix) <- rownames(atacAllsig)
  #Normalize your merged Matrix 
  normMatrix <- (mergedMatrix-rowMeans(mergedMatrix))/rowSds(mergedMatrix+0.5)
  assign("normMatrix", normMatrix, envir = globalenv())
}

# Cluster
cluster_data <- function(countmatrixcombo, num_centers){
  clusters <- kmeans(countmatrixcombo, centers = num_centers)
  clusters <<- clusters$cluster
}

## Run  ---------------------------------------------------
# Set thresholds
p <- 0.05
lfc <- 2
matrixPath <- "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/raw/atac/MEGA_ATAC_K562_WT_PMA_S_MANUAL_MACS2_counts.tsv"
time <- factor(c(0, 0, 360, 360, 4320, 4320))
biorep <- factor(c(1, 2, 1, 2, 1, 2))

# Create dds and set up experiment from path of count matrix
setupDeSeq(matrixPath)

# Compare each sample to 0
for (i in 1:length(samples)){
  result <- resultConvert(samples[i])
  coef <- coefConvert(samples[i])
  sig <- sigConvert(samples[i])
  resultVal <- simpleShrink(coef)
  sigVal <- resultVal[!is.na(resultVal$padj) & abs(resultVal$log2FoldChange)>lfc & resultVal$padj<p, ]
  
  assign(result, resultVal, envir = globalenv())
  assign(sig, sigVal, envir = globalenv())
}

# Summarize number of differential peaks per sample
lapply(samples, sigSummary)

# Compile all unique significant differential peaks into allsig file
list <- c()
for (i in 1:length(samples)){
  resultName <- resultConvert(samples[i])
  result <- get(resultName)
  matrix <- cbind(peakMatrix, result)
  matrix <- matrix[!is.na(matrix$padj) & abs(matrix$log2FoldChange)>lfc & matrix$padj<p, ]
  matrix$i <- i
  list[[i]] <- matrix
}
atacAllsig <- unique((do.call(rbind,list))[,1:3])

## Create LFC matrix 
atacLFC <- cbind(atacRes360$log2FoldChange, atacRes4320$log2FoldChange)
atacLFC <- data.frame(atacLFC)
colnames(atacLFC) <- samples
atacLFC$mean <- apply(atacLFC, 1, mean)
atacLFC$median <- apply(atacLFC, 1, median)
atacLFC <- cbind(peakMatrix[,1:3], atacLFC)
atacLFC$differential <- ""
atacLFC[rownames(atacLFC) %in% rownames(atacAllsig),]$differential <- TRUE
atacLFC[atacLFC$differential == "",]$differential <- FALSE

## Merged count matrix (for downstream filtering)
atacCounts <- countMatrix
newmatrix <- list()
for (i in c(1,3,5)){
  tmp <- matrix(data=c(atacCounts[,i], atacCounts[,i+1]), ncol=2, byrow=F)
  newcol <- apply(tmp, 1, mean)
  newmatrix[[i]] <- newcol
}
atacCounts <- do.call(cbind, newmatrix)
atacCounts <- data.frame(atacCounts)
rownames(atacCounts) <- rownames(peakMatrix)
colnames(atacCounts) <- c(0, samples)
atacCounts$mean <- apply(atacCounts, 1, mean)
atacCounts$median <- apply(atacCounts, 1, median)
atacCounts <- cbind(peakMatrix[,1:3], atacCounts)

# Write these files as txt files for downstream use 
saveRDS(atacLFC, "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/atac/atacLFC.rds")
saveRDS(atacCounts, "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/atac/atacCounts.rds")


## Clustering
# Normalize count matrix
normMatrix <- normalizeMatrix(countMatrix)
# Cluster based on which timepoint normalized counts reaches max or min: 
clusters <- data.frame(normMatrix)
clusters$whichmax <- apply(abs(clusters), 1, which.max)
clusters$cluster <- ""
clusters[clusters$X0>0 & clusters$whichmax==1,]$cluster <- 4
clusters[clusters$X360>0 & clusters$whichmax==2,]$cluster <- 2
clusters[clusters$X4320>0 & clusters$whichmax==3,]$cluster <- 3
clusters[clusters$X0<0 & clusters$whichmax==1,]$cluster <- 1
clusters[clusters$X360<0 & clusters$whichmax==2,]$cluster <- 5
clusters[clusters$X4320<0 & clusters$whichmax==3,]$cluster <- 6


clusters[clusters$X0<0 & clusters$whichmax==1,]$cluster <- 1
clusters[clusters$X360>0 & clusters$whichmax==2,]$cluster <- 2


clusters[clusters$X0<0 & clusters$whichmax==1 | clusters$X360>0 & clusters$whichmax==2,]$cluster <- 1
clusters[clusters$X4320>0 & clusters$whichmax==3,]$cluster <- 2
clusters[clusters$X0>0 & clusters$whichmax==1,]$cluster <- 3
clusters[clusters$X360<0 & clusters$whichmax==2 | clusters$X4320<0 & clusters$whichmax==3,]$cluster <- 4
c <- clusters$cluster
names(c) <- rownames(clusters)
order <- c(2,1,4,3)
# Set color palette
pal <- RColorBrewer::brewer.pal(8, "Spectral")
labs <- pal[c(1,2,7,8)]


# Lineplots to show cluster dynamics 
par(mfrow=c(1,1))
par(mar=c(3, 6, 3, 10))

pdf("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig2/atacUpClusters.pdf", width = 10, height=8.5)
plot(colMeans(normMatrix), type="n", ylim=c(-1.2,1.2), xaxt="n", xlab="PMA Treatment (hours)", ylab="Relative Expression")
axis(1, at=1:3, labels=c(0, 6, 72))
lines(colMeans(normMatrix[rownames(normMatrix) %in% rownames(clusters[clusters$cluster==1,]),]), col = labs[2], lwd = 3)
lines(colMeans(normMatrix[rownames(normMatrix) %in% rownames(clusters[clusters$cluster==2,]),]), col = labs[1], lwd = 3)
dev.off()

pdf("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig2/atacDownClusters.pdf", width = 10, height=8.5)
plot(colMeans(normMatrix), type="n", ylim=c(-1.2,1.2), xaxt="n", xlab="PMA Treatment (hours)", ylab="Relative Expression")
axis(1, at=1:3, labels=c(0, 6, 72))
lines(colMeans(normMatrix[rownames(normMatrix) %in% rownames(clusters[clusters$cluster==3,]),]), col = labs[4], lwd = 3)
lines(colMeans(normMatrix[rownames(normMatrix) %in% rownames(clusters[clusters$cluster==4,]),]), col = labs[3], lwd = 3)
dev.off()

# Heatmap
matrixOrdered <- normMatrix[order(rowMax(normMatrix), decreasing=TRUE),]
clusterOrder <- c[rownames(matrixOrdered)]
#check, proceed only if true
identical(rownames(matrixOrdered), names(clusterOrder))

matrixOrdered <- matrixOrdered[order(match(clusterOrder, order)),]
maxval <- 1.1
matrixOrdered[which(matrixOrdered>maxval)] <- maxval
matrixOrdered[which(matrixOrdered<(-maxval))] <- (-maxval)
clusterOrder <- clusterOrder[rownames(matrixOrdered)]
colnames(matrixOrdered) <- c(0, 6, 72)

rowColors <- c(rep(labs[1], table(clusters$cluster)[2]), 
               rep(labs[2], table(clusters$cluster)[1]), 
               rep(labs[3], table(clusters$cluster)[4]), 
               rep(labs[4], table(clusters$cluster)[3]))


pdf("~/Desktop/MEGAManuscript/Figure2/ATAC_heatmap.pdf", width = 8.5, height = 11)
par(mfrow=c(1,1))
par(mar=c(3, 6, 3, 10))
hotmap(matrixOrdered, labrow = FALSE, rowcolors = rowColors, colors = colorRampPalette(c("#0087B4", "black", "#C6A700"))(100))
addlegend(range=c(-maxval, maxval), palette = colorRampPalette(c("#0087B4", "black", "#C6A700")), bottominset = 0.2, xoffset=0.08)
dev.off()

## Combine all data: 
allData <- 
  cbind(atacRes360[,c(2,5)], 
      atacRes4320[,c(2,5)]) |> 
  as.data.frame() |> 
  setNames(c("log2FoldChange_360", "padj_360", "log2FoldChange_4320", "padj_4320"))
allData <- 
  cbind(allData, counts(dds))
allData <- 
  cbind(peakMatrix[,1:3], allData)
allData$cluster <- ""
allData[rownames(allData) %in% rownames(clusters[clusters$cluster==1,]),]$cluster <- "upEarly"
allData[rownames(allData) %in% rownames(clusters[clusters$cluster==2,]),]$cluster <- "upLate"
allData[rownames(allData) %in% rownames(clusters[clusters$cluster==3,]),]$cluster <- "downEarly"
allData[rownames(allData) %in% rownames(clusters[clusters$cluster==4,]),]$cluster <- "downLate"
allData[allData$cluster=="",]$cluster <- "static"

#Save
write.table(allData, 
            file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/atac/atacSummary.txt", quote = F, sep = "\t")

