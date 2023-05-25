## DIFFERENTIAL CUT&RUN ANALYIS 
## Identify number of differential CNR peaks across timecourse and cluster temporally 
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
  
  if (length(sampleNames)==6){
    colData <- data.frame(time, biorep)
    rownames(colData) <- sampleNames
    assign("colData", colData, envir = globalenv())
  } else if (length(sampleNames)==4){
    colData <- data.frame(time, biorep)
    colData <- colData[c(1,2,5,6),]
    rownames(colData) <- sampleNames
    assign("colData", colData, envir = globalenv())
  }
  
  dds <- DESeqDataSetFromMatrix(countMatrix, colData, ~biorep + time)
  dds <- DESeq(dds, test = "LRT", full = ~biorep + time, reduced = ~biorep)
  assign("dds", dds, envir = globalenv())
  
  if (length(sampleNames)==6){
    samples <- unique(time)[2:length(unique(time))]
    assign("samples", samples, envir = globalenv())
  } else if (length(sampleNames)==4){
    samples <- unique(time)[3:length(unique(time))]
    assign("samples", samples, envir = globalenv())
  }
  

}

# Convert samples to Coef objects for use with calling differential peaks
coefConvert <- function(samples){
  value <- paste0("time_", samples, "_vs_0")
  assign("coefs", value, envir = globalenv())
}

# Convert samples to Results objects for use with calling differential peaks
resultConvert <- function(samples){
  value <- paste0("res", samples)
  assign("results", value, envir = globalenv())
}

# Convert samples to Sig objects for use with calling differential peaks 
sigConvert <- function(samples){
  value <- paste0("res", samples, "sig")
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
  subMatrix <- countMatrix[rownames(countMatrix) %in% rownames(allsig), ]
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
  rownames(mergedMatrix) <- rownames(allsig)
  #Normalize your merged Matrix 
  normMatrix <- (mergedMatrix-rowMeans(mergedMatrix))/rowSds(mergedMatrix+0.5)
  assign("normMatrix", normMatrix, envir = globalenv())
}

## Run  ---------------------------------------------------
# Set thresholds
p <- 0.05
fc <- 2
time <- factor(c(0, 0, 360, 360, 4320, 4320))
biorep <- factor(c(1, 2, 1, 2, 1, 2))
matrixPaths <- list.files("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/raw/cnr/", full.names = TRUE)
assays <- unlist(strsplit(matrixPaths, "MEGA_"))[seq(2, length(matrixPaths)*2, by=2)]
assays <- unlist(strsplit(assays, "_K562_WT_PMA_S_counts.tsv"))
  
# Loop through each of the CNR datasets: 
for (x in 1:length(matrixPaths)){
  path <- matrixPaths[x]
  # Create dds and set up experiment from path of count matrix
  setupDeSeq(path)
  # Compare each sample to 0
  for (i in 1:length(samples)){
    result <- resultConvert(samples[i])
    coef <- coefConvert(samples[i])
    sig <- sigConvert(samples[i])
    resultVal <- simpleShrink(coef)
    sigVal <- resultVal[!is.na(resultVal$padj) & abs(resultVal$log2FoldChange)>fc & resultVal$padj<p, ]
    
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
    matrix <- matrix[!is.na(matrix$padj) & abs(matrix$log2FoldChange)>fc & matrix$padj<p, ]
    matrix$i <- i
    list[[i]] <- matrix
  }
  allsig <- unique((do.call(rbind,list))[,1:3])
  
  ## Create LFC matrix 
  if (x %in% c(1,2,3)){
    lfc <- cbind(res360$log2FoldChange, res4320$log2FoldChange)
  } else if (x ==4){
    lfc <- data.frame(res4320$log2FoldChange)
  }
  
  lfc <- data.frame(lfc)
  colnames(lfc) <- samples
  lfc$mean <- apply(lfc, 1, mean)
  lfc$median <- apply(lfc, 1, median)
  lfc <- cbind(peakMatrix[,1:3], lfc)
  lfc$differential <- ""
  lfc[rownames(lfc) %in% rownames(allsig),]$differential <- TRUE
  lfc[lfc$differential == "",]$differential <- FALSE
  
  ## Merged count matrix (for downstream filtering)
  counts <- countMatrix
  newmatrix <- list()
  if (x %in% c(1,2,3)){
    for (i in c(1,3,5)){
      tmp <- matrix(data=c(counts[,i], counts[,i+1]), ncol=2, byrow=F)
      newcol <- apply(tmp, 1, mean)
      newmatrix[[i]] <- newcol
    }
  } else if (x == 4){
    for (i in c(1,3)){
      tmp <- matrix(data=c(counts[,i], counts[,i+1]), ncol=2, byrow=F)
      newcol <- apply(tmp, 1, mean)
      newmatrix[[i]] <- newcol
    }
  }
  
  counts <- do.call(cbind, newmatrix)
  counts <- data.frame(counts)
  rownames(counts) <- rownames(peakMatrix)
  colnames(counts) <- c(0, samples)
  counts$mean <- apply(counts, 1, mean)
  counts$median <- apply(counts, 1, median)
  counts <- cbind(peakMatrix[,1:3], counts)
  
  # Write these files as txt files for downstream use 
  saveRDS(lfc, paste0("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/cnr/", assays[x], "LFC.rds"))
  saveRDS(counts, paste0("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/cnr/", assays[x], "Counts.rds"))

  ## Clustering
  # Normalize count matrix
  normMatrix <- normalizeMatrix(countMatrix)
  # Cluster based on which timepoint normalized counts reaches max or min: 
  clusters <- data.frame(normMatrix)
  clusters$whichmax <- apply(abs(clusters), 1, which.max)
  clusters$cluster <- ""
  
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
  
  pdf(paste0("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig3/", assays[x], "_heatmap.pdf"), width = 8.5, height = 11)
  par(mfrow=c(1,1))
  par(mar=c(3, 6, 3, 10))
  hotmap(matrixOrdered, labrow = FALSE, rowcolors = rowColors, colors = colorRampPalette(c("#0087B4", "black", "#C6A700"))(100))
  addlegend(range=c(-maxval, maxval), palette = colorRampPalette(c("#0087B4", "black", "#C6A700")), bottominset = 0.2, xoffset=0.08)
  dev.off()
  
  ## Save all data: 
  if (x %in% c(1,2,3)){
    allData <- 
      cbind(res360[,c(2,5)], 
            res4320[,c(2,5)]) |> 
      as.data.frame() |> 
      setNames(c("log2FoldChange_360", "padj_360", "log2FoldChange_4320", "padj_4320"))
    allData <- 
      cbind(peakMatrix[,1:3], allData)
    allData <- 
      cbind(allData, counts(dds))
    allData$cluster <- ""
    allData[rownames(allData) %in% rownames(clusters[clusters$cluster==1,]),]$cluster <- "upEarly"
    allData[rownames(allData) %in% rownames(clusters[clusters$cluster==2,]),]$cluster <- "upLate"
    allData[rownames(allData) %in% rownames(clusters[clusters$cluster==3,]),]$cluster <- "downEarly"
    allData[rownames(allData) %in% rownames(clusters[clusters$cluster==4,]),]$cluster <- "downLate"
    allData[allData$cluster=="",]$cluster <- "static"
    #Save: 
    write.table(allData, 
                file = paste0("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/cnr/", assays[x], "Summary.txt"), 
                quote = F, sep = "\t")
  } else if (x==4){
    allData <- res4320[,c(2,5)] |> 
      as.data.frame() |> 
      setNames(c("log2FoldChange_4320", "padj_4320"))
    allData <- 
      cbind(peakMatrix[,1:3], allData)
    allData <- 
      cbind(allData, counts(dds))
    allData$cluster <- ""
    allData[rownames(allData) %in% rownames(res4320sig[res4320sig$log2FoldChange>0,]),]$cluster <- "up"
    allData[rownames(allData) %in% rownames(res4320sig[res4320sig$log2FoldChange<0,]),]$cluster <- "down"
    allData[allData$cluster=="",]$cluster <- "static"
    #Save: 
    write.table(allData, 
                file = paste0("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/cnr/", assays[x], "Summary.txt"), 
                quote = F, sep = "\t")
  }

  
}

