## SUBSAMPLING POWER ANALYSIS 
## Determine number of differential loops at multiple replicate combinations and sub-sampled sequencing depths 

## Load libraries -----------------------------------------
library(DESeq2)
library(GenomicRanges)
library(RColorBrewer)
library(RNASeqPower)
library(ggplot2)
library(mariner)
library(hictoolsr)
library(magrittr)

## Process ------------------------------------------------

# Read in data
hicFiles <-
  list.files("~/Desktop/MEGA/subsample/", full.names = TRUE)
loopFiles <- 
  list.files("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/raw/hic/loops/", full.names = TRUE)

# Extract counts and convert to GInteractions
subCounts <-
  mergePairs(x = loopFiles,
             binSize = 10e03,
             radius = 2,
             column = "APScoreAvg") %>%
  extractCounts(hic = hicFiles, 
                chroms = c(1:22, "X"),
                res = 10e3,
                norm = "NONE",
                matrix = "observed") %>%
  as.data.frame()

subCounts <- as.data.frame(subCounts)

# Write to file
saveRDS(subCounts, "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/subsampleLoopCounts.rds")


## Differential analysis: 
subCounts <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/subsampleLoopCounts.rds")
subCounts <- subCounts[,c(1:3, 6:8, 20:ncol(subCounts))] #just get colums for loop coordinates and counts for all samples 
rownames(subCounts) <- paste0("loop", 1:nrow(subCounts))
sampleNames <- colnames(subCounts[7:(ncol(subCounts))])
sampleNames <- unlist(strsplit(sampleNames, ".hic"))
colnames(subCounts)[7:(ncol(subCounts))] <- sampleNames
countMatrix <- subCounts[,7:ncol(subCounts)]

## Read in full data to append to the countMatrix: 
fullData <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/loopCounts.rds") |> 
  as.data.frame()
fullData <- fullData[,c(1:3,6:8,21:ncol(fullData))]
rownames(fullData) <- rownames(subCounts)
sampleNames <- colnames(fullData[7:ncol(fullData)])
sampleNames <- unlist(strsplit(sampleNames, "_inter.hic"))
sampleNames <- paste0(sampleNames, "_1_full")
colnames(fullData)[7:(ncol(fullData))] <- sampleNames
fullData <- fullData[,7:ncol(fullData)]

countMatrix <- cbind(countMatrix, fullData)
countMatrix <- countMatrix[ , order(names(countMatrix))]


##Make ColData
colData <- data.frame(factor(c(rep(0, 50), rep(360, 50), rep(4320, 50))), 
                      factor(rep(c(rep(1, 15), rep(2, 15), rep(3, 10), rep(4, 10)), 3)), 
                      factor(rep(c(rep(c(rep(1, 5), rep(2, 5), rep(3, 5)), 2), rep(c(rep(1, 5), rep(2, 5)), 2)), 3)),
                      factor(rep(c("100M", "300M", "500M", "700M", "900M"), 30)))
colnames(colData) <- c("time", "biorep", "techrep", "depth")
rownames(colData) <- colnames(countMatrix)

#Set parameters
reps <- c(2, 3, 4)
depths <- c("100M", "300M", "500M", "700M", "900M")
depthNums <- c(100000000, 300000000, 500000000, 1500000000)
dispersion <- 0.001984929
count <- 38
p <- 0.05
lfc <- log2(1.5)
pal <- brewer.pal(8, "YlGnBu")

dataReps <- list() 
dataDepth <- list()
for (x in 1:length(reps)){
  #Set the correct replicate number
  rep <- reps[x]
  for (y in 1:length(depths)){
    #Set correct depth
    depth <- depths[y]
    #Select correct columns for these samples from the colData
    if (rep==2){
      keep <- 1:rep
      cd <- colData[colData$biorep %in% keep & colData$depth==depth,]
      cm <- countMatrix[colnames(countMatrix) %in% rownames(cd)]
      #take the sum of technical replicates: 
      merged <- data.frame(apply(cm[,1:3], 1, sum), 
                           apply(cm[,4:6], 1, sum), 
                           apply(cm[,7:9], 1, sum), 
                           apply(cm[,10:12], 1, sum), 
                           apply(cm[,13:15], 1, sum), 
                           apply(cm[,16:18], 1, sum))
      colnames(merged) <- c("MEGA_K562_WT_PMA_0_1", "MEGA_K562_WT_PMA_0_2",
                            "MEGA_K562_WT_PMA_360_1", "MEGA_K562_WT_PMA_360_2", 
                            "MEGA_K562_WT_PMA_4320_1", "MEGA_K562_WT_PMA_4320_2")
      cd <- cd[c(1, 4, 7, 10, 13, 16),]
      cd <- cd[,1:2]
      rownames(cd) <- colnames(merged)
      #filter for low counts: 
      filter <- data.frame(apply(cm, 1, median))
      colnames(filter) <- "median"
      filter$x <- ""
      filter <- rownames(filter[filter$median>5,])
      #Remove loops that fall above the filter threshold from the count matrix 
      merged <- merged[rownames(merged) %in% filter,]
      merged <- as.matrix(merged)
      #Run DESeq: 
      #Run DESeq
      dds <- DESeqDataSetFromMatrix(merged, cd, ~biorep + time)
      dds <- DESeq(dds, test = "LRT", full = ~biorep + time, reduced = ~biorep)
      #Compare 6h and 72h to 0h to determine LFC 
      loopRes360 <- lfcShrink(dds, coef = "time_360_vs_0", type="apeglm")
      loopRes360 <- data.frame(loopRes360)
      loopRes4320 <- lfcShrink(dds, coef= "time_4320_vs_0", type="apeglm")
      loopRes4320 <- data.frame(loopRes4320)
      #Filter for significant differential loops 
      loopRes360sig <- loopRes360[!is.na(loopRes360$padj) & abs(loopRes360$log2FoldChange)>lfc & loopRes360$padj<p, ]
      loopRes4320sig <- loopRes4320[!is.na(loopRes4320$padj) & abs(loopRes4320$log2FoldChange)>lfc & loopRes4320$padj<p, ]
      allsig <- c(rownames(loopRes360sig), rownames(loopRes4320sig)) |> unique()
      #Calculate power: 
      #count <- fpm(dds) |> median()
      power <- rnapower(alpha=0.5/33914, cv=sqrt(dispersion), depth=count, effect=2, n=rep)      
      #Determine median loop size: 
      distance <- subCounts[rownames(subCounts) %in% allsig,]$start2 - subCounts[rownames(subCounts) %in% allsig,]$start1
      if (length(distance)==0){
        distance <- 0
      }
      distance <-
        #Add to dataframe: 
        temp <- data.frame(depth, rep, nrow(loopRes360sig), nrow(loopRes4320sig), 
                           length(allsig), 
                           power, distance)
      colnames(temp) <- c("depth", "reps", "diff6h", "diff72", "diffTotal", "power", "distance")
      dataDepth[[y]] <- temp
      print(paste(rep, depth))
    }
    if (rep==3){
      keep <- 1:rep
      cd <- colData[colData$biorep %in% keep & colData$depth==depth,]
      cm <- countMatrix[colnames(countMatrix) %in% rownames(cd)]
      #take the sum of technical replicates: 
      merged <- data.frame(apply(cm[,1:3], 1, sum), 
                           apply(cm[,4:6], 1, sum), 
                           apply(cm[,7:8], 1, sum), 
                           apply(cm[,9:11], 1, sum), 
                           apply(cm[,12:14], 1, sum), 
                           apply(cm[,15:16], 1, sum), 
                           apply(cm[,17:19], 1, sum), 
                           apply(cm[,20:22], 1, sum), 
                           apply(cm[,23:24], 1, sum))
      colnames(merged) <- c("MEGA_K562_WT_PMA_0_1", "MEGA_K562_WT_PMA_0_2", "MEGA_K562_WT_PMA_0_3", 
                            "MEGA_K562_WT_PMA_360_1", "MEGA_K562_WT_PMA_360_2", "MEGA_K562_WT_PMA_360_3",
                            "MEGA_K562_WT_PMA_4320_1", "MEGA_K562_WT_PMA_4320_2", "MEGA_K562_WT_PMA_4320_3")
      cd <- cd[c(1, 4, 7, 9, 12, 15, 17, 20, 23),]
      cd <- cd[,1:2]
      rownames(cd) <- colnames(merged)
      #filter for low counts: 
      filter <- data.frame(apply(cm, 1, median))
      colnames(filter) <- "median"
      filter$x <- ""
      filter <- rownames(filter[filter$median>5,])
      #Remove loops that fall above the filter threshold from the count matrix 
      merged <- merged[rownames(merged) %in% filter,]
      merged <- as.matrix(merged)
      #Run DESeq: 
      #Run DESeq
      dds <- DESeqDataSetFromMatrix(merged, cd, ~biorep + time)
      dds <- DESeq(dds, test = "LRT", full = ~biorep + time, reduced = ~biorep)
      #Compare 6h and 72h to 0h to determine LFC 
      loopRes360 <- lfcShrink(dds, coef = "time_360_vs_0", type="apeglm")
      loopRes360 <- data.frame(loopRes360)
      loopRes4320 <- lfcShrink(dds, coef= "time_4320_vs_0", type="apeglm")
      loopRes4320 <- data.frame(loopRes4320)
      #Filter for significant differential loops 
      loopRes360sig <- loopRes360[!is.na(loopRes360$padj) & abs(loopRes360$log2FoldChange)>lfc & loopRes360$padj<p, ]
      loopRes4320sig <- loopRes4320[!is.na(loopRes4320$padj) & abs(loopRes4320$log2FoldChange)>lfc & loopRes4320$padj<p, ]
      allsig <- c(rownames(loopRes360sig), rownames(loopRes4320sig)) |> unique()
      #Calculate power: 
      #count <- fpm(dds) |> median()
      power <- rnapower(alpha=0.5/33914, cv=sqrt(dispersion), depth=count, effect=2, n=rep)      
      #Determine median loop size: 
      distance <- subCounts[rownames(subCounts) %in% allsig,]$start2 - subCounts[rownames(subCounts) %in% allsig,]$start1
      if (length(distance)==0){
        distance <- 0
      }
      distance <- median(distance)
      #Add to dataframe: 
      temp <- data.frame(depth, rep, nrow(loopRes360sig), nrow(loopRes4320sig), 
                         length(allsig), 
                         power, distance)
      colnames(temp) <- c("depth", "reps", "diff6h", "diff72", "diffTotal", "power", "distance")
      dataDepth[[y]] <- temp
      print(paste(rep, depth))
    }
    if (rep==4){
      keep <- 1:rep
      cd <- colData[colData$biorep %in% keep & colData$depth==depth,]
      cm <- countMatrix[colnames(countMatrix) %in% rownames(cd)]
      #take the sum of technical replicates: 
      merged <- data.frame(apply(cm[,1:3], 1, sum), 
                           apply(cm[,4:6], 1, sum), 
                           apply(cm[,7:8], 1, sum), 
                           apply(cm[,9:10], 1, sum), 
                           apply(cm[,11:13], 1, sum), 
                           apply(cm[,14:16], 1, sum), 
                           apply(cm[,17:18], 1, sum), 
                           apply(cm[,19:20], 1, sum), 
                           apply(cm[,21:23], 1, sum), 
                           apply(cm[,24:26], 1, sum), 
                           apply(cm[,27:28], 1, sum), 
                           apply(cm[,29:30], 1, sum))
      colnames(merged) <- c("MEGA_K562_WT_PMA_0_1", "MEGA_K562_WT_PMA_0_2", "MEGA_K562_WT_PMA_0_3", "MEGA_K562_WT_PMA_0_4", 
                            "MEGA_K562_WT_PMA_360_1", "MEGA_K562_WT_PMA_360_2", "MEGA_K562_WT_PMA_360_3", "MEGA_K562_WT_PMA_360_4",
                            "MEGA_K562_WT_PMA_4320_1", "MEGA_K562_WT_PMA_4320_2", "MEGA_K562_WT_PMA_4320_3", "MEGA_K562_WT_PMA_4320_4")
      cd <- cd[c(1, 4, 7, 9, 11, 14, 17, 19, 21, 24, 27, 29),]
      cd <- cd[,1:2]
      rownames(cd) <- colnames(merged)
      #filter for low counts: 
      filter <- data.frame(apply(cm, 1, median))
      colnames(filter) <- "median"
      filter$x <- ""
      filter <- rownames(filter[filter$median>5,])
      #Remove loops that fall above the filter threshold from the count matrix 
      merged <- merged[rownames(merged) %in% filter,]
      merged <- as.matrix(merged)
      #Run DESeq: 
      #Run DESeq
      dds <- DESeqDataSetFromMatrix(merged, cd, ~biorep + time)
      dds <- DESeq(dds, test = "LRT", full = ~biorep + time, reduced = ~biorep)
      #Compare 6h and 72h to 0h to determine LFC 
      loopRes360 <- lfcShrink(dds, coef = "time_360_vs_0", type="apeglm")
      loopRes360 <- data.frame(loopRes360)
      loopRes4320 <- lfcShrink(dds, coef= "time_4320_vs_0", type="apeglm")
      loopRes4320 <- data.frame(loopRes4320)
      #Filter for significant differential loops 
      loopRes360sig <- loopRes360[!is.na(loopRes360$padj) & abs(loopRes360$log2FoldChange)>lfc & loopRes360$padj<p, ]
      loopRes4320sig <- loopRes4320[!is.na(loopRes4320$padj) & abs(loopRes4320$log2FoldChange)>lfc & loopRes4320$padj<p, ]
      allsig <- c(rownames(loopRes360sig), rownames(loopRes4320sig)) |> unique()
      #Calculate power: 
      #count <- fpm(dds) |> median()
      power <- rnapower(alpha=0.5/33914, cv=sqrt(dispersion), depth=count, effect=2, n=rep)
      #Determine median loop size: 
      distance <- subCounts[rownames(subCounts) %in% allsig,]$start2 - subCounts[rownames(subCounts) %in% allsig,]$start1
      if (length(distance)==0){
        distance <- 0
      }
      distance <- median(distance)
      #Add to dataframe: 
      temp <- data.frame(depth, rep, nrow(loopRes360sig), nrow(loopRes4320sig), 
                         length(allsig), 
                         power, distance)
      colnames(temp) <- c("depth", "reps", "diff6h", "diff72", "diffTotal", "power", "distance")
      dataDepth[[y]] <- temp
      print(paste(rep, depth))
    }
    temp <- do.call(rbind, dataDepth)
    dataReps[[x]] <- temp
  }
}
data <- do.call(rbind, dataReps)
data <- data[,1:5]
data <- unique(data)

# Save data: 
saveRDS(data, file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/subsamplePowerData.rds")






