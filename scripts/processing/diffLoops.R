## DIFFERENTIAL LOOPING:
## This script is to be run after makeLoopCountMatrix.R to identify differential loops 

## Load libraries -----------------------------------------
library(DESeq2)
library(GenomicRanges)
library(RColorBrewer)
library(RNASeqPower)
library(ggplot2)
library(VennDiagram)

## Process ------------------------------------------------

# Set parameters
p <- 0.05
lfc <- log2(1.5)

# Read in loop matrix from makeLoopCountMatrix.R
loopMatrix <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/loopCounts.rds") |> 
  as.data.frame()
loopMatrix <- loopMatrix[,c(1:3,6:8,21:50)]

# Abbreviate names and replace column names
sampleNames <- colnames(loopMatrix[7:(ncol(loopMatrix))])
sampleNames <- unlist(strsplit(sampleNames, "_inter.hic"))
colnames(loopMatrix)[7:(ncol(loopMatrix))] <- sampleNames
# Set loop ID names for each loop to later access the loop coordinates 
rownames(loopMatrix) <- paste0("loop", 1:nrow(loopMatrix))
# For use with DESeq, get just the counts for each loop 
countMatrix <- loopMatrix[,7:36]

# Merge technical replicates
mergedMatrix <- list()
for (i in 1:12){
  x=c(1,4,7,9,11,14,17,19,21,24,27,29)[i]
  by=c(2,2,1,1,2,2,1,1,2,2,1,1)[i]
  sub <- matrix(apply((countMatrix[,x:(x+by)]), 1, sum))
  mergedMatrix[[i]] <- sub
}
mergedMatrix <- data.frame(do.call(cbind, mergedMatrix))
mergedMatrix <- (mergedMatrix)
# Set columnnames and rownames 
colnames(mergedMatrix) <- c("0_1", "0_2", "0_3", "0_4", "360_1", "360_2", "360_3", "360_4", "4320_1", "4320_2", "4320_3", "4320_4")
rownames(mergedMatrix) <- paste0("loop", 1:nrow(mergedMatrix))

# Manually create colData for DESeq
time <- factor(c(0,0,0,0,360,360,360,360,4320,4320,4320,4320))
rep <- factor(c(1,2,3,4,1,2,3,4,1,2,3,4))
colData <- data.frame(time, rep)
rownames(colData) <- colnames(mergedMatrix)

# Filter for low counts before DDS generation
filter <- data.frame(apply(countMatrix, 1, median))
colnames(filter) <- "median"
filter$x <- ""
filter <- rownames(filter[filter$median>5,])
# Remove loops that fall above the filter threshold from the count matrix 
counts <- mergedMatrix[rownames(mergedMatrix) %in% filter,]

# Run DESeq
dds <- DESeqDataSetFromMatrix(counts, colData, ~rep + time)
dds <- DESeq(dds, test = "LRT", full = ~rep + time, reduced = ~rep)
#Get Dispersion: 
d <- min(dispersions(dds))
# Determine median CPM for power analysis: 
median <- fpm(dds) |> median()

# Compare 6h and 72h to 0h to determine LFC 
loopRes360 <- lfcShrink(dds, coef = "time_360_vs_0", type="apeglm")
loopRes360 <- data.frame(loopRes360)
loopRes4320 <- lfcShrink(dds, coef= "time_4320_vs_0", type="apeglm")
loopRes4320 <- data.frame(loopRes4320)

#Filter for significant differential loops 
loopRes360sig <- loopRes360[!is.na(loopRes360$padj) & abs(loopRes360$log2FoldChange)>lfc & loopRes360$padj<p, ]
loopRes4320sig <- loopRes4320[!is.na(loopRes4320$padj) & abs(loopRes4320$log2FoldChange)>lfc & loopRes4320$padj<p, ]
#Compile significant diff loops from both timepoints
loopAllsig <- unique(c(rownames(loopRes360sig), rownames(loopRes4320sig)))
#Get coordinates of differential loops from the loopMatrix file 
diffLoops <- (loopMatrix[rownames(loopMatrix) %in% loopAllsig,])[,1:6]
diffLoops$seqnames1 <- unlist(strsplit(as.character(diffLoops$seqnames1), "chr"))[seq(2,(2*nrow(diffLoops)),by=2)]
diffLoops$seqnames2 <- unlist(strsplit(as.character(diffLoops$seqnames2), "chr"))[seq(2,(2*nrow(diffLoops)),by=2)]

# Clustering
# Determine clusters
diffLoopsCounts <- (counts[rownames(counts) %in% loopAllsig,])
# Get the treatment conditions to determine how to merge samples 
treatment <- c(0, 360, 4320)
x <- ncol(diffLoopsCounts)/length(treatment)
# Merge by the number determined above 
mergedMatrix <- c()
for (i in seq(1,ncol(diffLoopsCounts), by = x)){
  tempMatrix <- matrix(apply((diffLoopsCounts[,i:(i+3)]), 1, mean))
  mergedMatrix <- cbind(mergedMatrix,tempMatrix)
}
# Normalize 
normMatrix <- data.frame((mergedMatrix-rowMeans(mergedMatrix))/rowSds(mergedMatrix+0.5))
# Set correct row and col names
colnames(normMatrix) <- treatment
rownames(normMatrix) <- rownames(diffLoopsCounts)

# Using slopes to categorize samples: 
normMatrix$cluster <- ""
tmp <- normMatrix[,1:3]+abs(min(normMatrix[,1:3]))
tmp$cluster <- ""

for (i in 1:nrow(tmp)){
  val <- as.numeric(tmp[i,1:3])
  slope <- as.numeric(c((val[2]-val[1]), val[3]-val[2]))
  
  #Up early
  if((sum(slope)>0) & (abs(slope[1])>abs(slope[2]))){
    tmp[i,4] <- 1
  }
  
  #Up Late
  if((sum(slope)>0) & (abs(slope[1])<abs(slope[2]))){
    tmp[i,4] <- 2
  }
  
  #Down early
  if((sum(slope)<0) & (abs(slope[1])>abs(slope[2]))){
    tmp[i,4] <- 3
  }
  
  #Down late
  if((sum(slope)<0) & (abs(slope[1])<abs(slope[2]))){
    tmp[i,4] <- 4
  }
  
}
matrix <- data.frame(mergedMatrix)
rownames(matrix) <- rownames(tmp)
colnames(matrix) <- treatment
matrix$cluster <- tmp$cluster
normMatrix$cluster <- tmp$cluster

normMatrix[normMatrix$cluster=="" & normMatrix$`0`<0,]$cluster <- 2
normMatrix[normMatrix$cluster=="" & normMatrix$`0`>0,]$cluster <- 4

# Create gained/lost/static objects
gainedLoops <- (loopMatrix[rownames(loopMatrix) %in% (rownames(normMatrix[normMatrix$cluster==1 | normMatrix$cluster==2,])),])[,1:6]
lostLoops <- (loopMatrix[rownames(loopMatrix) %in% (rownames(normMatrix[normMatrix$cluster==3 | normMatrix$cluster==4,])),])[,1:6]
staticLoops <- (loopMatrix[!(rownames(loopMatrix) %in% loopAllsig),])[,1:6]

# Write 4320 FC to table
fc4320 <- cbind((loopMatrix[rownames(loopMatrix) %in% rownames(counts),])[,1:6], loopRes4320$log2FoldChange, loopRes4320$padj)
colnames(fc4320)[7:8] <- c("lfc", "padj")
fc4320$class <- ""
fc4320[rownames(fc4320) %in% rownames(gainedLoops),]$class <- "gained"
fc4320[rownames(fc4320) %in% rownames(lostLoops),]$class <- "lost"
fc4320[rownames(fc4320) %in% rownames(staticLoops),]$class <- "static"

saveRDS(fc4320, 
        file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/loopLFC.rds")

# Write 360 FC to table
fc360 <- cbind((loopMatrix[rownames(loopMatrix) %in% rownames(counts),])[,1:6], loopRes360$log2FoldChange, loopRes360$padj)
colnames(fc360)[7:8] <- c("lfc", "padj")
fc360$class <- ""
fc360[rownames(fc360) %in% rownames(gainedLoops),]$class <- "gained"
fc360[rownames(fc360) %in% rownames(lostLoops),]$class <- "lost"
fc360[rownames(fc360) %in% rownames(staticLoops),]$class <- "static"

# saveRDS(fc360, 
#         file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/loopLFC6h.rds")


# Loop counts with diff loop status: 
loops <- loopMatrix[,1:6]
loops$class <- ""
loops[rownames(loops) %in% rownames(staticLoops),]$class <- "static"
loops[rownames(loops) %in% rownames(gainedLoops),]$class <- "gained"
loops[rownames(loops) %in% rownames(lostLoops),]$class <- "lost"

saveRDS(loops, 
        file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/diffLoops.rds")

##PCA
vsd <- vst(dds, blind = FALSE)
pca <- plotPCA(vsd, intgroup = c("time", "rep"), returnData = TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))
ggplot(pca, aes(PC1, PC2, color = time, shape = rep)) + 
  geom_point(size = 3) + 
  xlab(paste0("PC1: ", percentVar[1], "%variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "%variance")) + 
  labs(title = "Loop Calls") +
  xlim(c(-2.5,2.5)) + ylim(c(-2.5,2.5)) +
  coord_fixed() + theme_classic()
ggsave("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/sup1/pca.pdf", last_plot())

## Save all data into summary file: 
allData <- 
  cbind(loopRes360[,c(2,5)], 
        loopRes4320[,c(2,5)]) |> 
  as.data.frame() |> 
  setNames(c("log2FoldChange_360", "padj_360", "log2FoldChange_4320", "padj_4320"))

allData <- merge(allData, loopMatrix, all = TRUE, by = 0)
allData <- allData[match(rownames(loopMatrix), allData$Row.names),]
allData$cluster <- ""
allData[allData$Row.names %in% rownames(gainedLoops),]$cluster <- "gained"
allData[allData$Row.names %in% rownames(lostLoops),]$cluster <- "lost"
allData[allData$cluster=="",]$cluster <- "static"

#Save
write.table(allData, 
            file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/hicSummary.txt", quote = F, sep = "\t", row.names = FALSE)


## Make heatmap showing differential loops: (everything from here out is for reviewers)
tmp <- normMatrix
order <- c(2,1,4,3)
# Set color palette
pal <- RColorBrewer::brewer.pal(8, "Spectral")
labs <- pal[c(1,2,7,8)]

normMatrix <- as.matrix(normMatrix[,1:3])
matrixOrdered <- normMatrix[order(rowMax(normMatrix), decreasing=TRUE),]
clusters <- tmp$cluster
names(clusters) <- rownames(tmp)
clusterOrder <- clusters[rownames(matrixOrdered)]
#check, proceed only if true
identical(rownames(matrixOrdered), names(clusterOrder))

matrixOrdered <- matrixOrdered[order(match(clusterOrder, order)),]
maxval <- 1.1
matrixOrdered[which(matrixOrdered>maxval)] <- maxval
matrixOrdered[which(matrixOrdered<(-maxval))] <- (-maxval)
clusterOrder <- clusterOrder[rownames(matrixOrdered)]
colnames(matrixOrdered) <- c(0, 360, 4320)

rowColors <- c(rep(labs[1], table(clusters)[2]), 
               rep(labs[2], table(clusters)[1]), 
               rep(labs[3], table(clusters)[4]), 
               rep(labs[4], table(clusters)[3]))

#pdf("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig2/rnaHeatmap.pdf", width = 8.5, height = 11, compress = FALSE)
par(mfrow=c(1,1))
par(mar=c(3, 6, 3, 10))
hotmap(matrixOrdered, labrow = FALSE, rowcolors = rowColors)
addlegend(range=c(-maxval, maxval), palette = colorRampPalette(c("deepskyblue2", "black", "gold")), bottominset = 0.2, xoffset=0.08)
#dev.off()


## Number of differential loops at 6h and 72h: 
data <- 
  data.frame(sample = c("diff_6h", "diff_72h"), 
             values = c(nrow(loopRes360sig), nrow(loopRes4320sig)))
ggplot(data, aes(x = sample, y = values, fill = sample)) + 
  geom_bar(stat = "identity") + 
  theme_classic() + 
  scale_fill_manual(values = c("gray80", "gray80")) + 
  theme(legend.position = "none")
ggsave(last_plot(), file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/reviewFigs/numDiffLoops.pdf")


## Numbers of loops identified after each timepoint vs 0
venn.diagram(list("4320" = rownames(loopRes4320sig), "360" = rownames(loopRes360sig)), fill = c("lightblue", "green"), 
             alpha = c(0.5, 0.5), lwd =0, filename = "~/Desktop/vennTest.png")

venn.diagram(list("360" = rownames(loopRes360sig), "4320" = rownames(loopRes4320sig)), fill = c("lightblue", "green"), 
             alpha = c(0.5, 0.5), lwd =0, filename = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/reviewFigs/diffLoopVenn.png")

## Loops ONLY differential at 6h: 

transient <- loopRes360sig[!rownames(loopRes360sig) %in% rownames(loopRes4320sig),]
transient$log2FoldChange |> range()
write(rownames(transient), "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/transient.txt")

## COMPARE VARIOUS ATTRIBUTES OF 6h loops and 72h loops: 

diff6 <- rownames(loopRes360sig)
diff72 <- rownames(loopRes4320sig)

#Size density: 
diff6Coord <- (loopMatrix[rownames(loopMatrix) %in% diff6,])[,1:6]
diff72Coord <- (loopMatrix[rownames(loopMatrix) %in% diff72,])[,1:6]

#make dataframe
sizes <- data.frame(class = c(rep("diff6", length(diff6)), 
                              rep("diff72", length(diff72))), 
                    values = c(abs(diff6Coord$start2 - diff6Coord$start1),
                               abs(diff72Coord$start2 - diff72Coord$start1)))
#test for significance
wilcox.test(abs(diff6Coord$start2 - diff6Coord$start1),
            abs(diff72Coord$start2 - diff72Coord$start1))
#Plot density:
ggplot(data = sizes, aes(x = values, fill = class)) + 
  geom_density(alpha = 0.75) + 
  theme_classic()
ggsave(last_plot(), filename = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/reviewFigs/diffLoopSizes.pdf")

#Fold Change Density: 
diff6FC <- (fc360[rownames(fc360) %in% diff6,])
diff72FC <- (fc4320[rownames(fc4320) %in% diff72,])

#make dataframe
fc <- data.frame(class = c(rep("diff6", nrow(diff6FC)), 
                              rep("diff72", nrow(diff72FC))), 
                    values = c(diff6FC$lfc, diff72FC$lfc))
#test for significance
wilcox.test(diff6FC$lfc, diff72FC$lfc)
#Plot density:
ggplot(data = fc, aes(x = values, fill = class)) + 
  geom_density(alpha = 0.75) + 
  theme_classic()

ggplot(data = fc, aes(x = abs(values), fill = class)) + 
  geom_density(alpha = 0.75) + 
  theme_classic()
ggsave(last_plot(), filename = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/reviewFigs/diffLoopLFC.pdf")

#PADJ Density: 

#make dataframe
padj <- data.frame(class = c(rep("diff6", nrow(diff6FC)), 
                           rep("diff72", nrow(diff72FC))), 
                 values = c(diff6FC$padj, diff72FC$padj))
#test for significance
wilcox.test(diff6FC$padj, diff72FC$padj)
#Plot density:
ggplot(data = padj, aes(x = values, fill = class)) + 
  geom_density(alpha = 0.75) + 
  theme_classic()
ggsave(last_plot(), filename = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/reviewFigs/diffLoopPadj.pdf")


## CTCF vs Non CTCF: 
ctcfLoops <- read.table("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/ctcfLoops.txt", header = F)$V1
nonCTCFloops <- read.table("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/nonCTCFLoops.txt", header = F)$V1

data <- data.frame(timepoint = c("6h", "6h", "72h", "72h"), 
                   category = c("CTCF", "nonCTCF", "CTCF", "nonCTCF"), 
                   values = c(length(diff6[diff6 %in% ctcfLoops]), 
                              length(diff6[diff6 %in% nonCTCFloops]), 
                              length(diff72[diff72 %in% ctcfLoops]), 
                              length(diff72[diff72 %in% nonCTCFloops])))

ggplot(data, aes(x = timepoint, y = values, fill = category)) + 
  geom_bar(stat = "identity", position = "stack") + 
  theme_classic()

#Frequency version: 

data$frequency <- ""
data$frequency[1:2] <- data[1:2,3]/(105+13)
data$frequency[3:4] <- data[3:4,3]/(1343+104)
data$frequency <- as.numeric(data$frequency)

ggplot(data, aes(x = timepoint, y = frequency, fill = category)) + 
  geom_bar(stat = "identity", position = "stack") + 
  theme_classic()
ggsave(last_plot(), filename = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/reviewFigs/diffLoopCTCF.pdf")

fisher.test(matrix(c(105,13,1343,104),nrow=2,ncol=2))$p.value #just gained


