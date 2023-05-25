## DIFFERENTIAL RNA-SEQ ANALYIS 
## Identify number of differential genes across timecourse and cluter temporally 
## This script (1) performs differential analysis and (2) generates heatmaps for normalized counts and lineplots showing clusters 

## Load libraries -----------------------------------------
library(DESeq2)
library(apeglm)
library(RColorBrewer)
library(ggplot2)
library(Sushi)
library(Sushi2)
library(gridExtra)

## Define functions  --------------------------------------
# Set up the dds (DESeq Data Sheet) and extract sample names 
setup_DeSeq <- function(txipath, configpath){
  #read in RDS and config
  txi <- readRDS(txipath)
  config <- read.delim(configpath, header = T, sep = '\t')
  #create colData for DDS
  col_data <- data.frame(time = factor(config$Time), rep = factor(config$Bio_Rep), row.names = config$Name)
  #transform txi into dds txi 
  dds_txi <- DESeqDataSetFromTximport(txi, colData = col_data, design = ~rep + time)
  #transform dds txi into dds
  dds <- DESeq(dds_txi, test = "LRT", full = ~rep + time, reduced = ~rep)
  assign("dds", dds, envir = globalenv())
  #extract sample names
  samples <- unique(as.vector(config$Time))
  samples <- samples[2:8]
  assign("samples", samples, envir = globalenv())
}

# Convert samples to Coefs
coef_convert <- function(samples){
  value <- paste0("time_", samples, "_vs_0")
  assign("coefs", value, envir = globalenv())
}

# Convert samples to Results
result_convert <- function(samples){
  value <- paste0("rnaRes_", samples)
  assign("results", value, envir = globalenv())
}

# Convert samples to Sig Results
sig_result_convert <- function(samples){
  value <- paste0("rnaRes", samples, "sig")
  assign("sigs", value, envir = globalenv())
}

# Simplify shrink to be used within loop
simple_shrink <- function(coef){
  lfcShrink(dds, coef = coef, type = "apeglm")
}

#Print out number of differential genes
sig_summary <- function(samples){
  input_name <- sig_result_convert(samples)
  input <- get(input_name)
  print(paste0("There are ", nrow(input), " differential genes after ", samples, " min."))
}

# Make normalized count matrix for clustering and plotting
make_norm_matrix <- function(dds){
  #normalize dds
  dds_norm <- vst(dds)
  #pull out only the genes in allsig
  count_matrix <- dds_norm[rownames(dds_norm) %in% rnaAllsig,]
  #make the data into a matrix that we can easily read
  count_matrix <- assay(count_matrix)
  #normalize to 0 and declare to global environment
  count_matrix_norm <<- (count_matrix-rowMeans(count_matrix))/rowSds(count_matrix+.5)
}

# Function to combine bioreps
combine_reps <- function(matrix, newcolnames){
  #first, establish an emptry matrix to build into
  newMatrix <- c()
  #this should equal half of your samples i.e. the number of cols you want in the final matrix
  for (i in seq(1, ncol(matrix)/2)){
    #build temp matrix with cols you want to combine
    #make sure that you're changing i to grab only the columns of interest
    tempMatrix <- matrix(data = c(matrix[,i], matrix[,i+8]), ncol = 2, byrow = F)
    #create a new column that's averaging the two previous columns by row 
    newCol <- apply(tempMatrix, 1, mean)
    #add this matrix on to the newMatrix
    newMatrix <- cbind(newMatrix, newCol)
  }
  #add rownames of old matrix to new matrix (we're not changing the number of rows or the order of rows)
  rownames(newMatrix) <- rownames(matrix)
  
  #how your new colnames are added to the old colnames 
  if(!is.null(newcolnames)){
    colnames(newMatrix) <- newcolnames
  }
  return(newMatrix)
}

# Kmeans clustering
cluster_data <- function(countmatrixcombo, num_centers){
  clusters <- kmeans(countmatrixcombo, centers = num_centers)
  clusters <<- clusters$cluster
}

## Run  ---------------------------------------------------
# Define locations of RDS and config###
txi_path <- "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/raw/rna/MEGA_K562_WT_PMA_S_txi.rds"
config_path <- "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/raw/rna/config_MEGA_K562_WT_PMA_S_210712_133920.tsv"

# Create DDS
setup_DeSeq(txi_path, config_path)
names(dds) <- unlist(strsplit(names(dds), "[.]"))[seq(1,(2*length(names(dds))),by=2)]

# Compare each sample to 0 with DESeq's LFC shrink function
for (i in samples){
  name <- result_convert(i)
  value <- coef_convert(i)
  value <- simple_shrink(value)
  assign(name, value, envir = globalenv())
}

# Pull out only the significant samples from each sample
for (i in samples){
  name <- sig_result_convert(i)
  input_name <- result_convert(i)
  input <- get(input_name)
  value <- input[!is.na(input$padj) & abs(input$log2FoldChange)>2 & input$padj<0.05, ]
  assign(name, value, envir = globalenv())
}

# List number of differential genes per sample
lapply(samples, sig_summary)

# Identify which samples have the most unique differential genes: 
sig <- list(rnaRes30sig, rnaRes90sig, rnaRes180sig, rnaRes360sig, rnaRes1440sig, rnaRes2880sig, rnaRes4320sig)
names(sig) <- samples
for (i in 1:length(sig)){
  x <- sig[i]
  x <- unlist(lapply(x, rownames))
  y <- sig[-i]
  y <- unlist(lapply(y, rownames))
  print(paste0("There are ", length((x[!x %in% y])), " unique genes after ", names(sig)[i], "min"))
}

## CLUSTERING
# Generate all significant, unique genes 
rnaAllsig <- c()
for (sample in seq(1:7)){
  input <- sig_result_convert(samples[sample])
  input <- get(input)
  list <- rownames(input)
  rnaAllsig <- c(rnaAllsig, list)
  rnaAllsig <- unique(rnaAllsig)
}

# Create LFC matrix 
rnaLFC <- cbind(rnaRes_30$log2FoldChange, rnaRes_90$log2FoldChange, rnaRes_180$log2FoldChange, rnaRes_360$log2FoldChange, 
                rnaRes_1440$log2FoldChange, rnaRes_2880$log2FoldChange, rnaRes_4320$log2FoldChange)
rnaLFC <- data.frame(rnaLFC)
colnames(rnaLFC) <- samples
rnaLFC$mean <- apply(rnaLFC, 1, mean)
rnaLFC$median <- apply(rnaLFC, 1, median)
rnaLFC$differential <- ""
rnaLFC[rownames(rnaLFC) %in% rnaAllsig,]$differential <- TRUE
rnaLFC[rnaLFC$differential=="",]$differential <- FALSE

# Merged count matrix (for downstream filtering)
rnaCounts <- counts(dds)
newmatrix <- list()
for (i in 1:8){
  tmp <- matrix(data=c(rnaCounts[,i], rnaCounts[,i+8]), ncol=2, byrow=F)
  newcol <- apply(tmp, 1, mean)
  newmatrix[[i]] <- newcol
}
rnaCounts <- do.call(cbind, newmatrix)
rnaCounts <- as.data.frame(newmatrix)
rownames(rnaCounts) <- rownames(dds)
colnames(rnaCounts) <- c(0, samples)
rnaCounts$mean <- apply(rnaCounts, 1, mean)
rnaCounts$median <- apply(rnaCounts, 1, median)

# Write these files as txt files for downstream use 
saveRDS(rnaLFC, "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/rna/rnaLFC.rds")
saveRDS(rnaCounts, "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/rna/rnaCounts.rds")

# Cluster and make heatmap: 
set.seed(1949)
k <- 6
normMatrix <- make_norm_matrix(dds)
normMatrix <- combine_reps(normMatrix, c(0, samples))
cluster_data(normMatrix, k)
order <- c(2, 6, 4, 1, 5, 3) #for this seed, this is the correct order for the clusters to be in
pal <- RColorBrewer::brewer.pal(8, "Spectral")
labs <- pal[c(8:6,3:1)]

# Line plots for up and down clusters: 
pdf("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig2/upClusters.pdf", width = 10, height=8.5)
par(mfrow=c(1,1))
par(mar=c(3, 6, 3, 10))
plot(colMeans(normMatrix), type="n", ylim=c(-1.8,1.8), xaxt="n", xlab="PMA Treatment (hours)", ylab="Relative Expression")
lines(colMeans(normMatrix[clusters==2,]), col=labs[1], lwd=3)
lines(colMeans(normMatrix[clusters==6,]), col=labs[2], lwd=3)
lines(colMeans(normMatrix[clusters==4,]), col=labs[3], lwd=3)
axis(1, at=1:8, labels=c(0, 0.5, 1.5, 3, 6, 24, 48, 72))
dev.off()

pdf("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig2/downClusters.pdf", width = 10, height=8.5)
par(mfrow=c(1,1))
par(mar=c(3, 6, 3, 10))
plot(colMeans(normMatrix), type="n", ylim=c(-1.8,1.8), xaxt="n", xlab="PMA Treatment (hours)", ylab="Relative Expression")
lines(colMeans(normMatrix[clusters==1,]), col=labs[4], lwd=3)
lines(colMeans(normMatrix[clusters==5,]), col=labs[5], lwd=3)
lines(colMeans(normMatrix[clusters==3,]), col=labs[6], lwd=3)
axis(1, at=1:8, labels=c(0, 0.5, 1.5, 3, 6, 24, 48, 72))
dev.off()

# Heatmap
matrixOrdered <- normMatrix[order(rowMax(normMatrix), decreasing=TRUE),]
clusterOrder <- clusters[rownames(matrixOrdered)]
#check, proceed only if true
identical(rownames(matrixOrdered), names(clusterOrder))

matrixOrdered <- matrixOrdered[order(match(clusterOrder, order)),]
maxval <- 2.5
matrixOrdered[which(matrixOrdered>maxval)] <- maxval
matrixOrdered[which(matrixOrdered<(-maxval))] <- (-maxval)
clusterOrder <- clusterOrder[rownames(matrixOrdered)]
colnames(matrixOrdered) <- c(0, 0.5, 1.5, 3, 6, 24, 48, 72)

rowColors <- c(rep(labs[1], table(clusters[clusters==2])), 
               rep(labs[2], table(clusters[clusters==6])), 
               rep(labs[3], table(clusters[clusters==4])), 
               rep(labs[4], table(clusters[clusters==1])), 
               rep(labs[5], table(clusters[clusters==5])), 
               rep(labs[6], table(clusters[clusters==3])))

pdf("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig2/rnaHeatmap.pdf", width = 8.5, height = 11, compress = FALSE)
par(mfrow=c(1,1))
par(mar=c(3, 6, 3, 10))
hotmap(matrixOrdered, labrow = FALSE, rowcolors = rowColors)
addlegend(range=c(-maxval, maxval), palette = colorRampPalette(c("deepskyblue2", "black", "gold")), bottominset = 0.2, xoffset=0.08)
dev.off()

## Write each of the genes in the cluster to text files: 
directory <- "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/rna/clusters/"
write(names(clusters[clusters==2]), paste0(directory, "upLate.txt"))
write(names(clusters[clusters==6]), paste0(directory, "upMid.txt"))
write(names(clusters[clusters==4]), paste0(directory, "upEarly.txt"))
write(names(clusters[clusters==1]), paste0(directory, "downLate.txt"))
write(names(clusters[clusters==5]), paste0(directory, "downMid.txt"))
write(names(clusters[clusters==3]), paste0(directory, "downEarly.txt"))

# Combine everything into a summary file 
tp <- list(rnaRes_30, rnaRes_90, rnaRes_180, rnaRes_360, 
           rnaRes_1440, rnaRes_2880, rnaRes_4320)
tp <- lapply(tp, as.data.frame)

foldchange <- 
  cbind(lapply(tp, '[[', 2) |> 
          as.data.frame() |> 
          setNames(paste0("log2FoldChange_", samples)), 
        lapply(tp, '[[', 5) |> 
          as.data.frame() |> 
          setNames(paste0("padj_", samples)))
rownames(foldchange) <- rownames(dds)
allData <- merge(foldchange, counts(dds), by = 0)
colnames(allData)[1] <- "ensembl"
allData$cluster <- ""
allData[allData$ensembl %in% names(clusters[clusters==2]),]$cluster <- "upLate"
allData[allData$ensembl %in% names(clusters[clusters==6]),]$cluster <- "upMid"
allData[allData$ensembl %in% names(clusters[clusters==4]),]$cluster <- "upEarly"
allData[allData$ensembl %in% names(clusters[clusters==1]),]$cluster <- "downLate"
allData[allData$ensembl %in% names(clusters[clusters==5]),]$cluster <- "downMid"
allData[allData$ensembl %in% names(clusters[clusters==3]),]$cluster <- "downEarly"
allData[allData$cluster=="",]$cluster <- "static"
#Save
write.table(allData, file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/rna/rnaSummary.txt", quote = F, sep = "\t")

## Comparing replicates for reviewer 2: 
data <- counts(dds)
data <- data[,c(9,16,1,8)]
data <- data.frame(data)
colnames(data) <- c("rep1_0", "rep1_72", "rep2_0", "rep2_72")
data$gene <- rownames(data)

plot1 <- 
  ggplot(data, aes(x = rep1_0, y = rep2_0)) + 
  geom_point() + 
  geom_smooth(method = lm, color = "black", linewidth = 0.5, se = FALSE) + 
  annotate("text", x = 1E5, y = 750000, label = paste0("R2 = ", (cor(data$rep1_0, data$rep2_0)^2 |> round(3)))) +
  theme_classic()

plot2 <- 
  ggplot(data, aes(x = rep1_72, y = rep2_72)) + 
  geom_point() + 
  geom_smooth(method = lm, color = "black", linewidth = 0.5, se = FALSE) + 
  annotate("text", x = 1E5, y = 750000, label = paste0("R2 = ", (cor(data$rep1_72, data$rep2_72)^2 |> round(3)))) +
  theme_classic()

plot3 <- 
  ggplot(data, aes(x = rep1_0, y = rep1_72)) + 
  geom_point() + 
  geom_smooth(method = lm, color = "black", linewidth = 0.5, se = FALSE) + 
  annotate("text", x = 1E5, y = 750000, label = paste0("R2 = ", (cor(data$rep1_0, data$rep1_72)^2 |> round(3)))) +
  theme_classic()

plot4 <- 
  ggplot(data, aes(x = rep2_0, y = rep2_72)) + 
  geom_point() + 
  geom_smooth(method = lm, color = "black", linewidth = 0.5, se = FALSE) + 
  annotate("text", x = 1E5, y = 750000, label = paste0("R2 = ", (cor(data$rep2_0, data$rep2_72)^2 |> round(3)))) +
  theme_classic()

grid.arrange(plot1, plot2, plot3, plot4, nrow=2)



plot(data$MEGA_RNA_K562_WT_PMA_0000_S_1.1.1, data$MEGA_RNA_K562_WT_PMA_0000_S_2.1.1)
plot(data$MEGA_RNA_K562_WT_PMA_4320_S_1.1.1, data$MEGA_RNA_K562_WT_PMA_4320_S_2.1.1)
plot(data$MEGA_RNA_K562_WT_PMA_0000_S_1.1.1, data$MEGA_RNA_K562_WT_PMA_4320_S_1.1.1)
plot(data$MEGA_RNA_K562_WT_PMA_0000_S_2.1.1, data$MEGA_RNA_K562_WT_PMA_4320_S_2.1.1)

cor(data$MEGA_RNA_K562_WT_PMA_0000_S_1.1.1, data$MEGA_RNA_K562_WT_PMA_0000_S_2.1.1)^2
cor(data$MEGA_RNA_K562_WT_PMA_4320_S_1.1.1, data$MEGA_RNA_K562_WT_PMA_4320_S_2.1.1)^2
cor(data$MEGA_RNA_K562_WT_PMA_0000_S_1.1.1, data$MEGA_RNA_K562_WT_PMA_4320_S_1.1.1)^2
cor(data$MEGA_RNA_K562_WT_PMA_0000_S_2.1.1, data$MEGA_RNA_K562_WT_PMA_4320_S_2.1.1)^2

## AP-1 family transcription: 
load("~/Desktop/geneAnno_gencode_v33_hg38.rda")
ref <- geneAnno1
# ap1 <- c("jun", "fos", "fra", "atf", "maf")
# #which ensembl ids are in AP-1 family 
# ap1Ensembl <- ref[lapply(ap1, grep, ref$hgnc_symbol, ignore.case = TRUE) |> unlist(),]$ensembl_gene_id

names <- c("JUN", "JUNB", "JUND", "FOSL1", "FOS", "FOSB", "FOSL2", "ATF3", "BATF3", "ATF6", "ATF7", "ATF1", "BATF", "ATF5", "ATF2", "ATF4", "MAF", "MAFG", "MAFB", "MAFF", "MAFK", "MAF1")
ap1Ensembl <- ref[ref$hgnc_symbol %in% names,]$ensembl_gene_id
data <- normMatrix[rownames(normMatrix) %in% ap1Ensembl,]
label <- ref[ref$ensembl_gene_id %in% rownames(data),]
label <- label[order(match(label$ensembl_gene_id, rownames(data))),]
label <- label$hgnc_symbol
rownames(data) <- label
colnames(data) <- 1:8

#Plot
par(mfrow=c(1,1))
par(mar=c(3, 6, 3, 10))
plot(colMeans(data), type="n", ylim=c(-3,3), xaxt="n", xlab="PMA Treatment (hours)", ylab="Relative Expression")
apply(data, 1, lines, lwd = 2, col=adjustcolor("gray", alpha.f=0.5))
lines(colMeans(data), col = "blue", lwd = 3)
axis(1, at=1:8, labels=c(0, 0.5, 1.5, 3, 6, 24, 48, 72))


##GGplot version
ggData <- melt(data)
median <- colMeans(data) |> 
  melt()
median$time <- as.integer(rownames(median)) 

ggplot(ggData, aes(x = Var2, y = value)) + 
  geom_line(aes(group = Var1), color = "gray", alpha = 50) + 
  geom_line(data = median, aes(x = time, y = value), lwd = 1, color = "blue")  + 
  geom_text_repel(data = ggData[ggData$Var2 == 8,], aes(label = Var1), nudge_x = 0.5) + 
  scale_x_continuous(labels = c("0", "30", "90", "180", "360", "1440", "2880", "4320"), breaks = 1:8) +
  xlab("PMA Treatment (h)") + ylab("Relative Expression") + 
  theme_classic()
ggsave(last_plot(), file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/reviewFigs/ap1Expression.pdf")






