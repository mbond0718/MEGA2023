## Script to generate unshrunken LFC from DESeq and make boxplots from Figure 3

##ATAC 
peakMatrix <- read.csv("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/raw/atac/MEGA_ATAC_K562_WT_PMA_S_MANUAL_MACS2_counts.tsv", header = T, sep = "\t")
rownames(peakMatrix) <- paste0("peak", 1:nrow(peakMatrix))
time <- factor(c(0, 0, 360, 360, 4320, 4320))
biorep <- factor(c(1, 2, 1, 2, 1, 2))
countMatrix <- peakMatrix[,4:ncol(peakMatrix)]
sampleNames <- colnames(countMatrix)
colData <- data.frame(time, biorep)
rownames(colData) <- sampleNames
dds <- DESeqDataSetFromMatrix(countMatrix, colData, ~biorep + time)
dds <- DESeq(dds, test = "LRT", full = ~biorep + time, reduced = ~biorep)
atacFC <- results(dds, contrast = c("time", "4320", "0"))
atacFC <- cbind(peakMatrix[,c(1:3)], atacFC$log2FoldChange) |> setNames(c("chr", "start", "stop", "lfc"))

##CNR 
#K27
peakMatrix <- read.csv("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/raw/cnr/MEGA_H3K27_K562_WT_PMA_S_counts.tsv", header = T, sep = "\t")
rownames(peakMatrix) <- paste0("peak", 1:nrow(peakMatrix))
time <- factor(c(0, 0, 360, 360, 4320, 4320))
biorep <- factor(c(1, 2, 1, 2, 1, 2))
countMatrix <- peakMatrix[,4:ncol(peakMatrix)]
sampleNames <- colnames(countMatrix)
colData <- data.frame(time, biorep)
rownames(colData) <- sampleNames
dds <- DESeqDataSetFromMatrix(countMatrix, colData, ~biorep + time)
dds <- DESeq(dds, test = "LRT", full = ~biorep + time, reduced = ~biorep)
k27FC <- results(dds, contrast = c("time", "4320", "0"))
k27FC <- cbind(peakMatrix[,c(1:3)], k27FC$log2FoldChange) |> setNames(c("chr", "start", "stop", "lfc"))

#CTCF
peakMatrix <- read.csv("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/raw/cnr/MEGA_CTCF_K562_WT_PMA_S_counts.tsv", header = T, sep = "\t")
rownames(peakMatrix) <- paste0("peak", 1:nrow(peakMatrix))
time <- factor(c(0, 0, 360, 360, 4320, 4320))
biorep <- factor(c(1, 2, 1, 2, 1, 2))
countMatrix <- peakMatrix[,4:ncol(peakMatrix)]
sampleNames <- colnames(countMatrix)
colData <- data.frame(time, biorep)
rownames(colData) <- sampleNames
dds <- DESeqDataSetFromMatrix(countMatrix, colData, ~biorep + time)
dds <- DESeq(dds, test = "LRT", full = ~biorep + time, reduced = ~biorep)
ctcfFC <- results(dds, contrast = c("time", "4320", "0"))
ctcfFC <- cbind(peakMatrix[,c(1:3)], ctcfFC$log2FoldChange) |> setNames(c("chr", "start", "stop", "lfc"))

#Jun
peakMatrix <- read.csv("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/raw/cnr/MEGA_Jun_K562_WT_PMA_S_counts.tsv", header = T, sep = "\t")
rownames(peakMatrix) <- paste0("peak", 1:nrow(peakMatrix))
time <- factor(c(0, 0, 360, 360, 4320, 4320))
biorep <- factor(c(1, 2, 1, 2, 1, 2))
countMatrix <- peakMatrix[,4:ncol(peakMatrix)]
sampleNames <- colnames(countMatrix)
colData <- data.frame(time, biorep)
rownames(colData) <- sampleNames
dds <- DESeqDataSetFromMatrix(countMatrix, colData, ~biorep + time)
dds <- DESeq(dds, test = "LRT", full = ~biorep + time, reduced = ~biorep)
junFC <- results(dds, contrast = c("time", "4320", "0"))
junFC <- cbind(peakMatrix[,c(1:3)], junFC$log2FoldChange) |> setNames(c("chr", "start", "stop", "lfc"))

#Rad21
peakMatrix <- read.csv("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/raw/cnr/MEGA_Rad21_K562_WT_PMA_S_counts.tsv", header = T, sep = "\t")
rownames(peakMatrix) <- paste0("peak", 1:nrow(peakMatrix))
time <- factor(c(0, 0, 4320, 4320))
biorep <- factor(c(1, 2, 1, 2))
countMatrix <- peakMatrix[,4:ncol(peakMatrix)]
sampleNames <- colnames(countMatrix)
colData <- data.frame(time, biorep)
rownames(colData) <- sampleNames
dds <- DESeqDataSetFromMatrix(countMatrix, colData, ~biorep + time)
dds <- DESeq(dds, test = "LRT", full = ~biorep + time, reduced = ~biorep)
radFC <- results(dds, contrast = c("time", "4320", "0"))
radFC <- cbind(peakMatrix[,c(1:3)], radFC$log2FoldChange) |> setNames(c("chr", "start", "stop", "lfc"))

##RNA
txi <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/raw/rna/MEGA_K562_WT_PMA_S_txi.rds")
config <- read.delim("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/raw/rna/config_MEGA_K562_WT_PMA_S_210712_133920.tsv", header = T, sep = '\t')
#create colData for DDS
col_data <- data.frame(time = factor(config$Time), rep = factor(config$Bio_Rep), row.names = config$Name)
#transform txi into dds txi 
dds_txi <- DESeqDataSetFromTximport(txi, colData = col_data, design = ~rep + time)
#transform dds txi into dds
dds <- DESeq(dds_txi, test = "LRT", full = ~rep + time, reduced = ~rep)
names(dds) <- unlist(strsplit(names(dds), "[.]"))[seq(1,(2*length(names(dds))),by=2)]
#unshrunken fc 
rnaFC <- results(dds, contrast = c("time", "4320", "0"))
rnaFC <- data.frame(rownames(rnaFC), rnaFC$log2FoldChange) |> setNames(c("gene_id", "lfc"))
genes <- loadDb("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/ref/hg38txdb.sqlite") |> genes()
genes <- data.frame(genes)
genes$seqnames <- paste0("chr", genes$seqnames)
rownames(genes) <- genes$gene_id
rnaFC <- merge(rnaFC, genes, by = "gene_id")
rnaFC <- rnaFC[,c(1,3:5,2)]

## Convert to GRanges
#ATAC 
atacGR <- GRanges(seqnames = Rle(atacFC$chr), 
                  ranges = IRanges(start = atacFC$start, end = atacFC$stop), 
                  lfc = atacFC$lfc,
                  name = rownames(atacFC))
#Enhancer
enhGR <- GRanges(seqnames = Rle(k27FC$chr), 
                 ranges = IRanges(start = k27FC$start, end = k27FC$stop), 
                 lfc = k27FC$lfc,
                 name = rownames(k27FC))
#Jun
junGR <- GRanges(seqnames = Rle(junFC$chr), 
                 ranges = IRanges(start = junFC$start, end = junFC$stop), 
                 lfc = junFC$lfc,
                 name = rownames(junFC))
#CTCF 
ctcfGR <- GRanges(seqnames = Rle(ctcfFC$chr), 
                  ranges = IRanges(start = ctcfFC$start, end = ctcfFC$stop), 
                  lfc = ctcfFC$lfc,
                  name = rownames(ctcfFC))
#Cohesin
radGR <- GRanges(seqnames = Rle(radFC$chr), 
                 ranges = IRanges(start = radFC$start, end = radFC$stop), 
                 lfc = radFC$lfc,
                 name = rownames(radFC))

#Cohesin
rnaGR <- GRanges(seqnames = Rle(rnaFC$seqnames), 
                 ranges = IRanges(start = rnaFC$start, end = rnaFC$stop), 
                 lfc = rnaFC$lfc,
                 name = rnaFC$gene_id)
rnaGR <- promoters(rnaGR)

## Loops
#Loops
loops <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/diffLoops.rds")
loops$loop <- rownames(loops)
gainedLoops <- loops[loops$class=="gained",]
lostLoops <- loops[loops$class=="lost",]
staticLoops <- loops[loops$class=="static",]

#Loop anchors
gainedAnc <- as_ginteractions(gainedLoops)
lostAnc <- as_ginteractions(lostLoops)
staticAnc <- as_ginteractions(staticLoops) |> 
  swapAnchors()
#Loop internal
gainedIn <- GRanges(seqnames = seqnames(anchors(gainedAnc)$first), 
                    ranges = IRanges(start = end(anchors(gainedAnc)$first), end = start(anchors(gainedAnc)$second)), 
                    loop = gainedAnc$loop)
lostIn <- GRanges(seqnames = seqnames(anchors(lostAnc)$first), 
                  ranges = IRanges(start = end(anchors(lostAnc)$first), end = start(anchors(lostAnc)$second)), 
                  loop = lostAnc$loop)
staticIn <- GRanges(seqnames = seqnames(anchors(staticAnc)$first), 
                    ranges = IRanges(start = end(anchors(staticAnc)$first), end = start(anchors(staticAnc)$second)), 
                    loop = staticAnc$loop)

intersect <- list(atacGR, enhGR, junGR, ctcfGR, radGR, rnaGR)
names <- c("ATAC", "H3K27ac", "Jun", "CTCF", "Rad21", "RNA")

#Palette
pal <- brewer.pal(8, "YlGnBu")

wilcox <- list()
medians <- list()
for (i in 1:length(intersect)){
  print(i)
  feat <- intersect[[i]]
  diff <- feat[feat$differential==TRUE,]
  ## Anchors - all
  gainedA <-   data.frame(subsetByOverlaps(feat, gainedAnc)$name, 
                          subsetByOverlaps(feat, gainedAnc)$lfc, 
                          "gained")
  lostA <-   data.frame(subsetByOverlaps(feat, lostAnc)$name, 
                        subsetByOverlaps(feat, lostAnc)$lfc, 
                        "lost")
  staticA <-   data.frame(subsetByOverlaps(feat, staticAnc)$name, 
                          subsetByOverlaps(feat, staticAnc)$lfc, 
                          "static")

  anc <- list(gainedA, lostA, staticA)
  for (a in 1:length(anc)){
    tmp <- anc[[a]]
    colnames(tmp) <- c("name", "lfc", "group")
    anc[[a]] <- tmp
  }
  anc <- do.call(rbind, anc)
  anc$intersection <- "anchor"
  
  ## Internal - all
  gainedI <-   data.frame(subsetByOverlaps(feat, gainedIn)$name, 
                          subsetByOverlaps(feat, gainedIn)$lfc, 
                          "gained")
  lostI <-   data.frame(subsetByOverlaps(feat, lostIn)$name, 
                        subsetByOverlaps(feat, lostIn)$lfc, 
                        "lost")
  staticI <-   data.frame(subsetByOverlaps(feat, staticIn)$name, 
                          subsetByOverlaps(feat, staticIn)$lfc, 
                          "static")
  
  int <- list(gainedI, lostI, staticI)
  for (a in 1:length(int)){
    tmp <- int[[a]]
    colnames(tmp) <- c("name", "lfc", "group")
    int[[a]] <- tmp
  }
  int <- do.call(rbind, int)
  int$intersection <- "internal"
  
  data <- rbind(anc, int)
  data$group <- factor(data$group, levels = c("gained", "static", "lost"))
  data <- na.omit(data)
  
  median <- 
    c(median(data[data$group=="gained" & data$intersection=="anchor",]$lfc), 
      median(data[data$group=="gained" & data$intersection=="internal",]$lfc), 
      median(data[data$group=="lost" & data$intersection=="anchor",]$lfc), 
      median(data[data$group=="lost" & data$intersection=="internal",]$lfc)) |> 
    as.data.frame() |> 
    cbind(c("gained", "gained", "lost", "lost")) |> 
    cbind(c("anchor", "interior", "anchor", "interior")) |> 
    cbind(rep(names[i], 4)) |> 
    setNames(c("lfc", "group", "intersection", "feature"))
  medians[[i]] <- median
  
  ggplot(data, aes(x = group, y = lfc, fill = intersection)) + 
    geom_boxplot(outlier.shape = NA) + 
    scale_fill_manual(values = rep(c(pal[8], pal[6]), 3)) + 
    ylab("LFC") + 
    coord_cartesian(ylim=c(-6.5,11)) + 
    #ylim(-6.5,11) + 
    ggtitle(paste0(names[i], " LFC")) + 
    geom_hline(yintercept = 0, lty=3) +
    theme_classic() + theme(legend.position = "none")
  #ggsave(last_plot(), filename = paste0("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig3/", names[i], "boxplotUnshrunken.pdf"))
  
  #Wilcox, comparing each to their corresponding static:
  stats <- data.frame(
    (wilcox.test(data[data$group=="gained" & data$intersection=="anchor",]$lfc, data[data$group=="static" & data$intersection=="anchor",]$lfc))$p.value,
    (wilcox.test(data[data$group=="gained" & data$intersection=="internal",]$lfc, data[data$group=="static" & data$intersection=="internal",]$lfc))$p.value,
    (wilcox.test(data[data$group=="lost" & data$intersection=="anchor",]$lfc, data[data$group=="static" & data$intersection=="anchor",]$lfc))$p.value,
    (wilcox.test(data[data$group=="lost" & data$intersection=="internal",]$lfc, data[data$group=="static" & data$intersection=="internal",]$lfc))$p.value
  )
  colnames(stats) <- c("gainedAnc", "gainedInt", "lostAnc", "lostInt")
  wilcox[[i]] <- stats
  
  #Wilcox, comapring each to 0 
  # stats <- data.frame(
  #   wilcox.test(data[data$group=="gained" & data$intersection=="anchor",]$lfc)$p.value, 
  #   wilcox.test(data[data$group=="gained" & data$intersection=="internal",]$lfc)$p.value, 
  #   wilcox.test(data[data$group=="static" & data$intersection=="anchor",]$lfc)$p.value, 
  #   wilcox.test(data[data$group=="static" & data$intersection=="internal",]$lfc)$p.value, 
  #   wilcox.test(data[data$group=="lost" & data$intersection=="anchor",]$lfc)$p.value, 
  #   wilcox.test(data[data$group=="lost" & data$intersection=="internal",]$lfc)$p.value
  # ) |> 
  #   setNames(c("gainedAnc", "gainedInt", "staticAnc", "staticInt", "lostAnc", "lostInt"))
  #   wilcox[[i]] <- stats
}

  


#Statistics 
wilcox <- do.call(rbind, wilcox)
rownames(wilcox) <- names
wilcox <- apply(wilcox, 1, p.adjust)

cutoffs <- c(1, 0.05, 0.01, 0.001, 0)

signif <- symnum(wilcox, corr = FALSE, na = FALSE,
                 cutpoints = cutoffs,
                 symbols = c("***", "**", "*", "ns"))

s <- signif
colnames(s) <- c("ATAC", "H3K27ac", "JUN", "CTCF", "RAD21", "RNA")
s[s=="***"] <- "*"


#Median plot
median <- do.call(rbind, medians)
median$feature <- factor(median$feature, levels = names)

ggplot(median[median$group=="gained",], aes(x = feature, y = lfc, color = intersection)) + 
  geom_point(shape = 8) + ylim(c(-2,4)) + scale_color_manual(values = c(pal[8], pal[6])) + 
  theme_classic() + theme(legend.position = "none") + geom_hline(yintercept = 0, lty=3) + 
  ylab("Feature Fold-Change")
ggsave(last_plot(), file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig3/summaryGainedUnshrunk.pdf")

ggplot(median[median$group=="lost",], aes(x = feature, y = lfc, color = intersection, shape = feature)) + 
  geom_point() + ylim(c(-2,4)) + scale_color_manual(values = c(pal[8], pal[6])) + 
  theme_classic() + theme(legend.position = "none") + geom_hline(yintercept = 0, lty=3) + 
  scale_shape_manual(values = c(8, 8, 8, 8, 8, 16)) + ylab("Feature Fold-Change")
ggsave(last_plot(), file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig3/summaryLostUnshrunk.pdf")

