## Script to intersect differential events from each dataset with differential loops:

## Load libraries -----------------------------------------
library(plotgardener)
library(purrr)
library(RColorBrewer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)
library(InteractionSet)
library(plyranges)
library(mariner)


## Read in data --------------------------------------------
#Loops
#loops <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/diffLoops.rds")
loops <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/loopLFC.rds")
loops$loop <- rownames(loops)
gainedLoops <- loops[loops$class=="gained",]
lostLoops <- loops[loops$class=="lost",]
staticLoops <- loops[loops$class=="static",]

#ATAC 
atac <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/atac/atacLFC.rds")
#Enhancer
enh <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/cnr/H3K27LFC.rds")
#Jun
jun <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/cnr/JunLFC.rds")
#CTCF
ctcf <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/cnr/CTCFLFC.rds")
#Rad21
rad <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/cnr/Rad21LFC.rds")
#RNA
rna <- read.delim("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/rna/rnaSummary.txt")

## Convert to GRanges: 
#Loop anchors
gainedAnc <- as_ginteractions(gainedLoops)
lostAnc <- as_ginteractions(lostLoops)
staticAnc <- as_ginteractions(staticLoops) |> 
  swapAnchors()
diffAnc <- c(gainedAnc, lostAnc)

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

#ATAC 
atacGR <- GRanges(seqnames = Rle(atac$chr), 
                  ranges = IRanges(start = atac$start, end = atac$stop), 
                  lfc = atac$`4320`,
                  name = rownames(atac), 
                  diff = atac$differential)
diffATAC <- atacGR[atacGR$diff==TRUE]
#Enhancer
enhGR <- GRanges(seqnames = Rle(enh$chr), 
                 ranges = IRanges(start = enh$start, end = enh$stop), 
                 lfc = enh$`4320`,
                 name = rownames(enh), 
                 diff = enh$differential)
#Jun
junGR <- GRanges(seqnames = Rle(jun$chr), 
                 ranges = IRanges(start = jun$start, end = jun$stop), 
                 lfc = jun$`4320`,
                 name = rownames(jun), 
                 diff = jun$differential)
#CTCF 
ctcfGR <- GRanges(seqnames = Rle(ctcf$chr), 
                  ranges = IRanges(start = ctcf$start, end = ctcf$stop), 
                  lfc = ctcf$`4320`,
                  name = rownames(ctcf), 
                  diff = ctcf$differential)
#Cohesin
radGR <- GRanges(seqnames = Rle(rad$chr), 
                 ranges = IRanges(start = rad$start, end = rad$stop), 
                 lfc = rad$`4320`,
                 name = rownames(rad), 
                 diff = rad$differential)
#Genes
genes <- loadDb("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/ref/hg38txdb.sqlite") |> genes()
genes <- data.frame(genes)
genes$seqnames <- paste0("chr", genes$seqnames)
rownames(genes) <- genes$gene_id

colnames(rna)[1] <- "gene_id"
#rna[!rna$cluster=="static",]
rna <- merge(rna, genes, by="gene_id")
rna <- rna[,c(33:35,1,32,8,23,31)]
rna$counts <- (rna$MEGA_RNA_K562_WT_PMA_4320_S_1.1.1+ rna$MEGA_RNA_K562_WT_PMA_4320_S_2.1.1)/2


#ABUNDANCE/TPM
txi <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/raw/rna/MEGA_K562_WT_PMA_S_txi.rds")
tpm <- data.frame(txi$abundance)
rownames(tpm) <- (unlist(strsplit(rownames(tpm), "[.]")))[seq(1,(2*nrow(tpm)),by=2)]
tpm <- tpm[,c(8,16)]
tpm$tpm <- (tpm$MEGA_RNA_K562_WT_PMA_4320_S_2.1.1+tpm$MEGA_RNA_K562_WT_PMA_4320_S_1.1.1)/2
tpm$gene_id <- rownames(tpm)
tpm <- tpm[,c(4,3)]

#merge
rna <- merge(rna, tpm, by = "gene_id")

rnaGR <- GRanges(seqnames = Rle(rna$seqnames), 
                 ranges = IRanges(start = rna$start, end = rna$end), 
                 cluster = rna$cluster,
                 name = rna$gene_id, 
                 fc = rna$log2FoldChange_4320, 
                 counts = rna$counts, 
                 tpm = rna$tpm)
rnaGR <- promoters(rnaGR)
diffRNA <- rnaGR[!rnaGR$cluster=="static"]



intersect <- list(atacGR, enhGR, junGR, ctcfGR, radGR, rnaGR)
names <- c("ATAC", "H3K27ac", "Jun", "CTCF", "Rad21", "RNA")

## Generate pie charts: 
#Loop
pie <- 
  data.frame(c(length(staticAnc), length(gainedAnc), length(lostAnc))) |> 
  cbind(c("static", "gained", "lost")) |> 
  setNames(c("num", "class"))
pie$class <- factor(pie$class, levels = c("static", "lost", "gained"))

ggplot(pie, aes(x = "", y = num, fill = class)) + 
  geom_bar(stat = "identity", width = 1) + 
  coord_polar("y", start=215) + 
  scale_fill_manual(values = c("gray60", "#931C1E", "#227B7F")) + 
  theme_void() + theme(legend.position="none") + 
  geom_text(aes(label = num), color = "white")
ggsave("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig2/loopPie.pdf", last_plot())

#Genes
pie <- 
  data.frame(c(length(rnaGR[rnaGR$cluster %in% c("static")]), 
               length(rnaGR[rnaGR$cluster %in% c("upEarly", "upLate", "upMid")]), 
               length(rnaGR[rnaGR$cluster %in% c("downEarly", "downMid", "downLate")]))) |> 
  cbind(c("static", "gained", "lost")) |> 
  setNames(c("num", "class"))
pie$class <- factor(pie$class, levels = c("static", "gained", "lost"))

ggplot(pie, aes(x = "", y = num, fill = class)) + 
  geom_bar(stat = "identity", width = 1) + 
  coord_polar("y", start=80) + 
  scale_fill_manual(values = c("gray60", "#227B7F", "#931C1E")) + 
  theme_void() + theme(legend.position="none") + 
  geom_text(aes(label = num), color = "white")
ggsave("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig2/rnaPie.pdf", last_plot())

#ATAC 
pie <- 
  data.frame(c(length(atacGR[atacGR$diff==FALSE]),
               length(atacGR[atacGR$diff==TRUE & atacGR$lfc>0]), 
               length(atacGR[atacGR$diff==TRUE & atacGR$lfc<0]))) |> 
  cbind(c("static", "gained", "lost")) |> 
  setNames(c("num", "class"))
pie$class <- factor(pie$class, levels = c("static", "gained", "lost"))

ggplot(pie, aes(x = "", y = num, fill = class)) + 
  geom_bar(stat = "identity", width = 1) + 
  coord_polar("y", start=80) + 
  scale_fill_manual(values = c("gray60", "#227B7F", "#931C1E")) + 
  theme_void() + theme(legend.position="none") + 
  geom_text(aes(label = num), color = "white")
ggsave("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig2/atacPie.pdf", last_plot())

## Intersect loops with rna
subsetByOverlaps(diffAnc, diffRNA)
subsetByOverlaps(gainedAnc, diffRNA)
subsetByOverlaps(lostAnc, diffRNA)

subsetByOverlaps(diffAnc, diffRNA)
length(subsetByOverlaps(gainedAnc, diffRNA[diffRNA$fc>0,])$loop)
length(subsetByOverlaps(gainedAnc, diffRNA[diffRNA$fc<0,])$loop)
length(subsetByOverlaps(lostAnc, diffRNA[diffRNA$fc>0,])$loop)
length(subsetByOverlaps(lostAnc, diffRNA[diffRNA$fc<0,])$loop)

data <- 
  data.frame(c(length(subsetByOverlaps(gainedAnc, diffRNA[diffRNA$fc>0,])$loop), 
               length(subsetByOverlaps(gainedAnc, diffRNA[diffRNA$fc<0,])$loop), 
               length(subsetByOverlaps(lostAnc, diffRNA[diffRNA$fc<0,])$loop), 
               length(subsetByOverlaps(lostAnc, diffRNA[diffRNA$fc>0,])$loop))) |> 
  cbind(c("gainedUp", "gainedDown", "lostDown", "lostUp")) |> 
  setNames(c("num", "class"))
data$class <- factor(data$class, levels = c("gainedUp", "gainedDown", "lostUp", "lostDown"))

ggplot(data, aes(x = class, y = num, fill = class)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_manual(values = c("#227B7F", "#658E8D", "#91575A", "#931C1E")) + 
  theme_classic() + theme(legend.position = "none") + xlab("") + ylab("Frequency")
ggsave("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig2/rnaConcordance.pdf", last_plot())

#RNA and interior: 
data <- 
  data.frame(c(length(subsetByOverlaps(gainedIn, diffRNA[diffRNA$fc>0,])$loop), 
               length(subsetByOverlaps(gainedIn, diffRNA[diffRNA$fc<0,])$loop), 
               length(subsetByOverlaps(lostIn, diffRNA[diffRNA$fc<0,])$loop), 
               length(subsetByOverlaps(lostIn, diffRNA[diffRNA$fc>0,])$loop))) |> 
  cbind(c("gainedUp", "gainedDown", "lostDown", "lostUp")) |> 
  setNames(c("num", "class"))
data$class <- factor(data$class, levels = c("gainedUp", "gainedDown", "lostUp", "lostDown"))

ggplot(data, aes(x = class, y = num, fill = class)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_manual(values = c("#227B7F", "#658E8D", "#91575A", "#931C1E")) + 
  theme_classic() + theme(legend.position = "none") + xlab("") + ylab("Frequency")
ggsave(last_plot(), file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/sup2/interiorConcordance.pdf")

## Counts of genes that are increasing at interiors of gained and lost loops
subsetByOverlaps(diffRNA[diffRNA$fc>0], gainedIn)$counts
subsetByOverlaps(diffRNA[diffRNA$fc>0], lostIn)$counts

#diff genes
data <- 
  rbind(data.frame(subsetByOverlaps(diffRNA[diffRNA$fc>0], gainedIn)$counts, "gained") |> setNames(c("counts", "loop")), 
      data.frame(subsetByOverlaps(diffRNA[diffRNA$fc>0], lostIn)$counts, "lost") |> setNames(c("counts", "loop")))

#all genes
rnaNaOm <- rnaGR[!is.na(rnaGR$fc)]

data <- 
  rbind(data.frame(subsetByOverlaps(rnaNaOm[rnaNaOm$fc>0], gainedIn)$tpm, "gained") |> setNames(c("tpm", "loop")), 
        data.frame(subsetByOverlaps(rnaNaOm[rnaNaOm$fc>0], lostIn)$tpm, "lost") |> setNames(c("tpm", "loop")))

#Log transform
log <- data
log$tpm <- log2(log$tpm)

#boxplot
ggplot(log, aes(x = loop, y = tpm, fill = loop)) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("#227B7F", "#931C1E")) + 
  theme_classic() + theme(legend.position = "none") + ylab("log2(tpm)")

#density 
ggplot(log, aes(x = tpm, fill = loop)) +
  geom_density(alpha = 0.8) + 
  scale_fill_manual(values = c("#227B7F", "#931C1E")) + 
  theme_classic() + xlab("log2(tpm)")
ggsave(last_plot(), file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/sup2/tpmDensity.pdf")

wilcox.test(log[log$loop=="gained",]$tpm, log[log$loop=="lost",]$tpm)

#barplot with errors
error <- 
  data.frame(mean = c(mean(data[data$loop=="gained",]$counts), 
                    mean(data[data$loop=="lost",]$counts)), 
           sd = c(sd(data[data$loop=="gained",]$counts), 
                  sd(data[data$loop=="lost",]$counts)), 
           loop = c("gained", "lost"))

ggplot(error, aes(x = loop, y = mean, fill = loop)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("#227B7F", "#931C1E")) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width = 0.2) + 
  theme_classic() + theme(legend.position = "none")

ggplot(data, aes(x = loop, y = counts, fill = loop)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("#227B7F", "#931C1E")) + 
  theme_classic() + theme(legend.position = "none")

#Log

error <- 
  data.frame(mean = c(mean(log[log$loop=="gained",]$counts), 
                      mean(log[log$loop=="lost",]$counts)), 
             sd = c(sd(log[log$loop=="gained",]$counts), 
                    sd(log[log$loop=="lost",]$counts)), 
             loop = c("gained", "lost"))
ggplot(error, aes(x = loop, y = mean, fill = loop)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("#227B7F", "#931C1E")) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width = 0.2) + 
  theme_classic() + theme(legend.position = "none") + ylab("log2(counts)")

#Stats: 
wilcox.test(data[data$loop=="gained",]$counts, data[data$loop=="lost",]$counts)
wilcox.test(log[log$loop=="gained",]$counts, log[log$loop=="lost",]$counts)


gainedGenes <- subsetByOverlaps(diffRNA[diffRNA$fc>0], gainedIn)$name
lostGenes <- subsetByOverlaps(diffRNA[diffRNA$fc>0], lostIn)$name
  

## Intersect loops with atac
gU <- subsetByOverlaps(gainedAnc, diffATAC[diffATAC$lfc>0,])$loop
gD <- subsetByOverlaps(gainedAnc, diffATAC[diffATAC$lfc<0,])$loop
lU <- subsetByOverlaps(lostAnc, diffATAC[diffATAC$lfc>0,])$loop
lD <- subsetByOverlaps(lostAnc, diffATAC[diffATAC$lfc<0,])$loop

remove <- c(gU[gU %in% gD], lU[lU %in% lD])
filtDiff <- diffAnc[!diffAnc$loop %in% remove]
filtGain <- filtDiff[filtDiff$lfc>0]
filtLoss <- filtDiff[filtDiff$lfc<0]

subsetByOverlaps(filtDiff, diffATAC)
subsetByOverlaps(filtGain, diffATAC)
subsetByOverlaps(filtLoss, diffATAC)

data <- 
  data.frame(c(length(subsetByOverlaps(filtGain, diffATAC[diffATAC$lfc>0,])$loop), 
               length(subsetByOverlaps(filtGain, diffATAC[diffATAC$lfc<0,])$loop), 
               length(subsetByOverlaps(filtLoss, diffATAC[diffATAC$lfc<0,])$loop), 
               length(subsetByOverlaps(filtLoss, diffATAC[diffATAC$lfc>0,])$loop))) |> 
  cbind(c("gainedUp", "gainedDown", "lostDown", "lostUp")) |> 
  setNames(c("num", "class"))
data$class <- factor(data$class, levels = c("gainedUp", "gainedDown", "lostUp", "lostDown"))

ggplot(data, aes(x = class, y = num, fill = class)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_manual(values = c("#227B7F", "#658E8D", "#91575A", "#931C1E")) + 
  theme_classic() + theme(legend.position = "none") + xlab("") + ylab("Frequency")
ggsave("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig2/atacConcordance.pdf", last_plot())

## Intersect loops interiors with atac
gU <- subsetByOverlaps(gainedIn, diffATAC[diffATAC$lfc>0,])$loop
gD <- subsetByOverlaps(gainedIn, diffATAC[diffATAC$lfc<0,])$loop
lU <- subsetByOverlaps(lostIn, diffATAC[diffATAC$lfc>0,])$loop
lD <- subsetByOverlaps(lostIn, diffATAC[diffATAC$lfc<0,])$loop

remove <- c(gU[gU %in% gD], lU[lU %in% lD])
filtDiff <- diffAnc[!diffAnc$loop %in% remove]
filtGain <- filtDiff[filtDiff$lfc>0]
filtLoss <- filtDiff[filtDiff$lfc<0]

subsetByOverlaps(filtDiff, diffATAC)
subsetByOverlaps(filtGain, diffATAC)
subsetByOverlaps(filtLoss, diffATAC)

data <- 
  data.frame(c(length(subsetByOverlaps(filtGain, diffATAC[diffATAC$lfc>0,])$loop), 
               length(subsetByOverlaps(filtGain, diffATAC[diffATAC$lfc<0,])$loop), 
               length(subsetByOverlaps(filtLoss, diffATAC[diffATAC$lfc<0,])$loop), 
               length(subsetByOverlaps(filtLoss, diffATAC[diffATAC$lfc>0,])$loop))) |> 
  cbind(c("gainedUp", "gainedDown", "lostDown", "lostUp")) |> 
  setNames(c("num", "class"))
data$class <- factor(data$class, levels = c("gainedUp", "gainedDown", "lostUp", "lostDown"))

ggplot(data, aes(x = class, y = num, fill = class)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_manual(values = c("#227B7F", "#658E8D", "#91575A", "#931C1E")) + 
  theme_classic() + theme(legend.position = "none") + xlab("") + ylab("Frequency")

subsetByOverlaps(diffAnc, diffATAC)
subsetByOverlaps(gainedAnc, diffATAC[diffATAC$lfc>0]) #gained concord
subsetByOverlaps(gainedAnc, diffATAC[diffATAC$lfc<0]) #gained not concord
subsetByOverlaps(lostAnc, diffATAC[diffATAC$lfc>0]) #lost not concord
subsetByOverlaps(lostAnc, diffATAC[diffATAC$lfc<0]) #lost concord

subsetByOverlaps(diffAnc, diffRNA)
subsetByOverlaps(gainedAnc, diffRNA[diffRNA$fc>0]) #gained concord
subsetByOverlaps(gainedAnc, diffRNA[diffRNA$fc<0]) #gained not concord
subsetByOverlaps(lostAnc, diffRNA[diffRNA$fc>0]) #lost not concord
subsetByOverlaps(lostAnc, diffRNA[diffRNA$fc<0]) #lost concord


subsetByOverlaps(diffAnc, diffRNA)
subsetByOverlaps(gainedAnc, diffRNA[diffRNA$fc>0,])$loop
subsetByOverlaps(gainedAnc, diffRNA[diffRNA$fc<0,])$loop
subsetByOverlaps(lostAnc, diffRNA[diffRNA$fc>0,])$loop
subsetByOverlaps(lostAnc, diffRNA[diffRNA$fc<0,])$loop

subsetByOverlaps(diffAnc, diffATAC)
a <- subsetByOverlaps(gainedAnc, diffATAC[diffATAC$lfc>0,])$loop
b <- subsetByOverlaps(gainedAnc, diffATAC[diffATAC$lfc<0,])$loop
c <- subsetByOverlaps(lostAnc, diffATAC[diffATAC$lfc>0,])$loop
d <- subsetByOverlaps(lostAnc, diffATAC[diffATAC$lfc<0,])$loop


## Ranked loops by padj, relationship with interior gene LFC

#Genes
genes <- loadDb("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/ref/hg38txdb.sqlite") |> genes()
genes <- data.frame(genes)
genes$seqnames <- paste0("chr", genes$seqnames)
rownames(genes) <- genes$gene_id

colnames(rna)[1] <- "gene_id"
rna <- merge(rna, genes, by="gene_id")
rna <- rna[,c(33:35,1,32,8)]
rnaGR <- GRanges(seqnames = Rle(rna$seqnames), 
                 ranges = IRanges(start = rna$start, end = rna$end), 
                 cluster = rna$cluster,
                 name = rna$gene_id, 
                 fc = rna$log2FoldChange_4320)
diffRNA <- rnaGR[!rnaGR$cluster=="static"]

loopPadj <- read.delim("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/hicSummary.txt")
loopPadj <- loopPadj[loopPadj$cluster=="lost",]
lostIn <-GRanges(seqnames = Rle(loopPadj$seqnames1), 
                 ranges = IRanges(start = loopPadj$end1, end = loopPadj$start2), 
                 loop = loopPadj$Row.names, 
                 padj = loopPadj$padj_4320, lfc = loopPadj$log2FoldChange_4320)
lostIn <- lostIn[order(lostIn$padj),]
lostIn <- lostIn[order(lostIn$lfc),]

data <- list()
for (i in 1:length(lostIn)){
  print(i)
  subset <- lostIn[1:i]
  olap <- subsetByOverlaps(diffRNA, subset)$fc |> 
    na.omit()
  median <- median(olap)
  data[[i]] <- median
}

data <- do.call(rbind, data) |> 
  as.data.frame()
colnames(data) <- "medianGeneFC"
data$rank <- 1:503

ggplot(data, aes(x = rank, y = medianGeneFC)) + 
  geom_line() + theme_minimal() + 
  ylab("Cumulative interior DEG LFC") + xlab("Lost loops ranked by LFC")




