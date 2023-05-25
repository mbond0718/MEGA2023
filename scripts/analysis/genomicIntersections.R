## INTERSECTION OF LOOP CLASSES WITH REGULATORY FEATURES
## This script generates the boxplots showing genomic intersections from figure 3

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

## Define function and parameters -------------------------- 
geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

#Palette
pal <- brewer.pal(8, "YlGnBu")

## Read in data --------------------------------------------
#Loops
loops <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/diffLoops.rds")
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
rna <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/rna/rnaLFC.rds")

## Convert to GRanges: 
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

#ATAC 
atacGR <- GRanges(seqnames = Rle(atac$chr), 
                     ranges = IRanges(start = atac$start, end = atac$stop), 
                     lfc = atac$`4320`,
                  name = rownames(atac), 
                  differential = atac$differential)
#Enhancer
enhGR <- GRanges(seqnames = Rle(enh$chr), 
                 ranges = IRanges(start = enh$start, end = enh$stop), 
                 lfc = enh$`4320`,
                 name = rownames(enh), 
                 differential = enh$differential)
#Jun
junGR <- GRanges(seqnames = Rle(jun$chr), 
                 ranges = IRanges(start = jun$start, end = jun$stop), 
                 lfc = jun$`4320`,
                 name = rownames(jun), 
                 differential = jun$differential)
#CTCF 
ctcfGR <- GRanges(seqnames = Rle(ctcf$chr), 
                  ranges = IRanges(start = ctcf$start, end = ctcf$stop), 
                  lfc = ctcf$`4320`,
                  name = rownames(ctcf), 
                  differential = ctcf$differential)
#Cohesin
radGR <- GRanges(seqnames = Rle(rad$chr), 
                 ranges = IRanges(start = rad$start, end = rad$stop), 
                 lfc = rad$`4320`,
                 name = rownames(rad), 
                 differential = rad$differential)
#Genes
genes <- loadDb("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/ref/hg38txdb.sqlite") |> genes()
genes <- data.frame(genes)
genes$seqnames <- paste0("chr", genes$seqnames)
rownames(genes) <- genes$gene_id

rna <- merge(rna, genes, by=0)
rna <- rna[,c(12:14,8,17,11)]
rnaGR <- GRanges(seqnames = Rle(rna$seqnames), 
                 ranges = IRanges(start = rna$start, end = rna$end), 
                 lfc = rna$`4320`,
                 name = rna$gene_id, 
                 differential = rna$differential)
rnaGR <- promoters(rnaGR)


intersect <- list(atacGR, enhGR, junGR, ctcfGR, radGR, rnaGR)
names <- c("ATAC", "H3K27ac", "Jun", "CTCF", "Rad21", "RNA")

wilcox <- list()
medians <- list()
for (i in 1:length(intersect)){
  print(i)
  feat <- intersect[[i]]
  diff <- feat[feat$differential==TRUE,]
  ## Anchors - all
  gainedA <-   data.frame(subsetByOverlaps(feat, gainedAnc)$name, 
                          subsetByOverlaps(feat, gainedAnc)$lfc, 
                          "gained", 
                          "all")
  lostA <-   data.frame(subsetByOverlaps(feat, lostAnc)$name, 
                          subsetByOverlaps(feat, lostAnc)$lfc, 
                        "lost", 
                        "all")
  staticA <-   data.frame(subsetByOverlaps(feat, staticAnc)$name, 
                          subsetByOverlaps(feat, staticAnc)$lfc, 
                          "static", 
                          "all")
  #Anchors - differential
  diffgainedA <-   data.frame(subsetByOverlaps(diff, gainedAnc)$name, 
                          subsetByOverlaps(diff, gainedAnc)$lfc, 
                          "gained", 
                          "diff")
  difflostA <-   data.frame(subsetByOverlaps(diff, lostAnc)$name, 
                        subsetByOverlaps(diff, lostAnc)$lfc, 
                        "lost", 
                        "diff")
  diffstaticA <-   data.frame(subsetByOverlaps(diff, staticAnc)$name, 
                          subsetByOverlaps(diff, staticAnc)$lfc, 
                          "static", 
                          "diff")
  
  anc <- list(gainedA, lostA, staticA, diffgainedA, difflostA, diffstaticA)
  for (a in 1:length(anc)){
    tmp <- anc[[a]]
    colnames(tmp) <- c("name", "lfc", "group", "diff")
    anc[[a]] <- tmp
  }
  anc <- do.call(rbind, anc)
  anc$intersection <- "anchor"
  
  ## Internal - all
  gainedI <-   data.frame(subsetByOverlaps(feat, gainedIn)$name, 
                          subsetByOverlaps(feat, gainedIn)$lfc, 
                          "gained",
                          "all")
  lostI <-   data.frame(subsetByOverlaps(feat, lostIn)$name, 
                        subsetByOverlaps(feat, lostIn)$lfc, 
                        "lost", 
                        "all")
  staticI <-   data.frame(subsetByOverlaps(feat, staticIn)$name, 
                          subsetByOverlaps(feat, staticIn)$lfc, 
                          "static", 
                          "all")
  ## Internal - differential 
  diffgainedI <-   data.frame(subsetByOverlaps(diff, gainedIn)$name, 
                          subsetByOverlaps(diff, gainedIn)$lfc, 
                          "gained",
                          "diff")
  difflostI <-   data.frame(subsetByOverlaps(diff, lostIn)$name, 
                        subsetByOverlaps(diff, lostIn)$lfc, 
                        "lost", 
                        "diff")
  diffstaticI <-   data.frame(subsetByOverlaps(diff, staticIn)$name, 
                          subsetByOverlaps(diff, staticIn)$lfc, 
                          "static", 
                          "diff")
  
  int <- list(gainedI, lostI, staticI, diffgainedI, difflostI, diffstaticI)
  for (a in 1:length(int)){
    tmp <- int[[a]]
    colnames(tmp) <- c("name", "lfc", "group", "diff")
    int[[a]] <- tmp
  }
  int <- do.call(rbind, int)
  int$intersection <- "internal"
  
  data <- rbind(anc, int)
  data$group <- factor(data$group, levels = c("gained", "static", "lost"))
  data <- na.omit(data)
  
  median <- 
    c(median(data[data$group=="gained" & data$intersection=="anchor" & data$diff=="all",]$lfc), 
    median(data[data$group=="gained" & data$intersection=="internal" & data$diff=="all",]$lfc), 
    median(data[data$group=="lost" & data$intersection=="anchor" & data$diff=="all",]$lfc), 
    median(data[data$group=="lost" & data$intersection=="internal" & data$diff=="all",]$lfc), 
    median(data[data$group=="gained" & data$intersection=="anchor" & data$diff=="diff",]$lfc), 
    median(data[data$group=="gained" & data$intersection=="internal" & data$diff=="diff",]$lfc), 
    median(data[data$group=="lost" & data$intersection=="anchor" & data$diff=="diff",]$lfc), 
    median(data[data$group=="lost" & data$intersection=="internal" & data$diff=="diff",]$lfc)) |> 
    as.data.frame() |> 
    cbind(c("gained", "gained", "lost", "lost", "gained", "gained", "lost", "lost")) |> 
    cbind(c("anchor", "interior", "anchor", "interior", "anchor", "interior", "anchor", "interior")) |> 
    cbind(c("all", "all", "all", "all", "diff", "diff", "diff", "diff")) |> 
    cbind(rep(names[i], 4)) |> 
    setNames(c("lfc", "group", "intersection", "diff", "feature"))
  medians[[i]] <- median
  #Plot:
  # ggplot(data, aes(x = group, y = lfc, fill = intersection)) + 
  #   geom_split_violin(alpha = 0.85) + 
  #   geom_boxplot(width = 0.1, outlier.shape = NA, coef = 0, aes(alpha = intersection), position = position_dodge(width = 0.5)) + 
  #   scale_fill_manual(values = rep(c(pal[8], pal[6]), 3)) + 
  #   scale_color_manual(values = rep(c(pal[8], pal[6]), 3)) + 
  #   ylab("LFC") + ylim(c(-5,12)) + 
  #   ggtitle(paste0(names[i], " LFC")) + 
  #   geom_hline(yintercept = 0, lty=3) +
  #   theme_classic()
  #ggsave(last_plot(), filename = paste0("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig3/", names[i], "violin.pdf"))
  
  # ggplot(data, aes(x = group, y = lfc, alpha = diff, fill = intersection)) + 
  #   geom_boxplot(outlier.shape = NA) + 
  #   scale_fill_manual(values = rep(c(pal[8], pal[6]), 3)) + 
  #   scale_alpha_manual(values = c(0.8, 0.5)) + 
  #   ylab("LFC") + ylim(-10, 15) + 
  #   ggtitle(paste0(names[i], " LFC")) + 
  #   geom_hline(yintercept = 0, lty=3) +
  #   theme_classic()
  #ggsave(last_plot(), filename = paste0("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig3/", names[i], "boxplot.pdf"))
  
  ggplot(data[data$diff=="all",], aes(x = group, y = lfc, fill = intersection)) + 
    geom_boxplot(outlier.shape = NA) + 
    scale_fill_manual(values = rep(c(pal[8], pal[6]), 3)) + 
    ylab("LFC") + ylim(-9, 13) + 
    ggtitle(paste0(names[i], " LFC")) + 
    geom_hline(yintercept = 0, lty=3) +
    theme_classic() + theme(legend.position = "none")
  ggsave(last_plot(), filename = paste0("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig3/", names[i], "boxplotAll.pdf"))
  
  #Wilcox: 
  # stats <- 
  #   data.frame(
  #   (wilcox.test(data[data$group=="gained" & data$intersection=="anchor",]$lfc)$p.value), 
  #   (wilcox.test(data[data$group=="gained" & data$intersection=="internal",]$lfc)$p.value), 
  #   (wilcox.test(data[data$group=="static" & data$intersection=="anchor",]$lfc)$p.value), 
  #   (wilcox.test(data[data$group=="static" & data$intersection=="internal",]$lfc)$p.value), 
  #   (wilcox.test(data[data$group=="lost" & data$intersection=="anchor",]$lfc)$p.value), 
  #   (wilcox.test(data[data$group=="lost" & data$intersection=="internal",]$lfc)$p.value)
  # )
  # colnames(stats) <- c("gainedAnchor", "gainedInternal", "staticAnchor", "staticInternal", "lostAnchor", "lostInternal")
  # wilcox[[i]] <- stats
  # 
  # #KS test: 
  # stats <- data.frame(
  #   (ks.test(data[data$group=="gained" & data$intersection=="anchor",]$lfc, data[data$group=="static" & data$intersection=="anchor",]$lfc))$p.value,
  #   (ks.test(data[data$group=="gained" & data$intersection=="internal",]$lfc, data[data$group=="static" & data$intersection=="internal",]$lfc))$p.value,
  #   (ks.test(data[data$group=="lost" & data$intersection=="anchor",]$lfc, data[data$group=="static" & data$intersection=="anchor",]$lfc))$p.value,
  #   (ks.test(data[data$group=="lost" & data$intersection=="internal",]$lfc, data[data$group=="static" & data$intersection=="internal",]$lfc))$p.value
  # )
  
  #T.test:
  stats <- data.frame(
    (wilcox.test(data[data$group=="gained" & data$intersection=="anchor",]$lfc, data[data$group=="static" & data$intersection=="anchor",]$lfc))$p.value,
    (wilcox.test(data[data$group=="gained" & data$intersection=="internal",]$lfc, data[data$group=="static" & data$intersection=="internal",]$lfc))$p.value,
    (wilcox.test(data[data$group=="lost" & data$intersection=="anchor",]$lfc, data[data$group=="static" & data$intersection=="anchor",]$lfc))$p.value,
    (wilcox.test(data[data$group=="lost" & data$intersection=="internal",]$lfc, data[data$group=="static" & data$intersection=="internal",]$lfc))$p.value
  )
  colnames(stats) <- c("gainedAnc", "gainedInt", "lostAnc", "lostInt")
  wilcox[[i]] <- stats
}

wilcox <- do.call(rbind, wilcox)
rownames(wilcox) <- names
wilcox <- apply(wilcox, 1, p.adjust)

cutoffs <- c(1, 0.05, 0.01, 0.001, 0)

signif <- symnum(wilcox, corr = FALSE, na = FALSE, 
                 cutpoints = cutoffs, 
                 symbols = c("***", "**", "*", "ns"))


median <- do.call(rbind, medians)
median$feature <- factor(median$feature, levels = names)

ggplot(median[median$group=="gained" & median$diff=="all",], aes(x = feature, y = lfc, color = intersection)) + 
  geom_point(shape = 8) + ylim(c(-2,4)) + scale_color_manual(values = c(pal[8], pal[6])) + 
  theme_classic() + theme(legend.position = "none") + geom_hline(yintercept = 0, lty=3) + 
  ylab("Feature Fold-Change")
ggsave(last_plot(), file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig3/summaryGained.pdf")

ggplot(median[median$group=="lost" & median$diff=="all",], aes(x = feature, y = lfc, color = intersection, shape = feature)) + 
  geom_point() + ylim(c(-2,4)) + scale_color_manual(values = c(pal[8], pal[6])) + 
  theme_classic() + theme(legend.position = "none") + geom_hline(yintercept = 0, lty=3) + 
  scale_shape_manual(values = c(8, 8, 8, 8, 8, 16)) + ylab("Feature Fold-Change")
ggsave(last_plot(), file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig3/summaryLost.pdf")

ggplot(median[median$group=="gained" & median$diff=="diff",], aes(x = feature, y = lfc, color = intersection)) + 
  geom_point(shape = 8) + ylim(c(-3,5.1)) + scale_color_manual(values = c(pal[8], pal[6])) + 
  theme_classic() + theme(legend.position = "none") + geom_hline(yintercept = 0, lty=3) + 
  ylab("Feature Fold-Change")
ggsave(last_plot(), file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig3/summaryGainedDiff.pdf")

ggplot(median[median$group=="lost" & median$diff=="diff",], aes(x = feature, y = lfc, color = intersection, shape = feature)) + 
  geom_point() + ylim(c(-3,5.1)) + scale_color_manual(values = c(pal[8], pal[6])) + 
  theme_classic() + theme(legend.position = "none") + geom_hline(yintercept = 0, lty=3) + 
  scale_shape_manual(values = c(8, 8, 8, 8, 8, 16)) + ylab("Feature Fold-Change")
ggsave(last_plot(), file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig3/summaryLostDiff.pdf")





