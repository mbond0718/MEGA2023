## Compare the p values from hic-dc and deseq

library(ggplot2)
library(mariner)
library(devtools)
install_github("EricSDavis/mariner")
library(cowplot)

## Read in data: 
data <- read.table("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/hicDC/diff180523.txt", header=T)
data <- na.omit(data)

padj <- 0.05
lfc <- log2(1.5)

deseq <- data[abs(data$lfc1)>lfc & data$padj1<0.05,] |> na.omit()
hicdc <- data[abs(data$lfc2)>lfc & data$padj2<0.05,] |> na.omit()

#Which loops are differential according to which differential loop caller: 
data$deseq <- ""
data$hicdc <- ""
data[data$loopname %in% deseq$loopname,]$deseq <- TRUE
data[data$loopname %in% hicdc$loopname,]$hicdc <- TRUE

## Assign loop identifiers: 
deseq$id <- paste0(deseq$coord1, ":", deseq$start1, "-", deseq$end1)
hicdc$id <- paste0(hicdc$coord1, ":", hicdc$start1, "-", hicdc$end1)

#Categorize: 
data$category <- ""
data[data$deseq==TRUE & data$hicdc=="",]$category <- "deseq"
data[data$deseq=="" & data$hicdc==TRUE,]$category <- "hicdc"
data[data$deseq==TRUE & data$hicdc==TRUE,]$category <- "shared"
data <- data[!data$category=="",]
table(data$category)
data$category <- factor(data$category, levels = c("deseq", "shared", "hicdc"))

# Make scatter plot comparing the p values in each dataset: 
scatter <- 
  ggplot(data, aes(x = -log10(padj1), y = -log10(padj2), color = category)) + 
  geom_point() + 
  scale_color_manual(values = c("#D6EAF1", "gray", "#A9D48A")) + 
  theme_classic() + 
  xlab("DESeq2 -log10(padj)") + ylab("HiC-DC+ -log10(padj)") + 
  geom_hline(yintercept = -log10(padj), lty = 3) + 
  geom_vline(xintercept = -log10(padj), lty = 3) + 
  theme(legend.position = "none") + 
  ylim(0,35) + xlim(0,15)
  

# Make density plot showing the density of p values for each of the classes: 
# Deseq: 
deseqDensity <- 
  ggplot(data[data$category=="deseq" | data$category=="shared",], aes(x = -log10(padj1), fill = category)) +
  geom_density(alpha = 0.85) + 
  scale_fill_manual(values = c("#D6EAF1", "gray")) + 
  theme_classic() + 
  xlab("") + 
  geom_vline(xintercept = -log10(median(data[data$category=="deseq",]$padj1)), color = "#D6EAF1") + 
  geom_vline(xintercept = -log10(median(data[data$category=="shared",]$padj1)), color = "gray") + 
  theme(legend.position = "none") + 
  scale_y_reverse() + xlim(0,15)
#statistical test: 
wilcox.test(data[data$category=="deseq",]$padj1, data[data$category=="shared",]$padj1) #difference is significant


# HicDC: 
hicdcDensity <- 
  ggplot(data[data$category=="hicdc" | data$category=="shared",], aes(x = -log10(padj2), fill = category)) +
  geom_density(alpha = 0.85) + 
  scale_fill_manual(values = c("gray", "#A9D48A")) + 
  theme_classic() + 
  xlab("") + 
  geom_vline(xintercept = -log10(median(data[data$category=="hicdc",]$padj2)), color = "#A9D48A") + 
  geom_vline(xintercept = -log10(median(data[data$category=="shared",]$padj2)), color = "gray") + 
  scale_y_reverse() + 
  coord_flip() + 
  theme(legend.position = "none") + xlim(0,35)
#statistical test: 
wilcox.test(data[data$category=="hicdc",]$padj2, data[data$category=="shared",]$padj2) #difference is significant


#Plot together" 
plot_grid(
  hicdcDensity, scatter, NULL, deseqDensity, 
  align = "h", 
  rel_heights = c(3,1), 
  rel_widths = c(1,3)
)
ggsave("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/reviewFigs/compareDiffLoops.pdf", last_plot(), 
       width = 10, height = 10)


