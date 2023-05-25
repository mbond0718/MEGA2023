## SUBSAMPLING POWER PLOT
## Plot the number of differential loops at multiple replicate combinations and sub-sampled sequencing depths 

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
data <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/subsamplePowerData.rds")

## PLOT: 
data$depth <- rep(depths, 3) |> as.factor()
data$reps <- data$reps |> as.factor()

practicalPowerPlot <- ggplot(data, aes(x=depth, y=diffTotal, color = reps, group = reps)) + 
  geom_point() + geom_line(lwd = 1) + scale_color_manual(values= pal[c(6:8)]) + 
  theme_classic() + theme(legend.position="none") + 
  theme(panel.grid.major.y = element_line(color = "gray",
                                          size = 0.25,
                                          linetype = 2)) + 
  xlab("Sequencing depth per replicate (Millions)") + 
  ylab("Total Differential Loops")


## Place power plot on a plotgardener page in the appropriate size: 
pdf("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig1/fig1e.pdf", width = 8.5, height=11)
pageCreate(width = 8.5, height = 11, showGuides = F, xgrid=0.0, ygrid=0.0)
plotGG(plot = practicalPowerPlot, x = 0, y = 3.5, width = 2.5, height = 2)
dev.off()
