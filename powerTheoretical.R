## THEORETICAL POWER ANALYSIS 
## Using the dispersion value from diffLoops.R, calculate theoretical power across various sequencing depths 

## Load libraries -----------------------------------------
library(DESeq2)
library(GenomicRanges)
library(RColorBrewer)
library(RNASeqPower)
library(ggplot2)

## Process ------------------------------------------------

# Palette
pal <- brewer.pal(8, "YlGnBu")

#Model dispersion over various sequencing depths: 
powerRep <- list()
power <- list()
for (x in 2:4){
  n <- x
  for (i in 1:251){
    z <- c(0:250)[i]
    powerRep[[i]] <- rnapower(alpha=0.5/33914, cv=sqrt(0.001984929), depth=z, effect=2, n=n)
  }
  
  powerData <- data.frame(do.call(rbind, powerRep))
  colnames(powerData) <- "Power"
  powerData$Depth <- c(0:250)
  powerData$Rep <- rep(n, nrow(powerData))
  
  power[[x]] <- powerData
  
}
power <- data.frame(do.call(rbind, power))
power$Rep <- as.character(power$Rep)

# Plot power
pdf("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig1/fig1d.pdf", width = 8.5, height=11)
pageCreate(width = 8.5, height = 11, showGuides = F, xgrid=0.0, ygrid=0.0)
theoreticalPower <- 
  ggplot(power, aes(x=Depth, y=Power, group=Rep)) + 
  geom_line(aes(color=Rep), size=1) + 
  coord_cartesian(xlim = c(0, 100)) + 
  scale_color_manual(values= pal[c(6:8)]) + 
  theme_classic() + theme(legend.position="none") + 
  theme(panel.grid.major.y = element_line(color = "gray",
                                          size = 0.25,
                                          linetype = 2)) + 
  geom_vline(xintercept = 38, color = "firebrick", size = 0.3) + 
  xlab("Sequencing Depth (CPM)")
plotGG(plot = theoreticalPower, x = 0, y = 3.5, width = 2.5, height = 2)
dev.off()
