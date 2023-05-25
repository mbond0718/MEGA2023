## SUPPLEMENTAL FIGURE 1: 
## Plotting the number of loops called in each time point, in the omega map, and total 

## Load libraries -----------------------------------------
library(ggplot2)
library(reshape2)

## Read in data -------------------------------------------
loopFiles <- 
  list.files("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/raw/hic/loops/", full.names = TRUE)
mergedLoops <- 
  readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/loopCounts.rds") |> 
  as.data.frame()
nMerged <- nrow(mergedLoops) |> 
  setNames("merged")

# Make list
loops <- lapply(loopFiles, read.table, header=T)
loops <- setNames(loops, c("0", "6", "72", "omega"))

# Get numbers of loops
nLoops <- lapply(loops, nrow) |> unlist()
nLoops <- c(nLoops, nMerged)
nLoops <- melt(nLoops)
nLoops$class <- rownames(nLoops)
nLoops$class <- factor(nLoops$class, levels = c(0, 6, 72, "omega", "merged"))

# Plot
ggplot(nLoops, aes(x = class, y = value)) + 
  geom_bar(stat = "identity") + 
  ggtitle("Number of Loops Called") + 
  theme_classic()
ggsave("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/sup1/loopNumBarplots.pdf")
