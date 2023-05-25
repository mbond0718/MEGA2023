## HIC VISUALIZATION: Figure 1A,1B,1C: 
## The purpose of this script is to generate Figure 1A-C, where hic data is visualized at: 
## (A) the compartment level
## (B) the tad level 
## (C) the loop level


## Load libraries -----------------------------------------
remotes::install_github("EricSDavis/mariner", ref = "dev", force = TRUE)
library(mariner)
library(plotgardener)
library(RColorBrewer)
library(purrr)

## Process ------------------------------------------------

# Read in data
hicPaths <- list.files("/Users/phanstiel11/Phanstiel Lab Dropbox/Marielle Bond/MEGA_shared/hic/hic_files/timepoints_normAfter/", full.names = TRUE)[c(1,3)]
comp <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/compartments/compartments.rds")
tad <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/tads/tadRanges.rds")
loops <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/uniqueLoops.rds")

# Format compartments: 
compartments <- list()
for (i in 1:3){
  compartments[[i]] <- comp[,c(1:3,i+3)]
}
#Assign column names
compartments <- lapply(compartments, setNames, nm = c("chr", "start", "end", "score"))

# Palette
pal <- brewer.pal(8, "YlGnBu")

# Create page
pdf("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig1/fig1ABC.pdf", width = 8.5, height = 11)
pageCreate(width = 8.5, height = 11, showGuides = F, xgrid=0.0, ygrid=0.0)

## Define params
#compartment params
pc <- pgParams(chrom = "7",
               chromstart = 30E6,
               chromend =50E6,
               assembly = "hg38",
               resolution = 100000,
               width = 2,
               x = 0.5,
               y = 0.5,
               buffer = 0.1,
               height = 2,
               zrange = c(0, 400),
               norm = "SCALE")

pcst <- pc
pcst$height <- 0.25
pcst$buffer <- 0.25

#tad params 
pt <- pgParams(chrom = "7",
               chromstart = 54990000-0.5E6,
               chromend = 55130000+1E6,
               assembly = "hg38",
               resolution = 25000,
               width = 2,
               x = 3,
               y = 0.5,
               buffer = 0.1,
               height = 2,
               zrange = c(0, 900),
               norm = "SCALE")

ptst <- pt
ptst$height <- 0.25
ptst$buffer <- 0.25
ptst$boxHeight <- 0.025

#loop params
pl <- pgParams(chrom = "7",
               chromstart = 54990000-330000,
               chromend = 55130000+220000,
               assembly = "hg38",
               resolution = 10000,
               width = 2,
               x = 5.5,
               y = 0.5,
               buffer = 0.1,
               height = 2,
               zrange = c(0, 275),
               norm = "SCALE")

plst <- pl
plst$chrom <- paste0('chr', plst$chrom)
plst$height <- 0.25
plst$buffer <- 0.25
plst$boxHeight <- 0.025

#Compartments: 
comp0data <- readHic(hicPaths[1], params = pc)
comp72data <- readHic(hicPaths[2], params = pc)
comp0 <- plotHicSquare(comp0data, params = pc, half = "top", matrix = "oe")
comp72 <- plotHicSquare(comp72data, params = pc, half = "bottom", matrix = "oe")
#Annotations: 
annoHeatmapLegend(
  plot = comp0, fontcolor = "black",
  x = 2.55, y = 0.5,
  width = 0.10, height = 0.75, fontsize = 6
)
#Eigenvectors: 
plotSignal(data = compartments[[1]], params = pcst, y=2.65, linecolor = c(pal[6], pal[8]), fill = c(pal[6], pal[8]), negData = T)
plotSignal(data = compartments[[3]], params = pcst, y=2.95, linecolor = c(pal[6], pal[8]), fill = c(pal[6], pal[8]), negData = T)
annoGenomeLabel(comp0, x=0.5, y=3.25, scale="bp", fontsize = 8)

#Tads
tad0data <- readHic(hicPaths[1], params = pt)
tad72data <- readHic(hicPaths[2], params = pt)
tad0 <- plotHicSquare(tad0data, params = pt, half = "top")
tad72 <- plotHicSquare(tad72data, params = pt, half = "bottom")
#Annotations: 
annoHeatmapLegend(
  plot = tad0, fontcolor = "black",
  x = 5.05, y = 0.5,
  width = 0.10, height = 0.75, fontsize = 6
)

#Tad ranges, may be replaced with insulation scores
pair0 <- plotPairs(tad[[1]], params = ptst, y=2.65, lwd = 1, spaceHeight = .5, 
                   fill = colorby("unique", palette = colorRampPalette(c("gray", pal[8]))))
pair72 <- plotPairs(tad[[3]], params = ptst, y=2.95, lwd = 1, spaceHeight = .5, 
                    fill = colorby("unique", palette = colorRampPalette(c("gray", pal[8]))))
annoGenomeLabel(tad0, x=3, y=3.25, scale="bp", fontsize = 8)

#Loops
loop0data <- readHic(hicPaths[1], params = pl)
loop72data <- readHic(hicPaths[2], params = pl)
loop0 <- plotHicSquare(loop0data, params = pl, half = "top")
loop72 <- plotHicSquare(loop72data, params = pl, half = "bottom")
#Annotations: 
annoHeatmapLegend(
  plot = loop0, fontcolor = "black",
  x = 7.55, y = 0.5,
  width = 0.10, height = 0.75, fontsize = 6
)
#Plot arches for loop, color by whether it's differential: 
#Plot
l0arc <- plotPairsArches(loops[[1]], params = plst, y = 2.65, x = 5.5, alpha = 1, fill = colorby("unique", palette = colorRampPalette(c("gray", pal[8]))))
l72arc <- plotPairsArches(loops[[3]], params = plst, y = 2.95, x = 5.5, alpha = 1, fill = colorby("unique", palette = colorRampPalette(c("gray", pal[8]))))
annoGenomeLabel(loop0, x=5.5, y=3.25, scale="bp", fontsize = 8)

#Add text labels: 
#Hic plots: 
plotText("0", x = pc$x+0.1, y = pc$y+0.1, fontsize = 8)
plotText("72", x = pc$x+pc$height-0.15, y = pc$y+pc$height-0.15, fontsize = 8)
plotText("0", x = pt$x+0.1, y = pt$y+0.1, fontsize = 8)
plotText("72", x = pt$x+pt$height-0.15, y = pt$y+pt$height-0.15, fontsize = 8)
plotText("0", x = pl$x+0.1, y = pl$y+0.1, fontsize = 8)
plotText("72", x = pl$x+pl$height-0.15, y = pl$y+pl$height-0.15, fontsize = 8)
#Signal tracks
plotText("0", x = pcst$x-0.15, y = 2.7, fontsize = 8)
plotText("72", x = pcst$x-0.15, y = 3.0, fontsize = 8)
plotText("0", x = ptst$x-0.15, y = 2.7, fontsize = 8)
plotText("72", x = ptst$x-0.15, y = 3.0, fontsize = 8)
plotText("0", x = plst$x-0.15, y = 2.7, fontsize = 8)
plotText("72", x = plst$x-0.15, y = 3.0, fontsize = 8)
#Resolution 
plotText("100 kb res", x = pc$x+pc$width, y = pc$y+pc$height+0.075, just = c("right", "center"), fontsize = 8)
plotText("25 kb res", x = pt$x+pt$width, y = pt$y+pt$height+0.075, just = c("right", "center"), fontsize = 8)
plotText("10 kb res", x = pl$x+pl$width, y = pl$y+pl$height+0.075, just = c("right", "center"), fontsize = 8)

dev.off()










