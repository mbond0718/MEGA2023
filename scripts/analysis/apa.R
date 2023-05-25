## AGGREGATE PEAK ANALYSIS: 
## This script is to be run after diffLoops.R generate aggregate peak plots of differential loop calls

## Load libraries -----------------------------------------
library(DESeq2)
library(GenomicRanges)
library(RColorBrewer)
library(RNASeqPower)

## Process ------------------------------------------------


##MEGA APA 
## Assemble all, WT, and FS loops into a list

allLoops <- 
  readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/loopCounts.rds") |> 
  as.data.frame()
allLoops <- 
  allLoops[,c(1:3,6:8)]
diffLoops <- 
  readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/diffLoops.rds") 
gainedLoops <- 
  diffLoops[diffLoops$class=="gained",]
gainedLoops <- 
  gainedLoops[,1:6]
lostLoops <- 
  diffLoops[diffLoops$class=="lost",]
lostLoops <- 
  lostLoops[,1:6]

hicPaths <- list.files("/Users/phanstiel11/Phanstiel Lab Dropbox/Marielle Bond/MEGA_shared/hic/hic_files/timepoints_normAfter/", full.names = TRUE)[c(1,3)]

loopList <- 
  list(allLoops = allLoops |> as.data.frame(),
       gainedLoops = gainedLoops,
       lostLoops = lostLoops)
loopList <- 
  lapply(loopList, as_ginteractions)


## Define resolution and buffer (pixels from center)
res <- 10e3
buffer <- 10

## Filter out short loop interactions
filteredLoops <- 
  lapply(X = loopList,
         FUN = filterBedpe,
         res = res,
         buffer = buffer) |>
  `names<-`(value = names(loopList))

lapply(filteredLoops, summary)

## Hi-C file paths
hic0Path <- hicPaths[1]
hic72Path <- hicPaths[2]

## Calculate APA matrices for loops from K562
loopApa0Hic <-
  lapply(X = filteredLoops,
         FUN = calcApa,
         hic = hic0Path,
         norm = "SCALE",
         buffer = buffer)

## Calculate APA matrices for loops from MK 
loopApa72Hic <-
  lapply(X = filteredLoops,
         FUN = calcApa,
         hic = hic72Path,
         norm = "SCALE",
         buffer = buffer)

## Get the number of loops for each condition
nLoops <- lapply(filteredLoops, length)

## Divide each matrix by nLoops
loopApa0Hic <- Map("/", loopApa0Hic, nLoops)
loopApa72Hic <- Map("/", loopApa72Hic, nLoops)


## Initiate plotgardener page
pdf("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig1/fig1f.pdf", width = 4.25, height = 3)
pageCreate(width = 4.25, height = 3, showGuides = F, xgrid=0.0, ygrid=0.0)

## Define shared parameters
p <- pgParams(x = 0.5,
              y = 0.5,
              width = 0.8,
              height = 0.8,
              space = 0.075,
              zrange = c(0, 125),
              palette = colorRampPalette(brewer.pal(9, 'YlGnBu')))

z <- list(c(0, 85), c(0, 175), c(0, 175))

## Define grid of coordinate positions
#xpos <- c(p$x, p$x + p$width + p$space, p$x + (p$width + p$space)*2)
xpos <- c(0.5, 1.5, 2.5)
ypos <- c(p$y, p$y + p$height + p$space, p$y + (p$height + p$space)*2)

s <- list(loopApa0Hic$allLoops, loopApa72Hic$allLoops)
g <- list(loopApa0Hic$gainedLoops, loopApa72Hic$gainedLoops)
l <- list(loopApa0Hic$lostLoops, loopApa72Hic$lostLoops)

## Plot column of static loops
static_plots <- 
  pmap(list(s, xpos[1], ypos[1:2], z[1]), \(a, x, y, z) {
    plotApa(params = p, apa = a, x = x, y = y, zrange = z)
  })
gained_plots <- 
  pmap(list(g, xpos[2], ypos[1:2], z[2]), \(a, x, y, z) {
    plotApa(params = p, apa = a, x = x, y = y, zrange = z)
  })
lost_plots <- 
  pmap(list(l, xpos[3], ypos[1:2], z[3]), \(a, x, y, z) {
    plotApa(params = p, apa = a, x = x, y = y, zrange = z)
  })

## Add legend
annoHeatmapLegend(plot = static_plots[[1]],
                  x = xpos[1] + 0.86,
                  y = ypos[1],
                  width = p$space,
                  height = p$height*0.75,
                  fontcolor = 'black')
annoHeatmapLegend(plot = gained_plots[[1]],
                  x = xpos[2] + 0.86,
                  y = ypos[1],
                  width = p$space,
                  height = p$height*0.75,
                  fontcolor = 'black')
annoHeatmapLegend(plot = lost_plots[[1]],
                  x = xpos[3] + 0.86,
                  y = ypos[1],
                  width = p$space,
                  height = p$height*0.75,
                  fontcolor = 'black')

## Add text labels
plotText(label = c("All loops", "Gained loops", "Lost loops"),
         x = xpos + p$width / 2,
         y = ypos[1] - p$space*3,
         just = c('center', 'bottom'))
plotText(label = paste0("n=", nLoops), 
         x = xpos + p$width / 2,
         y = ypos[1]-p$space, 
         fontsize = 10,
         just = c('center', 'bottom'))
plotText(label = c("K562", "MK"),
         x = xpos[1] - p$space,
         y = ypos[1:2] + p$height / 2,
         rot = 90,
         just = c('center', 'bottom'))

dev.off()
