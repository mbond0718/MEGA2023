## VISUALIZATION OF EXAMPLE GAINED AND LOST LOOP WITH REGULATORY DATA
## This script generates plotgardener visualizations of a gained and lost loop for Figure 3

## Load libraries -----------------------------------------
library(plotgardener)
library(purrr)
library(RColorBrewer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)
library(plyranges)

## Read in data --------------------------------------------
#Files
hicPaths <- list.files("/Users/phanstiel11/Phanstiel Lab Dropbox/Marielle Bond/MEGA_shared/hic/hic_files/timepoints_normAfter/", full.names = TRUE)[c(1,3)]
rnaPaths <- list.files("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/rna/", "bw", full.names = TRUE)[c(1,8)]
atacPaths <- list.files("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/atac/", "bw", full.names = TRUE)[c(1,3)]
h3k27Paths <- list.files("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/cnr/", "H3K27", full.names = TRUE)[c(4,6)]
ctcfPaths <- list.files("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/cnr/", "CTCF", full.names = TRUE)[c(4,6)]
junPaths <- list.files("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/cnr/", "Jun", full.names = TRUE)[c(4,6)]
radPaths <- list.files("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/cnr/", "Rad21", full.names = TRUE)[c(5,6)]

#Palette
pal <- brewer.pal(8, "YlGnBu")

#Read in loci and subset for specifically gained and lost loop of interest: 
loops <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/diffLoops.rds")
loops <- loops[loops$class=="gained" | loops$class=="lost",]
loops$seqnames1 <- unlist(strsplit(as.character(loops$seqnames1), "chr"))[seq(2, (2*nrow(loops)), by=2)]
loops$seqnames2 <- unlist(strsplit(as.character(loops$seqnames2), "chr"))[seq(2, (2*nrow(loops)), by=2)]

#Loop bedpe
loopanno <- loops
loopanno <- loopanno[rownames(loopanno) %in% c("loop24603", "loop21914"),]
loopanno <- loopanno[2:1,]

#loopanno <- loopanno[loopanno$class=="lost" & loopanno$seqnames1==11 & loopanno$start1>11.7E6 & loopanno$end2<12.7E6,]

#Continue formatting: 
loops <- loops[,c(1,2,5)]
colnames(loops) <- c("chr", "start", "stop")
loops <- loops[rownames(loops) %in% c("loop24603", "loop21914"),]
loops <- loops[2:1,]

loops <- loops[rownames(loops) %in% rownames(loopanno),]

## Calculate plotting regions
res <- 10e3
regions <- makeGRangesFromDataFrame(loops)

#Set x axis for panels a and b
x <- c(0.5, 4.5, 0.5, 0.5)
z <- c(300, 350, 300, 300)
width <- c(0.85E6, 0.6E6, 0.6E6, 0.6E6)
## Panels a and b, showing representative correlative gained and lost loop 
pdf("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/reviewFigs/fig3ab.pdf", width = 8.5, height = 11)
#create page 
pageCreate(width = 8.5, height = 11, showGuides = F, xgrid=0.0, ygrid=0.0)
#plot gained and lost loops for panels a and b 
for (i in 1:nrow(loops)){
  #Resize region
  reg <- 
    regions |>
    resize(width = 1, fix = 'center') |>
    shift(shift = res/2) |>
    resize(width = width[i], fix = 'center')
  
  ## Define params
  p <- pgParams(chrom = as.character(seqnames(reg)[i]),
                chromstart = start(reg)[i],
                chromend = end(reg)[i],
                assembly = "hg38",
                resolution = res,
                width = 3,
                x = x[[i]],
                y = 0.5,
                buffer = 0.05,
                height = 1,
                zrange = c(0, z[i]),
                norm = "SCALE")
  
  pst <- p
  pst$chrom <- paste0('chr', p$chrom)
  pst$height <- 0.2
  pst$buffer <- 0.01
  
  hic_ypos <- c(0.5, p$y + p$height + p$buffer)
  atac_ypos <- c(2.55, 2.55 + pst$height + pst$buffer)
  h3k_ypos <- c(2.97, 2.97 + pst$height + pst$buffer)
  jun_ypos <- c(3.39, 3.39 + pst$height + pst$buffer)
  ctcf_ypos <- c(3.81, 3.81 + pst$height + pst$buffer)
  rad_ypos <- c(4.23, 4.23 + pst$height + pst$buffer)
  rna_ypos <- c(4.65, 4.65 + pst$height + pst$buffer)
  tp <- c("0", "72")
  
  ## Plot hic
  hic <- pmap(list(hicPaths, hic_ypos), \(x, y) plotHicRectangle(data = x, params = p, y = y))
  
  ## ATAC 
  atacData <- lapply(atacPaths, \(x) readBigwig(file = x, params = pst)[c(1:3, 6)])
  atacRange <- c(0, lapply(atacData, `[[`, 4) |> unlist() |> max())
  ## Plot ATAC
  atac <- pmap(list(atacData, atac_ypos), \(x, y) plotSignal(data = x, params = pst, y = y, range = atacRange, linecolor = "#8ed368", fill = "#8ed368"))
  
  ### H3K27ac
  h3kData <- lapply(h3k27Paths, \(x) readBigwig(file = x, params = pst)[c(1:3, 6)])
  h3kRange <- c(0, lapply(h3kData, `[[`, 4) |> unlist() |> max())
  ## Plot H3K27ac
  h3k27 <- pmap(list(h3kData, h3k_ypos), \(x, y) plotSignal(data = x, params = pst, y = y, range = h3kRange, linecolor = pal[4], fill = pal[4]))
  
  ### Jun
  junData <- lapply(junPaths, \(x) readBigwig(file = x, params = pst)[c(1:3, 6)])
  junRange <- c(0, lapply(junData, `[[`, 4) |> unlist() |> max())
  ## Plot Jun
  jun <- pmap(list(junData, jun_ypos), \(x, y) plotSignal(data = x, params = pst, y = y, range = junRange, linecolor = pal[5], fill = pal[5]))
  
  ### CTCF
  ctcfData <- lapply(ctcfPaths, \(x) readBigwig(file = x, params = pst)[c(1:3, 6)])
  ctcfRange <- c(0, lapply(ctcfData, `[[`, 4) |> unlist() |> max())
  ## Plot CTCF
  ctcf <- pmap(list(ctcfData, ctcf_ypos), \(x, y) plotSignal(data = x, params = pst, y = y, range = ctcfRange, linecolor = pal[6], fill = pal[6]))
  
  ### RAD21
  radData <- lapply(radPaths, \(x) readBigwig(file = x, params = pst)[c(1:3, 6)])
  radRange <- c(0, lapply(radData, `[[`, 4) |> unlist() |> max())
  ## Plot CTCF
  rad <- pmap(list(radData, rad_ypos), \(x, y) plotSignal(data = x, params = pst, y = y, range = radRange, linecolor = pal[7], fill = pal[7]))
  
  ## RNA
  rnaData <- lapply(rnaPaths, \(x) readBigwig(file = x, params = pst)[c(1:3, 6)])
  rnaRange <- c(0, lapply(rnaData, `[[`, 4) |> unlist() |> max())
  ## Plot RNA
  rna <- pmap(list(rnaData, rna_ypos), \(x, y) plotSignal(data = x, params = pst, y = y, range = rnaRange, linecolor = pal[8], fill = pal[8]))
  
  ## Plot Genes
  genes <- plotGenes(params = pst, height = 0.25, y=5.10)
  
  ## Plot coords
  genCoord <- annoGenomeLabel(pst, x=x[[i]], y=5.4, scale="bp")
  
  #Highlight anchors
  anc1 <- annoHighlight(plot = hic[[1]], 
                        chrom = loops$chr[i], 
                        chromstart = loopanno[i,]$start1-2500, 
                        chromend = loopanno[i,]$end1+2500, 
                        y=2.55, height = 2.85, alpha = 0.25)
  
  anc2 <- annoHighlight(plot = hic[[1]], 
                        chrom = loops$chr[i], 
                        chromstart = loopanno[i,]$start2-2500, 
                        chromend = loopanno[i,]$end2+2500, 
                        y=2.55, height = 2.85, alpha = 0.25)
  
  #Annotate loop 
  anno1 <- annoPixels(hic[[1]], data = loopanno[i,], type = "arrow", shift = 3)
  anno2 <- annoPixels(hic[[2]], data = loopanno[i,], type = "arrow", shift = 3)
  
  annoHeatmapLegend(
    plot = hic[[1]],
    x = 3.6, y = 0.5, width = 0.1, height = 0.75,
    just = c("left", "top"), default.units = "inches"
  )
  
  ## Plot labels
  hictp <- plotText(label = tp, x = p$x + p$buffer, y = hic_ypos + p$buffer, just = c('left', 'top'))
  atactp <- plotText(label = tp, fontsize = 8, x = x[[i]]+0.05, y = atac_ypos + 0.15, just = c('right', 'bottom'))
  h3ktp <- plotText(label = tp, fontsize = 8, x = x[[i]]+0.05, y = h3k_ypos + 0.15, just = c('right', 'bottom'))
  juntp <- plotText(label = tp, fontsize = 8, x = x[[i]]+0.05, y = jun_ypos + 0.15, just = c('right', 'bottom'))
  ctcftp <- plotText(label = tp, fontsize = 8, x = x[[i]]+0.05, y = ctcf_ypos + 0.15, just = c('right', 'bottom'))
  radtp <- plotText(label = tp, fontsize = 8, x = x[[i]]+0.05, y = rad_ypos + 0.15, just = c('right', 'bottom'))
  rnatp <- plotText(label = tp, fontsize = 8, x = x[[i]]+0.05, y = rna_ypos + 0.15, just = c('right', 'bottom'))
  
  ## Plot labels
  ataclab <- plotText("ATAC", x=x[[i]]-0.2, y=mean(atac_ypos+0.2), fontface = "bold", fontsize = 8, rot = 90)
  h3klab <- plotText("H3K27ac", x=x[[i]]-0.2, y=mean(h3k_ypos+0.2), fontface = "bold", fontsize = 8, rot = 90)
  junlab <- plotText("Jun", x=x[[i]]-0.2, y=mean(jun_ypos+0.2), fontface = "bold", fontsize = 8, rot = 90)
  ctcflab <- plotText("CTCF", x=x[[i]]-0.2, y=mean(ctcf_ypos+0.2), fontface = "bold", fontsize = 8, rot = 90)
  radlab <- plotText("Rad21", x=x[[i]]-0.2, y=mean(rad_ypos+0.2), fontface = "bold", fontsize = 8, rot = 90)
  rnalab <- plotText("RNA", x=x[[i]]-0.2, y=mean(rna_ypos+0.2), fontface = "bold", fontsize = 8, rot = 90)
}
dev.off()
