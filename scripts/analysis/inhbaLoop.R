## VISUALIZATION OF EXAMPLE GAINED LOOP WITH GAINED GENE 
## This script generates plotgardener visualizations of the INHBA loop for Figure 2

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


#Palette
pal <- brewer.pal(8, "YlGnBu")

#Read in loci and subset for specifically gained and lost loop of interest: 
loops <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/diffLoops.rds")
loops <- loops[loops$class=="gained" | loops$class=="lost",]
loops$seqnames1 <- unlist(strsplit(as.character(loops$seqnames1), "chr"))[seq(2, (2*nrow(loops)), by=2)]
loops$seqnames2 <- unlist(strsplit(as.character(loops$seqnames2), "chr"))[seq(2, (2*nrow(loops)), by=2)]

#Loop bedpe
loopanno <- loops
loopanno <- loopanno[rownames(loopanno) %in% c("loop16018"),]

#Continue formatting: 
loops <- loops[,c(1,2,5)]
colnames(loops) <- c("chr", "start", "stop")
loops <- loops[rownames(loops) %in% c("loop16018"),]

## Calculate plotting regions
res <- 10e3
regions <- makeGRangesFromDataFrame(loops)


pdf("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/reviewFigs/inhba2D.pdf", width = 8.5, height = 11)
#create page 
pageCreate(width = 8.5, height = 11, showGuides = F, xgrid=0.0, ygrid=0.0)
#plot gained and lost loops for panels a and b 
  #Resize region
  reg <- 
    regions |>
    resize(width = 1, fix = 'center') |>
    shift(shift = res/2) |>
    resize(width = 0.5E6, fix = 'center')
  
  ## Define params
  p <- pgParams(chrom = as.character(seqnames(reg)[1]),
                chromstart = start(reg)[1],
                chromend = end(reg)[1],
                assembly = "hg38",
                resolution = res,
                width = 1.5,
                x = 0.5,
                y = 0.5,
                buffer = 0.05,
                height = 0.75,
                zrange = c(0, 300),
                norm = "SCALE")
  
  pst <- p
  pst$chrom <- paste0('chr', p$chrom)
  pst$height <- 0.2
  pst$buffer <- 0.01
  
  hic_ypos <- c(0.5, p$y + p$height + p$buffer)
  rna_ypos <- c(2.15, 2.15 + pst$height + pst$buffer)
  atac_ypos <- c(2.58, 2.58 + pst$height + pst$buffer)
  tp <- c("0", "72")
  
  ## Plot hic
  hic <- pmap(list(hicPaths, hic_ypos), \(x, y) plotHicRectangle(data = x, params = p, y = y))
  
  ## RNA
  rnaData <- lapply(rnaPaths, \(x) readBigwig(file = x, params = pst)[c(1:3, 6)])
  rnaRange <- c(0, lapply(rnaData, `[[`, 4) |> unlist() |> max())
  ## Plot RNA
  rna <- pmap(list(rnaData, rna_ypos), \(x, y) plotSignal(data = x, params = pst, y = y, range = rnaRange, linecolor = pal[8], fill = pal[8]))
  
  ## ATAC 
  atacData <- lapply(atacPaths, \(x) readBigwig(file = x, params = pst)[c(1:3, 6)])
  atacRange <- c(0, lapply(atacData, `[[`, 4) |> unlist() |> max())
  ## Plot ATAC
  atac <- pmap(list(atacData, atac_ypos), \(x, y) plotSignal(data = x, params = pst, y = y, range = atacRange, linecolor = "#8ed368", fill = "#8ed368"))
  
  ## Plot Genes
  genes <- plotGenes(params = pst, height = 0.25, y=3.1)
  
  ## Plot coords
  genCoord <- annoGenomeLabel(pst, x=0.5, y=3.37, scale="bp", fontsize = 7)
  
  #Highlight anchors
  anc1 <- annoHighlight(plot = hic[[1]], 
                        chrom = loops$chr, 
                        chromstart = loopanno$start1-2500, 
                        chromend = loopanno$end1+2500, 
                        y=2.05, height = 0.95, alpha = 0.25)
  
  anc2 <- annoHighlight(plot = hic[[1]], 
                        chrom = loops$chr[i], 
                        chromstart = loopanno[i,]$start2-2500, 
                        chromend = loopanno[i,]$end2+2500, 
                        y=2.05, height = 0.95, alpha = 0.25)
  
  #Annotate loop 
  anno1 <- annoPixels(hic[[1]], data = loopanno[i,], type = "circle", shift = 3)
  anno2 <- annoPixels(hic[[2]], data = loopanno[i,], type = "circle", shift = 3)
  
  ## Plot labels
  hictp <- plotText(label = tp, x = p$x + p$buffer, y = hic_ypos + p$buffer, just = c('left', 'top'))
  rnatp <- plotText(label = tp, fontsize = 7, x = p$x+0.05, y = rna_ypos + 0.15, just = c('right', 'bottom'))
  atactp <- plotText(label = tp, fontsize = 7, x = p$x+0.05, y = atac_ypos + 0.15, just = c('right', 'bottom'))
  
  ## Plot labels
  rnalab <- plotText("RNA", x=p$x-0.2, y=mean(rna_ypos+0.2), fontface = "bold", fontsize = 8, rot = 90)
  ataclab <- plotText("ATAC", x=p$x-0.2, y=mean(atac_ypos+0.2), fontface = "bold", fontsize = 8, rot = 90)
  
  ## Annotate heatmap legend
  annoHeatmapLegend(
    plot = hic[[1]],
    x = 2.1, y = 0.5, width = 0.1, height = 0.75,
    just = c("left", "top"), default.units = "inches"
  )
  
dev.off()
