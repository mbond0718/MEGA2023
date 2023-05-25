#Loops
loops <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/diffLoops.rds")
loops$loop <- rownames(loops)
gainedLoops <- loops[loops$class=="gained",]
lostLoops <- loops[loops$class=="lost",]
staticLoops <- loops[loops$class=="static",]

#RNA
rna <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/rna/rnaLFC.rds")

## Convert to GRanges: 
#Loop anchors
lostAnc <- as_ginteractions(lostLoops)

#Loop internal
lostIn <- GRanges(seqnames = seqnames(anchors(lostAnc)$first), 
                  ranges = IRanges(start = end(anchors(lostAnc)$first), end = start(anchors(lostAnc)$second)), 
                  loop = lostAnc$loop)

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
staticRNA <- rnaGR[rnaGR$differential==FALSE,]

subsetByOverlaps(lostIn, rnaGR)
olapAnc <- subsetByOverlaps(lostAnc, rnaGR)

loops <- lostLoops[!rownames(lostLoops) %in% olapAnc$loop & rownames(lostLoops) %in% subsetByOverlaps(lostIn, rnaGR)$loop,]
loopLFC <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/loopLFC.rds")

loops <- merge(loops, loopLFC, by = 0)
loops <- loops[order(loops$lfc),]
survey <- loops$Row.names

## Survey thru 
allLoops <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/diffLoops.rds")
loops <- allLoops[rownames(allLoops) %in% survey,]
loops <- loops[order(match(rownames(loops), survey)),]

loops$seqnames1 <- unlist(strsplit(as.character(loops$seqnames1), "chr"))[seq(2, (2*nrow(loops)), by=2)]
loops$seqnames2 <- unlist(strsplit(as.character(loops$seqnames2), "chr"))[seq(2, (2*nrow(loops)), by=2)]

loopanno <- loops

loops <- loops[,c(1,2,5)]
colnames(loops) <- c("chr", "start", "stop")

res <- 10e3
regions <- makeGRangesFromDataFrame(loops)

pdf("~/Desktop/MEGA_lostLoopSurvey_V2.pdf", width = 8.5, height = 11)
#plot gained and lost loops for panels a and b 
for (i in 1:100){
  #Resize region
  print(i)
  reg <- 
    regions |>
    resize(width = 1, fix = 'center') |>
    shift(shift = res/2) |>
    resize(width = (loops$stop[i]-loops$start[i])*2, fix = 'center')
  
  #create page 
  pageCreate(width = 8.5, height = 11, showGuides = F, xgrid=0.0, ygrid=0.0)
  
  ## Define params
  p <- pgParams(chrom = as.character(seqnames(reg)[i]),
                chromstart = start(reg)[i],
                chromend = end(reg)[i],
                assembly = "hg38",
                resolution = res,
                width = 3,
                x = 0.5,
                y = 0.5,
                buffer = 0.05,
                height = 1,
                zrange = c(0, 300),
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
  genCoord <- annoGenomeLabel(pst, x=0.5, y=5.4, scale="bp")
  
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
  
  ## Plot labels
  hictp <- plotText(label = tp, x = p$x + p$buffer, y = hic_ypos + p$buffer, just = c('left', 'top'))
  atactp <- plotText(label = tp, fontsize = 8, x = p$x+0.05, y = atac_ypos + 0.15, just = c('right', 'bottom'))
  h3ktp <- plotText(label = tp, fontsize = 8, x = p$x+0.05, y = h3k_ypos + 0.15, just = c('right', 'bottom'))
  juntp <- plotText(label = tp, fontsize = 8, x = p$x+0.05, y = jun_ypos + 0.15, just = c('right', 'bottom'))
  ctcftp <- plotText(label = tp, fontsize = 8, x = p$x+0.05, y = ctcf_ypos + 0.15, just = c('right', 'bottom'))
  radtp <- plotText(label = tp, fontsize = 8, x = p$x+0.05, y = rad_ypos + 0.15, just = c('right', 'bottom'))
  rnatp <- plotText(label = tp, fontsize = 8, x = p$x+0.05, y = rna_ypos + 0.15, just = c('right', 'bottom'))
  
  ## Plot labels
  ataclab <- plotText("ATAC", x=p$x-0.2, y=mean(atac_ypos+0.2), fontface = "bold", fontsize = 8, rot = 90)
  h3klab <- plotText("H3K27ac", x=p$x-0.2, y=mean(h3k_ypos+0.2), fontface = "bold", fontsize = 8, rot = 90)
  junlab <- plotText("Jun", x=p$x-0.2, y=mean(jun_ypos+0.2), fontface = "bold", fontsize = 8, rot = 90)
  ctcflab <- plotText("CTCF", x=p$x-0.2, y=mean(ctcf_ypos+0.2), fontface = "bold", fontsize = 8, rot = 90)
  radlab <- plotText("Rad21", x=p$x-0.2, y=mean(rad_ypos+0.2), fontface = "bold", fontsize = 8, rot = 90)
  rnalab <- plotText("RNA", x=p$x-0.2, y=mean(rna_ypos+0.2), fontface = "bold", fontsize = 8, rot = 90)
}
dev.off()




