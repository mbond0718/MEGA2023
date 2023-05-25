## Compare MEGA loop calls to 4DN loop calls

library(mariner)
library(VennDiagram)
library(venneuler)

megaK562 <- read.table("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/raw/hic/loops/0_loops.txt", header=T) |> 
  as_ginteractions() #all loops

megaK562 <- read.table("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/subsample100M_loops.txt", header=T) |> 
  as_ginteractions() #loops from subsampled data
GenomeInfoDb::seqlevelsStyle(megaK562) <- 'UCSC'

dekkerK562 <- read.table("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/4dn/buffer_dpn_4h.bedpe", header=F) |> 
  as_ginteractions()

overlap <- subsetByOverlaps(dekkerK562, megaK562, maxgap = 10E4) |> length() #overlap loop lists with 10kb wiggle room

##Venneuler
plot(venneuler(c(A = length(megaK562) - overlap,
                 B = length(dekkerK562) - overlap,
                 "A&B" = overlap)))

##Venn.Diagram
venn.diagram(list("MEGA" = 1:10204, "Dekker" = c(1:858, 10205:(10205+53))), fill = c("lightblue", "green"), 
             alpha = c(0.5, 0.5), lwd =0, filename = "~/Desktop/vennTest.png") #all loops 

venn.diagram(list("Dekker" = c(1:858, 10205:(10205+53)), "MEGA" = 1:10204), fill = c("lightblue", "green"), 
             alpha = c(0.5, 0.5), lwd =0, filename = "~/Desktop/vennTest.png", width = 5000, height = 5000) #all loops, order reversed


venn.diagram(list("Dekker" = c(1:898, 4392:(4392+13)), "MEGA" = 1:4392), fill = c("lightblue", "green"), 
             alpha = c(0.5, 0.5), lwd =0, filename = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/reviewFigs/dpnVenn.png", width = 5000, height = 5000) #all loops, order reversed



# Chart
venn.diagram(
  x = list(set1, set2, set3),
  category.names = c("Set 1" , "Set 2 " , "Set 3"),
  filename = '#14_venn_diagramm.png',
  output=TRUE
)

venn.diagram(
  x = list()
)



## Visualize loops: 
loops <- read.table("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/raw/hic/loops/0_loops.txt", header=T)
loops <- loops[order(loops$APScoreAvg, decreasing = T),] #order loops by strongest for visualization 
loops <- as_ginteractions(loops)
GenomeInfoDb::seqlevelsStyle(loops) <- 'ENSEMBL'

## Alternatively--loops ID'd in the Dekker samples: 
loops <- subsetByOverlaps(megaK562, dekkerK562, maxgap = 10E5)
loops <- loops[loops$APScoreAvg>10,]
GenomeInfoDb::seqlevelsStyle(loops) <- 'ENSEMBL'
loopsAlt <- loops
GenomeInfoDb::seqlevelsStyle(loopsAlt) <- 'UCSC'

#to annotate: 
anno <- read.table("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/raw/hic/loops/0_loops.txt", header=T)
anno <- as_ginteractions(anno)
GenomeInfoDb::seqlevelsStyle(anno) <- 'ENSEMBL'
annoAlt <- anno
GenomeInfoDb::seqlevelsStyle(annoAlt) <- 'UCSC'

## Read in hic paths: 
hic <- list("/Users/phanstiel11//Phanstiel Lab Dropbox/Marielle Bond/MEGA_shared/hic/hic_files/timepoints_normAfter/MEGA_K562_WT_0_inter.hic", 
            "/Users/phanstiel11//Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/4dn/MEGA_K562_WT_PMA_0_100M.hic", 
            "/Users/phanstiel11//Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/4dn/buffer_dpn_4h.hic",
            "/Users/phanstiel11//Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/4dn/dpn_dpn_5m.hic")


## Select which loops to plot
plot <- c(2,6,3)
x <- c(0.5, 3.5, 6.5)
w <- c(5,7,0.25)
res <- 10e3
z1 <- c(300,400,150)
z2 <- c(40,70,40)
z3 <- c(20,25,20)
z4 <- c(20,25,20)

pdf("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/reviewFigs/subsampledLoops.pdf", width = 8.5, height = 11)

## Create page
pageCreate(width = 8.5, height = 11, showGuides = F, xgrid=0.0, ygrid=0.0)

for (i in 1:length(plot)){
  
  #Resize region
  reg <- loops[plot[i]]
  #reg <- loops[i] 
  width <- start(anchors(reg)$second) - start(anchors(reg)$first)
  reg <- 
    reg |>
    resize(width = 1, fix = 'center') |>
    shift(shift = res/2) |>
    resize(width = width * w[i], fix = 'center')
  
  ## Define params
  p <- pgParams(chrom = as.character(seqnames(anchors(reg)$first)),
                chromstart = start(anchors(reg)$first),
                chromend = start(anchors(reg)$second),
                assembly = "hg38",
                resolution = res,
                width = 2.5,
                x = x[i],
                y = 0.5,
                buffer = 0.05,
                height = 1.5)
  
  p_alt <- pgParams(chrom = paste0("chr", as.character(seqnames(anchors(reg)$first))),
                chromstart = start(anchors(reg)$first),
                chromend = start(anchors(reg)$second),
                assembly = "hg38",
                resolution = res,
                width = 2.5,
                x = x[i],
                y = 0.5,
                buffer = 0.05,
                height = 1.5)

  ## Plot hic
  hic1 <- plotHicRectangle(hic[[1]], params = p, y = 0.5, norm = "SCALE", zrange = c(0,z1[i]))
  hic2 <- plotHicRectangle(hic[[2]], params = p_alt, y = 2.1, norm = "SCALE", zrange = c(0,z2[i]))
  hic3 <- plotHicRectangle(hic[[3]], params = p, y = 3.7, norm = "KR", zrange = c(0,z3[i]))
  #hic4 <- plotHicRectangle(hic[[4]], params = p, y = 5.3, norm = "KR", zrange = c(0,z4[i]))
  
  ## Annotate loops
  annoPixels(hic1, anno, type = "circle", shift = 3)
  annoPixels(hic2, annoAlt, type = "circle", shift = 3)
  annoPixels(hic3, anno, type = "circle", shift = 3)
  #annoPixels(hic4, anno, type = "circle", shift = 3)
  
  ## Anno heatmap legend 
  annoHeatmapLegend(
    plot = hic1,
    x = 6.1, y = 0.5, width = 0.1, height = 0.75,
    just = c("left", "top"), default.units = "inches"
  )
  annoHeatmapLegend(
    plot = hic2,
    x = 6.1, y = 2.1, width = 0.1, height = 0.75,
    just = c("left", "top"), default.units = "inches"
  )
  annoHeatmapLegend(
    plot = hic3,
    x = 6.1, y = 3.7, width = 0.1, height = 0.75,
    just = c("left", "top"), default.units = "inches"
  )
  
  ## Anno label genome
  annoGenomeLabel(plot = hic4, params = p, y = 6.9)
  
}

## Annotate each of the plots
plotText(label = "Bond et al", x = 2.05, y = 1)
plotText(label = "K562 full", x = 2.05, y = 1.15)

plotText(label = "Bond et al", x = 2.05, y = 2.85)
plotText(label = "K562 500M", x = 2.05, y = 3)

plotText(label = "Belaghzal", x = 2.05, y = 4.45)
plotText(label = "Buffer 4h", x = 2.05, y = 4.6)

#plotText(label = "Belaghzal", x = 0.05, y = 6.05)
#plotText(label = "DpnII 5min", x = 0.05, y = 6.2)

dev.off()














