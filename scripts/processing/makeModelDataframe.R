## PREP DATA FOR LINEAR MODEL IN FIGURE 4
## This script performs intersections between each feature and loops for either (1) count sum or (2) bigwig max

## Load libraries -----------------------------------------
library(sageseqr)
library(DESeq2)
library(GenomicInteractions)
library(purrr)
library(preprocessCore)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(mariner)

## Define functions --------------------------------------
as_grange <- function(df){
  GRanges(seqnames = Rle(df$chr), 
          ranges = IRanges(start = df$start, end = df$stop),
          name = rownames(df))
}

## Read in data -------------------------- 
#Loops
loops <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/loopCounts.rds") |> as.data.frame()
loops$loop <- paste0("loop", rownames(loops))
loops <- loops[,c(1:3,6:8,21:51)]
loopCounts <- loops

#ATAC 
atac <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/atac/atacCounts.rds")
atac <- atac[,c(1:4,6)]
atacBig <- (list.files("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/atac/", ".bw", full.names = TRUE))[c(1,3)]
#Enhancer
enh <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/cnr/H3K27Counts.rds")
enh <- enh[,c(1:4,6)]
enhBig <- (list.files("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/cnr/", ".bw", full.names = TRUE))[c(8,10)]
#Jun
jun <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/cnr/JunCounts.rds")
jun <- jun[,c(1:4,6)]
junBig <- (list.files("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/cnr/", ".bw", full.names = TRUE))[c(11,13)]
#CTCF
ctcf <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/cnr/CTCFCounts.rds")
ctcf <- ctcf[,c(1:4,6)]
ctcfBig <- (list.files("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/cnr/", ".bw", full.names = TRUE))[c(5,7)]
#Rad21
rad <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/cnr/Rad21Counts.rds")
rad <- rad[,1:5]
radBig <- (list.files("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/cnr/", ".bw", full.names = TRUE))[c(14,15)]
#RNA
rna <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/rna/binnedTPM.rds")
rna <- rna[,c(7:9,11,13)]
rnaBig <- (list.files("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/rna/", ".bw", full.names = TRUE))[c(1,8)]


## Convert to GRanges: 
#Loop anchors
loopAnc <- as_ginteractions(loops) |> 
  swapAnchors()
#Loop internal
loopIn <- GRanges(seqnames = seqnames(anchors(loopAnc)$first), 
                  ranges = IRanges(start = end(anchors(loopAnc)$first), end = start(anchors(loopAnc)$second)), 
                  loop = loopAnc$loop)

#ATAC 
atacGR <- as_grange(atac)
#Enhancer
enhGR <- as_grange(enh)
#Jun
junGR <- as_grange(jun)
#CTCF 
ctcfGR <- as_grange(ctcf)
#Cohesin
radGR <- as_grange(rad)
#Genes
genes <- loadDb("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/ref/hg38txdb.sqlite") |> 
  genes()
genes <- data.frame(genes)
genes$seqnames <- paste0("chr", genes$seqnames)
colnames(genes)[6] <- "gene"
rna <- merge(rna, genes, by="gene")
rna <- rna[,c(1,7:9,14:16)]
rnaGR <- GRanges(seqnames = Rle(rna$bin_seqnames), 
                 ranges = IRanges(start = rna$bin_start, end = rna$bin_end), 
                 gene = rna$gene)

#manual set test loop
loop <- loops$loop[2]

## Path to arguments for overlapping functions
granges <- list(atacGR, enhGR, junGR, ctcfGR, radGR, rnaGR)
bws <- list(atacBig, enhBig, junBig, ctcfBig, radBig, rnaBig)
dfs <- list(atac, enh, jun, ctcf, rad, rna)
letters <- list("A", "K", "J", "C", "R", "T")

## Function to get the sum of raw counts
getSum <- function(loop){
  
  print(loop)  
  
  #Get loop coords
  mcols <- (loops[loops$loop == loop,])
  coords <- mcols[,1:6]
  anc <- loopAnc[loopAnc$loop == loop,]
  int <- loopIn[loopIn$loop == loop,]
  
  #Delta loop strength: 
  loopStrength <- c(sum(mcols[7:16]), 
                    sum(mcols[27:36]))
  
  #Loop distance: 
  loopDist <- coords$start2 - coords$start1
  
  #Function to overlap each regulatory feature with loop anchors and internal 
  olap <- function(gr, df, letter){
    #anchor 1
    anc1 <- findOverlaps(gr, anchors(anc)$first) |> 
      queryHits()
    anc1 <- gr[anc1]$name
    anc1 <- (df[rownames(df) %in% anc1,])[,c(4,5)]
    if (nrow(anc1)>1){
      anc1sum <- c(sum(anc1[,1]), sum(anc1[,2]))
    } else if (nrow(anc1)==0){
      anc1sum <- c(0,0)
    } else if (nrow(anc1)==1){
      anc1sum <- anc1}
    #anchor 2
    anc2 <- findOverlaps(gr, anchors(anc)$second) |> 
      queryHits()
    anc2 <- gr[anc2]$name
    anc2 <- (df[rownames(df) %in% anc2,])[,c(4,5)]
    if (nrow(anc2)>1){
      anc2sum <- c(sum(anc2[,1]), sum(anc2[,2]))
    } else if (nrow(anc2)==0){
      anc2sum <- c(0,0)
    } else if (nrow(anc2)==1){
      anc2sum <- anc2}
    #internal
    inter <- findOverlaps(gr, int) |> 
      queryHits()
    inter <- gr[inter]$name
    inter <- (df[rownames(df) %in% inter,])[,c(4,5)]
    if (nrow(inter)>1){
      intsum <- c(sum(inter[,1]), sum(inter[,2]))
    } else if (nrow(inter)==0){
      intsum <- c(0,0)
    } else if (nrow(inter)==1){
      intsum <- inter}
    
    data <- data.frame(anc1sum[1], anc1sum[2], anc2sum[1], anc2sum[2], intsum[1], intsum[2])
    cols <- c("Anc1_Sum_0", "Anc1_Sum_72", "Anc2_Sum_0", "Anc2_Sum_72", "Int_Sum_0", "Int_Sum_72")
    colnames(data) <- paste0(letter, "_", cols)
    data
    
  }
  
  row <- pmap(list(granges, dfs, letters), \(x, y, z) olap(gr = x, df = y, letter = z))
  row <- do.call(cbind, row)
  
  loopDat <- data.frame(coords, loopStrength[1], loopStrength[2], loopDist)
  colnames(loopDat) <- c("L_chr1", "L_start1", "L_end1", "L_chr2", "L_start2", "L_end2", "L_counts_0", "L_counts_72", "L_Length")
  data <- cbind(loopDat, row)
  assign("data", data, envir = globalenv())
  
}
#Run Function
data <- lapply(loops$loop[1:100], getSum)
data <- do.call(rbind, data)
rownames(data) <- loops$loop
saveRDS(data, "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/model/modelSum.rds")

## Function to get the max from bigwig files: 
getMax <- function(loop){
  
  print(loop)  
  
  #Get loop coords
  mcols <- (loops[loops$loop == loop,])
  coords <- mcols[,1:6]

  #Delta loop strength: 
  loopStrength <- c(sum(mcols[7:16]), 
                    sum(mcols[27:36]))
  
  #Loop distance: 
  loopDist <- coords$start2 - coords$start1
  
  #Function to overlap each regulatory feature with loop anchors and internal 
  olap <- function(bw, letter){
    #anchor 1
    anc1max <- lapply(enhBig, \(x) readBigwig(file = x, chrom = coords[,1], chromstart = coords[,2], chromend = coords[,3])[c(6)])
    anc1max <- lapply(anc1max, max)
    #anchor 2
    anc2max <- lapply(enhBig, \(x) readBigwig(file = x, chrom = coords[,4], chromstart = coords[,5], chromend = coords[,6])[c(6)])
    anc2max <- lapply(anc2max, max)
    #internal
    intmax <- lapply(enhBig, \(x) readBigwig(file = x, chrom = coords[,1], chromstart = coords[,3], chromend = coords[,5])[6])
    intmax <- lapply(intmax, max)
    
    data <- data.frame(anc1max[[1]], anc1max[[2]], anc2max[[1]], anc2max[[2]], intmax[[1]], intmax[[2]])
    cols <- c("Anc1_Max_0", "Anc1_Max_72", "Anc2_Max_0", "Anc2_Max_72", "Int_Max_0", "Int_Max_72")
    colnames(data) <- paste0(letter, "_", cols)
    data
  }
  
  row <- pmap(list(bws, letters), \(x, y) olap(bw = x, letter = y))
  row <- do.call(cbind, row)
  
  loopDat <- data.frame(coords, loopStrength[1], loopStrength[2], loopDist)
  colnames(loopDat) <- c("L_chr1", "L_start1", "L_end1", "L_chr2", "L_start2", "L_end2", "L_counts_0", "L_counts_72", "L_Length")
  data <- cbind(loopDat, row)
  assign("data", data, envir = globalenv())
  
}
#Run Function
data <- lapply(loops$loop[1:10], getMax)
data <- do.call(rbind, data)
rownames(data) <- loops$loop
saveRDS(data, "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/model/modelMax.rds")


