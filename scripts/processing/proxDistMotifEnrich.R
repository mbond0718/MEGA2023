## PROXIMAL/DISTAL MOTIF ENRICHMENT
## Perform motif enrichment on (1) promoters of each gene cluster and (2) promoters looped to each gene cluster
## This script (1) prepared HOMER directories, performs intersections, generates matched backgrounds, (2) generates HOMER commands, (3) generate plots 

## Load libraries -----------------------------------------
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(InteractionSet)
library(nullranges)
library(RColorBrewer)
library(ggplot2)


## Read in data ----------------------------------------------------------------------

# Read in text files containing either all or differential elements 
#Loops
loops <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/loopCounts.rds") |> 
  as.data.frame()
loops <- loops[,c(1:3,6:8)]
rownames(loops) <- paste0("loop", 1:nrow(loops))

#ATAC
allATAC <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/atac/atacLFC.rds")

#Genes
clusters <- list.files("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/rna/clusters/", full.names = TRUE)
clusters <- clusters[c(4,6,5,1,3,2)]
clusters <- lapply(clusters, read.table)
names <- c("UpEarly", "UpMid", "UpLate", "DnEarly", "DnMid", "DnLate")

geneFC <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/rna/rnaLFC.rds")
geneCounts <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/rna/rnaCounts.rds")

## Convert to GRanges: 
#Loops
loopAnc <- GInteractions(anchor1 = GRanges(seqnames = Rle(loops$seqnames1), 
                                           ranges = IRanges(start = loops$start1, end = loops$end1)), 
                         anchor2 = GRanges(seqnames = Rle(loops$seqnames2), 
                                           ranges = IRanges(start = loops$start2, end = loops$end2)), 
                         name = rownames(loops))
#ATAC 
allATACGR <- GRanges(seqnames = Rle(allATAC$chr), 
                     ranges = IRanges(start = allATAC$start, end = allATAC$stop), 
                     peak = rownames(allATAC))

#Genes
txdb <- makeTxDbFromEnsembl(organism="Homo sapiens",
                            release=NA,
                            circ_seqs=NULL,
                            server="ensembldb.ensembl.org",
                            username="anonymous", password=NULL, port=0L,
                            tx_attrib=NULL)

genes <- genes(txdb)
genes <- data.frame(genes)
genes$seqnames <- paste0("chr", genes$seqnames)
genes <- GRanges(seqnames = Rle(genes$seqnames), 
                 ranges = IRanges(start = genes$start, end = genes$end), 
                 strand = genes$strand,
                 ensembl_gene_id = genes$gene_id)
promoters <- promoters(genes)
staticPromoters <- promoters[!promoters$ensembl_gene_id %in% diffGenes,]
#Get the coordinate information for each of the genes in each of the clusters: 
for (i in 1:length(clusters)){
  c <- clusters[[i]]
  clusters[[i]] <- promoters[promoters$ensembl_gene_id %in% c$V1,] 
}

## Define fucntions and run Homer  ----------------------------------------------------------------------
prepHomerMatch <- function(path, dirName, listOfCategories){
  #Define path to directory
  directoryPath <- paste0(path, dirName)
  #Create directories
  dir.create(directoryPath)
  dir.create(paste0(directoryPath, "/proximal"))
  dir.create(paste0(directoryPath, "/distal"))
  dir.create(paste0(directoryPath, "/proximal/commands"))
  dir.create(paste0(directoryPath, "/distal/commands"))
  for (i in 1:length(listOfCategories)){
    dir.create(paste0(directoryPath, "/proximal/input/", listOfCategories[[i]]), recursive = TRUE)
    dir.create(paste0(directoryPath, "/proximal/output/", listOfCategories[[i]]), recursive = TRUE)
    dir.create(paste0(directoryPath, "/distal/input/", listOfCategories[[i]]), recursive = TRUE)
    dir.create(paste0(directoryPath, "/distal/output/", listOfCategories[[i]]), recursive = TRUE)
  }
  
  ## Static genes to potentially be used in background: 
  proximalBackground <- subsetByOverlaps(staticPromoters, allATACGR)
  
  ## Intersections for each group of interest: 
  for (i in 1:length(clusters)){
    c <- clusters[[i]]
    focal <- geneCounts[rownames(geneCounts) %in% c$ensembl_gene_id,]
    background <- matchRanges(focal = focal, 
                              pool = geneCounts[rownames(geneCounts) %in% proximalBackground$ensembl_gene_id,], 
                              covar = ~median, method ="stratified")
    c <- subsetByOverlaps(allATACGR, c)$peak
    c <- (allATAC[rownames(allATAC) %in% c,])[,1:3]
    c$strand <- "+"
    c$peak <- rownames(c)
    c <- c[,c(5,1,2,3,4)]
    background <- subsetByOverlaps(allATACGR, promoters[promoters$ensembl_gene_id %in% rownames(background)])$peak
    background <- (allATAC[rownames(allATAC) %in% background,])[,1:3]
    background$strand <- "+"
    background$peak <- rownames(background)
    background <- background[,c(5,1,2,3,4)]
    #Write to a file:
    write.table(c, paste0(directoryPath, "/proximal/input/", listOfCategories[[i]], "/overlap.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
    write.table(background, paste0(directoryPath, "/proximal/input/", listOfCategories[[i]], "/background.txt"), quote = FALSE, sep = "\t", row.names = FALSE)}
  
  #Distal genes to be used in distal background
  distalBackground <- linkOverlaps(loopAnc, staticPromoters, allATACGR)$subject1
  distalBackground <- staticPromoters[distalBackground,]
  
  ## Intersections for each group of interest
  for (i in 1:length(clusters)){
    c <- clusters[[i]]
    focal <- geneCounts[rownames(geneCounts) %in% c$ensembl_gene_id,]
    background <- matchRanges(focal = focal, 
                              pool = geneCounts[rownames(geneCounts) %in% distalBackground$ensembl_gene_id,], 
                              covar = ~median, method ="stratified")
    link <- linkOverlaps(loopAnc, c, allATACGR)$subject2
    distal <- allATACGR[link] |> as.data.frame()
    distal$strand <- "+"
    distal <- distal[,c(6,1,2,3,5)]
    distal <- unique(distal)
    
    link <- linkOverlaps(loopAnc, staticPromoters[staticPromoters$ensembl_gene_id %in% rownames(background)], allATACGR)$subject2
    background <- allATACGR[link]
    background <- (allATAC[rownames(allATAC) %in% background$peak,])[,1:3]
    background$strand <- "+"
    background$peak <- rownames(background)
    background <- background[,c(5,1,2,3,4)]
    background <- unique(background)
    
    #Write to a file:
    write.table(distal, paste0(directoryPath, "/distal/input/", listOfCategories[[i]], "/overlap.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
    write.table(background, paste0(directoryPath, "/distal/input/", listOfCategories[[i]], "/background.txt"), quote = FALSE, sep = "\t", row.names = FALSE)}
  
  #Write commands
  #Proximal: 
  for (i in 1:length(clusters)){
    func <- "findMotifsGenome.pl "
    overlapPath <- paste0(directoryPath, "/proximal/input/", names[i], "/overlap.txt")
    g <- " hg38 "
    outPath <- paste0(directoryPath, "/proximal/output/", names[i], "/")
    size <- " -size 500 "
    backgroundPath <- paste0("-bg ", directoryPath, "/proximal/input/", names[i], "/background.txt")
    command <- paste0(func, overlapPath, g, outPath, size, backgroundPath, " -nomotif")
    write(command, paste0(directoryPath, "/proximal/commands/", names[i], ".sh"))}
  #Distal: 
  for (i in 1:length(clusters)){
    overlapPath <- paste0(directoryPath, "/distal/input/", names[i], "/overlap.txt")
    outPath <- paste0(directoryPath, "/distal/output/", names[i], "/")
    backgroundPath <- paste0("-bg ", directoryPath, "/distal/input/", names[i], "/background.txt")
    command <- paste0(func, overlapPath, g, outPath, size, backgroundPath, " -nomotif")
    write(command, paste0(directoryPath, "/distal/commands/", names[i], ".sh"))}
  
  #Make the commands executable
  setwd(directoryPath)
  system("chmod +x ./proximal/commands/*")
  system("chmod +x ./distal/commands/*")
  
}

#Run this function to create directories and overlap/background files: 
path <- "~/Desktop/"
dirName <- "220627_MotifHeatmap_matchedBack_500bp"
prepHomerMatch("~/Desktop/", "220627_MotifHeatmap_matchedBack_500bp")

## Plot results -----------------------------------------------------------------------------------------
pal <- brewer.pal(8, "YlGnBu")
path <- "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/atac/"
dirName <- "proxDistMotifEnrich"
proxDist <- "distal"
listOfCategories <- c("UpEarly", "UpMid", "UpLate", "DnEarly", "DnMid", "DnLate")
plotMotifs <- function(path, dirName, listOfCategories, proxDist){
  
  #Phase 1: read in each knownResults.txt file
  directoryPath <- paste0(path, dirName)
  motifs <- list()
  for (i in 1:length(listOfCategories)){
    motifs[[i]] <- read.delim(paste0(directoryPath, "/", proxDist, "/output/", names[i], "/knownResults.txt"))
    #motifs[[i]] <- read.delim(paste0(directoryPath, "/", proxDist, "/", names[i], "/output/knownResults.txt"))
  }
  names(motifs) <- listOfCategories
  
  #Phase 3: clean up each motif file
  newcolnames <- c("MotifName", "Consensus", "P.value", "Log P.value", "q value", 
                   "Number target sequences", "Percent_target_sequences", 
                   "Number background sequences", "Percent_background_sequences")
  abridge <- function(longname){
    temp <- strsplit(longname, '[(]')
    shortname <- unlist(temp)[1]
  }
  for (i in 1:length(motifs)){
    #Subset by the motif list
    motif <- motifs[[i]]
    #Assign new colnames
    colnames(motif) <- newcolnames
    #Assign abreviated names
    motif$MotifName <- sapply(1:nrow(motif), function(x) abridge((as.character(motif[,1]))[x]))
    #Remove % signs from percent target sequences and percent background sequences
    motif$Percent_target_sequences <- as.numeric(sub("%", "", motif$Percent_target_sequences))
    motif$Percent_background_sequences <- as.numeric(sub("%", "", motif$Percent_background_sequences))
    #Adjust p value: 
    motif$padj <- p.adjust(exp(motif$`Log P.value`))
    #Calculate log2 enrichment and log10pval
    motif$log2enrichment <- log2(motif$Percent_target_sequences/motif$Percent_background_sequences)
    motif$log10pval <- -log10(exp(motif$`Log P.value`))
    motif$log10padj <- -log10(motif$padj)
    motif$class <- listOfCategories[[i]]
    
    #Reassign to the list 
    motifs[[i]] <- motif
  }
  
  ##Alternate version, not by plotting specific motifs of interest: 
  data <- do.call(rbind, motifs) |> as.data.frame()
  data$class <- factor(data$class, names[length(listOfCategories):1])
  #Threshold for points with a significant log10p
  threshold <- 1.3
  temp <- data[data$log10pval>threshold,]
  order <- c()
  for (i in 1:length(listOfCategories)){
    #order[[i]] <- temp[temp$class==listOfCategories[i] & temp$log10pval>threshold,]$MotifName |> unique()
    x <- temp[temp$class==listOfCategories[i] & temp$log10pval>threshold,]$MotifName |> unique()
    order[[i]] <- x[1:10]
    
  }
  order <- unlist(order) |> unique()
  order <- order[!order %in% c("Pax7", "p53", "GATA3", "RAR:RXR", "PAX5", "GATA")]
  temp <- data[data$MotifName %in% order,]
  temp$MotifName <- factor(temp$MotifName, order)
  
  ggplot(temp, aes(x=MotifName, y=class, fill=log10pval)) + 
    geom_tile(color = "white") + scale_fill_gradientn(colors = brewer.pal(8, "Blues")) +
    coord_fixed() + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + 
    ggtitle(paste0(proxDist, " elements"))
  ggsave(last_plot(), filename = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig2/distHeatmap.pdf")

}
