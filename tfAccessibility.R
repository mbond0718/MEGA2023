## SPECIFIC TF ACCESSIBILITY
## This script performs intersections of gene clusters with promoters, but specifically gets ATAC peaks for motifs of interest and compares to ATAC counts

## Load libraries ----------------------------------------- 
library(reshape2)
library(GenomicRanges)
library(InteractionSet)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

## Read in data ---------------------------------------------
# Read in text files containing either all or differential elements 
#Loops
loops <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/loopCounts.rds") |> 
  as.data.frame()
loops <- loops[,c(1:3,6:8)]
rownames(loops) <- paste0("loop", 1:nrow(loops))

#ATAC
atacLFC <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/atac/atacLFC.rds")
atacCounts <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/atac/atacCounts.rds")
atacCounts <- atacCounts[,c(4:6)]

#Genes
clusters <- list.files("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/rna/clusters/", full.names = TRUE)
clusters <- clusters[c(4,6,5,1,3,2)]
clusters <- lapply(clusters, read.table)
names <- c("UpEarly", "UpMid", "UpLate", "DnEarly", "DnMid", "DnLate")
diffGenes <- unlist(clusters)

## Convert to GRanges: 
#Loops
loopAnc <- GInteractions(anchor1 = GRanges(seqnames = Rle(loops$seqnames1), 
                                           ranges = IRanges(start = loops$start1, end = loops$end1)), 
                         anchor2 = GRanges(seqnames = Rle(loops$seqnames2), 
                                           ranges = IRanges(start = loops$start2, end = loops$end2)), 
                         name = rownames(loops))
#ATAC 
atacGR <- GRanges(seqnames = Rle(atacLFC$chr), 
                  ranges = IRanges(start = atacLFC$start, end = atacLFC$stop), 
                  peak = rownames(atacLFC))

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
#Get the coordinate information for each of the differential genes
diffPromoters <- promoters[promoters$ensembl_gene_id %in% diffGenes,]


## Overlap promoters and distal anchors of looped promoters with ATAC peaks to get the inputs for HOMER  -----------------------------------------------------
path <- "~/Desktop/"
dirName <- "220708_FindMotifs_perCluster"

## Motifs to find
find <- c("hnf6", "cux2", "hnf1b", "ews-fli1", "nfkb", "p50", "hnf6b", "fosl2", "ews-erg", "xbox", "rfx1", "p65", "nf1-half", "nkx3.1", "foxo1", "pse", "p73", "znf165", "bcl6", "hlf", "nf1-fox", "tbet", "gre", "tr4", "zbtb12", "nfil3", "gfx", "tbx5", "atf1", "foxl2", "hnf4", "myog", "foxa1", "e2f", "dlx5", "cdx4", "reverb", "znf136", "eomes", "boris", "hoxc9", "zbtb18", "znf669", "nrf2", "hoxd11", "cre", "bach1", "hoxa11", "oct4-sox17")
find <- c("jun-ap1", "p65", "tbet", "e2a", "cdx4", "nrf2")
find <- paste0(find, ".motif")

prepHomer <- function(path, dirName){
  #Define path to directory
  directoryPath <- paste0(path, dirName)
  #Create directories
  dir.create(directoryPath)
  dir.create(paste0(directoryPath, "/proximal"))
  dir.create(paste0(directoryPath, "/distal"))
  dir.create(paste0(directoryPath, "/proximal/input"))
  dir.create(paste0(directoryPath, "/proximal/output"))
  dir.create(paste0(directoryPath, "/proximal/commands"))
  dir.create(paste0(directoryPath, "/distal/input"))
  dir.create(paste0(directoryPath, "/distal/output"))
  dir.create(paste0(directoryPath, "/distal/commands"))
  
  
  ## PROXIMAL ##
  
  ## Background using static genes only 
  proximalBackground <- subsetByOverlaps(atacGR, staticPromoters)$peak
  proximalBackground <- (atacCounts[rownames(atacCounts) %in% proximalBackground,])[,1:3]
  proximalBackground$strand <- "+"
  proximalBackground$peak <- rownames(proximalBackground)
  proximalBackground <- proximalBackground[,c(5,1,2,3,4)]
  
  ## Foreground using all differential promoters 
  proximalInput <- subsetByOverlaps(atacGR, diffPromoters)$peak
  proximalInput <- (atacCounts[rownames(atacCounts) %in% proximalInput,])[,1:3]
  proximalInput$strand <- "+"
  proximalInput$peak <- rownames(proximalInput)
  proximalInput <- proximalInput[,c(5,1,2,3,4)]
  
  ## Write files
  write.table(proximalInput, paste0(directoryPath, "/proximal/input/overlap.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
  write.table(proximalBackground, paste0(directoryPath, "/proximal/input/background.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
  
  
  ## DISTAL ##
  
  ## Background using static promoters 
  distalBackground <- linkOverlaps(loopAnc, staticPromoters, atacGR)$subject2
  distalBackground <- atacGR[distalBackground,] |> as.data.frame()
  distalBackground$strand <- "+"
  distalBackground$peak <- distalBackground$peak
  distalBackground <- distalBackground[,c(6,1,2,3,5)]
  distalBackground <- unique(distalBackground)
  
  ## Foreground using all differential promoters 
  distalInput <- linkOverlaps(loopAnc, diffPromoters, atacGR)$subject2
  distalInput <- atacGR[distalInput,] |> as.data.frame()
  distalInput$strand <- "+"
  distalInput$peak <- distalInput$peak
  distalInput <- distalInput[,c(6,1,2,3,5)]
  distalInput <- unique(distalInput)
  
  ## Foreground using only the cluster that the motif came from in the previous analysis: 
  for (i in 1:length(find)){
    dist <- linkOverlaps(loopAnc, promoters[promoters$ensembl_gene_id %in% clusters[[i]]$V1], atacGR)$subject2
    dist <- atacGR[dist,] |> as.data.frame()
    dist$strand <- "+"
    dist$peak <- dist$peak
    dist <- dist[,c(6,1,2,3,5)]
    dist <- unique(dist)
    write.table(dist, paste0(directoryPath, "/distal/input/", names[i], ".txt"), quote = F, sep = "\t", row.names = FALSE)
  }
  
  ## Write files
  write.table(distalBackground, paste0(directoryPath, "/distal/input/background.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
  
  ## COMMANDS ## 
  #parameters: 
  func <- "findMotifsGenome.pl "
  g <- " hg38 "
  size <- " -size 500 "
  motifPath <- paste0(" -find ~/Homer/motifs/")
  
  #Write distal using each cluster .txt file as input for specific motif 
  outPath <- paste0(directoryPath, "/distal/output/")
  backgroundPath <- paste0("-bg ", directoryPath, "/distal/input/background.txt")
  for (i in 1:length(find)){
    command <- 
      paste0(func, 
             directoryPath, "/distal/input/", names[i], ".txt", 
             g,
             outPath, 
             size, 
             backgroundPath, 
             motifPath, find[i], 
             " > ", outPath, find[i], ".txt")
    write(command, paste0(directoryPath, "/distal/commands/", find[i], ".sh"))
  }
  
  #Make the commands executable
  setwd(directoryPath)
  system("chmod +x ./proximal/commands/*")
  system("chmod +x ./distal/commands/*")
  
}
prepHomer("~/Desktop/", "220629_FindMotifs")

## Get accessibility for peaks that have these motifs and plot ---------------------------------------------------------------------------------------------

directoryPath <- "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/atac/individualMotifs/"
find <- c("nfkb", "nf1-fox", "nrf2")
find <- paste0(find, ".motif")
##Just distal FC and counts: 
for (i in 1:length(find)){
  data <- read.delim(paste0(directoryPath, find[i], ".txt"))
  counts <- (atacCounts[rownames(atacCounts) %in% data$PositionID,])
  colnames(counts) <- c("0", "6", "72")
  counts$reg <- "distal"
  counts <- melt(counts, id = "reg")
  scale <- (summary(counts$value))[5]*4
  #Plot
  ggplot(counts, aes(x=variable, y=value)) + 
    geom_boxplot(fill = "#6BAED6", alpha = 0.85, outlier.shape = NA) + ggtitle(unlist(strsplit(find[i], ".motif"))) +
    ylab("ATAC Log2(Fold Change)") + xlab("Time") + ylim(0,scale) +
    theme_classic()
  ggsave(last_plot(), filename = paste0("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig2/", unlist(strsplit(find[i], ".motif")), "AccessibilityCounts.pdf"))
}









