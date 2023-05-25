## HOMER MOTIF ENRICHMENT - ATAC PEAKS 
## Perform motif enrichment on gained/lost atac peaks, atac peaks at loop anchors, and atac peaks at gained/lost loops
## This script (1) prepares HOMER directories and performs genomic intersections, (2) generates commands to launch HOMER, and (3) generates plots

## Load libraries -----------------------------------------
library(GenomicRanges)
library(nullranges)
library(ggplot2)
library(InteractionSet)
library(dplyr)
library(gridExtra)

## Generate Homer Directories and Inputs  ------------------
#Function to create the directories that will be used with the primary function
#For the path argument, this will be where the directory will be created
#For the name argument, provide the name of the directory
makeOverlapDirs <- function(path, name){
  
  directoryPath <- paste0(path, name)
  dir.create(paste0(directoryPath, "/input/upATAC"), recursive = TRUE)
  dir.create(paste0(directoryPath, "/input/downATAC"), recursive = TRUE)
  dir.create(paste0(directoryPath, "/input/allLoopATAC"), recursive = TRUE)
  dir.create(paste0(directoryPath, "/input/gainedLoopATAC"), recursive = TRUE)
  dir.create(paste0(directoryPath, "/input/lostLoopATAC"), recursive = TRUE)
  
  dir.create(paste0(directoryPath, "/output/upATAC"), recursive = TRUE)
  dir.create(paste0(directoryPath, "/output/downATAC"), recursive = TRUE)
  dir.create(paste0(directoryPath, "/output/allLoopATAC"), recursive = TRUE)
  dir.create(paste0(directoryPath, "/output/gainedLoopATAC"), recursive = TRUE)
  dir.create(paste0(directoryPath, "/output/lostLoopATAC"), recursive = TRUE)
  
  dir.create(paste0(directoryPath, "/plots/upATAC"), recursive = TRUE)
  dir.create(paste0(directoryPath, "/plots/downATAC"), recursive = TRUE)
  dir.create(paste0(directoryPath, "/plots/allLoopATAC"), recursive = TRUE)
  dir.create(paste0(directoryPath, "/plots/gainedLoopATAC"), recursive = TRUE)
  dir.create(paste0(directoryPath, "/plots/lostLoopATAC"), recursive = TRUE)
  
  dir.create(paste0(directoryPath, "/commands"))
}

#Run this line (adjust arguments as necessary) to create the directories in the correct structure
makeOverlapDirs("~/Desktop/", "MEGA_newMotifEnrich")

#Provide paths to input files:
atacPeaksPath <- "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/atac/atacLFC.rds"
allLoopsPath <- "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/loopCounts.rds"
diffLoopsPath <- "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/diffLoops.rds"
#This should be the directory that you created with the previous command (makeOverlapDirs)
outDirPath <- "~/Desktop/MEGA_newMotifEnrich/input/"

#Function to create overlap and background files (input for Homer) into specified directory
atacLoopOverlap <- function(atacPeaksPath, allLoopsPath, diffLoopsPath, outDirPath){
  
  #Read in lists from input files: 
  atacPeaks <- readRDS(atacPeaksPath)
  diffPeaks <- atacPeaks[atacPeaks$differential==TRUE,]
  gainedPeaks <- diffPeaks[diffPeaks$`4320`>0,]
  lostPeaks <- diffPeaks[diffPeaks$`4320`<0,]
  
  allLoops <- readRDS(allLoopsPath) |> 
    as.data.frame()
  allLoops <- allLoops[,c(1:3,6:8)]
  diffLoops <- readRDS(diffLoopsPath)
  gainedLoops <- diffLoops[diffLoops$class=="gained",]
  lostLoops <- diffLoops[diffLoops$class=="lost",]

  #Adjust atacPeaks into correct format
  atacPeaks$name <- paste0("peak", 1:nrow(atacPeaks))
  atacPeaks <- atacPeaks[,c(1:3,9)]
  gainedPeaks$name <- rownames(gainedPeaks)
  gainedPeaks <- gainedPeaks[,c(1:3,9)]
  lostPeaks$name <- rownames(lostPeaks)
  lostPeaks <- lostPeaks[,c(1:3,9)]
  
  #Make GRanges
  allLoops <- GInteractions(anchor1 = GRanges(seqnames = Rle(allLoops$seqnames1), 
                                                    ranges = IRanges(start = allLoops$start1, end = allLoops$end1)), 
                                  anchor2 = GRanges(seqnames = Rle(allLoops$seqnames2), 
                                                    ranges = IRanges(start = allLoops$start2, end = allLoops$end2)), 
                                  name = rownames(allLoops))
  gainedLoops <- GInteractions(anchor1 = GRanges(seqnames = Rle(gainedLoops$seqnames1), 
                                                       ranges = IRanges(start = gainedLoops$start1, end = gainedLoops$end1)), 
                                     anchor2 = GRanges(seqnames = Rle(gainedLoops$seqnames2), 
                                                       ranges = IRanges(start = gainedLoops$start2, end = gainedLoops$end2)), 
                                     name = rownames(gainedLoops))
  lostLoops <- GInteractions(anchor1 = GRanges(seqnames = Rle(lostLoops$seqnames1), 
                                                     ranges = IRanges(start = lostLoops$start1, end = lostLoops$end1)), 
                                   anchor2 = GRanges(seqnames = Rle(lostLoops$seqnames2), 
                                                     ranges = IRanges(start = lostLoops$start2, end = lostLoops$end2)), 
                                   name = rownames(lostLoops))
  
  atacPeaks <- GRanges(seqnames = Rle(atacPeaks$chr), 
                       ranges = IRanges(start = atacPeaks$start, end = atacPeaks$stop), 
                       name = atacPeaks$name)
  
  #Create loopPeaks for downstream analysis: the atac peaks located within 5kb of ANY loop anchor
  loopPeaks <- subsetByOverlaps(atacPeaks, allLoops)
  
  ## Set 1: all gained ATAC peaks 
  overlap <- data.frame(gainedPeaks)
  background <- data.frame(atacPeaks)
  background <- background[,c(1:3,6)]
  overlap$strand <- "+"
  background$strand <- "+"
  overlap <- overlap[,c(4,1,2,3,5)]
  background <- background[,c(4,1,2,3,5)]
  colnames(overlap) <- c("peak", "seqnames", "start", "stop", "strand")
  colnames(background) <- c("peak", "seqnames", "start", "stop", "strand")
  
  write.table(overlap, paste0(outDirPath, "upATAC/overlap.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
  write.table(background, paste0(outDirPath, "upATAC/background.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
  
  ## Set 2: all lost ATAC peaks 
  overlap <- data.frame(lostPeaks)
  background <- data.frame(atacPeaks)
  background <- background[,c(1:3,6)]
  overlap$strand <- "+"
  background$strand <- "+"
  overlap <- overlap[,c(4,1,2,3,5)]
  background <- background[,c(4,1,2,3,5)]
  colnames(overlap) <- c("peak", "seqnames", "start", "stop", "strand")
  colnames(background) <- c("peak", "seqnames", "start", "stop", "strand")
  
  write.table(overlap, paste0(outDirPath, "downATAC/overlap.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
  write.table(background, paste0(outDirPath, "downATAC/background.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
  
  ## Set 3: all loops with an ATAC peak: 
  overlap <- subsetByOverlaps(atacPeaks, allLoops)
  background <- data.frame(atacPeaks)
  background <- background[,c(1:3,6)]
  overlap <- data.frame(overlap)
  overlap <- overlap[,c(1:3,6)]
  overlap$strand <- "+"
  background$strand <- "+"
  overlap <- overlap[,c(4,1,2,3,5)]
  background <- background[,c(4,1,2,3,5)]
  colnames(overlap) <- c("peak", "seqnames", "start", "stop", "strand")
  colnames(background) <- c("peak", "seqnames", "start", "stop", "strand")
  
  write.table(overlap, paste0(outDirPath, "allLoopATAC/overlap.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
  write.table(background, paste0(outDirPath, "allLoopATAC/background.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
  
  ## Set 4: all ganied Loops with an ATAC peak (with peaks at loops as background): 
  overlap <- subsetByOverlaps(loopPeaks, gainedLoops)
  background <- data.frame(loopPeaks)
  background <- background[,c(1:3,6)]
  overlap <- data.frame(overlap)
  overlap <- overlap[,c(1:3,6)]
  overlap$strand <- "+"
  background$strand <- "+"
  overlap <- overlap[,c(4,1,2,3,5)]
  background <- background[,c(4,1,2,3,5)]
  colnames(overlap) <- c("peak", "seqnames", "start", "stop", "strand")
  colnames(background) <- c("peak", "seqnames", "start", "stop", "strand")
  
  #alt background
  background <- sample_n(background, nrow(overlap))
  
  write.table(overlap, paste0(outDirPath, "gainedLoopATAC/overlap.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
  write.table(background, paste0(outDirPath, "gainedLoopATAC/background.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
  
  ## Set 5: all lost Loops with an ATAC peak (with peaks at loops as background): 
  overlap <- subsetByOverlaps(loopPeaks, lostLoops)
  background <- data.frame(loopPeaks)
  background <- background[,c(1:3,6)]
  overlap <- data.frame(overlap)
  overlap <- overlap[,c(1:3,6)]
  overlap$strand <- "+"
  background$strand <- "+"
  overlap <- overlap[,c(4,1,2,3,5)]
  background <- background[,c(4,1,2,3,5)]
  colnames(overlap) <- c("peak", "seqnames", "start", "stop", "strand")
  colnames(background) <- c("peak", "seqnames", "start", "stop", "strand")
  
  #alt background
  background <- sample_n(background, nrow(overlap))
  
  write.table(overlap, paste0(outDirPath, "lostLoopATAC/overlap.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
  write.table(background, paste0(outDirPath, "lostLoopATAC/background.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
  
}

#Run this line of code to create overlap and background .txt files into specified outDir 
atacLoopOverlap(atacPeaksPath, allLoopsPath, diffLoopsPath, outDirPath)


## Generate Homer Launch Commands  ------------------

#This is a command that will generate a command to be launched with Terminal for each of the groups 
makeHomerCommands <- function(path){
  
  #Define constants for final command: 
  func <- "findMotifsGenome.pl "
  overlapPath <- paste0(path, "/input/")
  g <- " hg38 "
  outPath <- paste0(path, "/output/")
  size <- " -size given "
  backgroundPath <- paste0("-bg ", path, "/input/")
  
  #Command: upATAC  with loop peaks as background 
  x <- "/upATAC/"
  command <- paste0(func, 
                    paste0(overlapPath, x, "overlap.txt"), 
                    g, 
                    paste0(outPath, x), 
                    size, 
                    paste0(backgroundPath, x, "background.txt"))
  
  write(command, paste0(path, "/commands/upATAC.sh"))
  
  #Command: downATAC  with loop peaks as background 
  x <- "/downATAC/"
  command <- paste0(func, 
                    paste0(overlapPath, x, "overlap.txt"), 
                    g, 
                    paste0(outPath, x), 
                    size, 
                    paste0(backgroundPath, x, "background.txt"))
  
  write(command, paste0(path, "/commands/downATAC.sh"))
  
  #Command: all loops with all peaks as background 
  x <- "/allLoopATAC/"
  command <- paste0(func, 
                    paste0(overlapPath, x, "overlap.txt"), 
                    g, 
                    paste0(outPath, x), 
                    size, 
                    paste0(backgroundPath, x, "background.txt"))
  
  write(command, paste0(path, "/commands/allLoopATAC.sh"))
  
  #Command: gained loops with loop peaks as background 
  x <- "/gainedLoopATAC/"
  command <- paste0(func, 
                    paste0(overlapPath, x, "overlap.txt"), 
                    g, 
                    paste0(outPath, x), 
                    size, 
                    paste0(backgroundPath, x, "background.txt"))
  
  write(command, paste0(path, "/commands/gainedLoopATAC.sh"))
  
  #Command: lost loops with loop peaks as background 
  x <- "/lostLoopATAC/"
  command <- paste0(func, 
                    paste0(overlapPath, x, "overlap.txt"), 
                    g, 
                    paste0(outPath, x), 
                    size, 
                    paste0(backgroundPath, x, "background.txt"))
  
  write(command, paste0(path, "/commands/lostLoopATAC.sh"))
  
}

#Set your path to the directory 
homerPath <- "~/Desktop/MEGA_newMotifEnrich/"

#Run this command
makeHomerCommands(homerPath)

#Within the specified directory, this should have created 6 scripts
#Go to the commands directory, and give execution permission with the following command: 
# <chmod +x *>
#Then, from this directory, launch each with <./command.sh> in 6 DIFFERENT TABS if running on local computer

## Generate Homer Motif Enrichment Plots  ------------------


###Load this Function to plot Known Motifs and save the plots to files
plotKnownMotifs <- function(dirPath){
  
  #PHASE 1: read in each knownResults.txt file: 
  upatac <- read.delim(paste0(dirPath, "output/upATAC/knownResults.txt"))
  downatac <- read.delim(paste0(dirPath, "output/downATAC/knownResults.txt"))
  all <- read.delim(paste0(dirPath, "output/allLoopATAC/knownResults.txt"))
  gained <- read.delim(paste0(dirPath, "output/gainedLoopATAC/knownResults.txt"))
  lost <- read.delim(paste0(dirPath, "output/lostLoopATAC/knownResults.txt"))
  
  #PHASE2: put all motif files into a list
  motifs <- list(upatac, downatac, all, gained, lost)
  names <- c("upATAC", "downATAC", "allLoops", "gainedLoops", "lostLoops")
  
  #PHASE3: clean up each motif file (new colnames, remove % signs)
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
    #Calculate log2 enrichment and log10pval
    motif$log2enrichment <- log2(motif$Percent_target_sequences/motif$Percent_background_sequences)
    motif$log10pval <- -log10(exp(motif$`Log P.value`))
    motif$class <- names[i]
    #Reassign to the list 
    motifs[[i]] <- motif
  }
  
  #PHASE 4: take care of any inf outliers in any of the motif sets 
  #Determine maximum for each log2enrichment column
  list <- list()
  for (i in 1:length(motifs)){
    motif <- motifs[[i]]
    enrichment <- motif$log2enrichment[!(motif$log2enrichment==Inf)]
    maxE <- as.list(max(enrichment))
    maxE$i <- i
    list[[i]] <- maxE
  }
  maxE <- unlist(do.call(rbind,list)[,1])
  #Determine maximum for each log10pval column
  list <- list()
  for (i in 1:length(motifs)){
    motif <- motifs[[i]]
    p <- motif$log10pval[!(motif$log10pval==Inf)]
    maxP <- as.list(max(p))
    maxP$i <- i
    list[[i]] <- maxP
  }
  maxP <- unlist(do.call(rbind,list)[,1])
  #Replace any Inf or -Inf with max + 1.1X max for enrichment and p value
  for (i in 1:length(motifs)){
    motif <- motifs[[i]]
    if(length(motif$log2enrichment[motif$log2enrichment==Inf | is.na(motif$log2enrichment)]>0)){
      motif$log2enrichment[motif$log2enrichment==Inf | is.na(motif$log2enrichment)] <- rep((maxE[i]+(maxE[i]*1.1)), length(motif$log2enrichment[motif$log2enrichment==Inf | is.na(motif$log2enrichment)]))
    } else {
      print(paste0("There are no values of Inf for ", i))
    }
    
    if(length(motif$log2enrichment[motif$log2enrichment==-Inf])>0){
      motif$log2enrichment[motif$log2enrichment==-Inf] <- rep((-maxE[i]+(-maxE[i]*1.01)), length(motif$log2enrichment[motif$log2enrichment==-Inf]))
    } else {
      print(paste0("There are no values of -Inf for ", i))
    }
    
    if(length(motif$log10pval[motif$log10pval==Inf | is.na(motif$log10pval)]>0)){
      motif$log10pval[motif$log10pval==Inf | is.na(motif$log10pval)] <- rep((maxP[i]+(maxP[i]*1.01)), length(motif$log10pval[motif$log10pval==Inf | is.na(motif$log10pval)]))
    } else {
      print(paste0("There are no values of Inf for ", i))
    }
    
    if(length(motif$log10pval[motif$log10pval==-Inf]>0)){
      motif$log10pval[motif$log10pval==-Inf] <- rep(0, length(motif$log10pval[motif$log10pval==-Inf]))
    } else {
      print(paste0("There are no values of -Inf for ", i))
    }
    motifs[[i]] <- motif
  }
  
  ##PHASE 5: alternate, top 10 motifs for each category
  moi <- list()
  for (i in 1:length(motifs)){
    motif <- motifs[[i]]
    moi[[i]] <- motif$MotifName[1:10]
  }
  moi <- unlist(moi) |> unique()
  #select from moi: 
  moi <- c("AP-1", "Fos", "JunB", "Jun-AP1", 
           "Gata1", "Gata4", "EKLF", "KLF1", 
           "CTCF", "BORIS", "HIC1",
           "GATA", "FOXA1:AR", "p53", "Etv2",
           "Foxo1", "LEF1", "DLX2", "Brn1")
  
  
  #PHASE 6: Make barplots 
  #cols <- brewer.pal(8, "YlGnBu")[4:8]
  cols <- c("#52A9B2", "#E0524C", "gray50", "#227B7F", "#931C1E")
  
  pdf("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig2/motifEnrichBarplots.pdf", width = 8.5, height = 11)
  plots <- list()
  for (i in 1:length(motifs)){
    motif <- motifs[[i]]
    motif <- motif[motif$MotifName %in% moi,]
    levels(motif$MotifName) <- moi
    plot <- ggplot(data = motif, aes(x = MotifName, y = log10pval)) + 
      geom_bar(stat = "identity", fill = cols[i], alpha = 0.75) + 
      xlab(names[i]) + scale_x_discrete(limits = moi) + 
      theme_classic() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
    plots[[i]] <- plot
    #pdf((paste0(dirPath, "plots", path[i], "knownMotifs_Enrichment.pdf")))
    #print(plot)
    #dev.off()
  }
  
  do.call("grid.arrange", c(plots, ncol=1))
  dev.off()
}

###Set your directory path (this should be the same directory that was created with the makeOverlapDirs command of the functionalPeakLoopOverlap Code)
dirPath <- "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA_shared/MEGAManuscript/Figure2/MotifsV4/"

###Run this line of code to visualize your motif data!
plotKnownMotifs(dirPath)

