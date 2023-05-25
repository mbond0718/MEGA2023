## Set 1: gained loops with gained atac peaks, all peaks as background: 
overlap <- subsetByOverlaps(atacPeaks[atacPeaks$name %in% rownames(gainedPeaks)], gainedLoops)
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

write.table(overlap, paste0(outDirPath, "gainedLoop/overlap.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(background, paste0(outDirPath, "gainedLoop/background.txt"), quote = FALSE, sep = "\t", row.names = FALSE)

## Set 2: lost loops with lost atac peaks, all peaks as background: 
overlap <- subsetByOverlaps(atacPeaks[atacPeaks$name %in% rownames(lostPeaks)], lostLoops)
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

write.table(overlap, paste0(outDirPath, "lostLoopCorr/overlap.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(background, paste0(outDirPath, "lostLoopCorr/background.txt"), quote = FALSE, sep = "\t", row.names = FALSE)

## Set 3: lost loops with gained atac peaks, all peaks as background: 
overlap <- subsetByOverlaps(atacPeaks[atacPeaks$name %in% rownames(gainedPeaks)], lostLoops)
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

write.table(overlap, paste0(outDirPath, "lostLoopAnti/overlap.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(background, paste0(outDirPath, "lostLoopAnti/background.txt"), quote = FALSE, sep = "\t", row.names = FALSE)

#Define constants for final command: 
func <- "findMotifsGenome.pl "
overlapPath <- paste0(path, "/input/")
g <- " hg38 "
outPath <- paste0(path, "/output/")
size <- " -size given "
backgroundPath <- paste0("-bg ", path, "/input/")

#Command: gained loop gained atac all peaks background 
x <- "/gainedLoop/"
command <- paste0(func, 
                  paste0(overlapPath, x, "overlap.txt"), 
                  g, 
                  paste0(outPath, x), 
                  size, 
                  paste0(backgroundPath, x, "background.txt"))

write(command, paste0(path, "/commands/gainedLoop.sh"))

#Command: lost loop lost atac all peaks background 
x <- "/lostLoopCorr/"
command <- paste0(func, 
                  paste0(overlapPath, x, "overlap.txt"), 
                  g, 
                  paste0(outPath, x), 
                  size, 
                  paste0(backgroundPath, x, "background.txt"))

write(command, paste0(path, "/commands/lostLoopCorr.sh"))

#Command: lost loop gained atac all peaks background 
x <- "/lostLoopAnti/"
command <- paste0(func, 
                  paste0(overlapPath, x, "overlap.txt"), 
                  g, 
                  paste0(outPath, x), 
                  size, 
                  paste0(backgroundPath, x, "background.txt"))

write(command, paste0(path, "/commands/lostLoopAnti.sh"))


## Visualize: 
all <- read.delim("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/atac/allMotifEnrichment/output/allLoopATAC/knownResults.txt")
gained <- read.delim("~/Desktop/MEGA_newMotifEnrich/output/gainedLoop/knownResults.txt")
lost <- read.delim("~/Desktop/MEGA_newMotifEnrich/output/lostLoopCorr/knownResults.txt")

motifs <- list(all, gained, lost)

names <- c("all", "gained", "lost")

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

#PHASE 6: Make barplots 
#cols <- brewer.pal(8, "YlGnBu")[4:8]
cols <- c("gray50", "#227B7F", "#931C1E")

pdf("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig2/motifEnrichBarplotsScale.pdf", width = 8.5, height = 11)
plots <- list()
for (i in 1:length(motifs)){
  motif <- motifs[[i]]
  motif <- motif[1:5,]
  motif$MotifName <- factor(motif$MotifName, levels = rev(motif$MotifName))
  plot <- 
    ggplot(data = motif, aes(x = MotifName, y = log10pval)) + 
    geom_bar(stat = "identity", fill = cols[i]) + coord_flip() +
    theme_classic() + ylim(0,60)
  plots[[i]] <- plot
  #pdf((paste0(dirPath, "plots", path[i], "knownMotifs_Enrichment.pdf")))
  #print(plot)
  #dev.off()
}

do.call("grid.arrange", c(plots, ncol=1))
dev.off()






