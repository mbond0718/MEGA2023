##SUPPLEMENTAL FIGURE 2: 
## Plot proximal according to distal TF enrichment and distal according to proximal TF enrichment 

pal <- brewer.pal(8, "YlGnBu")
path <- "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/atac/"
dirName <- "proxDistMotifEnrich"
proxDist <- "distal"
names <- c("UpEarly", "UpMid", "UpLate", "DnEarly", "DnMid", "DnLate")
plotMotifs <- function(path, dirName, listOfCategories, proxDist)

## Proximal:
directoryPath <- paste0(path, dirName)
proximal <- list()
for (i in 1:length(names)){
  proximal[[i]] <- read.delim(paste0(directoryPath, "/proximal/output/", names[i], "/knownResults.txt"))
}
  
names(proximal) <- names

#Phase 3: clean up each motif file
newcolnames <- c("MotifName", "Consensus", "P.value", "Log P.value", "q value", 
                 "Number target sequences", "Percent_target_sequences", 
                 "Number background sequences", "Percent_background_sequences")
abridge <- function(longname){
  temp <- strsplit(longname, '[(]')
  shortname <- unlist(temp)[1]
}
for (i in 1:length(proximal)){
  #Subset by the motif list
  motif <- proximal[[i]]
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
  motif$class <- names[[i]]
  
  #Reassign to the list 
  proximal[[i]] <- motif
}

##Alternate version, not by plotting specific motifs of interest: 
proximal <- do.call(rbind, proximal) |> as.data.frame()
proximal$class <- factor(proximal$class, names[length(names):1])
#Threshold for points with a significant log10p
threshold <- 1.3
proxOrder <- proximal[proximal$log10pval>threshold,]
order <- c()
for (i in 1:length(names)){
  #order[[i]] <- temp[temp$class==names[i] & temp$log10pval>threshold,]$MotifName |> unique()
  x <- proxOrder[proxOrder$class==names[i] & proxOrder$log10pval>threshold,]$MotifName |> unique()
  order[[i]] <- x[1:10]
  
}
proxOrder <- unlist(order) |> unique()
proxOrder <- proxOrder[!proxOrder %in% c("Pax7", "p53", "GATA3", "RAR:RXR", "PAX5", "GATA")]

## Distal:
distal <- list()
for (i in 1:length(names)){
  distal[[i]] <- read.delim(paste0(directoryPath, "/distal/output/", names[i], "/knownResults.txt"))
}

names(distal) <- names

#Phase 3: clean up each motif file
newcolnames <- c("MotifName", "Consensus", "P.value", "Log P.value", "q value", 
                 "Number target sequences", "Percent_target_sequences", 
                 "Number background sequences", "Percent_background_sequences")
abridge <- function(longname){
  temp <- strsplit(longname, '[(]')
  shortname <- unlist(temp)[1]
}
for (i in 1:length(distal)){
  #Subset by the motif list
  motif <- distal[[i]]
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
  motif$class <- names[[i]]
  
  #Reassign to the list 
  distal[[i]] <- motif
}

##Alternate version, not by plotting specific motifs of interest: 
distal <- do.call(rbind, distal) |> as.data.frame()
distal$class <- factor(distal$class, names[length(names):1])
#Threshold for points with a significant log10p
threshold <- 1.3
distalOrder <- distal[distal$log10pval>threshold,]
order <- c()
for (i in 1:length(names)){
  #order[[i]] <- temp[temp$class==names[i] & temp$log10pval>threshold,]$MotifName |> unique()
  x <- distalOrder[distalOrder$class==names[i] & distalOrder$log10pval>threshold,]$MotifName |> unique()
  order[[i]] <- x[1:10]
  
}
distalOrder <- unlist(order) |> unique()
distalOrder <- distalOrder[!distalOrder %in% c("Pax7", "p53", "GATA3", "RAR:RXR", "PAX5", "GATA")]

## Plot prox by dist: 
temp <- proximal[proximal$MotifName %in% distalOrder,]
temp$MotifName <- factor(temp$MotifName, distalOrder)

ggplot(temp, aes(x=MotifName, y=class, fill=log10pval)) + 
  geom_tile(color = "white") + scale_fill_gradientn(colors = brewer.pal(8, "Blues")) +
  coord_fixed() + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + 
  ggtitle("Proximal Elements")
ggsave(last_plot(), filename = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/sup2/proxByDist.pdf")

## Plot dist by prox: 
temp <- distal[distal$MotifName %in% proxOrder,]
temp$MotifName <- factor(temp$MotifName, proxOrder)

ggplot(temp, aes(x=MotifName, y=class, fill=log10pval)) + 
  geom_tile(color = "white") + scale_fill_gradientn(colors = brewer.pal(8, "Blues")) +
  coord_fixed() + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + 
  ggtitle("Distal Elements")
ggsave(last_plot(), filename = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/sup2/distByProx.pdf")













