## Manipulate the E-P pairs from Gasperini 2019

library(readxl)
library(GenomicRanges)
library(liftOver)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(mariner)


## Read in data
pairs <- read_xlsx("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/PerturbSeq/enh-genepairs_hg19.xlsx") |> 
  as.data.frame()

## Reformat 
enhancers <- GRanges(seqnames = Rle(pairs$chr.candidate_enhancer), 
                     ranges = IRanges(start = pairs$start.candidate_enhancer, end = pairs$stop.candidate_enhancer))
#add mcols
enhancers$ENSG <- pairs$ENSG
enhancers$target_gene_short <- pairs$target_gene_short
enhancers$confidence <- pairs$high_confidence_subset

## LiftOver
#Chain
path = system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")
ch = import.chain(path)
ch

## Lift Over: 
enh38 <- liftOver(enhancers, ch) |> unlist() |> unique()

## Read in gene information: 
genes <- loadDb("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/ref/hg38txdb.sqlite") |> genes() |> data.frame()
genes$seqnames <- paste0("chr", genes$seqnames)
rownames(genes) <- genes$gene_id
genes <- 
  GRanges(seqnames = Rle(genes$seqnames), 
        ranges = IRanges(start = genes$start, end = genes$end), 
        ENSG = genes$gene_id)
promoters <- promoters(genes)

## Get just the genes that are in the lifted over enhancers: 
promoters <- promoters[promoters$ENSG %in% enh38$ENSG,]

enh38 <- as.data.frame(enh38)
promoters <- as.data.frame(promoters)
merged <- merge(enh38, promoters, by = "ENSG")
hg38Pairs <- merged[,c(2:4,9:11,1,7,8)]
hg38Pairs <- as_ginteractions(hg38Pairs) #anchor 1 is enhancer, anchor 2 is promoter 

## Look at the sizes of these loops: 
hg38Pairs$dist <- pairdist(hg38Pairs)
plot(density(hg38Pairs$dist))

## Filter out the shorter range loops: 
hg38Pairs <- hg38Pairs[hg38Pairs$dist>40000]

## Read in loops: 
loops <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/diffLoops.rds")
loops$loop <- rownames(loops)
loops <- as_ginteractions(loops)

## Overlap enhancers from CRISPRi with our loops: 
subsetByOverlaps(hg38Pairs, loops, maxgap = 30E3) |> length()
subsetByOverlaps(hg38Pairs[hg38Pairs$confidence==TRUE], loops, maxgap = 30E3) |> length() #high conf only
subsetByOverlaps(hg38Pairs[hg38Pairs$confidence==FALSE], loops, maxgap = 30E3) |> length() #low conf only 

##Percentage of perturb seq pairs that overlap a MEGA loop: 
((subsetByOverlaps(hg38Pairs, loops, maxgap = 30E3) |> length()) / (length(hg38Pairs))) * 100 


## ALTERNATIVE APPROACH: TRY WITH FULCO ET AL 2019 AKA ABC: 
fulco <- read_xlsx("/Users/mariellebond/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/CRISPRi/41588_2019_538_MOESM3_ESM.xlsx", 
                   sheet = 16, skip = 1) |> 
  as.data.frame()
colnames(fulco)[23] <- "ABC"
fulco <- GRanges(seqnames = Rle(fulco$chr), 
                 ranges = IRanges(start = fulco$start, end = fulco$end), 
                 gene = fulco$Gene, 
                 class = fulco$class, 
                 significant = fulco$Significant, 
                 abc = fulco$ABC)
#lift over: 
fulco38 <- liftOver(fulco, ch) |> unlist() |> unique()

#Overlap: 
subsetByOverlaps(fulco38, loops, maxgap = 10e3) #all with 1 bin wiggle
subsetByOverlaps(fulco38[fulco38$significant==TRUE,], loops, maxgap = 10e3) #only significant hits: 

(subsetByOverlaps(fulco38[fulco38$significant==TRUE,], loops, maxgap = 10e3) |> length())/(length(fulco38[fulco38$significant==TRUE]))*100
