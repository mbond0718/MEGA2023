## BINNED TRANSCRIPTION (TPM())
## Instead of using counts per gene, look at transcription per bin of the genome 

## Load libraries -----------------------------------------
library(GenomicRanges)
library(DESeq2)
library(GenomicInteractions)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggplot2)

## Read in data -----------------------------------------
## Transcription data
txi <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/raw/rna/MEGA_K562_WT_PMA_S_txi.rds")
tpm <- data.frame(txi$abundance)
#Trim the transcript ID from the rownames of tpm
rownames(tpm) <- (unlist(strsplit(rownames(tpm), "[.]")))[seq(1,(2*nrow(tpm)),by=2)]
tpm$ensembl_gene_id <- rownames(tpm)
## Load txdb 
txdb <- loadDb("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/ref/hg38txdb.sqlite") |> genes() |> as.data.frame()
txdb$seqnames <- paste0("chr", txdb$seqnames)
txdb <- txdb[txdb$gene_id %in% tpm$ensembl_gene_id,]
tpm <- tpm[tpm$ensembl_gene_id %in% txdb$gene_id,]
colnames(txdb)[6] <- "ensembl_gene_id"
#Merge dataframes
tpm <- merge(tpm, txdb, by="ensembl_gene_id")
tpm <- tpm[order(tpm$seqnames),]
#Reorder columns
tpm <- tpm[,c(1,18:20,2:17)]

## Tile genome into bins with GRanges: 
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
binsize = 10E3
bins <- tileGenome(seqlengths = (seqlengths(txdb))[1:23], 
                   tilewidth = binsize+1, 
                   cut.last.tile.in.chrom = T)
bins$bin <- paste0("bin", 1:length(bins))


#Make GRanges for the genes that are present in the tpm dataframe 
genes <- GRanges(seqnames = Rle(tpm$seqnames),
                 ranges = IRanges(start = tpm$start, end = tpm$end), 
                 gene = tpm$ensembl_gene_id)
promoters <- promoters(genes)

#Intersect genes with bins: 
olap <- findOverlaps(genes, bins)
query <- queryHits(olap)
subject <- subjectHits(olap)

test <- pintersect(genes[query,], bins[subject,])
percentOverlap <- width(test)/width(genes[query,])

#Build a dataframe that contains the gene, its coordinates, the bin, and its coordinates: 
g <- (data.frame(genes[query,]))[,c(6,1:4)]
colnames(g) <- c("gene", paste0("gene_", colnames(g)[2:5]))
b <- (data.frame(bins[subject,]))[,c(6,1:3)]
colnames(b) <- c("bin", paste0("bin_", colnames(b)[2:4]))
d <- cbind(g,b)
d$adjust <- percentOverlap

#Transform tpm df to merge reps, just 0,6,72h
rna <- tpm[,c(5,13,9,17,12,20)]
rna$gene <- tpm$ensembl_gene_id
#Merge timepoints: 
rna$TP0 <- (rna$MEGA_RNA_K562_WT_PMA_0000_S_2.1.1+rna$MEGA_RNA_K562_WT_PMA_0000_S_1.1.1)/2
rna$TP6 <- (rna$MEGA_RNA_K562_WT_PMA_0360_S_2.1.1+rna$MEGA_RNA_K562_WT_PMA_0360_S_1.1.1)/2
rna$TP72 <- (rna$MEGA_RNA_K562_WT_PMA_4320_S_2.1.1+rna$MEGA_RNA_K562_WT_PMA_4320_S_1.1.1)/2
rna <- rna[,c(7:10)]

#Merge rna with the bin/gene dataframe
merged <- merge(d, rna, by="gene")
merged <- merged[order(merged$gene_seqnames),]

#Finally, normalize for the amount of the gene that's in each bin: 
merged$TP0_adj <- merged$TP0*merged$adjust
merged$TP6_adj <- merged$TP6*merged$adjust
merged$TP72_adj <- merged$TP72*merged$adjust
merged$mean_adj <- (merged$TP0_adj+merged$TP6_adj+merged$TP72_adj)/3

merged <- unique(merged)
saveRDS(merged, file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/rna/binnedTPM.rds")
