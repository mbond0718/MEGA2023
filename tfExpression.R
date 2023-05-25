## SPECIFIC TF EXPRESSION
## For selected motifs in Figure 2, this script looks at the expression of those genes from RNA-seq counts

## Load libraries -----------------------------------------
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(Homo.sapiens)
library(ggplot2)
library(gridExtra)

## Process -------------------------------------------------
#Generate txdb object for figuring out HGNC of genes as well as coords
load("~/Desktop/geneAnno_gencode_v33_hg38.rda")
genes <- geneAnno1
#Expression of genes: 
geneCounts <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/rna/rnaCounts.rds")
geneCounts$ensembl_gene_id <- rownames(geneCounts)
#Merge: 
data <- merge(genes, geneCounts, by="ensembl_gene_id")
#Genes of interest to plot: 
goi <- c("NFKB1", "FOXA1", "NRF1")

data <- data[data$hgnc_symbol %in% goi,]
data <- data[order(match(data$hgnc_symbol, goi)),]
#Generate list of ggplots
for (i in 1:length(goi)){
  x <- 
    matrix(data[i,9:16]) |> 
    data.frame()
  colnames(x) <- c("counts")
  x$counts <- as.numeric(x$counts)
  x$time <- factor(c(0, 0.5, 1.5, 3, 6, 24, 48, 72))
  x$gene <- data$hgnc_symbol[i]
  
  plot <- 
    ggplot(x, aes(x=time, y=counts, group = gene)) + 
    geom_line(color = "#9ECAE1") + 
    theme_classic() + 
    ggtitle(data$hgnc_symbol[i])
  print(plot)
  ggsave(last_plot(), filename = paste0("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig2/", data$hgnc_symbol[i], "Expression.pdf"))
}













